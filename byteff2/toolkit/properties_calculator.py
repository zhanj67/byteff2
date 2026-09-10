# Copyright (c) 2025 Bytedance Ltd. and/or its affiliates

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

from rdkit import Chem

from byteff2.toolkit.protocol import (
    DensityProtocol,
    DielectricProtocol,
    TransportProtocol,
    CompressibilityProtocol,
)
from bytemol.utils import setup_default_logging

logger = setup_default_logging()


class PropertiesLogger:
    """Enhanced logging for PropertiesCalculator with progress tracking."""

    def __init__(self, base_dir: Path, verbose: bool = True):
        self.base_dir = Path(base_dir)
        self.log_file = self.base_dir / "properties_calculator.log"
        self.verbose = verbose
        self.progress = 0

    def setup(self):
        """Setup file and console logging."""
        self.base_dir.mkdir(parents=True, exist_ok=True)

        # File handler
        file_handler = logging.FileHandler(str(self.log_file), mode='a')
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    def log_section(self, title: str, progress: int = 0):
        """Log a section header with progress."""
        self.progress = progress
        line = "=" * 80
        logger.info(line)
        logger.info(f"{title} [{progress}%]")
        logger.info(line)
        if self.verbose:
            print(f"\n{line}\n{title} [{progress}%]\n{line}")

    def log_step(self, step: str, details: str = ""):
        """Log a workflow step."""
        msg = f"  ✓ {step}"
        if details:
            msg += f": {details}"
        logger.info(msg)
        if self.verbose:
            print(msg)

    def log_result(self, name: str, value: float, unit: str = ""):
        """Log a calculated result."""
        msg = f"    {name:30s}: {value:12.6f} {unit}".rstrip()
        logger.info(msg)
        if self.verbose:
            print(msg)

    def log_error(self, msg: str):
        """Log an error."""
        logger.error(f"✗ {msg}")
        if self.verbose:
            print(f"✗ {msg}", file=sys.stderr)

    def log_summary(self, summary_dict: Dict):
        """Log a formatted summary."""
        line = "=" * 80
        logger.info(line)
        logger.info("CALCULATION SUMMARY")
        logger.info(line)
        for key, value in summary_dict.items():
            if isinstance(value, dict):
                logger.info(f"\n{key}:")
                for k, v in value.items():
                    logger.info(f"  {k:30s}: {v}")
            else:
                logger.info(f"{key:30s}: {value}")
        if self.verbose:
            print("\n" + line)
            print("CALCULATION SUMMARY")
            print(line)
            for key, value in summary_dict.items():
                if isinstance(value, dict):
                    print(f"\n{key}:")
                    for k, v in value.items():
                        print(f"  {k:30s}: {v}")
                else:
                    print(f"{key:30s}: {value}")


class PropertiesCalculator:
    """
    Unified interface for calculating molecular dynamics properties.

    Simplified input format:
    - solvent_name: Name of solvent (e.g., "DMC", "EC")
    - anion_name: Name of anion (e.g., "PF6", "TFSI")
    - ion_count: Number of ion pairs (default: auto-calculated from natoms)

    Default settings:
    - natoms: 5000
    - cation: LI (lithium)

    Supports:
    - Density calculation
    - Conductivity calculation
    - Viscosity calculation
    - Dielectric constant calculation
    - Compressibility calculation

    Can be called programmatically or via CLI.
    """

    # Fixed cation (always lithium)
    CATION = "LI"
    CATION_SMILES = "[Li+]"

    # Anion SMILES database
    ANION_DATABASE = {
        "PF6": "F[P-](F)(F)(F)(F)F",
        "BF4": "F[B-](F)(F)F",
        "ClO4": "[O-][Cl](=O)(=O)=O",
        "TFSI": "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F",
        "OTf": "[O-]S(=O)(=O)C(F)(F)F",
        "FSI": "FS(=O)(=O)[N-]S(=O)(=O)F",
    }

    # Solvent SMILES. Entries below are transcribed from Supplementary Table 1
    # ("Molecule Abbreviations") of the ByteFF2 paper SI, which is the naming
    # authority for this force field's training set.
    SOLVENT_DATABASE = {
        # carbonates
        "EC": "C1COC(=O)O1",                    # ethylene carbonate
        "DMC": "COC(=O)OC",                     # dimethyl carbonate
        "EMC": "CCOC(=O)OC",                    # ethyl methyl carbonate
        "PC": "CC1COC(=O)O1",                   # propylene carbonate
        "FEC": "C1C(OC(=O)O1)F",                # fluoroethylene carbonate
        "DFEC": "O1[C@H](F)[C@@H](F)OC1=O",     # difluoroethylene carbonate
        "TFPC": "C1C(OC(=O)O1)C(F)(F)F",        # 4-(trifluoromethyl)-1,3-dioxolan-2-one
        "FEMC": "COC(=O)OCC(F)(F)F",            # methyl 2,2,2-trifluoroethyl carbonate
        # esters / lactones
        "MA": "CC(=O)OC",                       # methyl acetate
        "EA": "CCOC(=O)C",                      # ethyl acetate
        "GBL": "C1CC(=O)OC1",                   # gamma-butyrolactone
        "HAC": "CC(=O)O",                       # acetic acid
        # ethers
        "TGDME": "COCCOCCOCCOC",                # triethylene glycol dimethyl ether
        "EMP": "COCCCOCC",                      # 1-methoxy-3-ethoxypropane
        "F3EMP": "COCCCOCC(F)(F)F",             # 1-methoxy-3-(2,2,2-trifluoroethoxy)propane
        # phosphate
        "TFP": "O=P(OCC(F)(F)F)(OCC(F)(F)F)OCC(F)(F)F",  # tris(2,2,2-trifluoroethyl) phosphate
        # others
        "AN": "CC#N",                           # acetonitrile
        "ACE": "CC(C)=O",                       # acetone
        "NOM": "C[N+](=O)[O-]",                 # nitromethane
        "BZ": "c1ccccc1",                       # benzene
        "Ani": "Nc1ccccc1",                     # aniline
        "EtCl": "CCCl",                         # ethyl chloride
        "EtSH": "CCS",                          # ethanethiol
        # not in the SI table, but standard and widely used
        "DEC": "CCOC(=O)OCC",                   # diethyl carbonate
        "H2O": "O",                             # water
    }

    def __init__(
        self,
        solvent: str,
        anion: str,
        li_count: int = 34,
        salt_to_solvent_ratio: Optional[float] = None,
        salt_to_solvent_ratio_str: Optional[str] = None,
        temperature: float = 298.0,
        base_dir: str = "./md_simulations",
        verbose: bool = True,
    ):
        """
        Initialize the properties calculator with simplified input.

        Parameters
        ----------
        solvent : str
            Name of solvent (e.g., "DMC", "EC")
        anion : str
            Name of anion (e.g., "PF6", "TFSI")
        li_count : int
            Number of Li ions (default: 34)
        salt_to_solvent_ratio : float, optional
            Ratio of salt pairs to solvent molecules (e.g., 1/10 or 0.1).
            If provided, solvent count = li_count * ratio (rounded down).
        salt_to_solvent_ratio_str : str, optional
            Alternative format: "1:10" means 1 salt pair : 10 solvent molecules.
        temperature : float
            Temperature in Kelvin (default: 298.0)
        base_dir : str
            Base directory for all simulations (default: "./md_simulations")
        """
        self.solvent = solvent
        self.anion = anion
        self.li_count = li_count
        self.temperature = temperature
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose

        # Setup logging
        self.logger = PropertiesLogger(self.base_dir, verbose=verbose)
        self.logger.setup()

        # Log initialization
        self.logger.log_section("PropertiesCalculator Initialization", progress=0)
        self.logger.log_step(f"Solvent: {solvent}, Anion: {anion}, Li count: {li_count}")

        # Validate and get SMILES
        self.solvent_smiles = self._get_solvent_smiles()
        self.anion_smiles = self._get_anion_smiles()

        # Parse salt:solvent ratio
        self.salt_to_solvent_ratio = self._parse_ratio(salt_to_solvent_ratio, salt_to_solvent_ratio_str)

        # Build components dictionary
        self.components = self._build_components()

        # Build SMILES dictionary
        self.smiles = self._build_smiles()

        self.logger.log_step(f"Components: {self.components}")
        self.logger.log_step(f"Working directory: {self.base_dir}")

        self.results = {}

    def _get_solvent_smiles(self) -> str:
        """Get SMILES for solvent."""
        if self.solvent not in self.SOLVENT_DATABASE:
            raise ValueError(
                f"Unknown solvent: {self.solvent}. "
                f"Known solvents: {list(self.SOLVENT_DATABASE.keys())}"
            )
        return self.SOLVENT_DATABASE[self.solvent]

    def _get_anion_smiles(self) -> str:
        """Get SMILES for anion."""
        if self.anion not in self.ANION_DATABASE:
            raise ValueError(
                f"Unknown anion: {self.anion}. "
                f"Known anions: {list(self.ANION_DATABASE.keys())}"
            )
        return self.ANION_DATABASE[self.anion]

    def _parse_ratio(self, ratio_float: Optional[float], ratio_str: Optional[str]) -> float:
        """Parse salt:solvent ratio from float or string format."""
        if ratio_str:
            # Parse "1:10" format
            parts = ratio_str.split(":")
            if len(parts) == 2:
                try:
                    return float(parts[0]) / float(parts[1])
                except (ValueError, ZeroDivisionError):
                    raise ValueError(f"Invalid ratio format: {ratio_str}")

        if ratio_float is not None:
            return ratio_float

        # Default: 1:10 (1 salt pair : 10 solvent molecules)
        return 1.0 / 10.0

    def _build_components(self) -> Dict[str, int]:
        """Build components dictionary from salt:solvent ratio."""
        # Calculate solvent count from ratio
        # salt_to_solvent_ratio = salt_pairs / solvent_molecules
        # solvent_count = li_count / ratio (rounded down)
        solvent_count = int(self.li_count / self.salt_to_solvent_ratio)

        components = {
            self.solvent: solvent_count,
            self.CATION: self.li_count,
            self.anion: self.li_count,
        }
        logger.info(f"Components: {components} (ratio {self.salt_to_solvent_ratio:.4f})")
        return components

    def _build_smiles(self) -> Dict[str, str]:
        """Build SMILES dictionary."""
        return {
            self.solvent: self.solvent_smiles,
            self.CATION: self.CATION_SMILES,
            self.anion: self.anion_smiles,
        }

    def _count_atoms(self) -> int:
        """Exact atom count of the requested composition, hydrogens included."""
        total = 0
        for name, count in self.components.items():
            mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles[name]))
            total += count * mol.GetNumAtoms()
        return total

    def _create_config(self, protocol_name: str, params_dir: str, output_dir: str) -> dict:
        """Create configuration for a protocol."""
        # protocol.search_mixture treats `components` as a ratio and scales it up
        # until it reaches `natoms`. Passing the composition's exact atom count
        # makes that scale factor land on 1, preserving the requested counts.
        config = {
            "protocol": protocol_name,
            "params_dir": params_dir,
            "output_dir": output_dir,
            "working_dir": str(self.base_dir / f"{protocol_name.lower()}_working_dir"),
            "temperature": int(self.temperature),
            "natoms": self._count_atoms(),
            "components": self.components,
            "smiles": self.smiles,
        }
        return config

    def _save_config(self, config: dict, config_file: str) -> None:
        """Save configuration to JSON file."""
        os.makedirs(os.path.dirname(config_file), exist_ok=True)
        with open(config_file, "w") as f:
            json.dump(config, f, indent=2)
        logger.info(f"Saved config to {config_file}")

    def _run_protocol(self, protocol_class, config: dict) -> dict:
        """Run a protocol and return results."""
        try:
            protocol = protocol_class(config)
            logger.info(f"Running {protocol_class.__name__}")
            protocol.run_protocol()
            protocol.post_process()

            # Read results from output file
            output_dir = config["output_dir"]
            results_file = os.path.join(output_dir, "results.json")

            # Try different result file names based on protocol
            if not os.path.exists(results_file):
                if protocol_class == DensityProtocol:
                    results_file = os.path.join(output_dir, "density_results.json")
                elif protocol_class == DielectricProtocol:
                    results_file = os.path.join(output_dir, "dielectric_results.json")

            if os.path.exists(results_file):
                with open(results_file, "r") as f:
                    return json.load(f)
            else:
                logger.warning(f"Results file not found for {protocol_class.__name__}")
                return {}
        except Exception as e:
            logger.error(f"Error running {protocol_class.__name__}: {e}")
            raise

    def calculate_density(self) -> dict:
        """Calculate density."""
        protocol_name = "Density"
        params_dir = str(self.base_dir / "density_params")
        output_dir = str(self.base_dir / "density_results")

        config = self._create_config(protocol_name, params_dir, output_dir)
        config_file = str(self.base_dir / "density_config.json")
        self._save_config(config, config_file)

        logger.info("Calculating density...")
        result = self._run_protocol(DensityProtocol, config)
        self.results["density"] = result
        return result

    def calculate_transport(self) -> dict:
        """Calculate conductivity and viscosity."""
        protocol_name = "Transport"
        params_dir = str(self.base_dir / "transport_params")
        output_dir = str(self.base_dir / "transport_results")

        config = self._create_config(protocol_name, params_dir, output_dir)
        config_file = str(self.base_dir / "transport_config.json")
        self._save_config(config, config_file)

        logger.info("Calculating conductivity and viscosity...")
        result = self._run_protocol(TransportProtocol, config)
        self.results["transport"] = result
        return result

    def calculate_dielectric(self) -> dict:
        """Calculate dielectric constant."""
        protocol_name = "Dielectric"
        params_dir = str(self.base_dir / "dielectric_params")
        output_dir = str(self.base_dir / "dielectric_results")

        config = self._create_config(protocol_name, params_dir, output_dir)
        config_file = str(self.base_dir / "dielectric_config.json")
        self._save_config(config, config_file)

        # Add dielectric-specific parameters
        config["npt_steps"] = 2000000
        config["nvt_steps"] = 6000000
        config["dipole_interval"] = 500

        logger.info("Calculating dielectric constant...")
        result = self._run_protocol(DielectricProtocol, config)
        self.results["dielectric"] = result
        return result

    def calculate_compressibility(self) -> dict:
        """Calculate compressibility."""
        protocol_name = "Compressibility"
        params_dir = str(self.base_dir / "compressibility_params")
        output_dir = str(self.base_dir / "compressibility_results")

        config = self._create_config(protocol_name, params_dir, output_dir)
        config_file = str(self.base_dir / "compressibility_config.json")
        self._save_config(config, config_file)

        # Add compressibility-specific parameters
        config["npt_steps"] = 5000000

        logger.info("Calculating compressibility...")
        result = self._run_protocol(CompressibilityProtocol, config)
        self.results["compressibility"] = result
        return result

    def calculate(
        self,
        properties: Optional[List[str]] = None,
    ) -> dict:
        """
        Calculate specified properties.

        Parameters
        ----------
        properties : List[str], optional
            List of properties to calculate. Options:
            ["density", "conductivity", "viscosity", "dielectric", "compressibility"]
            If None, calculates all.

        Returns
        -------
        dict
            Dictionary containing all calculated results.
        """
        if properties is None:
            properties = ["density", "conductivity", "viscosity", "dielectric"]

        properties_lower = [p.lower() for p in properties]
        all_results = {}

        self.logger.log_section(
            f"Starting Property Calculations: {', '.join(properties_lower)}",
            progress=10
        )

        # Density
        if "density" in properties_lower:
            try:
                self.logger.log_step("Calculating density...")
                all_results["density"] = self.calculate_density()
                if all_results["density"]:
                    density_val = all_results["density"].get("density", 0)
                    self.logger.log_result("Density", density_val, "g/mL")
            except Exception as e:
                self.logger.log_error(f"Failed to calculate density: {e}")

        # Transport (conductivity + viscosity)
        if any(p in properties_lower for p in ["conductivity", "viscosity"]):
            try:
                self.logger.log_step("Calculating transport properties...")
                transport_result = self.calculate_transport()
                all_results["transport"] = transport_result
                if "conductivity" in properties_lower:
                    cond = transport_result.get("conductivity_onsager", 0)
                    all_results["conductivity"] = cond
                    self.logger.log_result("Conductivity (Onsager)", cond, "mS/cm")
                if "viscosity" in properties_lower:
                    visc = transport_result.get("viscosity", 0)
                    all_results["viscosity"] = visc
                    self.logger.log_result("Viscosity", visc, "cP")
            except Exception as e:
                self.logger.log_error(f"Failed to calculate transport properties: {e}")

        # Dielectric
        if "dielectric" in properties_lower:
            try:
                self.logger.log_step("Calculating dielectric constant...")
                all_results["dielectric"] = self.calculate_dielectric()
                if all_results["dielectric"]:
                    diel_val = all_results["dielectric"].get("dielectric", 0)
                    self.logger.log_result("Dielectric Constant", diel_val)
            except Exception as e:
                self.logger.log_error(f"Failed to calculate dielectric constant: {e}")

        # Compressibility
        if "compressibility" in properties_lower:
            try:
                self.logger.log_step("Calculating compressibility...")
                all_results["compressibility"] = self.calculate_compressibility()
            except Exception as e:
                self.logger.log_error(f"Failed to calculate compressibility: {e}")

        # Save and display summary
        self.logger.log_section("Calculation Complete", progress=100)
        self._save_summary(all_results)
        self._display_summary(all_results)

        return all_results

    def _save_summary(self, results: dict) -> None:
        """Save results summary to JSON."""
        summary = {
            "solvent": self.solvent,
            "anion": self.anion,
            "cation": self.CATION,
            "li_count": self.li_count,
            "salt_to_solvent_ratio": self.salt_to_solvent_ratio,
            "components": self.components,
            "temperature": self.temperature,
            "results": results,
        }
        summary_file = self.base_dir / "summary.json"
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2)
        self.logger.log_step(f"Saved summary to {summary_file}")

    def _display_summary(self, results: dict) -> None:
        """Display formatted results summary."""
        summary_data = {
            "System Composition": {
                "Solvent": self.solvent,
                "Anion": self.anion,
                "Cation": self.CATION,
                "Li Count": self.li_count,
                "Salt:Solvent Ratio": f"1:{int(1/self.salt_to_solvent_ratio)}",
            },
            "Simulation Parameters": {
                "Temperature": f"{self.temperature} K",
                "Components": str(self.components),
            },
        }

        # Add calculated properties
        if results:
            props_data = {}
            if "density" in results and results["density"]:
                props_data["Density (g/mL)"] = f"{results['density'].get('density', 0):.6f}"
            if "conductivity" in results:
                props_data["Conductivity (mS/cm)"] = f"{results['conductivity']:.6f}"
            if "viscosity" in results:
                props_data["Viscosity (cP)"] = f"{results['viscosity']:.6f}"
            if "dielectric" in results and results["dielectric"]:
                props_data["Dielectric Constant"] = f"{results['dielectric'].get('dielectric', 0):.6f}"

            if props_data:
                summary_data["Calculated Properties"] = props_data

        self.logger.log_summary(summary_data)

    @classmethod
    def from_simple_config(
        cls,
        solvent: str,
        anion: str,
        li_count: int = 34,
        salt_to_solvent_ratio: Optional[float] = None,
        salt_to_solvent_ratio_str: Optional[str] = None,
        **kwargs
    ) -> "PropertiesCalculator":
        """
        Create calculator from simple parameters.

        Parameters
        ----------
        solvent : str
            Solvent name
        anion : str
            Anion name
        li_count : int
            Number of Li ions (default: 34)
        salt_to_solvent_ratio : float, optional
            Ratio of salt pairs to solvent molecules
        salt_to_solvent_ratio_str : str, optional
            Ratio in "1:10" format
        **kwargs : dict
            Additional parameters (temperature, base_dir)
        """
        return cls(
            solvent=solvent,
            anion=anion,
            li_count=li_count,
            salt_to_solvent_ratio=salt_to_solvent_ratio,
            salt_to_solvent_ratio_str=salt_to_solvent_ratio_str,
            **kwargs
        )


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Calculate MD properties (density, conductivity, viscosity, dielectric)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 1:10 LiFSI:DMC (34 Li, 340 DMC)
  python -m byteff2.toolkit.properties_calculator \\
      --solvent DMC --anion FSI --ratio 1:10

  # Custom Li count
  python -m byteff2.toolkit.properties_calculator \\
      --solvent EC --anion TFSI --li-count 50 --ratio 1:8

  # Only density
  python -m byteff2.toolkit.properties_calculator \\
      --solvent DMC --anion PF6 --ratio 1:10 \\
      --calculate density --work-dir ./density_only
        """
    )

    parser.add_argument(
        "--solvent",
        type=str,
        required=True,
        help="Solvent name (e.g., DMC, EC, H2O)",
    )
    parser.add_argument(
        "--anion",
        type=str,
        required=True,
        help="Anion name (e.g., PF6, TFSI, FSI)",
    )
    parser.add_argument(
        "--li-count",
        type=int,
        default=34,
        help="Number of Li ions (default: 34)",
    )
    parser.add_argument(
        "--ratio",
        type=str,
        default="1:10",
        help="Salt:solvent ratio (default: 1:10)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=298.0,
        help="Temperature in Kelvin (default: 298.0)",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default="./md_simulations",
        help="Working directory (default: ./md_simulations)",
    )
    parser.add_argument(
        "--calculate",
        type=str,
        default="density,conductivity,viscosity,dielectric",
        help="Comma-separated properties to calculate",
    )
    parser.add_argument(
        "--log",
        type=str,
        default="properties_calculator.log",
        help="Log file (default: properties_calculator.log)",
    )

    args = parser.parse_args()

    try:
        # Create calculator
        calculator = PropertiesCalculator(
            solvent=args.solvent,
            anion=args.anion,
            li_count=args.li_count,
            salt_to_solvent_ratio_str=args.ratio,
            temperature=args.temperature,
            base_dir=args.work_dir,
        )

        # Parse properties to calculate
        properties = [p.strip() for p in args.calculate.split(",")]

        # Run calculations
        results = calculator.calculate(properties=properties)

        # Print summary
        print("\n" + "="*60)
        print("CALCULATION SUMMARY")
        print("="*60)
        print(f"Solvent:  {args.solvent}")
        print(f"Anion:    {args.anion}")
        print(f"Cation:   LI (fixed)")
        print(f"Li count: {calculator.li_count}")
        print(f"Ratio:    1:{int(1/calculator.salt_to_solvent_ratio)}")
        print(f"Components: {calculator.components}")
        print(f"Temperature: {args.temperature} K")
        print("\nResults:")
        print(json.dumps(results, indent=2))
        print("="*60 + "\n")

        return 0

    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    exit(main())
