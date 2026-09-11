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
"""EXPERIMENTAL: property calculations for solvent-free ionic liquids.

The salt-in-solvent workflow lives in
:mod:`byteff2.toolkit.properties_calculator`; this module mirrors it for a
neat ionic liquid -- one or more cations, one or more anions, zero solvents --
and runs the very same protocols on the result.

Two deliberate differences from PropertiesCalculator:

* Every species is declared by SMILES (``smiles={"EMIM": "CCn1cc[n+](C)c1"}``),
  the way ``example/8_encapsulation/run_full_smiles.py`` declares its species.
  Names are labels only -- there is no name lookup against the ByteFF2
  inventory, because that inventory is a salt-in-carbonate inventory and
  contains almost no IL species. The inventory is consulted by structure, for
  one purpose: an ion whose SMILES is not one the force field was trained on is
  warned about, since its parameters are extrapolated and unvalidated. LI, PF6,
  TFSI and FSI are in the inventory and pass without a warning.
* Every ion is assumed to carry a charge of exactly +1 or -1, which is checked
  against the SMILES. Charge neutrality then means the cations and the anions
  each total ``ion_pairs`` molecules.

Because these species sit outside the training inventory, their force field
parameters are extrapolated and unvalidated -- hence "experimental".
"""

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

from rdkit import Chem

from byteff2.toolkit.properties_calculator import (
    _ANION_DATABASE,
    _CATION_DATABASE,
    PropertiesLogger,
)
from byteff2.toolkit.protocol import (
    CompressibilityProtocol,
    DensityProtocol,
    DielectricProtocol,
    TransportProtocol,
)
from bytemol.utils import setup_default_logging

logger = setup_default_logging()


def _inventory_ions() -> Dict[str, str]:
    """Canonical SMILES -> inventory name, for every parameterised ion.

    Matching is by structure rather than by name: this module takes names purely
    as labels, so an ion counts as parameterised when its SMILES is one the
    force field was trained on, whatever the user chose to call it.
    """
    table = {}
    for name, smi in {**_CATION_DATABASE, **_ANION_DATABASE}.items():
        table[Chem.MolToSmiles(Chem.MolFromSmiles(smi))] = name
    return table


_INVENTORY_IONS = _inventory_ions()


class IonicLiquidPropertiesCalculator:
    """
    Unified interface for calculating MD properties of a neat ionic liquid.

    Input format:
    - cation: one cation name ("EMIM"), several names (["EMIM", "BMIM"]) with an
      optional cation_ratio ("5:5"), or a dict {"EMIM": 5, "BMIM": 5}
    - anion: same three forms, with anion_ratio
    - smiles: {name: SMILES} for every cation and anion named above
    - ion_pairs: number of cation-anion pairs; the cations share that total
      according to cation_ratio and the anions according to anion_ratio, so the
      box is charge neutral for any mixture of +1 cations and -1 anions.

    Supports:
    - Density calculation
    - Conductivity calculation
    - Viscosity calculation
    - Dielectric constant calculation
    - Compressibility calculation

    Can be called programmatically or via CLI.
    """

    def __init__(
        self,
        cation: Union[str, Sequence[str], Dict[str, float]],
        anion: Union[str, Sequence[str], Dict[str, float]],
        smiles: Dict[str, str],
        ion_pairs: int = 34,
        cation_ratio: Optional[Union[str, Sequence[float]]] = None,
        anion_ratio: Optional[Union[str, Sequence[float]]] = None,
        temperature: float = 298.0,
        base_dir: str = "./md_simulations",
        verbose: bool = True,
    ):
        """
        Initialize the ionic liquid properties calculator.

        Parameters
        ----------
        cation : str, list of str, or dict
            One cation name ("EMIM"), several names (["EMIM", "BMIM"]), or a
            dict mapping names to their molar proportions ({"EMIM": 5, "BMIM": 5}).
        anion : str, list of str, or dict
            Same, for the anions.
        smiles : dict
            {name: SMILES} covering every cation and anion. ByteFF2 predicts
            parameters from the molecular graph, so arbitrary ions do run -- but
            an ion outside the training inventory has extrapolated, unvalidated
            parameters, and is warned about when the calculator is built.
        ion_pairs : int
            Number of cation-anion pairs (default: 34). Both the cation total and
            the anion total equal this number.
        cation_ratio : str or list of float, optional
            Molar proportions between the cations, in the order given by
            `cation`: "5:5" or [5, 5]. Defaults to equal parts. Ignored when
            `cation` is a dict (the dict values are the proportions).
        anion_ratio : str or list of float, optional
            Same, for the anions.
        temperature : float
            Temperature in Kelvin (default: 298.0)
        base_dir : str
            Base directory for all simulations (default: "./md_simulations")
        """
        self.cations, self.cation_ratio = self._parse_species(cation, cation_ratio, "cation")
        self.anions, self.anion_ratio = self._parse_species(anion, anion_ratio, "anion")
        overlap = set(self.cations) & set(self.anions)
        if overlap:
            raise ValueError(f"Species named both as a cation and as an anion: {sorted(overlap)}")

        # Human-readable labels, e.g. "EMIM" or "EMIM/BMIM"
        self.cation = "/".join(self.cations)
        self.anion = "/".join(self.anions)
        self.ion_pairs = ion_pairs
        self.temperature = temperature
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose

        # Setup logging
        self.logger = PropertiesLogger(self.base_dir, verbose=verbose)
        self.logger.setup()

        # Log initialization
        self.logger.log_section("IonicLiquidPropertiesCalculator Initialization", progress=0)
        self.logger.log_step(
            f"Cation: {self._describe(self.cations, self.cation_ratio)}, "
            f"Anion: {self._describe(self.anions, self.anion_ratio)}, "
            f"Ion pairs: {ion_pairs}"
        )

        # Validate SMILES and ionic charges
        self.smiles = self._parse_smiles(smiles)

        # Build components dictionary
        self.components = self._build_components()

        self.logger.log_step(f"Components: {self.components}")
        self.logger.log_step(f"Working directory: {self.base_dir}")

        self.results = {}

    @staticmethod
    def _describe(names: List[str], ratio: List[float]) -> str:
        """Render "EMIM/BMIM (5:5)" for the log line."""
        desc = "/".join(names)
        if len(names) > 1:
            desc += " (" + ":".join(f"{r:g}" for r in ratio) + ")"
        return desc

    @staticmethod
    def _parse_species(
        species: Union[str, Sequence[str], Dict[str, float]],
        ratio: Optional[Union[str, Sequence[float]]],
        kind: str,
    ) -> tuple[List[str], List[float]]:
        """Normalise an ion spec into (names, molar proportions)."""
        if isinstance(species, dict):
            names = [str(s).strip() for s in species.keys()]
            proportions = [float(v) for v in species.values()]
        else:
            if isinstance(species, str):
                names = [s.strip() for s in species.split(",")]
            else:
                names = [str(s).strip() for s in species]
            if ratio is None:
                proportions = [1.0] * len(names)
            elif isinstance(ratio, str):
                try:
                    proportions = [float(x) for x in ratio.split(":")]
                except ValueError:
                    raise ValueError(f"Invalid {kind}_ratio format: {ratio!r} (expected e.g. '5:5')")
            else:
                proportions = [float(x) for x in ratio]

        if not names or any(not n for n in names):
            raise ValueError(f"Empty {kind} name in {species!r}")
        if len(set(names)) != len(names):
            raise ValueError(f"Duplicate {kind} names in {names}")
        if len(proportions) != len(names):
            raise ValueError(
                f"{kind}_ratio has {len(proportions)} entries but {len(names)} {kind}s were given: {names}"
            )
        if any(r <= 0 for r in proportions):
            raise ValueError(f"{kind}_ratio entries must be positive, got {proportions}")
        return names, proportions

    def _parse_smiles(self, smiles: Optional[Dict[str, str]]) -> Dict[str, str]:
        """Canonicalise the SMILES and check every ion carries the right charge."""
        smiles = smiles or {}
        expected = {name: +1 for name in self.cations}
        expected.update({name: -1 for name in self.anions})

        missing = [name for name in expected if name not in smiles]
        if missing:
            raise ValueError(
                f"No SMILES given for {missing}. Ionic liquids are submitted by SMILES: "
                f"pass smiles={{{missing[0]!r}: '<SMILES>', ...}} covering every cation and anion."
            )
        extra = [name for name in smiles if name not in expected]
        if extra:
            raise ValueError(f"smiles lists {extra}, which are neither cations nor anions of this system")

        parsed = {}
        for name, charge in expected.items():
            smi = smiles[name]
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError(f"smiles[{name!r}] is not valid SMILES: {smi!r}")
            net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            if net_charge != charge:
                kind = "cation" if charge > 0 else "anion"
                raise ValueError(
                    f"{name} was given as a {kind}, so it must carry a net charge of {charge:+d}, "
                    f"but its SMILES {smi!r} has {net_charge:+d}. Only singly charged ions are supported."
                )
            parsed[name] = Chem.MolToSmiles(mol)
            inventory_name = _INVENTORY_IONS.get(parsed[name])
            if inventory_name is None:
                logger.warning(
                    f"{name} is outside the ByteFF2 training inventory; its force field "
                    f"parameters will be extrapolated and are unvalidated."
                )
            else:
                logger.info(
                    f"{name} is the inventory species {inventory_name}; "
                    f"its force field parameters are the validated ByteFF2 ones."
                )
        return parsed

    @staticmethod
    def _apportion(total: int, weights: Sequence[float]) -> List[int]:
        """Split `total` into integers proportional to `weights` (largest remainder)."""
        wsum = float(sum(weights))
        exact = [total * w / wsum for w in weights]
        counts = [int(x) for x in exact]
        remainder = total - sum(counts)
        for i in sorted(range(len(weights)), key=lambda i: exact[i] - counts[i], reverse=True)[:remainder]:
            counts[i] += 1
        return counts

    def _build_components(self) -> Dict[str, int]:
        """Split ion_pairs over the cations and (separately) over the anions.

        Every ion is singly charged, so giving both sides the same total keeps
        the box neutral whatever the within-side proportions are.
        """
        if self.ion_pairs < max(len(self.cations), len(self.anions)):
            raise ValueError(
                f"ion_pairs={self.ion_pairs} is too small for {len(self.cations)} cations "
                f"and {len(self.anions)} anions"
            )
        cation_counts = self._apportion(self.ion_pairs, self.cation_ratio)
        anion_counts = self._apportion(self.ion_pairs, self.anion_ratio)
        for names, counts, ratio, kind in (
            (self.cations, cation_counts, self.cation_ratio, "cation"),
            (self.anions, anion_counts, self.anion_ratio, "anion"),
        ):
            if any(c < 1 for c in counts):
                raise ValueError(
                    f"{kind} ratio {ratio} leaves {kind}s with zero molecules out of "
                    f"{self.ion_pairs} ion pairs; increase ion_pairs or adjust the ratio"
                )

        components = dict(zip(self.cations, cation_counts))
        components.update(zip(self.anions, anion_counts))
        logger.info(
            f"Components: {components} (ion pairs {self.ion_pairs}, "
            f"cation ratio {':'.join(f'{r:g}' for r in self.cation_ratio)}, "
            f"anion ratio {':'.join(f'{r:g}' for r in self.anion_ratio)})"
        )
        return components

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
            "cation": self.cation,
            "cations": self.cations,
            "cation_ratio": self.cation_ratio,
            "anion": self.anion,
            "anions": self.anions,
            "anion_ratio": self.anion_ratio,
            "ion_pairs": self.ion_pairs,
            "components": self.components,
            "temperature": self.temperature,
            "smiles": self.smiles,
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
                "Cation": self.cation,
                "Cation Ratio": ":".join(f"{r:g}" for r in self.cation_ratio),
                "Anion": self.anion,
                "Anion Ratio": ":".join(f"{r:g}" for r in self.anion_ratio),
                "Ion Pairs": self.ion_pairs,
                "Solvent": "none (neat ionic liquid)",
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


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Calculate MD properties of a neat ionic liquid "
                    "(density, conductivity, viscosity, dielectric)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 34 pairs of neat EMIM TFSI
  python -m byteff2.toolkit.ionic_liquid_properties_calculator \\
      --cation EMIM --anion TFSI \\
      --smiles "EMIM=CCn1cc[n+](C)c1" \\
      --smiles "TFSI=[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F"

  # Mixed cations, mixed anions (17 EMIM + 17 BMIM, 17 TFSI + 17 BF4)
  python -m byteff2.toolkit.ionic_liquid_properties_calculator \\
      --cation EMIM,BMIM --cation-ratio 5:5 \\
      --anion TFSI,BF4 --anion-ratio 5:5 --ion-pairs 34 \\
      --smiles "EMIM=CCn1cc[n+](C)c1" --smiles "BMIM=CCCCn1cc[n+](C)c1" \\
      --smiles "TFSI=[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F" \\
      --smiles "BF4=[B-](F)(F)(F)F"
        """
    )

    parser.add_argument(
        "--cation",
        type=str,
        required=True,
        help="Cation name, or comma-separated names for a mixture (e.g., EMIM or EMIM,BMIM)",
    )
    parser.add_argument(
        "--cation-ratio",
        type=str,
        default=None,
        help="Molar proportions between cations, e.g. 5:5 (default: equal parts)",
    )
    parser.add_argument(
        "--anion",
        type=str,
        required=True,
        help="Anion name, or comma-separated names for a mixture (e.g., TFSI or TFSI,BF4)",
    )
    parser.add_argument(
        "--anion-ratio",
        type=str,
        default=None,
        help="Molar proportions between anions, e.g. 5:5 (default: equal parts)",
    )
    parser.add_argument(
        "--ion-pairs",
        type=int,
        default=34,
        help="Number of cation-anion pairs (default: 34)",
    )
    parser.add_argument(
        "--smiles",
        type=str,
        action="append",
        default=None,
        metavar="NAME=SMILES",
        required=True,
        help="SMILES for one ion, e.g. EMIM=CCn1cc[n+](C)c1. "
             "Repeat the flag for every cation and anion.",
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

    args = parser.parse_args()

    smiles = {}
    for entry in args.smiles or []:
        name, sep, smi = entry.partition("=")
        if not sep or not name.strip() or not smi.strip():
            parser.error(f"--smiles expects NAME=SMILES, got {entry!r}")
        smiles[name.strip()] = smi.strip()

    try:
        # Create calculator
        calculator = IonicLiquidPropertiesCalculator(
            cation=args.cation,
            anion=args.anion,
            smiles=smiles,
            ion_pairs=args.ion_pairs,
            cation_ratio=args.cation_ratio,
            anion_ratio=args.anion_ratio,
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
        print(f"Cation:   {calculator.cation}")
        if len(calculator.cations) > 1:
            print(f"Cation ratio: {':'.join(f'{r:g}' for r in calculator.cation_ratio)}")
        print(f"Anion:    {calculator.anion}")
        if len(calculator.anions) > 1:
            print(f"Anion ratio: {':'.join(f'{r:g}' for r in calculator.anion_ratio)}")
        print(f"Solvent:  none (neat ionic liquid)")
        print(f"Ion pairs: {calculator.ion_pairs}")
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
