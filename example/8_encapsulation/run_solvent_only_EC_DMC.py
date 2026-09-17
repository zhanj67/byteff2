#!/usr/bin/env python
"""

Pure-solvent (zero-salt) run

"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.common import save
from byteff2.toolkit.properties_calculator import PropertiesCalculator

BASE_DIR = "./run_solvent_only_EC_DMC" 

calc = PropertiesCalculator(
    solvent=["EC", "DMC"],
    solvent_ratio="5:5",
    li_count=0,
    total_solvent=340,  # 170 EC + 170 DMC, matching run_smiles_EC_DMC's solvent counts exactly
    temperature=298.0,
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    # Conductivity is skipped from the requested/reported properties (meaningless
    # for a neutral system anyway).
    save(calc.calculate(properties=["density", "viscosity", "dielectric"]), BASE_DIR)
