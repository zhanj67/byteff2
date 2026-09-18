#!/usr/bin/env python
"""Sulfolane solvent and LiTfNCN anion in SMILES: 1:3 LiTfNCN in SL.

Both Sulfolane and TfNCN anion are not in the ByteFF2 inventory, so it is declared through custom_smiles.
ByteFF2 predicts parameters from the molecular graph and will run it, but the
molecule sits outside the training set, so its parameters are extrapolated.

Interface
---------
Inputs  : solvent="SL", custom_smiles={"SL": "C1CCS(=O)(=O)C1"}  (by SMILES)
          anion="N#C[N-]S(=O)(=O)C(F)(F)F"                       (by SMILES)
          li_count=60, salt_to_solvent_ratio_str="1:3"
Output  : ./run_full_electrolyte_by_smiles/  (results.json + per-property subdirectories)
Properties calculated (PROPERTIES):
    density, viscosity, conductivity
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.common import PROPERTIES, save
from byteff2.toolkit.properties_calculator import PropertiesCalculator

BASE_DIR = "./run_full_electrolyte_by_smiles"

calc = PropertiesCalculator(
    solvent="SL",
    anion="TfN",
    li_count=60,
    salt_to_solvent_ratio_str="1:3",
    custom_smiles={"SL": "O=S1(=O)CCCC1","TfN":"N#C[N-]S(=O)(=O)C(F)(F)F"},
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    save(calc.calculate(properties=["density", "viscosity", "conductivity"]), BASE_DIR)

