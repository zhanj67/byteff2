#!/usr/bin/env python
"""Solvent by SMILES: 1:10 LiPF6 in DMSO.

DMSO is not in the ByteFF2 inventory, so it is declared through custom_smiles.
ByteFF2 predicts parameters from the molecular graph and will run it, but the
molecule sits outside the training set, so its parameters are extrapolated.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.common import PROPERTIES, save
from byteff2.toolkit.properties_calculator import PropertiesCalculator

BASE_DIR = "./run_by_smiles"

calc = PropertiesCalculator(
    solvent="DMSO",
    anion="PF6",
    li_count=34,
    salt_to_solvent_ratio_str="1:10",
    custom_smiles={"DMSO": "CS(C)=O"},
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    save(calc.calculate(properties=PROPERTIES), BASE_DIR)
