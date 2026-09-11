#!/usr/bin/env python
"""Run by inventory names: 1:10 LiPF6 in EC:DMC = 3:7. 
PF6, EC, DMC are in the ByteFF2 inventory, therefore can use their name directly, no need to manually input SMILES
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.common import PROPERTIES, save
from byteff2.toolkit.properties_calculator import PropertiesCalculator

BASE_DIR = "./run_by_name"

calc = PropertiesCalculator(
    solvent=["EC", "DMC"],
    solvent_ratio="3:7",
    anion="PF6",
    li_count=34,
    salt_to_solvent_ratio_str="1:10",
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    save(calc.calculate(properties=PROPERTIES), BASE_DIR)
