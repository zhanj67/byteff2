#!/usr/bin/env python
"""

Pure-solvent run for sulfolane solvent, 
Sulfolane is not in the ByteFF2 inventory, so it is declared through custom_smiles.

"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.common import save
from byteff2.toolkit.properties_calculator import PropertiesCalculator

BASE_DIR = "./run_solvent_only_SL"  

calc = PropertiesCalculator(
    solvent=["SL"],
    li_count=0,
    total_solvent=340,  
    temperature=298.0,
    custom_smiles={"SL": "O=S1(=O)CCCC1"},
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    # Conductivity is skipped from the requested/reported properties (meaningless
    # for a neutral system anyway). Note this doesn't skip any compute: "viscosity"
    # alone still runs the same single TransportProtocol/NEMD job that would also
    # produce conductivity_onsager -- it's just not requested/promoted here.
    save(calc.calculate(properties=["density", "viscosity", "dielectric"]), BASE_DIR)
