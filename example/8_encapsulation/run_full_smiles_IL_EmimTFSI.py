#!/usr/bin/env python
"""

To run an ionic liquid: one type cation, one type of anion, zero solvents.
Every ionic species has to be declared by SMILES, because name of ionic species are not available in utils/mol_inventory.py inventory.

Element coverage still matters: ByteFF2 was trained on C, F, Li, N, O, P, S
(and H), so anions such as BF4 or DFOB are out of reach entirely, while TFSI,
FSI and PF6 at least stay inside the trained element set. The imidazolium
cations below are C/N/H, but no imidazolium sits in the training inventory, so
their parameters are extrapolated and unvalidated.

Every ion carries +1 or -1, so the box is neutral when the cations and the
anions each total ion_pairs molecules: 34 EMIM to 34 TFSI

Interface
---------
Calculator : IonicLiquidPropertiesCalculator (default choice for system with no solvent)
Inputs  : cation=["EMIM"], 
          anion=["TFSI"],  
          ion_pairs=34, smiles={name: SMILES for every ion}  (all by SMILES)
Output  : ./run_full_smiles_IL_EmimTFSI/  (results.json + per-property subdirectories)
Properties calculated:
    density, conductivity, viscosity
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


from byteff2.toolkit.common import PROPERTIES, save
from byteff2.toolkit.ionic_liquid_properties_calculator import IonicLiquidPropertiesCalculator

BASE_DIR = "./run_full_smiles_EmimTFSI_IL"

calc = IonicLiquidPropertiesCalculator(
    cation=["EMIM"],
    anion=["TFSI"],
    ion_pairs=34,
    smiles={
        "EMIM": "CCn1cc[n+](C)c1",
        "TFSI": "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F"
    },
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    save(calc.calculate(properties=["density", "viscosity" ,"conductivity"]), BASE_DIR)
