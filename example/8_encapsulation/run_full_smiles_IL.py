#!/usr/bin/env python
"""Run all by SMILES -- EXPERIMENTAL ionic liquid variant

A neat ionic liquid: one or more cations, one or more anions, zero solvents.
Every ion is declared by SMILES, the way run_full_smiles.py declares its
species -- there is no name lookup against the ByteFF2 inventory here.

Element coverage still matters: ByteFF2 was trained on C, F, Li, N, O, P, S
(and H), so anions such as BF4 or DFOB are out of reach entirely, while TFSI,
FSI and PF6 at least stay inside the trained element set. The imidazolium
cations below are C/N/H, but no imidazolium sits in the training inventory, so
their parameters are extrapolated and unvalidated.

Every ion carries +1 or -1, so the box is neutral when the cations and the
anions each total ion_pairs molecules: 17 EMIM + 17 BMIM against 17 TFSI +
17 FSI below.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.common import PROPERTIES, save
from byteff2.toolkit.ionic_liquid_properties_calculator import IonicLiquidPropertiesCalculator

BASE_DIR = "./run_full_smiles_IL"

calc = IonicLiquidPropertiesCalculator(
    cation=["EMIM", "Li","BMIM"],
    cation_ratio="5:4:1",
    anion=["TFSI", "FSI","PF6"],
    anion_ratio="5:4:1",
    ion_pairs=34,
    smiles={
        "EMIM": "CCn1cc[n+](C)c1",
        "BMIM": "CCCCn1cc[n+](C)c1",
        "Li":   "[Li+]",
        "TFSI": "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F",
        "FSI": "[N-](S(=O)(=O)F)S(=O)(=O)F",
        'PF6': "[P-:1]([F:2])([F:3])([F:4])([F:5])([F:6])[F:7]"
    },
    base_dir=BASE_DIR,
)

if __name__ == "__main__":
    save(calc.calculate(properties=["density", "conductivity"]), BASE_DIR)
