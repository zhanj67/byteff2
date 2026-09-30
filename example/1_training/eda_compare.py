import json
import os

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "1")

import h5py
import pyscf
from rdkit import Chem
from gpu4pyscf.properties import eda

DATA = "/home/jinyi/local/src/byteff2/example/1_training/example_data"
SYSTEM = "AN_AN_0"
SNAPSHOT = 0

smiles = json.load(open(f"{DATA}/example.json"))[SYSTEM]
with h5py.File(f"{DATA}/example.h5", "r") as f:
    g = f[SYSTEM]
    coords = g["coords"][SNAPSHOT]
    ref = {k: float(g[k][SNAPSHOT]) for k in g.keys() if k != "coords"}

mols, offset = [], 0
for smi in smiles:
    rd = Chem.MolFromSmiles(smi, sanitize=False)
    atoms = sorted(rd.GetAtoms(), key=lambda a: a.GetAtomMapNum())
    n = len(atoms)
    xyz = coords[offset:offset + n]
    charge = sum(a.GetFormalCharge() for a in atoms)
    mols.append(pyscf.M(atom=[(a.GetSymbol(), tuple(x)) for a, x in zip(atoms, xyz)],
                        basis="def2-TZVPD", charge=charge, unit="Angstrom"))
    offset += n
assert offset == len(coords)

eda_result, dft_result = eda.eval_ALMO_EDA_2_energies(mols, xc="wB97M-V")
print("EDA_RESULT", eda_result)
print("DFT_RESULT", dft_result)
print("REFERENCE", ref)
