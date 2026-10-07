import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_conf_gpu4pyscf import main

MONOMERS = {
    "DMS": "[C:1]([S:2]([C:3]([H:9])([H:10])[H:11])(=[O:4])=[O:5])([H:6])([H:7])[H:8]",
}

MONO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "monomers")
os.makedirs(MONO_DIR, exist_ok=True)
os.chdir(MONO_DIR)
for name, smi in MONOMERS.items():
    if not os.path.exists(f"{name}.xyz"):
        main(smi, name)
        print("DONE", name, flush=True)
