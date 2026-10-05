import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_conf_gpu4pyscf import main

MONOMERS = {
    "EMIM": "[C:1]([C:2]([n+:3]1[c:4]([H:14])[c:5]([H:15])[n:6]([C:7]([H:16])([H:17])[H:18])[c:8]1[H:19])([H:12])[H:13])([H:9])([H:10])[H:11]",
    "LI": "[Li+:1]",
    "CL": "[Cl-:1]",
}

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "monomers"))
for name, smi in MONOMERS.items():
    if not os.path.exists(f"{name}.xyz"):
        main(smi, name)
        print("DONE", name, flush=True)
