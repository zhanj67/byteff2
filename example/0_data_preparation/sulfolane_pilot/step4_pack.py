import glob
import json
import os

import h5py
import numpy as np

from bytemol.core import Molecule

HERE = os.path.dirname(os.path.abspath(__file__))
DATASET = "sulfolane_pilot"
PAIRS = ["DMS_DMS"]
KEYMAP = {
    "elec_int_energy": "ELEC",
    "cls_elec_int_energy": "CLS ELEC",
    "elec_pauli_int_energy": "ELEC_PAULI",
    "disp_int_energy": "DISP",
    "frozen_int_energy": "FROZEN",
    "polarization_int_energy": "POLARIZATION",
    "charge_transfer_int_energy": "CHARGE TRANSFER",
    "preparation_int_energy": "PREPARATION",
    "total_int_energy": "TOTAL",
}

os.chdir(HERE)
os.makedirs("packed", exist_ok=True)
smiles_map, meta = {}, []
with h5py.File(f"packed/{DATASET}.h5", "w") as h5:
    for pair in PAIRS:
        entry = f"{pair}_0"
        coords, mind, energies, smiles = [], [], {k: [] for k in KEYMAP}, None
        confs = sorted(glob.glob(f"dimers/{pair}/conf_*"), key=lambda p: int(p.rsplit("_", 1)[-1]))
        for conf in confs:
            res = json.load(open(f"{conf}/eda.json")) if os.path.exists(f"{conf}/eda.json") else {"error": "missing"}
            if "error" in res:
                print("skip", conf, res["error"][:80])
                continue
            mols = [Molecule.from_xyz(x) for x in sorted(glob.glob(f"{conf}/*.xyz"))]
            c1, c2 = (m.conformers[0].coords for m in mols)
            smiles = [m.get_mapped_smiles(isomeric=True) for m in mols]
            coords.append(np.concatenate([c1, c2]))
            mind.append(np.linalg.norm(c1[:, None] - c2[None], axis=-1).min())
            for k, v in KEYMAP.items():
                energies[k].append(res[v])
        if not coords:
            continue
        g = h5.create_group(entry)
        g["coords"] = np.array(coords)
        g["min_dists"] = np.array(mind)
        for k, v in energies.items():
            g[k] = np.array(v)
        smiles_map[entry] = smiles
        meta.append(f"{DATASET},{entry}")
        print(entry, len(coords), "confs")

json.dump(smiles_map, open(f"packed/{DATASET}.json", "w"), indent=2)
open("packed/meta.txt", "w").write("\n".join(meta) + "\n")
