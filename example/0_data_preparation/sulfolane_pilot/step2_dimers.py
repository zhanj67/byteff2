import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(os.path.dirname(HERE), "generate_dimer.py")
NCONFS = 20
PAIRS = [("DMS", "DMS")]

os.chdir(HERE)
for m1, m2 in PAIRS:
    save_dir = f"dimers/{m1}_{m2}"
    if os.path.isdir(save_dir):
        continue
    xyz2 = f"monomers/{m2}.xyz"
    if m1 == m2:  # distinct name so generate_dimer writes two files, not one file with two frames
        xyz2 = f"monomers/{m2}_B.xyz"
        text = open(f"monomers/{m2}.xyz").read().replace(f"name={m2}", f"name={m2}_B")
        open(xyz2, "w").write(text)
    subprocess.run([sys.executable, GEN, "--mol1", f"monomers/{m1}.xyz", "--mol2", xyz2,
                    "--save_dir", save_dir, "--nconfs", str(NCONFS)], check=True)
    print("DONE", save_dir, flush=True)
