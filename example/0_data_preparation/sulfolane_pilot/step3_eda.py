import glob
import json
import multiprocessing as mp
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GPUS = ["0", "1", "2"]
PAIRS = ["DMS_DMS"]


def worker(gpu, conf_dirs):
    os.environ["CUDA_VISIBLE_DEVICES"] = gpu
    sys.path.insert(0, os.path.dirname(HERE))
    from bytemol.core import Molecule
    from calculate_eda_gpu4pyscf import run_eda

    for conf_dir in conf_dirs:
        out = f"{conf_dir}/eda.json"
        if os.path.exists(out):
            continue
        try:
            mols = [Molecule.from_xyz(x) for x in sorted(glob.glob(f"{conf_dir}/*.xyz"))]
            assert len(mols) == 2
            res = run_eda([m.conformers[0].to_ase_atoms() for m in mols],
                          [int(sum(m.formal_charges)) for m in mols])
        except Exception as e:
            res = {"error": repr(e)}
        json.dump(res, open(out, "w"), indent=1)
        print("FAIL" if "error" in res else "DONE", conf_dir, flush=True)


if __name__ == "__main__":
    os.chdir(HERE)
    jobs = []
    for pair in PAIRS:
        jobs += sorted(glob.glob(f"dimers/{pair}/conf_*"), key=lambda p: int(p.rsplit("_", 1)[-1]))
    ctx = mp.get_context("spawn")
    procs = [ctx.Process(target=worker, args=(g, jobs[i::len(GPUS)])) for i, g in enumerate(GPUS)]
    for p in procs:
        p.start()
    for p in procs:
        p.join()
