# Data Preparation Example

This directory contains example scripts for preparing molecular data for training `ByteFF-Pol`.

## Usage

### 1. Generate Conformers

To generate optimized molecular conformers:

```bash
PYTHONPATH=$(git rev-parse --show-toplevel):${PYTHONPATH} python generate_conf.py --mol_name ACT --mapped_smiles "[C:1]([C:2](=[O:3])[C:4]([H:8])([H:9])[H:10])([H:5])([H:6])[H:7]"
```

This script uses geomeTRIC optimizer with B3LYP/def2-SVPD to optimize the molecular structure and saves it as an XYZ file.

### 2. Generate Dimers

To create molecular dimer configurations:

```bash
PYTHONPATH=$(git rev-parse --show-toplevel):${PYTHONPATH} python generate_dimer.py --mol1 ACT.xyz --mol2 EC.xyz --save_dir ACT_EC_dimer --nconfs 100
```

This script using optimized structure of step 1 to generates multiple dimer configurations by randomly rotating and translating molecules while maintaining minimum/maximum distance constraints.

### 3. Calculate EDA

To perform Energy Decomposition Analysis on dimers:

```bash
PYTHONPATH=$(git rev-parse --show-toplevel):${PYTHONPATH} python calculate_eda.py --input_dir ACT_EC_dimer/conf_0
```

This script performs ALMO-EDA calculations using wB97M-V/def2-TZVPD to extract interaction energy components (electrostatic, exchange, polarization, charge transfer, etc.).

### 3b. Calculate EDA without Q-Chem (gpu4pyscf)

`calculate_eda.py` needs a Q-Chem license. `calculate_eda_gpu4pyscf.py` is a drop-in replacement: the same `--input_dir` and the same output JSON keys (`ELEC`, `CLS ELEC`, `ELEC_PAULI`, `DISP`, `FROZEN`, `POLARIZATION`, `CHARGE TRANSFER`, `PREPARATION`, `TOTAL`) in kcal/mol. It runs gpu4pyscf `eda.eval_ALMO_EDA_2_energies` (ALMO-EDA2, wB97M-V/def2-TZVPD) on a GPU.

```bash
PYTHONPATH=$(git rev-parse --show-toplevel):${PYTHONPATH} python calculate_eda_gpu4pyscf.py --input_dir ACT_EC_dimer/conf_0
```

Mapping from gpu4pyscf output (kJ/mol) to Q-Chem EDA2 labels: `electrostatic`→ELEC, `classical electrostatic`→CLS ELEC, `electrostatic + pauli`→ELEC_PAULI (Q-Chem's "PAULI" line already includes electrostatics), `dispersion`→DISP, `frozen`→FROZEN, `polarization`, `charge transfer`, `total`. PREPARATION is 0 because the monomers are rigid.

**Validation against the Q-Chem data in `example/1_training/example_data/example.h5`** (geometry 0 of each entry, kcal/mol, gpu4pyscf / Q-Chem). Full numbers are in `gpu4pyscf_validation/validation_results.json`.

| dimer | atoms | elec | elec_pauli | disp | polarization | charge_transfer | total | max abs diff |
|---|---|---|---|---|---|---|---|---|
| AN_AN_0 | 12 | -17.930 / -17.931 | 11.125 / 11.125 | -6.037 / -6.038 | -1.151 / -1.151 | -1.194 / -1.194 | 2.742 / 2.741 | 0.0009 |
| DMC_PF6_0 | 19 | -2.023 / -2.023 | 0.533 / 0.534 | -1.787 / -1.788 | -1.802 / -1.803 | -0.177 / -0.177 | -3.233 / -3.234 | 0.0017 |
| DOL_FSI_0 | 20 | -3.266 / -3.266 | 3.438 / 3.438 | -3.872 / -3.870 | -1.276 / -1.276 | -0.519 / -0.520 | -2.229 / -2.227 | 0.0019 |
| EC_EC_0 | 20 | -0.374 / -0.374 | 5.354 / 5.355 | -2.994 / -2.993 | -0.742 / -0.742 | -0.571 / -0.571 | 1.047 / 1.048 | 0.0008 |
| EC_LI_0 | 11 | 22.774 / 22.774 | 73.387 / 73.388 | -2.866 / -2.866 | -35.397 / -35.398 | -6.425 / -6.424 | 28.699 / 28.699 | 0.0010 |
| FEC_FEC_0 | 20 | -53.397 / -53.397 | 44.698 / 44.698 | -14.691 / -14.691 | -4.777 / -4.777 | -5.057 / -5.057 | 20.173 / 20.173 | 0.0006 |
| FSI_LI_0 | 10 | -86.040 / -86.042 | -85.109 / -85.112 | -0.260 / -0.259 | -10.794 / -10.794 | -0.578 / -0.579 | -96.741 / -96.743 | 0.0023 |
| PF6_LI_0 | 8 | -82.414 / -82.415 | -82.228 / -82.228 | -0.085 / -0.085 | -4.636 / -4.635 | -0.276 / -0.277 | -87.225 / -87.225 | 0.0011 |
| PF6_PF6_0 | 14 | 49.749 / 49.750 | 77.411 / 77.412 | -4.543 / -4.545 | -3.564 / -3.563 | -0.336 / -0.337 | 68.968 / 68.967 | 0.0021 |

Environment notes:
- Tested with `gpu4pyscf-cuda11x==1.4.3`, `cupy-cuda11x`, and `pyscf==2.9.0` on NVIDIA driver 515 (CUDA 11.7). Newer gpu4pyscf builds fail on this driver ("PTX was compiled with an unsupported toolchain"). pyscf 2.14 is incompatible with gpu4pyscf 1.4.3.
- Cost: about 9 min for a 12-atom dimer on an RTX 3080 that other jobs were also using.

### 1b. Generate Conformers without Q-Chem (gpu4pyscf)

`generate_conf_gpu4pyscf.py` is a drop-in replacement for `generate_conf.py`: same arguments, same method (B3LYP-D3BJ/def2-SVPD, geomeTRIC, GAU convergence), and the same output `<mol_name>.xyz`. Single atoms (e.g. `[Li+:1]`, `[Cl-:1]`) are written without optimization.

```bash
PYTHONPATH=$(git rev-parse --show-toplevel):${PYTHONPATH} python generate_conf_gpu4pyscf.py --mol_name ACT --mapped_smiles "[C:1]([C:2](=[O:3])[C:4]([H:8])([H:9])[H:10])([H:5])([H:6])[H:7]"
```

## Tutorial: new training data without Q-Chem (Cl⁻ pilot)

`cl_pilot/` generates EDA training data for Cl⁻ paired with Emim⁺, Li⁺ and Cl⁻ (20 geometries per pair, 60 clusters) and packs it in the same format as `example/1_training/example_data/example.h5`. Each step is a script with its settings written inline; edit the constants at the top to make your own dataset. Run all steps from `cl_pilot/` with:

```bash
export PYTHONPATH=$(git rev-parse --show-toplevel):${PYTHONPATH}
```

**Step 1: monomers** (`step1_monomers.py` → `monomers/<NAME>.xyz`). List each monomer's atom-mapped SMILES; every atom, including H, needs a map number:

```python
MONOMERS = {
    "EMIM": "[C:1]([C:2]([n+:3]1[c:4]([H:14])[c:5]([H:15])[n:6]([C:7]([H:16])([H:17])[H:18])[c:8]1[H:19])([H:12])[H:13])([H:9])([H:10])[H:11]",
    "LI": "[Li+:1]",
    "CL": "[Cl-:1]",
}
```

To get a mapped SMILES from a plain SMILES:

```python
from rdkit import Chem
m = Chem.AddHs(Chem.MolFromSmiles("CC[n+]1ccn(C)c1"))
for a in m.GetAtoms():
    a.SetAtomMapNum(a.GetIdx() + 1)
print(Chem.MolToSmiles(m))
```

```bash
CUDA_VISIBLE_DEVICES=0 python step1_monomers.py
```

**Step 2: dimer geometries** (`step2_dimers.py` → `dimers/<M1>_<M2>/conf_*/`). Set `PAIRS` and `NCONFS`; it calls `generate_dimer.py` for each pair. For self-pairs (e.g. `("CL", "CL")`) the script writes a renamed copy `CL_B.xyz`, because `generate_dimer.py` would otherwise append both molecules into one file. `generate_dimer.py` always uses seed 42, so two pairs of single atoms get identical separations; use a different seed per pair for production data.

```python
NCONFS = 20
PAIRS = [("CL", "EMIM"), ("CL", "LI"), ("CL", "CL")]
```

```bash
python step2_dimers.py
```

**Step 3: EDA** (`step3_eda.py` → `conf_*/eda.json`). Spreads all `conf_*` folders across the GPUs listed in `GPUS` and calls `run_eda` from `calculate_eda_gpu4pyscf.py`. A failed geometry writes `{"error": ...}` instead of stopping the batch, and rerunning skips finished geometries. Run it detached, since it takes hours:

```python
GPUS = ["0", "1", "2"]
PAIRS = ["CL_EMIM", "CL_LI", "CL_CL"]
```

```bash
nohup python step3_eda.py > step3.log 2>&1 &
grep -c '^DONE' step3.log; grep '^FAIL' step3.log
```

**Step 4: pack** (`step4_pack.py` → `packed/cl_pilot.{h5,json}`, `packed/meta.txt`). Writes one entry per pair (`<M1>_<M2>_0`) with the same keys as `example.h5`: `coords` (mol1 atoms then mol2), `min_dists`, and the nine `*_int_energy` arrays in kcal/mol. Failed geometries are skipped.

```bash
python step4_pack.py
python ../../1_training/preprocess.py --conf preprocess_cl_pilot.yaml   # → processed_data/
```

**Pilot results** (60/60 clusters succeeded):

| pair | geometries | total interaction energy (kcal/mol) | min distance (Å) |
|---|---|---|---|
| Cl⁻–Emim⁺ | 20 | −90.8 to +5.6 | 1.55–7.08 |
| Cl⁻–Li⁺ | 20 | −147.2 to −34.6 | 1.56–9.61 |
| Cl⁻–Cl⁻ | 20 | +34.4 to +345.3 | 1.56–9.61 |

Geometries with very large forces (e.g. Cl⁻–Cl⁻ below 2 Å) are kept; the training loss down-weights them through `force_cutoff: [50, 80]`.

**Throughput** (3 shared RTX 3080s): the EDA step dominates. A 20-atom cluster (Cl⁻–Emim⁺) takes about 33 min per GPU (29–43 min), so 1,000 such clusters take about 7.6 days on 3 GPUs. Single-atom pairs take about 2.5 min each. Cost rises steeply with cluster size.

**Fine-tuning:** `example/1_training/train_cl_pilot.yaml` and `finetune_cl_pilot.py` fine-tune `optimal.pt` on `example.h5` plus this set. This is a pipeline demonstration only. With 60 clusters and equal dataset weights, the Cl⁻ loss dominated and the example-data validation loss rose from 4.49 to 10.64 over 40 epochs. A real run needs more Cl⁻ data, a smaller `loss_weight` for the Cl⁻ dataset, and held-out Cl⁻–cation entries for validation.
