# Training example: setup notes

## Run
```bash
PYTHONPATH=$(git rev-parse --show-toplevel) python preprocess.py --conf preprocess_example.yaml
PYTHONPATH=$(git rev-parse --show-toplevel) python train.py --conf train.yaml
```
- preprocess → `processed_data/{dataset_config.yaml, processed_data_shard0.pkl}` (~25 MB).
- train → `training_logs/` (`optimal.pt`, `history.json`, `ckpt/`, logs). The script refuses to start
  if `training_logs/` exists; delete it or pass `--restart`.

## Fixes needed on this machine
1. Git LFS: `conda install -c conda-forge git-lfs && git lfs install && git lfs pull`
   (otherwise h5py fails with "file signature not found").
2. PyTorch was CPU-only (`2.6.0+cpu`). Installed `torch==2.0.1+cu117` from the PyTorch cu117
   index to match driver 515 / CUDA 11.7.
3. OOM on a 10 GB RTX 3080 at `batch_size: 96`. `train.yaml` is now set to 16 (local change).

## What training does
Fine-tunes the released checkpoint `byteff2/trained_models/optimal.pt`. The learning rate is 0 for
Graph and MMBondedConj (frozen), 2e-5 for ChargeVolume, and 2e-4 for Exp6Pol. Validation runs
every 4 epochs, a checkpoint every 12, and early stop after 20 validations without improvement.
Observed run: 83 epochs in about 1 h, best validation loss 4.38 at epoch 20, stopped manually.

## Data generation
- Recipe (`example/0_data_preparation`): conformers at B3LYP/def2-SVPD with geomeTRIC → random
  dimers → ALMO-EDA2 at wB97M-V/def2-TZVPD.
- Verified reproducible with gpu4pyscf; see `PYSCF_ALMO_EDA_GUIDE.md`.
- ALMO-EDA availability: Q-Chem (commercial, the original implementation) and gpu4pyscf (free).
  Gaussian, ORCA, and Psi4 do not provide ALMO-EDA (ORCA has LED, Psi4 has SAPT; both give
  different decompositions).
- Cost: about 9 GPU-min per 12-atom dimer here. Larger dimers scale steeply. Full paper-scale
  data (the user reports ~6M points) is HPC-scale.
- Adding a new element such as boron needs new dimer data covering it. The frozen GNN has not
  seen boron, so consider unfreezing Graph (lr > 0) when fine-tuning.
