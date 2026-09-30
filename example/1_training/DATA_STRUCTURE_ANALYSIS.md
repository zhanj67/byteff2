# example.h5 / example.json structure

## Files
- `example_data/example.h5` — Git LFS object. If `file` reports ASCII text, run `git lfs install && git lfs pull`.
- `example_data/example.json` — `{entry: [mapped_smiles_mol1, mapped_smiles_mol2]}`.
- `example_data/meta.txt` — `example,<entry>` lines consumed by `preprocess.py`.

## Layout
- 950 entries named `<MOL1>_<MOL2>_<k>`, k = 0–4 → 190 unique pairs of 19 species
  (AN DEC DEGDM DMC DMET DOL EA EC EFA EMC FEA FEC FSI GBL LI MA PC PF6 TFSI), self-pairs included.
- Each entry has 20 dimer geometries → 19,000 data points total.

Per entry (N = n_atoms_mol1 + n_atoms_mol2):

| key | shape | meaning (kcal/mol unless noted) |
|---|---|---|
| `coords` | (20, N, 3) | Å; mol1 atoms first, then mol2; within each molecule ordered by SMILES atom-map number |
| `elec_int_energy` | (20,) | ALMO-EDA electrostatic |
| `cls_elec_int_energy` | (20,) | classical electrostatic |
| `elec_pauli_int_energy` | (20,) | electrostatic + Pauli (Pauli itself not stored) |
| `disp_int_energy` | (20,) | dispersion |
| `frozen_int_energy` | (20,) | electrostatic + Pauli + dispersion |
| `polarization_int_energy` | (20,) | polarization |
| `charge_transfer_int_energy` | (20,) | charge transfer |
| `total_int_energy` | (20,) | frozen + polarization + charge transfer |
| `preparation_int_energy` | (20,) | always 0 (rigid monomers) |
| `min_dists` | (20,) | Å, closest intermolecular atom pair |

Molecular charge = sum of formal charges in the SMILES (e.g. LI = +1, TFSI/FSI/PF6 = −1).
Sanity check: LI_LI_0 geometry 0 has total = 113.6 at 2.92 Å ≈ 332/2.92 (pure Li⁺–Li⁺ Coulomb).

## Provenance (verified)
Values are ALMO-EDA version 2 at wB97M-V/def2-TZVPD. Reproduced with gpu4pyscf
`eda.eval_ALMO_EDA_2_energies` to within 0.001 kcal/mol on AN_AN_0 geometry 0.
See `PYSCF_ALMO_EDA_GUIDE.md`.
