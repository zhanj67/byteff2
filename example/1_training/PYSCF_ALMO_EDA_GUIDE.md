# Reproducing example.h5 energies with gpu4pyscf ALMO-EDA

## Method
`gpu4pyscf.properties.eda.eval_ALMO_EDA_2_energies(mol_list, xc="wB97M-V")`,
basis `def2-TZVPD`, one `pyscf.M` per monomer with its formal charge.
Reference example: gpu4pyscf `examples/properties/36-amlo_eda.py`.

Output is kJ/mol; divide by 4.184. Key mapping:

| gpu4pyscf key | example.h5 key |
|---|---|
| `electrostatic` | `elec_int_energy` |
| `classical electrostatic` | `cls_elec_int_energy` |
| `electrostatic` + `pauli` | `elec_pauli_int_energy` |
| `dispersion` | `disp_int_energy` |
| `frozen` | `frozen_int_energy` |
| `polarization` | `polarization_int_energy` |
| `charge transfer` | `charge_transfer_int_energy` |
| `total` | `total_int_energy` |

## Validation (AN_AN_0, geometry 0, kcal/mol)
| | gpu4pyscf | reference |
|---|---|---|
| electrostatic | -17.930 | -17.931 |
| classical elec | -8.365 | -8.365 |
| elec + Pauli | 11.125 | 11.125 |
| dispersion | -6.038 | -6.038 |
| frozen | 5.087 | 5.086 |
| polarization | -1.151 | -1.151 |
| charge transfer | -1.194 | -1.194 |
| total | 2.742 | 2.741 |

Cost: ~9 min for this 12-atom dimer on a shared RTX 3080.

## Environment on this machine (conda env `byt_min`)
NVIDIA driver 515 / CUDA 11.7 is too old for gpu4pyscf ≥ 1.5 (error: "PTX was compiled with an
unsupported toolchain", then `failed in block_diag kernel`). Working combination:
- `gpu4pyscf-cuda11x==1.4.3` (oldest release containing `properties/eda`)
- `cupy-cuda11x`
- `pyscf==2.9.0` (2.14 raises `eigh() takes from 2 to 3 positional arguments`)

The "Failed to set CUDA shm size" warning is harmless. Select a GPU with `CUDA_VISIBLE_DEVICES`.

## Script
`eda_compare.py` (this directory): set `SYSTEM` / `SNAPSHOT` at the top, run `python eda_compare.py`.
It reads one entry/geometry from example.h5, builds monomers from the mapped SMILES in
example.json (atoms sorted by map number, charge from formal charges), runs the EDA, and prints
results (kJ/mol) next to the reference (kcal/mol).

## Not usable for this data
- Gaussian 16: no ALMO-EDA.
- Plain supermolecular E(dimer) − E(A) − E(B) gives only the total, no components.
