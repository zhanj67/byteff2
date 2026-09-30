# EDA components of 10 random entries (seed 42)

Means over the 20 geometries of each entry, kcal/mol. Not representative of all 950 entries.

| entry | total (mean ± std) | frozen | pol | CT | disp | mean min_dist (Å) |
|---|---|---|---|---|---|---|
| EFA_LI_4 | -4.45 ± 8.02 | 2.71 | -6.48 | -0.67 | -0.34 | 4.76 |
| DEC_DMET_4 | 0.22 ± 1.66 | 0.79 | -0.23 | -0.35 | -1.51 | 3.83 |
| AN_DOL_0 | 8.12 ± 27.10 | 12.81 | -1.85 | -2.83 | -2.71 | 4.21 |
| FEA_PC_4 | -0.06 ± 2.56 | 0.56 | -0.34 | -0.28 | -1.27 | 4.11 |
| DMC_DOL_1 | 2.50 ± 12.45 | 3.50 | -0.49 | -0.51 | -1.37 | 4.10 |
| DEGDM_MA_0 | 0.68 ± 3.35 | 1.40 | -0.30 | -0.42 | -1.46 | 3.84 |
| DEGDM_FEA_3 | -0.75 ± 0.75 | -0.38 | -0.16 | -0.22 | -1.56 | 3.68 |
| DEC_FEA_2 | -0.40 ± 0.53 | 0.01 | -0.17 | -0.24 | -1.37 | 3.78 |
| FEA_MA_4 | 0.51 ± 2.23 | 0.96 | -0.15 | -0.30 | -0.91 | 4.30 |
| DEC_DEGDM_4 | 2.36 ± 8.07 | 4.26 | -0.82 | -1.08 | -4.00 | 2.82 |

Frozen already includes dispersion (see DATA_STRUCTURE_ANALYSIS.md), so total = frozen + pol + CT.

## Observations
- Each entry spans short-range repulsive to near-equilibrium geometries; a few close contacts
  (e.g. AN_DOL_0 max 122 kcal/mol) dominate the means and std.
- Neutral–neutral pairs: dispersion is the largest attractive term; pol and CT are small.
- Ion–solvent (EFA_LI_4): polarization dominates the attraction (-6.5 vs -0.3 dispersion).
- train.yaml losses: one MSE each for Pol, Disp, ElecPauli, CT (weight 1), plus total (weight 10).
