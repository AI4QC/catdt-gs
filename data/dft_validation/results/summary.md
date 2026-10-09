# DFT validation results

Job status: ok 176

## Barriers, tier `oc20` (n = 18 steps; reported apart: ['S19_H2_to_2H'])
- MAE(Ea) = 0.087 eV; max |err| = 0.475 eV
- within 0.1 eV of DFT: 14/18 (78%)
- MAE(reaction energy) = 0.067 eV

## Barriers, tier `spin` (n = 7 steps; reported apart: none)
- MAE(Ea) = 0.117 eV; max |err| = 0.198 eV
- within 0.1 eV of DFT: 4/7 (57%)
- MAE(reaction energy) = 0.170 eV

## Barriers, tier `u_spin` (n = 1 steps; reported apart: none)
- MAE(Ea) = 0.191 eV; max |err| = 0.191 eV
- within 0.1 eV of DFT: 0/1 (0%)
- MAE(reaction energy) = 0.344 eV

## Sibling gaps (n = 8 pairs)
- MAE(gap) = 0.086 eV; max |err| = 0.292 eV
- MAE(absolute adsorption energy, same species) = 0.126 eV (n = 16)
- pruning threshold 0.30 eV / MAE(gap) = 3.5

## Adsorption sites (7 adsorbate/surface sets)
- UMA-selected site is the DFT global minimum, or degenerate with it to 0.01 eV, in 6/7 sets
  - mos2_sv2__CO: 3 configurations, UMA pick site_00, DFT minimum site_00, UMA pick lies 0.0 eV above the DFT minimum
  - mos2_sv2__H: 22 configurations, UMA pick site_00, DFT minimum site_02, UMA pick lies 0.0 eV above the DFT minimum
  - ni_tiox__H: 17 configurations, UMA pick site_00, DFT minimum site_00, UMA pick lies 0.0 eV above the DFT minimum
  - pdmo0001__CO: 4 configurations, UMA pick site_00, DFT minimum site_00, UMA pick lies 0.0 eV above the DFT minimum
  - pdmo0001__H: 4 configurations, UMA pick site_00, DFT minimum site_02, UMA pick lies 0.196 eV above the DFT minimum
  - tio2_110__CH3: 8 configurations, UMA pick site_00, DFT minimum site_00, UMA pick lies 0.0 eV above the DFT minimum
  - tio2_110__H: 5 configurations, UMA pick site_00, DFT minimum site_00, UMA pick lies 0.0 eV above the DFT minimum
- otherwise within: pdmo0001__H 0.196 eV
