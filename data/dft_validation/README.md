# DFT validation of the UMA energetics: VASP inputs

Inputs for the DFT calculations behind Supplementary Table S20 (barriers and reaction energies),
Supplementary Section 9 (energy gaps between sibling intermediates) and Supplementary Table S26
(adsorption sites).
Every structure was generated with UMA (`uma-s-1p1`, task `oc20`), the potential used by CatDT, and
the DFT calculations are run at those geometries, so that the comparison measures the potential.

## Systems

Compact models of the benchmark chemistries, at most 40 atoms.

| Model | Slab atoms | Elementary steps |
|---|---|---|
| Ni(111) p(3x3), 3 layers | 27 | S01 CO+H->CHO, S02 CO+H->COH, S03 CO->C+O, S04 CHO+H->CHOH, S05 COH+H->CHOH, S06 C3H7->C3H6+H, S07 C3H6->C3H5+H |
| Cu(111) p(3x3), 3 layers | 27 | S08 CO2+H->HCOO, S09 CO2+H->COOH, S10 HCOO+H->H2COO |
| Ru(0001) p(3x3), 3 layers | 27 | S11 N2->2N, S12 N+H->NH |
| hcp-PdMo(0001) 3x3, 4 layers | 36 | S13 CO+H->CHO, S14 CO2->CO+O |
| MoS2 monolayer 3x3, in-plane double S vacancy | 25 | S15 CO+H->CHO, S16 CH3O->CH3+O |
| 2D Mo2C (Mo-terminated) 3x3 | 27 | S17 CO2->CO+O, S18 CO+H->CHO |
| rutile TiO2(110) p(2x1), 3 trilayers | 36 | S19 H2->2H (DFT+U check) |
| Ti3O4 on Ni(111) p(3x3), from the VSSR-reconstructed TiOx overlayer | 34 | adsorption sites only |

Sibling pairs (same parent, same composition): CHO/COH, CH2O/CHOH and 1-/2-propyl on Ni(111);
HCOO/COOH on Cu(111); CHO/COH on PdMo, MoS2 and Mo2C; CH3O/CH2OH on MoS2.

Adsorption-site sets (every symmetry-distinct site, relaxed): H on the Ni/TiOx interface; CO and H on
PdMo(0001); CO and H on MoS2 with the double S vacancy; H and CH3 on TiO2(110).

## How the structures were obtained

- Adsorption geometry: every symmetry-distinct site (top, bridge, fcc/hcp hollow; sites over a vacancy
  on the defective surface), two orientations each, relaxed with UMA; the lowest is kept. Heights are
  measured from the atoms under the site, so that sites in the troughs of a corrugated surface are
  reached.
- Elementary steps: initial and final states built by CatDT's Agent 4/5 loop from the lowest-energy
  placement of the reactant; transition states from CatDT's two-phase CI-NEB tool on UMA (no climbing
  image to 0.2 eV/A, then climbing to 0.05 eV/A, IDPP start), all images of a band evaluated in one
  batched UMA call. A step is accepted only after the Agent 5 checks: both endpoints carry exactly the
  species of the step, the reacting bond is activated at the transition state, and the recorded path
  is continuous.
- Fixed atoms are the bottom layer of the three-layer metals, the bottom two of
  PdMo, bottom S of MoS2, bottom Mo of Mo2C, bottom trilayer of TiO2).

UMA energies of every structure are in the extxyz comment lines under `uma_reference/` and in
`manifest.csv`; the paths are in `uma_reference/steps/*/neb_path.extxyz`.

## DFT settings

Main tier (`01`-`04`): the level of theory of the OC20 data on which the UMA `oc20` task is
trained, taken from `fairchem/data/oc/utils/vasp_flags.py` and used in the same way by CatTSunami for
transition-state single points.

| Setting | Value |
|---|---|
| Functional | RPBE (`GGA = RP`) |
| PAW data sets | standard PBE sets, ASE default choice |
| Cut-off | `ENCUT = 350` |
| Spin | `ISPIN = 1` |
| Smearing | `ISMEAR = 1, SIGMA = 0.2` |
| Symmetry | `ISYM = 0` |
| k-mesh | `round(40/a) x round(40/b) x 1`, Monkhorst-Pack, with a and b as in `calculate_surface_k_points` |
| Electronic convergence | `EDIFF = 1E-6` (single points), `1E-5` (relaxations) |
| Relaxation | `IBRION = 2, EDIFFG = -0.03, NSW = 300` |

Energies are the `energy(sigma->0)` of the last ionic step.

Physics tiers, for what the OC20 level leaves out:

- `05_spin_polarised/`: every Ni(111) calculation again with `ISPIN = 2` (initial moment 1 muB per Ni).
- `06_dftU_spin/`: H2 dissociation on TiO2(110), which creates two Ti3+ centres, with DFT+U
  (Dudarev, U(Ti 3d) = 4.2 eV) and spin polarisation (`ISMEAR = 0, SIGMA = 0.05`).

## Layout

```
01_barrier_single_points/<step>/{IS,TS,FS}/       POSCAR INCAR KPOINTS POTCAR
02_sibling_pairs/<pair>/{A_<species>,B_<species>}/
03_clean_slabs/<model>/
04_site_relaxations/<model>__<adsorbate>/site_NN/  site_00 = site selected by UMA
05_spin_polarised/...                              Ni models, ISPIN = 2
06_dftU_spin/...                                   TiO2, DFT+U
manifest.csv                                       one row per calculation
uma_reference/                                     UMA structures, energies, paths, summaries
```

## What is here, and what is not

Every calculation directory holds the VASP input (`POSCAR`, `INCAR`, `KPOINTS`), the relaxed
structure (`CONTCAR`) and the convergence record (`OSZICAR`, `vasp.out`). `manifest.csv` has one
row per calculation with its system, role, tier, k-mesh, composition and UMA energy, and
`results/` holds the comparison tables. `uma_reference/` holds the UMA structures and paths the
DFT calculations were run at.

`POTCAR` files are not included: the VASP PAW data sets are licensed and may not be
redistributed. `OUTCAR` files are omitted for size; the energies they contain are in
`manifest.csv`, in `results/` and, step by step, in Supplementary Table S20.
