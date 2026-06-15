# AdsorbDiff Predictor - Usage Guide

## Overview
This module wraps AdsorbDiff to predict adsorption sites and orientations on
surfaces. It can handle a separate adsorbate or an adsorbate already embedded
in the surface structure.

## Requirements
- AdsorbDiff repo (default: `deps/AdsorbDiff`)
- Python deps: `ase`, `numpy`, and AdsorbDiff dependencies

## Quick start (Python)
```python
from adsorbdiff_predictor import AdsorbDiffPredictor

predictor = AdsorbDiffPredictor(
    adsorbdiff_root="deps/AdsorbDiff",
    checkpoint_path=None,
    use_gpu=True,
    num_sites=10,
)

result = predictor.predict(
    surface="slab.vasp",
    adsorbate="*CO",
    has_adsorbate=False,
)
print(result.summary())
```

## CLI usage
```bash
python adsorbdiff_predictor.py slab.vasp \
  --adsorbate "*CO" \
  --adsorbdiff-root deps/AdsorbDiff \
  --num-samples 5
```

If the surface already contains the adsorbate:
```bash
python adsorbdiff_predictor.py ads_slab.vasp \
  --has-adsorbate
```

## Inputs
- `surface`: file path, ASE Atoms, or pymatgen Structure
- `adsorbate`: SMILES string, database id, ASE Atoms, or file path
- `has_adsorbate`: set to `True` if the adsorbate is already in the structure
- `adsorbate_tag`: tag value used to identify adsorbate atoms (default: 2)

## Outputs
- `best_config.vasp`: best relaxed configuration
- `trajectories/`: diffusion and relaxation trajectories (if enabled)
- `PredictionOutput` includes a ranked list of results

## Notes
- If no relaxation checkpoint is provided, the predictor ranks by a proxy energy
  and sets `use_proxy_energy=True` in the output.
- You can query the internal adsorbate database with
  `AdsorbDiffPredictor.list_available_adsorbates()`.

## Convenience helper
- `predict_adsorption_site(...)`

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

AdsorbDiff itself maintains its own license. Please refer to the AdsorbDiff repository for details.
