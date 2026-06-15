# VSSR-MC Predictor - Usage Guide

## Overview
This module wraps the surface-sampling VSSR-MC algorithm to sample surface
reconstructions using virtual adsorption sites and Monte Carlo sampling.

## Requirements
- surface-sampling repo (default: `deps/surface-sampling`)
- Python deps: `ase`, `numpy`, and surface-sampling dependencies

## Quick start (Python)
```python
from vssr_mc_predictor import VSSRMCPredictor

predictor = VSSRMCPredictor(
    surface_sampling_root="deps/surface-sampling",
    model_type="CHGNetNFF",
    device="cuda",
)

result = predictor.sample(
    surface="slab.vasp",
    adsorbates=["O", "H"],
    total_sweeps=100,
    temperature=1.0,
)
print(result.summary())
```

## CLI usage
```bash
python vssr_mc_predictor.py slab.vasp \
  --surface-sampling-root deps/surface-sampling \
  --model-type CHGNetNFF \
  --adsorbates O H \
  --total-sweeps 50 \
  --temperature 1.0
```

## Key parameters
- `canonical`: use fixed composition (NVT) sampling
- `num_adsorbates`: number of adsorbates when `canonical=True`
- `auto_offset`: auto-generate `offset_data` if not provided
- `chem_pots`: chemical potentials for semi-grand canonical mode

## Outputs
Files written to `output_dir` include:
- `mc.log`
- `all_virtual_ads.cif`
- `*_mcmc_structures.pkl`
- `offset_data.json` (if auto generated)
- `settings.json`

## Convenience helper
- `sample_surface_reconstruction(...)`

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

The surface-sampling package itself maintains its own license. Please refer to the surface-sampling repository for details.
