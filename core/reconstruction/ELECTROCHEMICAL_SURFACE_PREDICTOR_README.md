# Electrochemical Surface Predictor - Usage Guide

## Overview
This module integrates CP-MACE with VSSR-MC to sample surface structures under
fixed potential and pH conditions.

## Requirements
- CP-MACE repo (default: `deps/CP-MACE`)
- surface-sampling repo (default: `deps/surface-sampling`)
- Python deps: `ase`, `numpy`, `torch`, and surface-sampling dependencies

## Quick start (Python)
```python
from electrochemical_surface_predictor import ElectrochemicalSurfacePredictor

predictor = ElectrochemicalSurfacePredictor(
    cp_mace_root="deps/CP-MACE",
    surface_sampling_root="deps/surface-sampling",
    model_paths=["model1.model", "model2.model"],
    device="cuda",
)

result = predictor.sample_surface(
    surface="slab.vasp",
    adsorbates=["O", "OH"],
    potential_she=1.0,
    ph=12.0,
    temperature=300.0,
    total_sweeps=100,
)
print(result.summary())
```

## CLI usage
```bash
python electrochemical_surface_predictor.py slab.vasp \
  --adsorbates O OH \
  --potential 1.0 \
  --ph 12 \
  --total-sweeps 100 \
  --cp-mace-root deps/CP-MACE \
  --surface-sampling-root deps/surface-sampling
```

## Key parameters
- `potential_she`: electrode potential vs SHE (V)
- `ph`: solution pH
- `total_sweeps`, `sweep_size`: MC sampling controls
- `force_threshold`, `fermi_threshold`: uncertainty thresholds
- `model_paths`: CP-MACE models for ensemble calculations (optional)
- `initial_electron_number`: required if not stored in the structure

## Outputs
The run writes the following files to `output_dir`:
- `final_structure.xyz`
- `lowest_energy_structure.xyz`
- `all_structures.xyz`
- `high_uncertainty_structures.xyz` (if any)
- `settings.json`

## Convenience helper
- `sample_electrochemical_surface(...)`

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

CP-MACE and surface-sampling packages maintain their own licenses. Please refer to their respective repositories for details.
