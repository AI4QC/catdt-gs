# SurFF Predictor - Usage Guide

## Overview
This module wraps SurFF to predict exposed surfaces for a bulk crystal using
slab generation, relaxation, and Wulff construction.

## Requirements
- SurFF repo (default: `deps/SurFF`)
- Python deps: `ase`, `numpy`, `pymatgen`, `pandas`, and SurFF dependencies

## Quick start (Python)
```python
from surff_predictor import SurFFPredictor

predictor = SurFFPredictor(surff_root="deps/SurFF", use_gpu=True)
result = predictor.predict(structure="POSCAR", crystal_id="crystal", top_n=5)
print(result.summary())
```

## CLI usage
```bash
python surff_predictor.py POSCAR \
  --surff-root deps/SurFF \
  --top-n 5 \
  --output-dir output/surff
```

## Outputs
When `output_dir` is set, the predictor writes a structured output tree:
- `crystal/`: input crystal POSCAR
- `slabs/`: generated slabs
- `bulks/`: oriented bulk structures
- `lmdb/`: relaxation datasets
- `traj/`: relaxation trajectories
- `wulff/`: Wulff plot (if enabled)

CLI mode also writes a CSV summary to `output_dir/<crystal_id>_results.csv`.

## Notes
- `save_wulff_shape=True` enables the Wulff plot export.
- Use `--cpu` if GPU is not available.

## Convenience helper
- `predict_surface_exposure(...)`

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

SurFF itself maintains its own license. Please refer to the SurFF repository for details.
