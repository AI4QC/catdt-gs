# CatMAP Predictor - Usage Guide

## Overview
This module provides a thin wrapper around CatMAP for microkinetic modeling.
It helps define species and reactions, generate CatMAP input files, and run
descriptor scans or single point calculations.

## Requirements
- CatMAP repo (default: `deps/catmap`)
- Python deps: `numpy` and CatMAP dependencies

## Quick start (Python)
```python
from catmap_predictor import CatMAPPredictor

predictor = CatMAPPredictor()

predictor.add_reaction('*_s + CO_g -> CO*')
predictor.add_reaction('CO* + O* -> CO2_g + 2*')

predictor.add_gas('CO', formation_energy=0.0, frequencies=[2170])
predictor.add_gas('O2', formation_energy=0.0, frequencies=[1580])
predictor.add_gas('CO2', formation_energy=-2.45, frequencies=[1333, 2349, 667, 667])

predictor.add_adsorbate('CO', {'Pt': 1.7, 'Pd': 1.55})
predictor.add_adsorbate('O', {'Pt': 1.62, 'Pd': 1.55})

predictor.set_conditions(temperature=500, pressures={'CO_g': 1.0, 'O2_g': 0.33})
result = predictor.run(output_dir='output/catmap')

result.plot_rate()
result.plot_coverage()
```

## Single point calculation
```python
single = predictor.run_single_point(surface='Pt', output_dir='output/catmap_pt')
print(single.summary())
```

## Output files
The run generates these inputs in `output_dir`:
- `energies.txt`
- `model.mkm` (or `single_point.mkm` for single point)

CatMAP writes its own result files alongside these inputs.

## Built-in model builders
The module provides helper constructors for common mechanisms:
- `create_co_oxidation_model(...)`
- `create_her_model(...)`
- `create_orr_model(...)`
- `create_co2rr_model(...)`
- `create_nrr_model(...)`

## Running the built-in demo
The script includes a CO oxidation example in `__main__`:
```bash
python catmap_predictor.py
```
It runs a descriptor scan and writes results to `output/catmap_test`.

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

CatMAP itself maintains its own license. Please refer to the CatMAP repository for details.
