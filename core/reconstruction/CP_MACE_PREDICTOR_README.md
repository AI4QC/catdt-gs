# CP-MACE Predictor - Usage Guide

## Overview
This module wraps CP-MACE to train models and run constant-potential molecular dynamics.
It exposes a Python API and a small CLI for training and simulation.

## Requirements
- CP-MACE repo (default: `deps/CP-MACE`)
- Python deps: `ase`, `numpy`, `torch`, `yaml`, and CP-MACE dependencies
- GPU is optional; pass `--cpu` to force CPU

## Quick start (Python)
```python
from cp_mace_predictor import CPMACEPredictor

predictor = CPMACEPredictor(cp_mace_root="deps/CP-MACE", use_gpu=True)

# Train
train_result = predictor.train(
    train_file="train.extxyz",
    model_name="my_model",
    max_num_epochs=300,
)

# Simulate
sim_result = predictor.simulate(
    structure="init.xyz",
    model_paths=[train_result.model_path],
    target_potential=-3.36,
    temperature=300.0,
    steps=1000,
)
print(sim_result.summary())
```

## CLI usage
Train:
```bash
python cp_mace_predictor.py train train.extxyz \
  --cp-mace-root deps/CP-MACE \
  --model-name my_model \
  --max-epochs 300
```

Simulate:
```bash
python cp_mace_predictor.py simulate init.xyz \
  --model-paths my_model.model \
  --target-potential -3.36 \
  --temperature 300 \
  --steps 1000
```

## Input data format
Training data must be extended XYZ with these keys:
- `REF_energy` (eV)
- `REF_forces` (eV/Angstrom)
- `electron` (electron number)
- `potential` (V)

You can generate this file with `CPMACEPredictor.prepare_training_data()`.

## Key parameters
- `train()`: `model_type`, `r_max`, `batch_size`, `max_num_epochs`,
  `energy_weight`, `forces_weight`, `potential_weight`
- `simulate()`: `initial_electron_number`, `timestep`, `integrator`,
  `force_threshold`, `fermi_threshold`, `save_frequency`

## Outputs
Training:
- model file: `<model_name>.model` (or compiled model if available)
- log file: `train.log`

Simulation:
- trajectory: `atoms.traj`
- final structure: `final.xyz`
- optional: `high_uncertainty.xyz` (if uncertainty thresholds are triggered)

## Convenience helpers
- `train_cp_mace(...)`
- `run_cp_md(...)`
- `run_slow_growth(...)`

## License

This wrapper script is part of CatDT and is licensed under the Apache-2.0 License.

CP-MACE itself maintains its own license. Please refer to the CP-MACE repository for details.
