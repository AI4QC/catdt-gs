# <span style="font-size:larger;">CP-MACE</span>

<img width="90%" alt="image" src="https://github.com/user-attachments/assets/fc7589c4-030e-4281-9207-bb0bb2998347" />

This repository contains the MACE reference implementation developed by Ruoyu Wang and Shaoheng Fang.

- [About CP-MACE](#about-mace)
- [Installation](#installation)
- [Dataset](#dataset)
- [Usage](#usage)
- [References](#references)
- [Contact](#contact)

## About CP-MACE

CP-MACE (Constant-Potential MACE) is an extension of the MACE framework that enables constant-potential (CP) molecular simulations via machine learning force field (MLFF). It incorporates the number of electrons as an input and learns to predict the Fermi level as an additional output, allowing accurate modeling of electrochemical interfaces under grand canonical ensemble conditions.

CP-MACE is built on the [MACE](https://github.com/ACEsuit/mace) architecture, which combines equivariant message passing neural networks with the Atomic Cluster Expansion (ACE) formalism for high-accuracy interatomic potential modeling.


## Installation

The easiest way to install it is to first install the regular MACE framework and then replace the original MACE directory with our modified version.

Alternatively, if you prefer installing from source:
1. Clone the original MACE repository:
```sh
git clone https://github.com/ACEsuit/mace.git
```
2. Replace the ./mace directory with the CP-MACE version provided.
3. Install it via pip:
```sh
pip install ./mace
```

## Dataset
CP-MACE uses the same xyz file format as the regular MACE framework, with one important addition: for each structure, you must include two extra tags after the atom count (i.e., the second line of each structure block) with the following format:
```sh
electron=XX potential=YY 
```
`electron` refers to the net charge in the system and `potential` refers to the Fermi level of the system. You can also supply the electron number as `electron`, which might be more convenient. For datasets containing only a single system type, both options are equivalent. However, if your dataset includes different systems with varying atom counts or chemical compositions, we strongly recommend using net charge instead of electron number, as the absolute electron numbers may differ significantly across systems and impair model training. An example for the data is:

```sh
202
Lattice="11.798185 0.0 0.0 0.0 11.798185 0.0 0.0 0.0 33.342577" Properties=species:S:1:pos:R:3:REF_forces:R:3 REF_energy=-873.98192526 pbc="T T T" potential=-3.407347 electron=661.7
H        8.11218000       7.31211000      11.63614000      -0.35109000      -0.57477900      -0.54562900
H        7.33506000       6.74093000      12.76225000      -0.53156400      -0.57238000       0.33610800
H        4.18877000       6.81587000      18.30567000       1.48010700      -0.05893200       0.44728300
……
```

## Usage

### Training

The training script for CP-MACE closely follows that of MACE. Below is a simple example:
```sh
mace_run_train \
    --name="MACE_model" \
    --train_file="train.xyz" \
    --valid_fraction=0.05 \
    --config_type_weights='{"Default":1.0}' \
    --energy_key='REF_energy' \
    --forces_key='REF_forces' \
    --E0s='average' \
    --model="FermiMACE" \
    --loss="fermi_weighted" \
    --error_table="Fermi_PerAtomRMSE" \
    --forces_weight=100.0 \
    --energy_weight=1.0 \
    --potential_weight=10.0 \
    --hidden_irreps='128x0e + 128x1o' \
    --r_max=5.0 \
    --batch_size=10 \
    --max_num_epochs=300 \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --restart_latest \
    --device=cuda \
    --seed=1
```

Several additional arguments are introduced in CP-MACE to support constant-potential training:

Use `--model=FermiMACE` to enable the node augmentation method (recommended). Alternatively, the global state method is also available by setting `--model=FermiMACE_2`.

The loss function should be adjusted for Fermi level prediction by setting `--loss=fermi_weighted`.

Use `--error_table=Fermi_PerAtomRMSE` to customize the error reporting.

Tune the weight of the potential term via `--potential_weight=XX`.

### Simulation

We provide some scripts in `./simulation` to run free MD, slow growth or metadynamics simulations. You can modify or extend the scripts based on your own needs.

Taking slow growth as an example, we provide an initial structure `init.xyz` and two pre-trained models in `./simulation/slow_growth`. Parameters are configured in `inputs.yml`:

```sh
save_dir: result
save_freq: 1
steps: 100
read_velocity: True
force_threshold: 0.15
fermi_threshold: 0.04
T_init: 300.
integrator: NoseHoover
integrator_config: {"timestep": 1, "temperature": 300., "ttime": 40., "constraints": [[0,137,204,1.36091]], "increm": 0.001, "Mne": 660.74, "eta_length":2, "targetmu": -3.36}
```
`steps` defines the total number of simulation steps. 

If `read_velocity` is True, initial velocities will be read from init.xyz; otherwise, velocities will be initialized according to a Boltzmann distribution at temperature `T_init`. 

In `integrator_config`, in addition to thermostat-related parameters, `timestep` is the simulation time step and `targetmu` is the target electrode potential. Currently, the slow growth calculation only supports using the distance between two atoms as the collective variable (CV). In `constraints`, the second and third parameters specify the indexes of the two atoms (as in the xyz file), and the fourth parameter is their distance. `increm` defines the change in distance at each time step. To perform a free MD simulation, simply set `constraints` to an empty list [].

You can run the simulation using:
```sh
python simulate.py > log
```
The log file will record the standard deviation of forces and Fermi levels, as well as the distance and gradient of the target CV. If the standard deviation of forces or Fermi levels exceeds the specified threshold, the structure will be saved to `save_dir/structures_force.xyz` or `save_dir/structures_fermi.xyz`, respectively.
## References
Please cite the paper below for using CP-MACE:

[Constant-Potential Machine Learning Force Field for the Electrochemical Interface](https://doi.org/10.1021/acs.jctc.5c00784)

Ruoyu Wang, Shaoheng Fang, Qixing Huang, Yuanyue Liu*, Journal of Chemical Theory and Computation, 2025, DOI: 10.1021/acs.jctc.5c00784

## Contact

If you encounter any issues or have questions, feel free to open an issue on this repository.
