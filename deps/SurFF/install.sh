#!/bin/bash

# Ensure the script exits on any error
set -e

# Create a Conda environment named 'surff' with Python 3.9
conda create -n surff python=3.9 -y

# Activate the Conda environment
source activate surff

# Upgrade pip to the latest version
pip install --upgrade pip

# Install PyTorch with CUDA 11.8
pip install torch==2.1.0+cu118 -f https://download.pytorch.org/whl/torch_stable.html

# Install PyTorch Geometric and related libraries
pip install torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu118.html

# Install other dependencies
pip install ase==3.22.1 black==22.3.0 e3nn==0.4.4 matplotlib numba orjson \
             pre-commit==2.10.* pytest lmdb pyyaml submitit \
             syrupy==3.0.6 tensorboard tqdm wandb pymatgen==2023.5.10 requests \
             flask

echo "Environment 'surff' created and required packages installed successfully."
