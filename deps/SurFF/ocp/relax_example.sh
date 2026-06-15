#!/bin/bash

# Get the current working directory
current_dir=$(pwd)

# Print the current working directory
echo "The current working directory is: $current_dir"

traj_save_dir=traj/example
checkpoint_pth=checkpoints/2024-03-03-23-57-52/best_checkpoint.pt

relax_dataset_dir=data/example
python main.py \
--mode run-relaxations \
--config-yml configs/equiformer_v2_002_relax.yml \
--checkpoint $checkpoint_pth \
--amp \
--cpu \
--task.relax_opt.traj_dir=$traj_save_dir \
--task.relax_opt.maxstep=0.03 \
--task.relax_dataset.src=$relax_dataset_dir \
--task.relaxation_steps=200 \
--task.relaxation_fmax=0.05 \




