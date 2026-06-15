#!/bin/bash

traj_save_dir=traj/2024-04-12-12-11-44/exp
checkpoint_pth=checkpoints/2024-04-12-12-11-44/best_checkpoint.pt
relax_dataset_dir=data/is2re

python main.py \
--mode run-relaxations \
--config-yml configs/equiformer_v2_002_relax.yml \
--checkpoint $checkpoint_pth \
--cpu \
--task.relax_opt.traj_dir=$traj_save_dir \
--task.relax_opt.maxstep=0.03 \
--task.relax_dataset.src=$relax_dataset_dir \
--task.relaxation_steps=300 \
--task.relaxation_fmax=0.03 \




