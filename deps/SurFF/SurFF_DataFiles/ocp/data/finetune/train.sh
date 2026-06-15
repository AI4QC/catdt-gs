#!/bin/bash

python main.py \
    --mode train \
    --config-yml configs/ef/equiformer_v2_002_finetune.yml \
    --amp \
    --print-every 100 \
    --local_rank 3 \
    --optim.max_epochs=5 \
    --optim.lr_initial=1e-5 \
    --checkpoint checkpoints/2024-03-03-23-57-52/best_checkpoint.pt \











