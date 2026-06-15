#!/bin/bash

python main.py \
    --mode train \
    --config-yml configs/ef/equiformer_v2_002.yml \
    --amp \
    --print-every 100 \
    --local_rank 6 \
    --model.num_layers=6 \
    --optim.max_epochs=20 \







