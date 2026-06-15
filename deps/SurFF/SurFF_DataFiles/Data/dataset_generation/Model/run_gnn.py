import numpy as np
import pandas as pd

from CFG import CFG
from main import Model_Trainer

cfg = CFG()

cfg.trainset = [
    r"Data/dataset_generation/LMDB/element_2_1.lmdb",
    r"Data/dataset_generation/LMDB/element_2_2.lmdb",
    r"Data/Data/dataset_generation/LMDB/element_2_3.lmdb",
    r"Data/dataset_generation/LMDB/element_2_4.lmdb",
    r"Data/dataset_generation/LMDB/element_3_1.lmdb",
    r"Data/dataset_generation/LMDB/element_3_2.lmdb",
    r"Data/dataset_generation/LMDB/element_3_3.lmdb",
    r"Data/dataset_generation/LMDB/element_3_4.lmdb",
    r"Data/dataset_generation/LMDB/element_3_5.lmdb",
    r"Data/dataset_generation/LMDB/element_3_6.lmdb",
    r"Data/dataset_generation/LMDB/element_3_7.lmdb",
    ]

cfg.testset = [
    r"Data/dataset_generation/LMDB/element_2_0.lmdb",
    r"Data/dataset_generation/LMDB/element_3_0.lmdb",
    ]

cfg.valset = [
    r"Data/dataset_generation/LMDB/element_2_all.lmdb",
    # r"Data/dataset_generation/LMDB/element_3_all.lmdb",
    ]

cfg.atom_feature = 96
cfg.num_edge_fea = 128
cfg.conv_feature = 256
cfg.conv_channels = 8
cfg.fc_feature = 256
cfg.out_feature = 1
cfg.num_conv_layers = 8
cfg.num_fc_layers = 1
cfg.conv_block = 1
cfg.layer_fea = False

cfg.batch_size = 16
cfg.accumulation_batch = 1
cfg.epoch = 500
cfg.lr = 1e-4
cfg.weight_decay = 1e-3
cfg.lr_patience = 5
cfg.dropout = 0
cfg.gamma = 0.1

cfg.print_freq = 500
cfg.test_freq = 500
cfg.save_freq = 100000

cfg.device = "cuda"

N = 5
cfg.save_dir = "element_3_7"
test_results = []
val_results = []

for i in range(N):
    cfg.save_name = f"{cfg.save_dir}_{i}.pth"

    trainer = Model_Trainer(cfg)
    print(f"model parameters: {sum(p.numel() for p in trainer.model.parameters())}")
    test_result = trainer.train()
    test_results.append(test_result[2])
    true_results = test_result[1]

    pred_mean, _ = trainer.predict()
    val_results.append(pred_mean)

test_results = np.array(test_results).squeeze().T
test_save = np.concatenate([true_results, test_results], axis=-1)
test_save = pd.DataFrame(test_save)
test_save.to_csv(f"Model/save_dir/{cfg.save_dir}/test.csv", index=False)

test_mean = np.mean(test_results, axis=-1)
test_mae = np.mean(np.abs(test_mean - true_results.squeeze()))

val_results = np.array(val_results).T
val_mean = np.mean(val_results, axis=-1).reshape(-1, 1)
val_std = np.std(val_results, axis=-1).reshape(-1, 1)

val_save = np.concatenate([val_mean, val_std], axis=-1)
val_save = pd.DataFrame(val_save)
val_save.to_csv(f"Model/save_dir/{cfg.save_dir}/val.csv", index=True)
