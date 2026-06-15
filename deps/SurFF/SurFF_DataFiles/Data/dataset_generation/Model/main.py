import os
import time
from random import sample
from typing import Dict, List, Tuple
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
from torch import nn
from torch.optim.lr_scheduler import MultiStepLR, ReduceLROnPlateau, LambdaLR
from torch_geometric.loader import DataLoader

from data import GCN_Dataset, Normalizer, TrainingLogger
from model import GCN

import torch
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from CFG import CFG

from torch.cuda.amp import autocast, GradScaler


class Model_Trainer:
    """default model trainer"""

    def __init__(self, cfg: CFG):
        """provide configurations in cfg object"""

        self.cfg = cfg

        self.data(pth_train=self.cfg.trainset, pth_test=self.cfg.testset, pth_val=self.cfg.valset)
        self.init_model()

    def data(self, pth_train: [str], pth_test: [str], pth_val=None) -> None:
        """Generate dataloader for training & testing"""

        train_set = GCN_Dataset(pth_train)
        self.train_loader = DataLoader(train_set, batch_size=self.cfg.batch_size, shuffle=True)

        test_set = GCN_Dataset(pth_test)
        self.test_loader = DataLoader(test_set, batch_size=self.cfg.batch_size // 2, shuffle=False)

        if self.cfg.valset is not None:
            val_set = GCN_Dataset(pth_val)
            self.val_loader = DataLoader(val_set, batch_size=self.cfg.batch_size // 2, shuffle=False)

        sample_data_list = torch.tensor(
            [train_set[i].y for i in sample(range(len(train_set)), min(1000, len(train_set)))])
        self.norm_output = Normalizer(sample_data_list)

        print(self.norm_output.mean, self.norm_output.std)

    def init_model(self) -> None:
        if self.cfg.model_path is None:
            self.model = GCN(atom_feature=self.cfg.atom_feature,
                             edge_dim=self.cfg.num_edge_fea,
                             conv_feature=self.cfg.conv_feature,
                             fc_feature=self.cfg.fc_feature,
                             out_feature=self.cfg.out_feature,
                             num_fc_layers=self.cfg.num_fc_layers,
                             num_conv_layers=self.cfg.num_conv_layers,
                             conv_channels=self.cfg.conv_channels,
                             conv_block=self.cfg.conv_block,
                             dropout=self.cfg.dropout,
                             ).to(self.cfg.device)

            initialize_weights(self.model)

        else:  # load an existing model if path is provided
            print("Loading model from: ", self.cfg.model_path)
            save_dict = torch.load(self.cfg.model_path)
            cfg = save_dict["cfg"]

            # cfg.trainset = self.cfg.trainset
            # cfg.testset = self.cfg.testset
            # cfg.batch_size = self.cfg.batch_size
            # cfg.model_path = self.cfg.model_path
            # self.cfg = cfg

            self.model = GCN(atom_feature=cfg.atom_feature,
                             edge_dim=cfg.num_edge_fea,
                             conv_feature=cfg.conv_feature,
                             fc_feature=cfg.fc_feature,
                             out_feature=cfg.out_feature,
                             num_fc_layers=cfg.num_fc_layers,
                             num_conv_layers=cfg.num_conv_layers,
                             conv_channels=cfg.conv_channels,
                             conv_block=cfg.conv_block,
                             dropout=cfg.dropout,
                             layer_fea=cfg.layer_fea).to(self.cfg.device)

            self.model.load_state_dict(save_dict["model"])
            self.norm_output = save_dict["norm_output"]
            # self.optimizer.load_state_dict(save_dict['optimizer'])

        # param_groups = [
        #     {'params': [p for n, p in self.model.named_parameters() if 'fc' not in n], 'lr': self.cfg.lr},
        #     {'params': [p for n, p in self.model.named_parameters() if 'fc' in n], 'lr': self.cfg.lr * 0.1}
        # ]

        # self.scheduler = WarmupThenDecaySchedule(base_lr=self.cfg.lr, warmup_steps=5000)

        # gaussian_params = list(map(id, self.model.edge_embed.parameters()))
        # rest_params = filter(lambda p: id(p) not in gaussian_params, self.model.parameters())

        # Use weight decay on all parameters EXCEPT the gaussian parameters
        # self.optimizer = torch.optim.AdamW([
        #     {'params': self.model.edge_embed.parameters(), 'weight_decay': 0.0},
        #     {'params': rest_params, 'weight_decay': self.cfg.weight_decay}],
        #     lr=self.cfg.lr, betas=(0.9, 0.999), eps=1e-6, )

        self.optimizer = torch.optim.AdamW(self.model.parameters(),
                                           betas=(0.9, 0.998),
                                           eps=1e-6,
                                           weight_decay=self.cfg.weight_decay,
                                           lr=self.cfg.lr, )

        if self.cfg.model_path is not None:
            self.optimizer.load_state_dict(save_dict['optimizer'])

        self.scheduler = ReduceLROnPlateau(self.optimizer,
                                           factor=self.cfg.gamma,
                                           verbose=True,
                                           patience=self.cfg.lr_patience,
                                           )

        self.warmup_steps = 5000
        self.warm_up = LambdaLR(self.optimizer, self.warmup_lr_schedule)

        self.logger = TrainingLogger()

    def train(self):

        self.scaler = GradScaler()

        # reset_all_weights(self.model)
        start_time = time.time()

        self.warmup_steps = 2000
        loss_f = torch.nn.L1Loss()

        step = 0
        weight_update_step = self.cfg.step
        accumulation_steps = self.cfg.accumulation_batch  # You can adjust this value

        total_loss = []
        num_data = 0

        for epoch in range(self.cfg.epoch):

            for data in self.train_loader:
                # if data.batch.shape[0]> 1000:
                #     print(f"WARNING: large_batch_find {weight_update_step}, {data.batch.shape[0]}, skipping batch")
                #     continue

                try:
                    batch_loss = self.train_data(data, loss_f)

                    total_loss.append(batch_loss * len(data))  # Normalize the loss
                    num_data += len(data)

                    step += 1

                except RuntimeError as e:
                    print(f"WARNING:{weight_update_step} {step}, skipping batch, {e}")
                    print(data)
                    torch.cuda.empty_cache()
                    self.optimizer.zero_grad()
                    step = 0
                    continue

                if step % accumulation_steps == 0:

                    # lr = self.scheduler(weight_update_step)
                    # for param_group in self.optimizer.param_groups:
                    #     param_group['lr'] = lr

                    self.scaler.unscale_(self.optimizer)
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), 5.0)
                    self.scaler.step(self.optimizer)
                    self.scaler.update()
                    self.optimizer.zero_grad()

                    if weight_update_step < self.warmup_steps:
                        self.warm_up.step()

                    if weight_update_step % self.cfg.save_freq == self.cfg.save_freq - 1:
                        self.save_model(
                            save_name=f"epoch_{epoch + 1:03d}_weight_update_{weight_update_step + 1:08d}.pth")

                    if weight_update_step % self.cfg.test_freq == self.cfg.test_freq - 1 and self.cfg.test_freq != 0:
                        torch.cuda.empty_cache()
                        test_mae, _, _ = self.test(self.test_loader, info='test')
                        self.logger.log_testing_loss(test_mae)
                        if weight_update_step > self.warmup_steps:
                            self.scheduler.step(test_mae)
                        # print(self.model.edge_embed.means.weight.max(), self.model.edge_embed.means.weight.min())

                    if weight_update_step % self.cfg.print_freq == self.cfg.print_freq - 1 and self.cfg.print_freq != 0:
                        train_loss = sum(total_loss) / num_data
                        lr = self.optimizer.param_groups[0]['lr']

                        current_memory_allocated = torch.cuda.memory_reserved(self.cfg.device) / (
                                1024 ** 3)  # Convert to GB
                        max_memory = torch.cuda.get_device_properties(self.cfg.device).total_memory / (
                                1024 ** 3)  # Convert to GB

                        print(
                            f"{format_elapsed_time(time.time() - start_time)}: {epoch + 1} epoch {weight_update_step + 1} weight update, "
                            f"train loss = {train_loss:.4f}, lr = {lr:.6f}, cuda: {current_memory_allocated:.1f}/{max_memory:.1f}")

                        self.logger.log_training_loss(train_loss)
                        self.logger.log_other_info("learning_rate", lr)
                        self.logger.record_epoch(epoch + 1, weight_update_step + 1)
                        self.logger.save_logs(filename=f"log_{self.cfg.save_dir}.json")

                        total_loss = []
                        num_data = 0
                        # print(f"kernel {self.model.edge_embed.means.weight.max()}, {self.model.edge_embed.means.weight.min()}")
                        # print(f"mul {self.model.edge_embed.mul.weight.max()}, {self.model.edge_embed.mul.weight.min()}")
                        # print(f"bias {self.model.edge_embed.bias.weight.max()}, {self.model.edge_embed.bias.weight.min()}")

                    weight_update_step += 1

                    lr = self.optimizer.param_groups[0]['lr']
                    if lr < 1e-5 and weight_update_step > self.warmup_steps:
                        print("Reach minimum learning rate.")
                        print(f"Test starting: save_dir {self.cfg.save_dir}")
                        mae, true_dict, pred_dict = self.test(self.test_loader)

                        self.save_model(save_name=self.cfg.save_name)

                        root_dir = os.path.dirname(os.path.realpath(__file__))
                        save_dir = os.path.join(root_dir, "save_dir", self.cfg.save_dir)
                        self.show_performance(true_dict, pred_dict, save_dir=save_dir)

                        self.train_loader.dataset.close_env()
                        self.test_loader.dataset.close_env()
                        return mae, true_dict, pred_dict

        print(f"Test starting: save_dir {self.cfg.save_dir}")
        mae, true_dict, pred_dict = self.test(self.test_loader)

        self.save_model(save_name=self.cfg.save_name)

        root_dir = os.path.dirname(os.path.realpath(__file__))
        save_dir = os.path.join(root_dir, "save_dir", self.cfg.save_dir)
        self.show_performance(true_dict, pred_dict, save_dir=save_dir)

        self.train_loader.dataset.close_env()
        self.test_loader.dataset.close_env()

        return mae, true_dict, pred_dict

    def train_data(self, data, loss_f):
        self.model.train()

        with autocast():
            num = [torch.sum(data.tags[data.ptr[i]:data.ptr[i + 1]]) for i in range(len(data.ptr) - 1)]
            num = torch.stack(num).unsqueeze(-1).to(self.cfg.device)

            train_out = data.y.to(self.cfg.device)
            pred_out = self.model(data.to(self.cfg.device))
            pred_out = pred_out / num

        loss = loss_f(pred_out, train_out)
        # loss = torch.abs(pred_out / (train_out + 1e-2) - 1).mean()

        self.scaler.scale(loss / self.cfg.accumulation_batch).backward()  # Normalize the loss for gradient accumulation

        # loss.backward()

        return loss.item()

    def test(self, test_loader=None, info='test', eval_num: int = 1):

        if test_loader is None:
            test_loader = self.test_loader

        self.model.eval()
        # for m in self.model.modules():
        #     if m.__class__.__name__.startswith('Dropout'):
        #         m.train()

        true_dict = []
        pred_dict = []
        std_dict = []

        for data in tqdm(test_loader):
            test_out = data.y
            num = [torch.sum(data.tags[data.ptr[i]:data.ptr[i + 1]]) for i in range(len(data.ptr) - 1)]
            num = torch.stack(num).unsqueeze(-1).to(self.cfg.device)

            batch_list = []
            for _ in range(eval_num):
                pred_out = self.model(data.to(self.cfg.device))  # [N, output_size]
                pred_out = (pred_out/num).detach().cpu().numpy()
                batch_list.append(pred_out)

            batch_list = np.hstack(batch_list)  # [N, eval_num, output_size]
            pred_mean = np.mean(batch_list, axis=1, keepdims=True)  # [N, output_size]
            pred_std = np.std(batch_list, axis=1, keepdims=True)

            true_dict.append(test_out.detach().cpu().numpy())
            pred_dict.append(pred_mean)
            std_dict.append(pred_std)

        true_dict = np.vstack(true_dict)  # [N, output_size] where N is total test size
        pred_dict = np.vstack(pred_dict)  # [N, output_size]
        std_dict = np.vstack(std_dict)

        lb = pred_dict - std_dict
        ub = pred_dict + std_dict
        percent_in_range = np.sum((lb < true_dict) & (true_dict < ub)) / len(true_dict)

        abs_ev = 0.01
        percent_in_abs = np.sum(np.abs(true_dict - pred_dict) < abs_ev) / len(true_dict)

        mae = np.mean(np.abs(true_dict - pred_dict))
        print(f"{info} mae = {mae:.4f}, percent in 1 sigma: {percent_in_range * 100:.2f}%, "
              f"percent in {abs_ev}ev: {percent_in_abs * 100:.2f}%")
        return mae, true_dict, pred_dict

    def predict(self, eval_num: int = 1):

        self.model.eval()
        for m in self.model.modules():
            if m.__class__.__name__.startswith('Dropout'):
                m.train()

        pred_mean = []
        pred_std = []

        for data in tqdm(self.val_loader):

            data.to(self.cfg.device)
            batch_list = []
            for _ in range(eval_num):
                pred_out = self.model(data)  # [N, output_size]
                pred_out = pred_out.detach().cpu().numpy()
                batch_list.append(pred_out)

            batch_list = np.hstack(batch_list)  # [N, eval_num, output_size]
            batch_mean = np.mean(batch_list, axis=1)  # [N, output_size]
            batch_std = np.std(batch_list, axis=1)  # [N, output_size]

            pred_mean.append(batch_mean)
            pred_std.append(batch_std)

        pred_mean = np.concatenate(pred_mean, axis=0)  # [N, output_size] where N is total test size
        pred_std = np.concatenate(pred_std, axis=0)  # [N, output_size]

        self.val_loader.dataset.close_env()

        return pred_mean, pred_std

    def show_performance(self, true_dict, pred_dict, save_dir):
        """ calculate mea and r2 values based on model prediction
            also generate parity plot for each output"""
        normalized_mae = np.mean(np.abs(true_dict - pred_dict), axis=0)

        r2_list = []

        for i in range(self.cfg.out_feature):
            r2 = self.rsquared(true_dict[:, i], pred_dict[:, i])
            r2_list.append(r2)

        data = {'var_name': 'model_prediction', 'r2': r2_list, "normalized_mae": normalized_mae}
        resutls = pd.DataFrame(data)
        pd_save_path = os.path.join(save_dir, "class_performance.csv")
        resutls.to_csv(pd_save_path)
        print(resutls)

        for i in range(self.cfg.out_feature):
            x = true_dict[:, i]
            y = pred_dict[:, i]
            xy = np.vstack([x, y])
            d = gaussian_kde(xy)(xy)

            idx = d.argsort()
            x, y, d = x[idx], y[idx], d[idx]

            plt.scatter(x, y, c=d, s=1, cmap="rainbow")
            plt.xlabel("True values")
            plt.ylabel("Model Prediction")
            plt.title(f"MAE={normalized_mae.tolist()[0]:.4f}")
            plt.axis('square')
            plt.axline((0, 0), (1, 1), linewidth=1, color='black')

            img_save_path = os.path.join(save_dir, "{}.png".format('model_prediction'))

            plt.savefig(img_save_path, dpi=600)
            plt.clf()

    @staticmethod
    def rsquared(x: np.ndarray, y: np.ndarray) -> float:
        """ Return R^2 where x and y are array-like."""

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        return r_value ** 2

    def save_model(self, save_name: str):
        save_dict = {"norm_output": self.norm_output,
                     "model": self.model.state_dict(),
                     'optimizer': self.optimizer.state_dict(),
                     "cfg": self.cfg,
                     }

        root_dir = os.path.dirname(os.path.realpath(__file__))
        save_dir = os.path.join(root_dir, "save_dir", self.cfg.save_dir)
        save_path = os.path.join(save_dir, save_name)

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        torch.save(save_dict, save_path)

    def warmup_lr_schedule(self, step):
        if step < self.warmup_steps:
            return float(step) / float(max(1, self.warmup_steps))
        return 1.0


def format_elapsed_time(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02d}h {int(minutes):02d}m {int(seconds):02d}s"


def reset_all_weights(model: nn.Module) -> None:
    @torch.no_grad()
    def weight_reset(m: nn.Module):
        reset_parameters = getattr(m, "reset_parameters", None)
        if callable(reset_parameters):
            m.reset_parameters()

    model.apply(fn=weight_reset)


def initialize_weights(model):
    for module in model.modules():
        if isinstance(module, nn.Linear):
            nn.init.xavier_uniform_(module.weight)
            # if module.bias is not None:
            #     nn.init.uniform_(module.bias, -0.01, 0.01)


class WarmupThenDecaySchedule:
    def __init__(self, base_lr, warmup_steps=5000):
        self.base_lr = base_lr
        self.warmup_steps = warmup_steps

    def __call__(self, step_num):

        if step_num <= self.warmup_steps:
            factor = step_num / self.warmup_steps
        else:
            factor = max(0.05, 1 - (step_num - self.warmup_steps) / 1000000)

        return self.base_lr * factor
