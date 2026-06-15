import pickle
import random

import lmdb
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from torch.nn.functional import one_hot
# from torch.utils.data import Dataset
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import DataLoader
import torch
import random

import time
import json
import os

import time
import json
import os


class TrainingLogger:
    def __init__(self, log_dir="logs"):
        self.log_dir = log_dir
        self.training_loss = []
        self.testing_loss = []
        self.other_info = {}
        self.start_time = time.time()
        self.epoch = 0
        self.iter = 0

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

    def log_training_loss(self, loss):
        self.training_loss.append(float(loss))

    def log_testing_loss(self, loss):
        self.testing_loss.append(float(loss))

    def log_other_info(self, key, value):
        if key not in self.other_info:
            self.other_info[key] = []
        self.other_info[key].append(float(value))

    def record_epoch(self, epoch=None, iter=None):
        self.epoch = epoch
        self.iter = iter

    def save_logs(self, filename="logs.json"):
        log_data = {
            "training_loss": self.training_loss,
            "testing_loss": self.testing_loss,
            "other_info": self.other_info,
            "epochs": self.epoch,
            "iterations": self.iter,
            "elapsed_time": time.time() - self.start_time,
        }

        # Overwrite existing logs by opening the file in write mode ('w')
        with open(os.path.join(self.log_dir, filename), "w") as f:
            json.dump(log_data, f, indent=2)

    def __str__(self):
        return f"TrainingLogger(training_loss={self.training_loss}, testing_loss={self.testing_loss}, other_info={self.other_info}, current_epoch={self.current_epoch})"


def plot_losses_from_file(log_file, save_plot=False, plot_filename_prefix="losses_plot"):
    with open(log_file, "r") as f:
        log_data = json.load(f)

    training_loss = log_data["training_loss"]
    testing_loss = log_data["testing_loss"]

    lr = log_data["other_info"]["learning_rate"]

    total_epochs = log_data["epochs"]
    total_iterations = log_data["iterations"]

    epochs = list(range(1, total_epochs + 1))

    # Plot training loss in a separate figure

    plt.figure(figsize=(12, 10))

    plt.subplot(2, 2, 1)
    plt.plot(np.linspace(1, len(training_loss), len(training_loss)) * total_iterations / len(training_loss),
             training_loss,
             label="Training Loss")
    # plt.xlabel("Epochs")
    plt.ylabel("Loss")
    plt.yscale("log")

    if min(training_loss) < 1:
        plt.ylim(top=1)
    plt.title("Training Loss")
    plt.legend()
    plt.title(f"(Epochs: {total_epochs}, Iterations: {total_iterations})")

    plt.subplot(2, 2, 2)
    plt.plot(np.linspace(1, len(testing_loss), len(testing_loss)) * total_iterations / len(testing_loss),
             testing_loss,
             label="Testing MAE")

    min_loss_value = min(testing_loss)
    min_loss_index = testing_loss.index(min_loss_value)
    plt.scatter((1 + min_loss_index) * total_iterations / len(testing_loss), min_loss_value, label="Min Loss",
                color="red")
    # plt.xlabel("Epochs")
    plt.ylabel("MAE")
    plt.yscale("log")
    plt.ylim(top=1)
    plt.xlim(0)
    plt.title("Testing MAE")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(np.linspace(1, len(lr), len(lr)) * total_iterations / len(lr),
             lr,
             label="Learning Rate")
    # plt.xlabel("Epochs")
    plt.ylabel("Learning Rate")
    # set y axis to log scale
    plt.yscale("log")
    plt.title("Learning Rate")
    plt.legend()
    plt.title(f"(Epochs: {total_epochs}, Iterations: {total_iterations})")

    plt.tight_layout()

    if save_plot:
        plt.savefig(f"{plot_filename_prefix}_testing.png")
    else:
        plt.show()

    plt.clf()


class Normalizer:
    """ variable normalizer"""

    def __init__(self, var):
        self.mean = torch.mean(var, dim=0)
        self.std = torch.std(var, dim=0)
        self.max, _ = torch.max(var, dim=0)
        self.min, _ = torch.min(var, dim=0)
        self.dif = self.max - self.min

        # self.std = 1.0

    def norm(self, var):
        norm_var = var / self.std

        # norm_var = (var - self.mean) / self.std
        # norm_var = (var - self.min)/self.dif
        return norm_var

    def denorm(self, var):
        denorm_var = var * self.std
        # denorm_var = var * self.std + self.mean
        # denorm_var = var * self.dif + self.min
        return denorm_var


class GaussianDistance(object):
    """
    Expands the distance by Gaussian basis.
    Unit: angstrom
    """

    def __init__(self, dmin, dmax, step, var=None):
        """
        Parameters
        ----------
        dmin: float
          Minimum interatomic distance
        dmax: float
          Maximum interatomic distance
        step: float
          Step size for the Gaussian filter
        """
        assert dmin < dmax
        assert dmax - dmin > step
        self.filter = np.arange(dmin, dmax + step, step)
        if var is None:
            var = step
        self.var = var

    def expand(self, distances):
        """
        Apply Gaussian disntance filter to a numpy distance array
        Parameters
        ----------
        distance: np.array shape n-d array
          A distance matrix of any shape
        Returns
        -------
        expanded_distance: shape (n+1)-d array
          Expanded distance matrix with the last dimension of length
          len(self.filter)
        """
        return np.exp(-(distances[..., np.newaxis] - self.filter) ** 2 /
                      self.var ** 2)


class GCN_Dataset(Dataset):
    def __init__(self, path: [str]):
        super().__init__()

        self.env_list = []
        self.l_list = []

        for p in path:
            env = lmdb.open(p, subdir=False, meminit=False, map_async=True,)
            L = pickle.loads(env.begin().get(key='length'.encode('ascii')))

            self.env_list.append(env)
            self.l_list.append(L)

        self.length = sum(self.l_list)

        # self.gdf = GaussianDistance(0, 6, 6 / (num_edge_fea - 1))

    def len(self):
        return self.length

    def get(self, idx):

        for i, l in enumerate(self.l_list):
            if idx < sum(self.l_list[:(i+1)]):
                env_i = i
                env_idx = idx - sum(self.l_list[:i])
                break

        data = pickle.loads(self.env_list[env_i].begin().get(key=f'{env_idx}'.encode('ascii')))
        # data.x = one_hot(torch.atleast_1d(data.x.type(torch.long).squeeze()), num_classes=96).type(torch.float)

        # edge_attr = self.gdf.expand(data.edge_attr.numpy()).squeeze()
        # data.edge_attr = torch.from_numpy(edge_attr).type(torch.float)

        return data

    def close_env(self):
        for env in self.env_list:
            env.close()
