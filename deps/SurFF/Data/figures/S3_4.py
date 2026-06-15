# visualize data the stat for additional data added for each active learning iteration
import pdb

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, gridspec
from scipy.stats import gaussian_kde
from tqdm import tqdm


# a figure of 1x3 subplots
# 1. the distribution of all uncertainty values and sampled uncertainty values
# 2. the model prediction vs. true value for the sampled data

def plot(result_folder):
    all_prediction = pd.read_csv(f"{result_folder}/val.csv")
    test_prediction = pd.read_csv(f"{result_folder}/test.csv")
    sampled_data = pd.read_csv(f"{result_folder}/sampled_slab_data.csv")

    fig = plt.figure(figsize=(8, 4), dpi=600)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])  # Adjust the ratio here

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # 1. the distribution of all uncertainty values and sampled uncertainty values
    all_std = all_prediction.iloc[:, 2].values * 1000  # y1
    sampled_std = sampled_data["weight"].values * 1000  # y2

    ax1 = plt.subplot(1, 2, 1, )
    ax1.boxplot([all_std, sampled_std], labels=["All", "Sampled"],
                # set flier properties, transparency
                flierprops=dict(marker='o', color='black', markersize=3, alpha=0.1),
                )
    # set y lim
    ax1.set_ylim([0, 1000])
    ax1.set_title("Uncertainty Distribution")
    ax1.set_ylabel("Uncertainty (meV/$\AA^2$)")

    # 2. the model prediction vs. true value for the sampled data
    y = test_prediction.iloc[:, 0].values * 1000
    x = np.mean(test_prediction.iloc[:, 1:].values, axis=1) * 1000

    ax2 = plt.subplot(1, 2, 2)

    # density x-y parity plot
    # Calculate MAE
    mae = np.mean(np.abs(x - y))

    # Calculate the point density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    scatter = ax2.scatter(x, y, c=z, s=10, cmap='viridis', alpha=0.5)
    # plt.colorbar(scatter, ax=ax2, label='Density')

    # Add the diagonal line
    lims = [
        np.min([ax2.get_xlim(), ax2.get_ylim()]),  # find the lower limit of x and y
        np.max([ax2.get_xlim(), ax2.get_ylim()]),  # find the upper limit of x and y
    ]
    ax2.plot(lims, lims, 'k-', alpha=0.75, zorder=0, color='black')
    ax2.set_aspect('equal')
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)

    # Annotate MAE
    ax2.text(np.max(lims) * 0.65, np.min(lims) * 1.3 + 0.01, f'MAE: {mae:.1f}', fontsize=12, color='black')

    # Labels and titles
    ax2.set_xlabel('True Values (meV/$\AA^2$)')
    ax2.set_ylabel('Predicted Values (meV/$\AA^2$)')
    ax2.set_title('Accuracy of Model Prediction')
    save_name = result_folder.split("/")[-1]
    # global title
    plt.suptitle(f"Iteration: {save_name}")

    plt.tight_layout()
    plt.savefig(f"Data/figures/fig_save/3_4_{save_name}.png")


result_folders = [r"Data/dataset_generation/Model/save_dir/element_2_3",
                  ]

for result_folder in tqdm(result_folders):
    plot(result_folder)
    print(f"Plotting {result_folder} done!")
