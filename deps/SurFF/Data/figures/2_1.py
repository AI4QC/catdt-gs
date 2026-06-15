import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline

results_dir_list = [
                    r"Data/dataset_generation/Model/save_dir/element_2_1",
                    r"Data/dataset_generation/Model/save_dir/element_2_2",
                    r"Data/dataset_generation/Model/save_dir/element_2_3",
                    r"Data/dataset_generation/Model/save_dir/element_2_4",
                    r"Data/dataset_generation/Model/save_dir/element_3_1",
                    r"Data/dataset_generation/Model/save_dir/element_3_2",
                    r"Data/dataset_generation/Model/save_dir/element_3_3",
                    r"Data/dataset_generation/Model/save_dir/element_3_4",
                    r"Data/dataset_generation/Model/save_dir/element_3_5",
                    r"Data/dataset_generation/Model/save_dir/element_3_6",
                    r"Data/dataset_generation/Model/save_dir/element_3_7",
]
# get all sub-folder in dir

def process_data(dir):

    data_pth = os.path.join(dir, "test.csv")

    data = pd.read_csv(data_pth)

    data_true = data.iloc[:, 0].values
    data_pred = data.iloc[:, 1:].values

    mean_pred = data_pred.mean(axis=1)
    std_pred = data_pred.std(axis=1)

    error = abs(mean_pred - data_true)
    mae = error.mean()

    num_pth = os.path.join(dir, "slab_data.csv")

    num_data = pd.read_csv(num_pth)
    # drop unconverged data
    num_data = num_data[num_data['converge'] == True]
    length = len(num_data)

    return mae, error, length


def plot_error_data_boxplot(error_data):
    # Create the box plot
    plt.boxplot(error_data)

    # Adding titles and labels
    plt.title('Error Data Distribution')
    plt.xlabel('Sublist Number')
    plt.ylabel('Error Values')
    plt.grid(True)

    # Show the plot
    plt.show()


def plot_error_data_violin(error_data, num_data):
    plt.figure(figsize=(10, 7), dpi=600)

    # change the font size
    plt.rcParams.update({'font.size': 12})

    # Create the violin plot
    # parts = plt.violinplot(error_data, showmeans=True, showmedians=False, showextrema=False, widths=0.8)

    # create the box plot
    parts = plt.boxplot(error_data, showmeans=True, widths=0.5, patch_artist=True,
                        showfliers=False,  # Don't show outliers
                        # change box line width
                        boxprops=dict(color='steelblue', facecolor='aliceblue', linewidth=2),
                        whiskerprops=dict(color='steelblue'),
                        capprops=dict(color='steelblue'),
                        medianprops=dict(color='aliceblue'),  # don't show median line
                        meanprops=dict(marker='o', markeredgecolor='black',
                                       markerfacecolor='steelblue'),
                                                              )  # Customize mean marker
    # Calculate means
    means = [np.mean(sublist) for sublist in error_data]

    # plot means as text
    for i in range(len(means)):
        plt.text(i + 1, means[i] + 0.0005, round(means[i], 4), horizontalalignment='center', fontweight='bold',
                 fontsize=12, color='steelblue')

    # X-axis values (1, 2, 3, ..., N)
    x_values = np.arange(1, len(error_data) + 1)
    # Interpolation for a smooth line
    xnew = np.linspace(x_values.min(), x_values.max(), 300)
    spl = make_interp_spline(x_values, means, k=2)
    smooth_means = spl(xnew)

    # Plot the smooth line
    plt.plot(xnew, smooth_means, color='steelblue', label='Error mean', linewidth=5, alpha=0.4)

    # Customize the violin plot appearance
    # for pc in parts['bodies']:
    #     pc.set_facecolor('deepskyblue')
        # pc.set_alpha(0.4)

    # Adding titles and labels
    plt.title('Active Learning Progess', fontsize=18, fontweight='bold')
    plt.xlabel('Batch', fontsize=18, fontweight='bold')
    plt.ylabel('Prediction Error\n (eV/Angstrom squared)', fontsize=18, fontweight='bold')
    # change x axis label
    plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
               ['2-1', '2-2', '2-3', '2-4', '3-1', '3-2', '3-3', '3-4', '3-5', '3-6', '3-7'])
    plt.ylim(0, 0.02)
    plt.tick_params(axis='both', which='major', labelsize=16)

    # add a horizontal grid to the plot,
    plt.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
    # add a vertical line between data 3 and 4
    plt.axvline(x=4.5, color='grey', linestyle='--', alpha=.5, lw=3)
    # add test to position
    plt.text(3, 0.0005, 'Bi-metallic', horizontalalignment='center', fontsize=18, color='grey', fontweight='bold')
    plt.text(6, 0.0005, 'Tri-metallic', horizontalalignment='center', fontsize=18, color='grey', fontweight='bold')

    # Create a proxy artist for the violin plot
    violin_proxy = Line2D([0], [0], color='deepskyblue', lw=6, label='Error Distribution', alpha=0.3)

    # SECONDARY Y-AXIS FOR ACCUMULATIVE DATA POINTS
    ax2 = plt.gca().twinx()

    # Generate some accumulative data points for demonstration
    # In a real scenario, this would be your actual accumulative data
    accumulative_data = np.cumsum(num_data)

    # Plotting the accumulative data points as an inverted bar plot
    bars = ax2.bar(range(1, len(num_data) + 1), accumulative_data, color='gray', alpha=0.5, width=0.7)
    ax2.set_ylim(0, 20000)
    # ax2 y-axis ticks, show only 0, 2500, 5000, 7500, 10000
    ax2.set_yticks([0, 2500, 5000, 7500, 10000,])

    # Label each bar with its corresponding accumulative data value
    j = 0
    for bar, value in zip(bars, accumulative_data):
        label_x_pos = j+0.65  # Adjust +0.1 or as needed for label positioning
        ax2.text(label_x_pos, bar.get_y() + bar.get_height() - 500, str(value), va='center', fontsize=12, fontweight='bold')
        j = j + 1

    # Inverting the secondary y-axis
    ax2.invert_yaxis()

    # Setting the label for the secondary y-axis
    ax2.set_ylabel('Accumulative Data Size', fontsize=18, fontweight='bold')
    ax2.tick_params(axis='y', labelsize=16)

    # Update legend to include all elements
    plt.legend(handles=[violin_proxy, Line2D([0], [0], lw=3, color='steelblue', label='Error Mean', alpha=0.6),
                        Line2D([0], [0], color='gray', label='Accumulative Data', linewidth=7,alpha=0.5)], loc='upper right')

    plt.tight_layout()
    # Show the plot
    plt.show()


mae_list = []
error_list = []
num_data = []

for dir in results_dir_list:

    mae, error, num = process_data(dir)
    mae_list.append(mae)
    error_list.append(error)
    num_data.append(num)

# plot error
plot_error_data_violin(error_list, num_data)
