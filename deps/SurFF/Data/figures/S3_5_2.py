# visualize surface energy in csv file

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os

# fontsize 20, dpi 300
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.dpi': 300})


def plot_energy_distribution(csv_pth_dir, name):
    # plot max force distribution

    all_data = pd.DataFrame()

    # get all csv files in the directory
    files = [os.path.join(csv_pth_dir, file) for file in os.listdir(csv_pth_dir) if file.endswith('.csv')]
    for file in files:
        print(f"Processing {file}")
        data = pd.read_csv(file)
        all_data = all_data.append(data, ignore_index=True)

    # drop unconverged data
    all_data = all_data[all_data['converge'] == True]

    print(f"Total number of data points: {len(all_data)}")

    surface_energies = all_data['surface_energy']

    plt.figure(figsize=(7, 5))
    # remove 3-sigma outliers
    sigma = np.std(surface_energies)
    mean = np.mean(surface_energies)
    surface_energies = surface_energies[(surface_energies > mean - 3 * sigma) & (surface_energies < mean + 3 * sigma)]
    weights = np.ones_like(surface_energies) / len(surface_energies)

    plt.hist(surface_energies, bins=100, weights=weights, color='gray')

    plt.xlabel('Surface Energy (eV/$\AA^2$)')
    plt.ylabel('Data Fraction')
    plt.xlim(0, 0.25)
    plt.ylim(0, 0.03)
    plt.title(f'{name}: Surface Energy Distribution')

    plt.text(0.9, 0.9, f'Total Surface Energy Data\n{len(all_data)}',
             ha='right', va='top', transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig(f'Data/figures/fig_save/S3_5_{name}_surface_energy.png')


# plot max force distribution: train
csv_pth_dir = r"Data/Surface/Surface_Energy/data/trainset"

plot_energy_distribution(csv_pth_dir, name='Train')

# plot max force distribution: test AL
csv_pth_dir = r"Data/Surface/Surface_Energy/data/AL_testset"

plot_energy_distribution(csv_pth_dir, name='Test_AL')

# plot max force distribution: test ID
csv_pth_dir = r"Data/Surface/Surface_Energy/data/ID_testset"

plot_energy_distribution(csv_pth_dir, name='Test_ID')


# plot max force distribution: test OOD
csv_pth_dir = r"Data/Surface/Surface_Energy/data/OOD_testset"

plot_energy_distribution(csv_pth_dir, name='Test_OOD')

