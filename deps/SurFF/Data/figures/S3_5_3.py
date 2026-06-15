# visualize graph data SP energy in json file

import json
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# fontsize 20, dpi 300
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.dpi': 300})


def plot_force_distribution(json_pth, name, fraction=1.0):
    # plot max force distribution

    all_images = []

    for json_pth in json_pth:
        with open(json_pth, 'r') as f:
            print(f"Processing {json_pth}")
            images = json.load(f)
            all_images.extend(images)

    print(f"Total number of images: {len(all_images)}")

    expanded_images = [image for images in tqdm(all_images) for image in images[:int(len(images) * fraction)]]

    plt.figure(figsize=(7, 5))
    energy = np.array([image['energy'] for image in tqdm(expanded_images)])
    ref_energy = np.array([image['ref_energy'] for image in tqdm(expanded_images)])
    lattice = np.array([image['lattice'] for image in tqdm(expanded_images)])

    area = np.linalg.norm(np.cross(lattice[:, 0], lattice[:, 1]), axis=1)

    surface_energy = (energy - ref_energy) / (area * 2)  # eV/A^2

    # remove 3-sigma outliers
    sigma = np.std(surface_energy)
    mean = np.mean(surface_energy)
    surface_energy = surface_energy[(surface_energy > mean - 3 * sigma) & (surface_energy < mean + 3 * sigma)]
    weights = np.ones_like(surface_energy) / len(surface_energy)

    plt.hist(surface_energy, bins=100, weights=weights, color='gray')
    if fraction < 1.0:
        plt.xlabel(f'Energy (eV/$\AA^2$)\nFirst {int(fraction * 100)}% of Relaxation Trajectory')
    else:
        plt.xlabel('Energy (eV/$\AA^2$)')
    plt.ylabel('Data Fraction')
    plt.xlim(0, 0.3)
    plt.ylim(0, 0.03)
    plt.title(f'{name}: Energy Distribution')

    # add number of data points to upper right corner
    plt.text(0.9, 0.9, f'Total SP Data\n{len(expanded_images)}', ha='right', va='top', transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig(f'Data/figures/fig_save/S3_5_{name}_energy_{fraction}.png')


# plot max force distribution: train
json_pth = [r"Data/Surface/Traj/all_images_element_1_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_2_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_ref_energy.json", ]

plot_force_distribution(json_pth, name='Train', fraction=1)

# plot max force distribution: test AL
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_al_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_al_ref_energy.json", ]

plot_force_distribution(json_pth, name='Test_AL', fraction=1)

# plot max force distribution: test ID
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_id_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_id_ref_energy.json",
            ]

plot_force_distribution(json_pth, name='Test_ID', fraction=1)

# plot max force distribution: test OOD
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_ood_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_ood_ref_energy.json", ]

plot_force_distribution(json_pth, name='Test_OOD', fraction=1)
