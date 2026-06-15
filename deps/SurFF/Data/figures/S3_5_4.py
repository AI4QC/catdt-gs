# visualize graph data trajectory energy difference between relaxed and unrelaxed structures

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

    first_image = [images[0] for images in all_images]
    last_image = [images[-1] for images in all_images]

    plt.figure(figsize=(7, 5))

    first_energy = np.array([image['energy'] for image in tqdm(first_image)])
    last_energy = np.array([image['energy'] for image in tqdm(last_image)])

    lattice = np.array([image['lattice'] for image in tqdm(first_image)])

    area = np.linalg.norm(np.cross(lattice[:, 0], lattice[:, 1]), axis=1)

    surface_energy = (first_energy - last_energy) / (area * 2)  # eV/A^2

    # remove 3-sigma outliers
    sigma = np.std(surface_energy)
    mean = np.mean(surface_energy)
    surface_energy = surface_energy[(surface_energy > 0) & (surface_energy < mean + 3 * sigma)]
    weights = np.ones_like(surface_energy) / len(surface_energy)

    plt.hist(surface_energy, bins=100, color='gray', weights=weights)
    if fraction < 1.0:
        plt.xlabel(f'Energy (eV/$\AA^2$)\nFirst {int(fraction*100)}% of Relaxation Trajectory')
    else:
        plt.xlabel('Energy (eV/$\AA^2$)')
    plt.ylabel('Data Fraction')
    plt.xlim(0, 0.05)
    plt.ylim(0, 0.08)
    plt.title(f'{name}: $\Delta$Energy Distribution')

    # add number of data points to upper right corner
    plt.text(0.9, 0.9, f'Total Trajectory Data\n{len(first_image)}', ha='right', va='top', transform=plt.gca().transAxes)



    plt.tight_layout()
    plt.savefig(f'Data/figures/fig_save/S3_5_{name}_energy_dif_{fraction}.png')


# plot max force distribution: train
json_pth = [r"Data/Surface/Traj/all_images_element_1_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_2_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_ref_energy.json", ]

plot_force_distribution(json_pth, name='Train', fraction=1)

# plot max force distribution: test AL
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_al_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_al_ref_energy.json", ]

plot_force_distribution(json_pth, name='Test_AL', fraction=1)

# plot max force distribution: test AL
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_id_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_id_ref_energy.json",
            ]

plot_force_distribution(json_pth, name='Test_ID', fraction=1)

# plot max force distribution: test OOD
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_ood_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_ood_ref_energy.json", ]

plot_force_distribution(json_pth, name='Test_OOD', fraction=1)
