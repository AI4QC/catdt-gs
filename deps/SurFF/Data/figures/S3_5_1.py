# visualize graph data max force in json file

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
    max_forces = np.array([np.max(image['forces']) for image in tqdm(expanded_images)])
    # remove 3-sigma outliers
    sigma = np.std(max_forces)
    mean = np.mean(max_forces)
    max_forces = max_forces[(max_forces > mean - 3 * sigma) & (max_forces < mean + 3 * sigma)]
    weights = np.ones_like(max_forces) / len(max_forces)

    plt.hist(max_forces, bins=100, weights=weights, color='gray')
    if fraction < 1.0:
        plt.xlabel(f'Max Force (eV/$\AA$)\nFirst {int(fraction*100)}% of Relaxation Trajectory')
    else:
        plt.xlabel('Max Force (eV/$\AA$)')
    plt.ylabel('Data Fraction')
    plt.xlim(0, 1.5)
    plt.ylim(0, 0.05)
    plt.title(f'{name}: Max Force Distribution')

    # # add accumulated distribution to right y-axis
    # ax2 = plt.gca().twinx()
    # ax2.hist(max_forces, bins=100, histtype='step', cumulative=True, density=True, color='black')
    # # set second y-axis to show 0-1
    # ax2.set_ylim(0, 1)

    # add number of data points to upper right corner
    plt.text(0.9, 0.9, f'Total SP Data\n{len(expanded_images)}', ha='right', va='top', transform=plt.gca().transAxes)



    plt.tight_layout()
    plt.savefig(f'Data/figures/fig_save/S3_5_{name}_max_force_{fraction}.png')


# plot max force distribution: train
json_pth = [r"Data/Surface/Traj/all_images_element_1_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_2_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_ref_energy.json", ]

plot_force_distribution(json_pth, name='Train', fraction=1)

# plot max force distribution: test AL
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_al_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_al_ref_energy.json", ]

plot_force_distribution(json_pth, name='Test_AL', fraction=1)

# plot max force distribution: test OOD
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_ood_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_ood_ref_energy.json", ]

plot_force_distribution(json_pth, name='Test_OOD', fraction=1)

# plot max force distribution: test ID
json_pth = [r"Data/Surface/Traj/all_images_element_2_test_id_ref_energy.json",
            r"Data/Surface/Traj/all_images_element_3_test_id_ref_energy.json",
            ]

plot_force_distribution(json_pth, name='Test_ID', fraction=1)