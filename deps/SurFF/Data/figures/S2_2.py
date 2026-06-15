# crystal cells structure before and after optimization

import os
from pymatgen.io.vasp import Poscar
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# Directories containing the POSCAR files before and after optimization
dir_before = r'Data/Crystal/Raw_crystal/ID'
dir_after = r'Data/Crystal/Opt_crystal/ID'
# Lists to store the maximum lattice length change and angle change
length_changes = []
angle_changes = []

# Loop through the files in the directory before optimization
for filename in tqdm(os.listdir(dir_after)):
    path_before = os.path.join(dir_before, filename)
    path_after = os.path.join(dir_after, filename)

    # Load the POSCAR files
    structure_before = Poscar.from_file(path_before).structure
    structure_after = Poscar.from_file(path_after).structure

    # Calculate maximum lattice length change
    length_before = np.array(structure_before.lattice.abc)
    length_after = np.array(structure_after.lattice.abc)
    max_length_change = np.max(np.abs(length_after - length_before))
    length_changes.append(max_length_change)

    # Calculate lattice angle change
    angle_before = np.array(structure_before.lattice.angles)
    angle_after = np.array(structure_after.lattice.angles)
    max_angle_change = np.max(np.abs(angle_after - angle_before))
    angle_changes.append(max_angle_change)

# Plotting
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

plt.figure(figsize=(10, 5))
font_size = 12

# Histogram of Maximum Lattice Length Changes
plt.subplot(1, 2, 1)

# Bins with explicit underflow and overflow
bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]  # Adding np.inf for overflow
bins_label = ['<0.01', '0.01-0.05', '0.05-0.1', '0.1-0.2', '>0.2']
# convert into bars data
bar_data = np.histogram(length_changes, bins=bins)
# plot as bar
plt.bar(bins_label, bar_data[0], color='Gray')
plt.title('Histogram of Maximum Lattice Length Change', fontdict={'fontsize': font_size})
plt.xlabel('Length Change (Å)', fontdict={'fontsize': font_size})
plt.ylabel('Frequency', fontdict={'fontsize': font_size})
plt.ylim(0, 2600)

# Histogram of Maximum Lattice Angle Changes
plt.subplot(1, 2, 2)

# Bins with explicit underflow and overflow
bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]  # Adding np.inf for overflow
bins_label = ['<0.01', '0.01-0.05', '0.05-0.1', '0.1-0.2', '>0.2']
# convert into bars data
bar_data = np.histogram(angle_changes, bins=bins)
# plot as bar
plt.bar(bins_label, bar_data[0], color='Gray')
plt.title('Histogram of Maximum Lattice Angle Change', fontdict={'fontsize': font_size})
plt.xlabel('Angle Change (Degrees)', fontdict={'fontsize': font_size})
plt.ylabel('Frequency', fontdict={'fontsize': font_size})
plt.ylim(0, 2600)

plt.tight_layout()
# save the figure
plt.savefig('Data/figures/fig_save/S2_2.png')
