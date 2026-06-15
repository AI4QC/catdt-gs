import itertools

import pandas as pd
from collections import defaultdict


# Function to parse formula and return elements and their counts
def parse_formula(formula):
    import re
    pattern = r'([A-Z][a-z]*)(\d*)'
    pairs = re.findall(pattern, formula)
    return {element: int(count) if count else 1 for element, count in pairs}


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# Example DataFrame
csv_pth = [r"Data/dataset_generation/Model/save_dir/element_2_0/sampled_slab_data.csv",
             #  r"Data/dataset_generation/Model/save_dir/element_2_1/sampled_slab_data.csv",
             #    r"Data/dataset_generation/Model/save_dir/element_2_2/sampled_slab_data.csv",
             #    r"Data/dataset_generation/Model/save_dir/element_2_3/sampled_slab_data.csv",
             # r"Data/dataset_generation/Model/save_dir/element_2_4/sampled_slab_data.csv",
           r"Data/dataset_generation/Model/save_dir/element_3_0/sampled_slab_data.csv",
              # r"Data/dataset_generation/Model/save_dir/element_3_1/sampled_slab_data.csv",
              #   r"Data/dataset_generation/Model/save_dir/element_3_2/sampled_slab_data.csv",
              #   r"Data/dataset_generation/Model/save_dir/element_3_3/sampled_slab_data.csv",
              #   r"Data/dataset_generation/Model/save_dir/element_3_4/sampled_slab_data.csv",
              #   r"Data/dataset_generation/Model/save_dir/element_3_5/sampled_slab_data.csv",
              #   r"Data/dataset_generation/Model/save_dir/element_3_6/sampled_slab_data.csv",
                # r"Data/dataset_generation/Model/save_dir/element_3_7/sampled_slab_data.csv",
           ]

# read all csv and concat
df_lst = []
for pth in csv_pth:
    df_lst.append(pd.read_csv(pth))

df = pd.concat(df_lst, axis=0, ignore_index=True)

# Count pairs
pair_counts = defaultdict(int)
for formula in df['formula']:
    elements = parse_formula(formula)
    # Generate all unique pairs
    for element1, element2 in itertools.combinations(elements.keys(), 2):
        pair = tuple(sorted([element1, element2]))
        pair_counts[pair] += 1

# Convert to DataFrame for easier manipulation
pair_df = pd.DataFrame(list(pair_counts.items()), columns=['Pair', 'Count'])

# Get unique elements
elements = ['Ag', 'Al', 'As', 'Au', 'Ba', 'Be', 'Bi', 'Ca', 'Cd', 'Co', 'Cr', 'Cu', 'Fe', 'Ga', 'Ge', 'Hf', 'In', 'Ir',
            'K', 'Li', 'Mg', 'Mn', 'Mo', 'Na', 'Nb', 'Ni', 'Os', 'Pb', 'Pd', 'Pt', 'Rb', 'Re', 'Rh', 'Ru', 'Sb', 'Sc',
            'Si', 'Sn', 'Sr', 'Ta', 'Te', 'Ti', 'Tl', 'V', 'W', 'Y', 'Zn', 'Zr', 'Cs']

# Initialize a matrix for the plot
matrix = np.zeros((len(elements), len(elements)))

# Fill the matrix with counts
for (element1, element2), count in pair_counts.items():
    i, j = elements.index(element1), elements.index(element2)
    matrix[i, j] = matrix[j, i] = count

# Define a colormap, blues
cmap = plt.cm.Blues

min_count = 1
max_count = 20

# Normalize color scale based on the count range
norm = mcolors.Normalize(vmin=min_count, vmax=max_count)

# size font size
plt.rcParams.update({'font.size': 12})
# Plot
plt.figure(figsize=(10, 10), dpi=600)
ax = plt.gca()

# Set the aspect of the plot to equal, so each cell will be square-shaped
ax.set_aspect('equal')

for i in range(len(elements)):
    for j in range(len(elements)):
        count = matrix[i, j]
        if count > 0:
            # Calculate circle size (ensure it does not exceed grid size)
            # Adjust 'min_size' and 'max_size' as needed
            min_size, max_size = 0.1, 1
            size = np.interp(count, [min_count, max_count], [min_size, max_size])
            # Draw the circle
            circle = plt.Circle((j, i), size / 2, color=cmap(norm(count)), alpha=0.7)
            ax.add_artist(circle)

# Set the limits of the plot
ax.set_xlim(-0.5, len(elements) - 0.5)
ax.set_ylim(-0.5, len(elements) - 0.5)

# Add horizontal and vertical dashed lines every 5 elements
for index in range(len(elements)):
    if index % 5 == 0:
        ax.axhline(y=index - 0.5, color='grey', linestyle='--', linewidth=0.7)
        ax.axvline(x=index - 0.5, color='grey', linestyle='--', linewidth=0.7)

# Set ticks
ax.set_xticks(np.arange(len(elements)))
ax.set_xticklabels(labels=elements, rotation=90)  # Horizontal labels
ax.set_yticks(np.arange(len(elements)))
ax.set_yticklabels(labels=elements)

# Title and labels with increased font size
plt.title('Data Distribution: Initial Round', fontsize=32)  # Adjust fontsize as needed
plt.xlabel('Element 1', fontsize=32)  # Adjust fontsize as needed
plt.ylabel('Element 2', fontsize=32)  # Adjust fontsize as needed
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', fraction=0.046,
                    pad=0.04)

# Modify a specific tick label
tick = np.linspace(0, max_count, 5)
cbar.set_ticks(tick)
cbar.set_ticklabels([int(i) for i in tick[:-1]] + [f'{int(tick[-1])}+'])
cbar.ax.tick_params(labelsize=24)  # Adjust font size as needed

plt.tight_layout()
# plt.savefig(r"Data/figures/fig_save/2_2.png")
plt.show()