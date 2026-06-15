# compare predicted surface energies of model before and after fine-tuning

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({'font.size': 16,})

data = pd.read_csv("results/element_exp2&3/wulff_results.csv")


pred = data["surface_energy_pred"].values * 1000
true = data["surface_energy_true"].values * 1000

dif = pred - true

# plot the scatter plot of dif, show horizontal line at y=0

plt.figure(figsize=(8, 6),dpi=300)
plt.scatter(np.arange(len(dif)), dif, color='gray', alpha=0.3)
plt.axhline(0, color='gray', linestyle='--')
plt.xlabel('Number of surfaces')

plt.ylabel('Difference in Surface Energy (meV/Å²)')
plt.title('Difference in Predicted & True Surface Energy')
plt.tight_layout()

plt.show()


# show details of one crystal
crystal_id = 3024
crystal = data[data["crystal_id"] == crystal_id]
crystal_pred = crystal["surface_energy_pred"].values * 1000
crystal_true = crystal["surface_energy_true"].values * 1000
crystal_dif = crystal_pred - crystal_true

# plot the scatter plot of dif, show horizontal line at y=0

plt.figure(figsize=(8, 6),dpi=300)
plt.scatter(np.arange(len(crystal_dif)), crystal_dif, color='gray', alpha=0.5)

plt.axhline(0, color='gray', linestyle='--')
plt.xlabel('Number of surfaces')

plt.ylabel('Difference in Surface Energy (meV/Å²)')
plt.ylim(-30, 30)
plt.title(f'Difference in Predicted & True Surface Energy \n Crystal ID: {crystal_id}')
plt.tight_layout()

plt.show()