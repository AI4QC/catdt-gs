# PLOT statistics of the predicted surface energy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

prediction_csv = r"D:\OneDrive\OneDrive - National University of Singapore\Pyprojects\ocp\results\all\wulff_results.csv"
prediction = pd.read_csv(prediction_csv)

save_dir = r"Data/figures/fig_save"

# histogram of the predicted surface energy
predicted_surface_energy = prediction['surface_energy_pred'].values
plt.hist(predicted_surface_energy, bins=50, color='gray')
plt.xlabel('Predicted surface energy (eV/Å²)')
plt.ylabel('Quantity')
plt.title('Histogram of the predicted surface energy')
plt.savefig(f'{save_dir}/S5_2_01.png')

# statistics of the predicted surface area
predictied_surface_area = prediction['area_pred'].values
h = plt.hist(predictied_surface_area, bins=50, color='gray')
plt.xlabel('Predicted surface area fraction')
plt.ylabel('Quantity')
# log scale y
plt.yscale('log')
plt.ylim(1e2, 1.5e5)

plt.title('Histogram of the predicted surface area fraction')
plt.savefig(f'{save_dir}/S5_2_02.png')

