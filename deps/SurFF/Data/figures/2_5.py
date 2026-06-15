from ase.data import covalent_radii
from ase.io import read
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# Read the OUTCAR file
outcar_pth =r"D:\Temp_Data\Intermatallic_Crystal_Nanoparticle\surface_model\surface_trajactory\sample_local_0421_element_3_test\element_3_test\slab\1077_15_set\vasprun.xml"

atoms = read(outcar_pth, index=":")  # '-1' means the last frame, adjust accordingly

energy = np.array([atoms[i].get_potential_energy() for i in range(len(atoms))])
forces = [atoms[i].get_forces() for i in range(len(atoms))]
max_forces = [np.max(np.linalg.norm(forces[i], axis=1)) for i in range(len(atoms))]

energy = energy - energy[0]  # set the initial energy to zero

plt.rcParams.update({'font.size': 18,})

# plot the energy
plt.figure(figsize=(15, 2))
plt.plot(energy, linewidth=5)
plt.show()

# plot the max forces
plt.figure(figsize=(15, 2), linewidth=5)
plt.plot(max_forces, linewidth=5)
plt.show()


