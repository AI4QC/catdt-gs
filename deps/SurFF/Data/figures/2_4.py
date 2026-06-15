from ase.data import covalent_radii
from ase.io import read
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# Read the OUTCAR file
outcar_pth =r"D:\Temp_Data\Intermatallic_Crystal_Nanoparticle\surface_model\surface_trajactory\sample_local_0421_element_3_test\element_3_test\slab\1077_15_set\vasprun.xml"

atoms = read(outcar_pth, index="-1")  # '-1' means the last frame, adjust accordingly

# Extract positions, forces, and atomic numbers
positions = atoms.get_positions()
forces = atoms.get_forces()
atomic_numbers = atoms.get_atomic_numbers()

# Get atomic radii based on atomic numbers
radii = [covalent_radii[number] for number in atomic_numbers] * 5

# Create a 3D scatter plot for atoms
fig = go.Figure()

# Add atoms to the plot with sizes based on their radii
fig.add_trace(go.Scatter3d(
    x=positions[:, 0],
    y=positions[:, 1],
    z=positions[:, 2],
    mode='markers',
    marker=dict(
        size=[radius * 50 for radius in radii],  # Scale radii for better visualization
        color=atomic_numbers,
        colorscale='Viridis',
        colorbar=dict(title='Atomic Number')
    )
))

# Add lines to visualize forces as line sticks
for pos, force in zip(positions, forces):
    fig.add_trace(go.Scatter3d(
        x=[pos[0], pos[0] + force[0]],
        y=[pos[1], pos[1] + force[1]],
        z=[pos[2], pos[2] + force[2]],
        mode='lines',
        line=dict(color='gray', width=8)
    ))



# Add coordinate axes
axis_length = 3  # Adjust the length of the axes as needed
fig.add_trace(go.Scatter3d(
    x=[0, 2.5],
    y=[0, 0],
    z=[0, 0],
    mode='lines',
    line=dict(color='blue', width=2),
    name='X axis'
))

fig.add_trace(go.Scatter3d(
    x=[0, 0],
    y=[0, 25],
    z=[0, 0],
    mode='lines',
    line=dict(color='green', width=2),
    name='Y axis'
))

fig.add_trace(go.Scatter3d(
    x=[0, 0],
    y=[0, 0],
    z=[0, 25],
    mode='lines',
    line=dict(color='purple', width=2),
    name='Z axis'
))

## Set the layout, disable the grid, hide ticks and axis titles, and set background color to blank
fig.update_layout(
    scene=dict(
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=False,
            title='',
            backgroundcolor="rgba(0, 0, 0, 0)"
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=False,
            title='',
            backgroundcolor="rgba(0, 0, 0, 0)"
        ),
        zaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=False,
            title='',
            backgroundcolor="rgba(0, 0, 0, 0)"
        ),
        bgcolor="rgba(0, 0, 0, 0)"
    )
)
# Show the plot
fig.show()