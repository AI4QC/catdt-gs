"""CatDT Template: Standalone Visualization

Render structures, trajectories, and energy diagrams from existing data.
Useful for post-processing results from other computations.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

OUTPUT_DIR = os.path.join("output", "viz", datetime.now().strftime("%Y%m%d_%H%M%S"))
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── 1. Render a single structure ──────────────────────────────
from ase.io import read
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

viz = CatalystSurfaceVisualizer(
    quality="high",           # low | medium | high | ultra
    renderer="tachyon",       # tachyon | ospray | opengl
    elevation=30,
    azimuth=45,
    auto_expand=True,         # expand small structures to supercell
    supercell=(2, 2, 1),
)

# Single structure
# atoms = read("structure.vasp")
# viz.visualize_structure(atoms, os.path.join(OUTPUT_DIR, "structure.png"), title="My Structure")

# ── 2. Render a trajectory as GIF ─────────────────────────────
# trajectory = read("trajectory.traj", index=":")
# viz.visualize_trajectory(
#     trajectory,
#     output_file=os.path.join(OUTPUT_DIR, "trajectory.gif"),
#     fps=5,
#     titles=[f"Frame {i}" for i in range(len(trajectory))],
# )

# ── 3. Energy diagram ─────────────────────────────────────────
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

plotter = EnergyDiagramPlotter(style=DiagramStyle(
    color_scheme="material",   # material | nature | science | elegant
    figsize=(12, 6),
    dpi=300,
    interp_method="spline",
    show_energies=True,
    use_shadow=True,
))

# Example: reaction pathway energy profile
energies = [0.0, 0.85, -0.30, 1.20, -0.50, 0.45, -1.80]
labels = ["*CO", "TS1", "*CHO", "TS2", "*CH2O", "TS3", "CH3OH(g)"]
is_ts = [False, True, False, True, False, True, False]

plotter.plot(energies=energies, labels=labels, is_ts_list=is_ts)
plotter.save(os.path.join(OUTPUT_DIR, "energy_diagram.png"))
plotter.export_data(os.path.join(OUTPUT_DIR, "energy_diagram.csv"))
print(f"Saved: energy_diagram.png")

# ── 4. Quick one-liner energy diagram ─────────────────────────
from core.viz.energy_diagram_plotter import quick_plot

quick_plot(
    energies=[0.0, 0.85, -0.30],
    labels=["IS", "TS", "FS"],
    is_ts_list=[False, True, False],
    output=os.path.join(OUTPUT_DIR, "quick_diagram.png"),
    color_scheme="nature",
)

# ── 5. Complete workflow visualization ─────────────────────────
# from core.viz.visualization_manager import CatalysisVisualizationManager
# mgr = CatalysisVisualizationManager(output_dir=OUTPUT_DIR)
#
# # Surfaces
# mgr.visualize_generated_surfaces(
#     surface_files={"111": "slab_111.vasp", "100": "slab_100.vasp"},
#     surface_energies={"111": 0.082, "100": 0.095},
# )
#
# # MC trajectory
# mgr.visualize_mc_reconstruction("mc.traj", "Cu111", fps=5, max_frames=50)
#
# # NEB trajectory
# mgr.visualize_neb_from_file("neb.traj", reaction_name="CO_to_CHO", n_images=11)
#
# # HTML summary
# mgr.generate_summary_html(viz_outputs, "summary.html")

print(f"\nAll outputs in: {OUTPUT_DIR}")
