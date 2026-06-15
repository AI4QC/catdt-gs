# Visualization Module Reference

CatDT provides three visualization subsystems that should be used together to produce publication-quality outputs.

## Auto-Visualization Policy

**MANDATORY**: Every computation producing structures or trajectories MUST generate visualization automatically. Do not return results without visual output.

| Computation | Must Produce | How |
|---|---|---|
| Structure relaxation | Before/after PNG + trajectory GIF | `visualize_trajectory(traj, "relax.gif")` + `visualize_structure(initial/final)` |
| NEB barrier | Energy diagram PNG + NEB path GIF | `EnergyDiagramPlotter.plot()` + `visualize_trajectory(neb_frames)` |
| MD simulation | Trajectory GIF + T/E plot | `visualize_trajectory(md_traj)` + matplotlib T(t), E(t) |
| MC reconstruction | MC trajectory GIF + energy history | `visualize_mc_reconstruction(traj)` + energy/acceptance plot |
| Adsorption sites | Top-N site PNGs | `visualize_adsorption_sites(structures)` |
| Surface prediction | Top-N slab PNGs | `visualize_generated_surfaces(files)` |
| Mechanism search | Best pathway energy diagram | `EnergyDiagramPlotter.plot(pathway_energies)` |

## 1. CatalystSurfaceVisualizer

OVITO-based 3D structure rendering with ray-tracing.

```python
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

viz = CatalystSurfaceVisualizer(
    elevation=30, azimuth=45,
    quality='high',          # low|medium|high|ultra
    renderer='tachyon',      # tachyon|ospray|opengl
    background_color='white',
    color_scheme='jmol',
    show_cell=False,
    auto_expand=True,        # expand small structures to supercell
    expand_threshold=15,     # atoms below this → expand
    supercell=(2, 2, 1),
)

# Single structure → PNG
viz.visualize_structure(atoms, output_file="structure.png", title="Cu(111)")

# Trajectory → animated GIF
viz.visualize_trajectory(
    trajectory=[atoms1, atoms2, ...],  # List[Atoms]
    output_file="trajectory.gif",
    fps=5,
    titles=["Step 1", "Step 2", ...],
)
```

**Quality settings:**
| Quality | Resolution | Samples/pixel |
|---------|-----------|---------------|
| low | 800x800 | 4 |
| medium | 1200x1200 | 8 |
| high | 2048x2048 | 16 |
| ultra | 4096x4096 | 32 |

## 2. EnergyDiagramPlotter

Publication-quality free energy profile diagrams.

```python
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

plotter = EnergyDiagramPlotter(style=DiagramStyle(
    color_scheme='material',     # material|nature|science|elegant
    figsize=(12, 6),
    dpi=300,
    interp_method='spline',      # spline|quadratic
    show_energies=True,
    use_shadow=True,
    ylabel='Free Energy (eV)',
))

# From lists
plotter.plot(
    energies=[0.0, 0.85, -0.3, 1.2, -1.5],
    labels=['*CO', 'TS1', '*CHO', 'TS2', '*CH2O'],
    is_ts_list=[False, True, False, True, False],
    reference_index=0,
)
plotter.save("energy_diagram.png")

# From dicts
plotter.plot([
    {'label': '*CO', 'energy': 0.0, 'is_ts': False},
    {'label': 'TS1', 'energy': 0.85, 'is_ts': True},
    {'label': '*CHO', 'energy': -0.3, 'is_ts': False},
])

# Export data
plotter.export_data("energies.csv")
```

**Convenience function:**
```python
from core.viz.energy_diagram_plotter import quick_plot

quick_plot(
    energies=[0.0, 0.85, -0.3],
    labels=['IS', 'TS', 'FS'],
    is_ts_list=[False, True, False],
    output="diagram.png",
    color_scheme="nature",
)
```

Chemical formula auto-formatting: `CH3OH` → `CH₃OH`, `*CO` renders with adsorption marker.

## 3. CatalysisVisualizationManager

High-level orchestrator that generates a complete 7-type visualization suite.

```python
from core.viz.visualization_manager import CatalysisVisualizationManager

mgr = CatalysisVisualizationManager(output_dir="output/viz")

# Individual visualization methods:

# 1. Surfaces
mgr.visualize_generated_surfaces(
    surface_files={"111": "slab_111.vasp", "100": "slab_100.vasp"},
    surface_energies={"111": 0.082, "100": 0.095},
) # → Dict[str, str] (miller_index → PNG path)

# 2. Adsorption sites
mgr.visualize_adsorption_sites(
    adsorbate_structures={"site_0": "ads_0.vasp", "site_1": "ads_1.vasp"},
    adsorption_energies={"site_0": -1.5, "site_1": -1.3},
    surface_name="Cu111",
)

# 3. MC reconstruction trajectory
mgr.visualize_mc_reconstruction(
    trajectory_file="mc_trajectory.traj",  # or List[Atoms]
    surface_name="Cu111",
    fps=5,
    max_frames=50,
) # → GIF path

# 4. Reaction pathway animation
mgr.visualize_reaction_pathway(
    pathway_structures={"*CO": "co.vasp", "*CHO": "cho.vasp"},
    sequence=["*CO", "*CHO", "*CH2O"],
    surface_name="Cu111",
    fps=2,
    hold_frames=3,  # hold each intermediate for N frames
)

# 5. Energy diagram
mgr.visualize_energy_diagram(
    energies=[0.0, 0.85, -0.3],
    labels=["*CO", "TS1", "*CHO"],
    is_ts_list=[False, True, False],
    barriers={"step1": 0.85},
    title="CO Hydrogenation on Cu(111)",
)

# 6. NEB trajectory
mgr.visualize_neb_trajectory(
    trajectory=neb_frames,         # List[Atoms]
    reaction_name="CO_to_CHO",
    fps=2,
    is_optimized=True,
) # → GIF path

# Also from file:
mgr.visualize_neb_from_file(
    traj_file="neb.traj",
    n_images=11,
)

# 7. Complete workflow (all 7 types at once)
outputs = mgr.visualize_complete_workflow(workflow_result, "Cu111")

# 8. HTML summary page
mgr.generate_summary_html(outputs, "summary.html")
```

**Output directory structure:**
```
output_dir/visualizations/
├── 01_surfaces/          # Surface slab PNGs
├── 02_adsorption/        # Adsorption site PNGs
├── 03_mc_reconstruction/ # MC trajectory GIF
├── 04_reaction_pathway/  # Pathway animation GIF
├── 05_energy_diagram/    # Energy diagram PNG + CSV
├── 06_neb_trajectories/  # NEB path GIFs
├── 07_kmc_dynamics/      # KMC surface dynamics GIF
└── visualization_summary.html
```

## Trajectory File Handling

```python
from ase.io import read, Trajectory

# Read trajectory
frames = read("trajectory.traj", index=":")  # all frames

# Write trajectory
traj = Trajectory("output.traj", "w")
for atoms in frame_list:
    traj.write(atoms)
traj.close()

# Supported formats: .traj (ASE), .xyz, XDATCAR (VASP), .pdb
```

## Visualization Recipe: NEB with Full Output

```python
import os, sys
CATDT_ROOT = os.environ.get("CATDT_ROOT", ".")
sys.path.insert(0, CATDT_ROOT)

from ase.io import read
from core.pathway.barrier_predictor import BarrierPredictor
from core.pathway.fairchem_predictor import FairchemPredictor
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter

# Setup
output_dir = "output/neb_demo"
os.makedirs(output_dir, exist_ok=True)

# NEB calculation
fairchem = FairchemPredictor(fairchem_root=os.path.join(CATDT_ROOT, "deps/fairchem"))
bp = BarrierPredictor(fairchem_predictor=fairchem)
result = bp.predict_from_structures(
    reactant=read("reactant.vasp"),
    product=read("product.vasp"),
    n_frames=11,
)

# Auto-visualization: energy diagram
plotter = EnergyDiagramPlotter()
plotter.plot(
    energies=result.energies,
    labels=[f"Frame {i}" for i in range(len(result.energies))],
)
plotter.save(os.path.join(output_dir, "energy_profile.png"))

# Auto-visualization: NEB trajectory GIF
viz = CatalystSurfaceVisualizer(quality='high')
viz.visualize_trajectory(
    result.neb_frames,
    output_file=os.path.join(output_dir, "neb_path.gif"),
    fps=3,
)

# Auto-visualization: transition state structure
viz.visualize_structure(
    result.transition_state,
    output_file=os.path.join(output_dir, "transition_state.png"),
    title=f"TS (Ea = {result.activation_energy_forward:.2f} eV)",
)

print(f"Ea_fwd = {result.activation_energy_forward:.3f} eV")
print(f"Ea_rev = {result.activation_energy_reverse:.3f} eV")
print(f"E_rxn  = {result.reaction_energy:.3f} eV")
print(f"Converged: {result.converged}")
```

## Key Files

| File | Class | Lines |
|------|-------|-------|
| `core/viz/catalyst_surface_visualizer.py` | `CatalystSurfaceVisualizer` | 714 |
| `core/viz/energy_diagram_plotter.py` | `EnergyDiagramPlotter`, `DiagramStyle` | 871 |
| `core/viz/visualization_manager.py` | `CatalysisVisualizationManager` | 922 |
| `camel_agents/visualized_digital_twin.py` | `VisualizedGasSolidDigitalTwin` | 366 |
