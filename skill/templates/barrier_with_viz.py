"""CatDT Template: NEB Barrier Calculation with Visualization

Computes activation barriers via Nudged Elastic Band (two-phase CI-NEB).
Auto-generates: energy diagram, NEB trajectory GIF, transition state structure.

NOTE: Staged-atom NEB for hydrogenation steps may give unreliable barriers.
      H on surface hollow is NOT a true local minimum.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

REACTANT_PATH = "reactant.vasp"     # Initial state structure
PRODUCT_PATH = "product.vasp"       # Final state structure
N_FRAMES = 11                       # NEB images (odd number recommended for CI-NEB)
FMAX = 0.05                         # Force convergence criterion (eV/A)
MAX_STEPS = 400                     # Max optimization steps
MODEL = "uma-s-1p1"                 # UMA model: "uma-s-1p1" or "uma-m-1p1"
OUTPUT_DIR = os.path.join("output", "neb", datetime.now().strftime("%Y%m%d_%H%M%S"))

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from ase.io import read
from core.pathway.fairchem_predictor import FairchemPredictor
from core.pathway.barrier_predictor import BarrierPredictor

# Initialize
fairchem = FairchemPredictor(
    fairchem_root=os.path.join(CATDT_ROOT, "deps/fairchem"),
    model_name=MODEL,
    use_gpu=True,
)
bp = BarrierPredictor(fairchem_predictor=fairchem, work_dir=OUTPUT_DIR)

# Run NEB
reactant = read(REACTANT_PATH)
product = read(PRODUCT_PATH)

result = bp.predict_from_structures(
    reactant=reactant,
    product=product,
    n_frames=N_FRAMES,
    fmax=FMAX,
    max_steps=MAX_STEPS,
    output_dir=OUTPUT_DIR,
)

# Print results
print(f"Reaction: {REACTANT_PATH} → {PRODUCT_PATH}")
print(f"Ea_fwd  = {result.activation_energy_forward:.3f} eV")
print(f"Ea_rev  = {result.activation_energy_reverse:.3f} eV")
print(f"E_rxn   = {result.reaction_energy:.3f} eV")
print(f"TS index: {result.transition_state_index}")
print(f"Converged: {result.converged} (fmax={result.fmax_final:.4f})")

# ── Auto-Visualization ────────────────────────────────────────
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter

# 1. Energy diagram
plotter = EnergyDiagramPlotter()
plotter.plot(
    energies=result.energies,
    labels=[f"F{i}" for i in range(len(result.energies))],
)
energy_png = os.path.join(OUTPUT_DIR, "energy_profile.png")
plotter.save(energy_png)
plotter.export_data(os.path.join(OUTPUT_DIR, "energy_profile.csv"))
print(f"Saved: {energy_png}")

# 2. NEB trajectory GIF
viz = CatalystSurfaceVisualizer(quality="high")
viz.visualize_trajectory(
    result.neb_frames,
    output_file=os.path.join(OUTPUT_DIR, "neb_path.gif"),
    fps=3,
    titles=[f"Frame {i} ({result.energies[i]:.2f} eV)" for i in range(len(result.neb_frames))],
)
print(f"Saved: neb_path.gif")

# 3. Transition state structure
viz.visualize_structure(
    result.transition_state,
    output_file=os.path.join(OUTPUT_DIR, "transition_state.png"),
    title=f"TS (Ea = {result.activation_energy_forward:.2f} eV)",
)
print(f"Saved: transition_state.png")

# 4. Reactant and product
viz.visualize_structure(reactant, os.path.join(OUTPUT_DIR, "reactant.png"), title="Reactant")
viz.visualize_structure(product, os.path.join(OUTPUT_DIR, "product.png"), title="Product")

print(f"\nAll outputs in: {OUTPUT_DIR}")
