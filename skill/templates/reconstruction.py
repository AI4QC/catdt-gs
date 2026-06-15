"""CatDT Template: Surface Reconstruction (VSSR-MC)

Simulates surface reconstruction under operating conditions using Monte Carlo.
Supports both thermal (gas-solid) and electrochemical (liquid-solid) modes.
Auto-generates MC trajectory GIF and energy history plot.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

SURFACE_PATH = "slab.vasp"           # Surface with adsorbate
ADSORBATE_ELEMENTS = ["Cu"]          # Elements to add/remove in MC moves
TEMPERATURE_K = 500.0                # Temperature in Kelvin
TOTAL_SWEEPS = 100                   # MC sweep count
SWEEP_SIZE = 20                      # Moves per sweep
MODEL_TYPE = "CHGNetNFF"             # "CHGNetNFF" or "UMA"
CANONICAL = False                    # True = fixed adsorbate count
NUM_ADSORBATES = 0                   # Only for canonical mode
OUTPUT_DIR = os.path.join("output", "vssr_mc", datetime.now().strftime("%Y%m%d_%H%M%S"))

# Electrochemical mode (set both to enable Pourbaix)
POTENTIAL_SHE = None  # e.g., 1.23 (V vs SHE)
PH = None             # e.g., 14.0

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from core.reconstruction.vssr_mc_predictor import VSSRMCPredictor

predictor = VSSRMCPredictor(
    surface_sampling_root=os.path.join(CATDT_ROOT, "deps/surface-sampling"),
    model_type=MODEL_TYPE,
    device="cuda",
    potential_she=POTENTIAL_SHE,
    ph=PH,
)

result = predictor.sample(
    surface=SURFACE_PATH,
    adsorbates=ADSORBATE_ELEMENTS,
    temperature=TEMPERATURE_K,
    total_sweeps=TOTAL_SWEEPS,
    sweep_size=SWEEP_SIZE,
    canonical=CANONICAL,
    num_adsorbates=NUM_ADSORBATES,
    output_dir=OUTPUT_DIR,
)

print(f"Best structure energy: {result.lowest_energy_structure.energy:.3f} eV")
print(f"Acceptance rate: {sum(result.acceptance_history) / len(result.acceptance_history):.2%}")

# ── Auto-Visualization ────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

# 1. MC trajectory GIF
viz = CatalystSurfaceVisualizer(quality="high")
trajectory = [s.atoms for s in result.structures if s.atoms is not None]
if trajectory:
    # Sample every N frames for manageable GIF
    step = max(1, len(trajectory) // 30)
    sampled = trajectory[::step]
    viz.visualize_trajectory(
        sampled,
        output_file=os.path.join(OUTPUT_DIR, "mc_trajectory.gif"),
        fps=5,
    )
    print(f"Saved: mc_trajectory.gif")

# 2. Energy history plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.plot(result.energy_history, "b-", linewidth=0.5)
ax1.set_ylabel("Energy (eV)")
ax1.set_title("MC Energy History")
ax2.plot(result.acceptance_history, "r-", linewidth=0.5)
ax2.set_ylabel("Acceptance")
ax2.set_xlabel("MC Step")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "mc_energy_history.png"), dpi=200)
plt.close(fig)
print(f"Saved: mc_energy_history.png")

# 3. Best structure PNG
viz.visualize_structure(
    result.lowest_energy_structure.atoms,
    output_file=os.path.join(OUTPUT_DIR, "best_structure.png"),
    title=f"Lowest E = {result.lowest_energy_structure.energy:.3f} eV"
)

print(f"\nAll outputs in: {OUTPUT_DIR}")
