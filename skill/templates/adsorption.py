"""CatDT Template: Adsorption Site Prediction (AdsorbDiff)

Predicts optimal adsorption sites for an adsorbate on a surface.
Auto-generates visualization of top adsorption configurations.

IMPORTANT: Adsorbate names MUST use * prefix: "*CO", "*OH", "*H", NOT "CO".
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

SURFACE_PATH = "slab.vasp"        # Path to surface structure
ADSORBATE = "*CO"                  # MUST use * prefix: *CO, *OH, *H, *O, *OOH
NUM_SITES = 5                      # Number of sites to sample
OUTPUT_DIR = os.path.join("output", "adsorption", datetime.now().strftime("%Y%m%d_%H%M%S"))

# If you have a Miller index instead of a file:
BUILD_SLAB = False
ELEMENT = "Cu"
MILLER = (1, 1, 1)
LAYERS = 4
VACUUM = 15.0

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Build slab if needed
if BUILD_SLAB:
    from ase.build import fcc111, fcc100, bcc110
    from ase.io import write
    from ase.constraints import FixAtoms
    import numpy as np

    slab = fcc111(ELEMENT, size=(3, 3, LAYERS), vacuum=VACUUM)
    # Fix bottom 2 layers
    z = slab.positions[:, 2]
    mask = z < z.min() + (z.max() - z.min()) * 0.5
    slab.set_constraint(FixAtoms(mask=mask))
    SURFACE_PATH = os.path.join(OUTPUT_DIR, f"{ELEMENT}{''.join(map(str, MILLER))}.vasp")
    write(SURFACE_PATH, slab, format="vasp")
    print(f"Built slab: {SURFACE_PATH}")

# Ensure * prefix
if not ADSORBATE.startswith("*"):
    ADSORBATE = f"*{ADSORBATE}"
    print(f"Auto-added * prefix: {ADSORBATE}")

# Adsorption prediction
from core.reconstruction.adsorbdiff_predictor import AdsorbDiffPredictor

predictor = AdsorbDiffPredictor(
    adsorbdiff_root=os.path.join(CATDT_ROOT, "deps/AdsorbDiff"),
    use_gpu=True,
)
result = predictor.predict(
    surface=SURFACE_PATH,
    adsorbate=ADSORBATE,
    num_sites=NUM_SITES,
    output_dir=OUTPUT_DIR,
)

# Print results
print(f"\n{ADSORBATE} adsorption on {SURFACE_PATH}:")
for i, r in enumerate(result.results):
    print(f"  Site {i}: E = {r.energy:.3f} eV")

# ── Auto-Visualization ────────────────────────────────────────
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

viz = CatalystSurfaceVisualizer(quality="high")

for i, r in enumerate(result.results[:3]):  # Top 3
    png_path = os.path.join(OUTPUT_DIR, f"adsorption_site_{i}.png")
    viz.visualize_structure(
        r.final_structure, output_file=png_path,
        title=f"{ADSORBATE} site {i} (E={r.energy:.3f} eV)"
    )
    print(f"Saved: {png_path}")

print(f"\nAll outputs in: {OUTPUT_DIR}")
