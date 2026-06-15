"""CatDT Template: Surface Prediction (SurFF)

Predicts exposed crystal facets and generates slab structures from bulk.
Auto-generates visualization of all predicted surfaces.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

BULK_STRUCTURE = "POSCAR"          # Path to bulk structure (VASP, CIF, etc.)
TOP_N = 3                          # Number of top facets to predict
OUTPUT_DIR = os.path.join("output", "surff", datetime.now().strftime("%Y%m%d_%H%M%S"))

# If no bulk file, generate from element:
ELEMENT = None  # e.g., "Cu", "Pt", "Ni" — set to generate bulk automatically
CRYSTAL = "fcc"  # fcc, bcc, hcp
LATTICE_CONST = None  # None = use ASE default

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Generate bulk if needed
if ELEMENT and not os.path.exists(BULK_STRUCTURE):
    from ase.build import bulk
    from ase.io import write
    atoms = bulk(ELEMENT, CRYSTAL, a=LATTICE_CONST)
    BULK_STRUCTURE = os.path.join(OUTPUT_DIR, f"bulk_{ELEMENT}.vasp")
    write(BULK_STRUCTURE, atoms, format="vasp")
    print(f"Generated bulk: {BULK_STRUCTURE}")

# Surface prediction
from core.surface.surff_predictor import SurFFPredictor

predictor = SurFFPredictor(
    surff_root=os.path.join(CATDT_ROOT, "deps/SurFF"),
    use_gpu=True,
)
result = predictor.predict(structure=BULK_STRUCTURE, top_n=TOP_N, output_dir=OUTPUT_DIR)

# Print results
for i, surface in enumerate(result.surfaces[:TOP_N]):
    print(f"Surface {i+1}: {surface.miller_index}")
    print(f"  Energy: {surface.surface_energy:.4f} eV/A^2")
    print(f"  Area fraction: {surface.area_fraction:.2%}")

# ── Auto-Visualization ────────────────────────────────────────
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from ase.io import read

viz = CatalystSurfaceVisualizer(quality="high", renderer="tachyon")

for i, surface in enumerate(result.surfaces[:TOP_N]):
    slab_file = os.path.join(OUTPUT_DIR, f"slab_{i}.vasp")
    if os.path.exists(slab_file):
        atoms = read(slab_file)
        png_path = os.path.join(OUTPUT_DIR, f"surface_{surface.miller_index}.png")
        viz.visualize_structure(
            atoms, output_file=png_path,
            title=f"{surface.miller_index} (E={surface.surface_energy:.4f})"
        )
        print(f"  Saved: {png_path}")

print(f"\nAll outputs in: {OUTPUT_DIR}")
