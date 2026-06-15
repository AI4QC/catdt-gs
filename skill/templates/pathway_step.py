"""CatDT Template: Single-Step Endpoint Design + Validation + NEB

Standalone reproduction of the Agent4 + Agent5 + NEB-Tool slice for ONE
elementary step. Does NOT require an LLM or the full CAMEL pipeline.

Three modes (set MODE):
  - "edit"     : you supply a reactant.vasp and an EDIT plan (atoms_to_add /
                 atoms_to_remove); template builds the product, validates,
                 and runs NEB. This mimics the Agent4 → Agent5 → NEB chain.
  - "endpoints": you already have reactant.vasp and product.vasp; template
                 only runs Agent5-equivalent validation and NEB.
  - "auto"     : template auto-detects: if PRODUCT_PATH exists → "endpoints",
                 else → "edit" (default).

WARNING (critical, see SKILL.md §7):
  Staged-atom NEB (H placed on a surface hollow near an organic adsorbate)
  CANNOT reliably compute hydrogenation barriers. The H site is not a true
  local minimum. Only dissociation / bond-breaking steps give trustworthy Ea.
"""

import os, sys, json
from datetime import datetime
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

REACTANT_PATH = "reactant.vasp"          # initial state (must contain the adsorbate already on the surface)
PRODUCT_PATH = "product.vasp"            # final state; if missing, template tries "edit" mode
MODE = "auto"                            # "edit" | "endpoints" | "auto"

# Edit plan (only used in "edit" mode) — Agent4-style declarative spec.
# Each entry in ATOMS_TO_ADD: {"species": "H", "position": [x, y, z]}.
# ATOMS_TO_REMOVE: list of indices into REACTANT to drop in the product.
ATOMS_TO_ADD = [
    # Example: stage an H near the *CO carbon for *CO + *H → *CHO
    # {"species": "H", "position": [4.52, 3.18, 13.40]},
]
ATOMS_TO_REMOVE = []                     # indices into reactant atoms list

# Adsorbate / surface index hints (required for validation).
# Indices in REACTANT atoms list; surface = everything that is NOT adsorbate.
REACTANT_ADSORBATE_INDICES = []          # e.g. [-2, -1] — last two atoms are *CO
PRODUCT_ADSORBATE_INDICES = None         # auto-derived in "edit" mode

REACTION_TYPE = "hydrogenation"          # "hydrogenation" | "dissociation" | "association" | ...
STEP_NAME = "step"                       # used in filenames

# NEB settings
N_FRAMES = 11
FMAX = 0.05
MAX_STEPS = 400
MODEL = "uma-s-1p1"                      # "uma-s-1p1" or "uma-m-1p1"
USE_TWO_PHASE = True                     # phase-1 NEB → phase-2 CI-NEB

OUTPUT_DIR = os.path.join(
    "output", "pathway_step", datetime.now().strftime("%Y%m%d_%H%M%S")
)

# ── Helpers ────────────────────────────────────────────────────
def _resolve_adsorbate_indices(atoms, raw):
    """Normalise negative indices, default to []."""
    if not raw:
        return []
    n = len(atoms)
    return [(idx + n) % n for idx in raw]


def _build_product_via_edit(reactant, atoms_to_add, atoms_to_remove,
                            reactant_adsorbate_indices):
    """Apply the Agent4-style edit plan to construct the product.

    Element-multiset constraint: NEB requires reactant and product to have
    matching element counts. We don't enforce that here — Agent5 validation
    will catch it.
    """
    from ase import Atoms
    from ase.atom import Atom

    product = reactant.copy()

    # Remove atoms (sort descending so indices stay valid).
    for idx in sorted(set(int(i) for i in atoms_to_remove), reverse=True):
        if 0 <= idx < len(product):
            del product[idx]
        else:
            print(f"[WARN] atoms_to_remove index {idx} out of range (skipped)")

    # Add atoms.
    new_indices = []
    for spec in atoms_to_add:
        symbol = str(spec["species"])
        pos = [float(c) for c in spec["position"]]
        product.append(Atom(symbol=symbol, position=pos))
        new_indices.append(len(product) - 1)

    # Derive product adsorbate indices: original adsorbate indices that survived
    # the removal, plus all newly added atoms.
    removed = set(int(i) for i in atoms_to_remove)
    survived_ads = [i for i in (reactant_adsorbate_indices or []) if i not in removed]
    # Re-index survived adsorbate after deletions.
    reindexed = []
    for orig_idx in survived_ads:
        shift = sum(1 for r in removed if r < orig_idx)
        reindexed.append(orig_idx - shift)
    product_ads = sorted(set(reindexed + new_indices))

    return product, product_ads


# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from ase.io import read, write

# 1. Load reactant.
print(f"[Step 1] Reading reactant from {REACTANT_PATH}")
reactant = read(REACTANT_PATH)
reactant_ads = _resolve_adsorbate_indices(reactant, REACTANT_ADSORBATE_INDICES)
if not reactant_ads:
    raise ValueError(
        "REACTANT_ADSORBATE_INDICES must be set — the validator and NEB tool "
        "both need to know which atoms are adsorbate vs. surface."
    )
reactant.info["adsorbate_indices"] = reactant_ads
reactant.info["surface_indices"] = [i for i in range(len(reactant)) if i not in set(reactant_ads)]

# 2. Resolve product.
mode = MODE
if mode == "auto":
    mode = "endpoints" if Path(PRODUCT_PATH).is_file() else "edit"
print(f"[Step 2] Mode = {mode}")

if mode == "endpoints":
    product = read(PRODUCT_PATH)
    if PRODUCT_ADSORBATE_INDICES is None:
        # Fall back: assume same set of indices as reactant (overlap by symbol)
        # works only when atom ordering is preserved; user should set explicitly.
        product_ads = list(reactant_ads)
        print("[WARN] PRODUCT_ADSORBATE_INDICES not set — falling back to reactant's. "
              "Set it explicitly for correctness.")
    else:
        product_ads = _resolve_adsorbate_indices(product, PRODUCT_ADSORBATE_INDICES)
elif mode == "edit":
    print(f"[Step 2a] Building product via edit plan: "
          f"+{len(ATOMS_TO_ADD)} atoms, -{len(ATOMS_TO_REMOVE)} atoms")
    product, product_ads = _build_product_via_edit(
        reactant, ATOMS_TO_ADD, ATOMS_TO_REMOVE, reactant_ads,
    )
else:
    raise ValueError(f"Unknown MODE={MODE!r}; expected 'edit' | 'endpoints' | 'auto'")

product.info["adsorbate_indices"] = product_ads
product.info["surface_indices"] = [i for i in range(len(product)) if i not in set(product_ads)]

# Persist endpoints for audit.
reactant_path = Path(OUTPUT_DIR) / f"{STEP_NAME}_reactant.vasp"
product_path = Path(OUTPUT_DIR) / f"{STEP_NAME}_product.vasp"
write(reactant_path, reactant)
write(product_path, product)
print(f"  Reactant: {len(reactant)} atoms, adsorbate {reactant_ads}")
print(f"  Product : {len(product)} atoms, adsorbate {product_ads}")

# 3. Agent5-equivalent validation.
print("[Step 3] Validating endpoints (EnhancedNEBValidator = programmatic gate)")
from core.pathway.enhanced_neb_validator import EnhancedNEBValidator

validator = EnhancedNEBValidator(auto_fix=True, fix_attempts=3)
report = validator.validate_endpoint_pair(
    initial=reactant,
    final=product,
    expected_reaction=REACTION_TYPE,
)
print(f"  Status: {report.status}")
for issue in getattr(report, "issues", []):
    sev = getattr(issue, "severity", "?")
    msg = getattr(issue, "message", str(issue))
    print(f"    [{sev}] {msg}")

with open(Path(OUTPUT_DIR) / f"{STEP_NAME}_validation.json", "w") as fh:
    json.dump({
        "status": str(report.status),
        "reaction_type": str(getattr(report, "reaction_type", REACTION_TYPE)),
        "n_issues": len(getattr(report, "issues", [])),
    }, fh, indent=2)

status_str = str(report.status).lower()
if "fail" in status_str:
    print("[ABORT] Validation FAILED — fix endpoints before NEB. "
          "Common fixes: align element multisets, separate atoms < 0.8 A, "
          "lift staged atoms above the surface top z + 1.5 A.")
    sys.exit(1)

# 4. Reaction-type guardrail (Skill §7 known limitation).
if REACTION_TYPE.lower() == "hydrogenation":
    print("[NOTE] Hydrogenation step detected. Staged-atom NEB barriers for "
          "*X + *H → *XH are typically NOT reliable on metal surfaces — H on "
          "the surface hollow is not a true local minimum. Treat the resulting "
          "Ea as an upper-bound estimate only.")

# 5. Run NEB (deterministic Tool, not an Agent).
print(f"[Step 4] Running CI-NEB ({N_FRAMES} frames, fmax={FMAX}, max_steps={MAX_STEPS})")
from core.pathway.fairchem_predictor import FairchemPredictor
from core.pathway.barrier_predictor import BarrierPredictor

fairchem = FairchemPredictor(
    fairchem_root=os.path.join(CATDT_ROOT, "deps/fairchem"),
    model_name=MODEL,
    use_gpu=True,
)
bp = BarrierPredictor(fairchem_predictor=fairchem, work_dir=OUTPUT_DIR, verbose=True)

result = bp.predict_from_structures(
    reactant=reactant,
    product=product,
    n_frames=N_FRAMES,
    fmax=FMAX,
    max_steps=MAX_STEPS,
    relax_endpoints=False,           # endpoints already vetted
    use_two_phase=USE_TWO_PHASE,
    reaction_name=STEP_NAME,
    reactant_adsorbate_indices=reactant_ads,
    product_adsorbate_indices=product_ads,
)

print(f"  Ea_fwd  = {result.activation_energy_forward:.3f} eV")
print(f"  Ea_rev  = {result.activation_energy_reverse:.3f} eV")
print(f"  E_rxn   = {result.reaction_energy:.3f} eV")
print(f"  TS at frame {result.transition_state_index}/{result.n_frames}")
print(f"  Converged: {result.converged}")

with open(Path(OUTPUT_DIR) / f"{STEP_NAME}_neb_result.json", "w") as fh:
    json.dump({
        "step_name": STEP_NAME,
        "reaction_type": REACTION_TYPE,
        "Ea_fwd_eV": float(result.activation_energy_forward),
        "Ea_rev_eV": float(result.activation_energy_reverse),
        "E_rxn_eV": float(result.reaction_energy),
        "ts_index": int(result.transition_state_index),
        "n_frames": int(result.n_frames),
        "converged": bool(result.converged),
        "fmax_final": float(getattr(result, "fmax_final", float("nan"))),
        "model": MODEL,
        "validation_status": str(report.status),
    }, fh, indent=2)

# 6. Auto-visualization (mandatory per SKILL.md §5).
print("[Step 5] Visualizing endpoints, NEB path, energy diagram, and TS")
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter

viz = CatalystSurfaceVisualizer(quality="high")
viz.visualize_structure(reactant, str(Path(OUTPUT_DIR) / f"{STEP_NAME}_reactant.png"),
                        title=f"Reactant ({STEP_NAME})")
viz.visualize_structure(product, str(Path(OUTPUT_DIR) / f"{STEP_NAME}_product.png"),
                        title=f"Product ({STEP_NAME})")
viz.visualize_structure(result.transition_state, str(Path(OUTPUT_DIR) / f"{STEP_NAME}_TS.png"),
                        title=f"TS ({STEP_NAME}, Ea={result.activation_energy_forward:.2f} eV)")
viz.visualize_trajectory(
    result.neb_frames,
    output_file=str(Path(OUTPUT_DIR) / f"{STEP_NAME}_neb_path.gif"),
    fps=3,
    titles=[f"F{i} ({result.energies[i]:.2f} eV)" for i in range(len(result.neb_frames))],
)

plotter = EnergyDiagramPlotter()
plotter.plot(
    energies=result.energies,
    labels=[f"F{i}" for i in range(len(result.energies))],
)
plotter.save(str(Path(OUTPUT_DIR) / f"{STEP_NAME}_energy_profile.png"))
plotter.export_data(str(Path(OUTPUT_DIR) / f"{STEP_NAME}_energy_profile.csv"))

print(f"\nAll outputs in: {OUTPUT_DIR}")
