"""CatDT Template: UniMech Mechanism Search

Discovers competing reaction pathways from initial to target state.
Three modes: agent_guided (fast), systematic (thorough), both (comprehensive).
Auto-generates best pathway energy diagram.

See references/unimech.md for full API documentation.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

INITIAL_STATE = "*CO"              # Starting adsorbate
TARGET_STATE = "CH4(g)"            # Target product
SURFACE_ID = "Cu111"               # Surface identifier
SURFACE_PATH = "slab_with_co.vasp" # Surface structure with initial adsorbate
CLEAN_SLAB_PATH = "clean_slab.vasp"# Clean surface (no adsorbate)
TEMPERATURE_K = 523.0              # Reaction temperature
PRESSURE_PA = 1e6                  # Total pressure

# Search parameters
EXPLORATION_MODE = "agent_guided"  # agent_guided | systematic | both_sequential | both_parallel
MAX_DEPTH = 8                      # Max pathway length
BEAM_WIDTH = 8                     # Frontier size
MAX_EVALUATIONS = 40               # Max UMA energy calls
DELTA_KEEP = 0.10                  # eV: siblings within this gap → keep both
DELTA_PRUNE = 0.30                 # eV: siblings beyond this gap → prune higher
ENERGY_BACKEND = "thermal_uma"     # thermal_uma | electro_che

OUTPUT_DIR = os.path.join("output", "mechanism", datetime.now().strftime("%Y%m%d_%H%M%S"))

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from camel_agents.tools import CatDTTools

tools = CatDTTools(output_base_dir=OUTPUT_DIR)

# Step 1: Initialize mechanism context
print("Step 1: Initializing mechanism context...")
ctx = tools.initialize_mechanism_context_tool(
    initial_state=INITIAL_STATE,
    target_state=TARGET_STATE,
    surface_id=SURFACE_ID,
    surface_structure_path=SURFACE_PATH,
    clean_slab_path=CLEAN_SLAB_PATH,
    environment={"T": TEMPERATURE_K, "P": PRESSURE_PA},
)
print(f"  Context: {ctx.get('summary', 'initialized')}")

# Step 2: Generate candidate steps
print("Step 2: Generating candidate elementary steps...")
candidates = tools.generate_candidate_steps_tool(
    use_generic_ops=True,
    electro=(ENERGY_BACKEND == "electro_che"),
)
print(f"  Generated {candidates.get('total_candidates', '?')} candidates")

# Step 3: Run mechanism search
print(f"Step 3: Running {EXPLORATION_MODE} search...")
result = tools.run_mechanism_search_tool(
    max_depth=MAX_DEPTH,
    beam_width=BEAM_WIDTH,
    delta_keep=DELTA_KEEP,
    delta_prune=DELTA_PRUNE,
    max_state_evaluations=MAX_EVALUATIONS,
    energy_backend=ENERGY_BACKEND,
    exploration_mode=EXPLORATION_MODE,
)

# Print results
retained = result.get("retained_pathways", [])
pruned = result.get("pruned_pathways", [])
stats = result.get("search_stats", {})
print(f"\nSearch complete:")
print(f"  Retained pathways: {len(retained)}")
print(f"  Pruned pathways: {len(pruned)}")
print(f"  Depth reached: {stats.get('depth_reached', '?')}")
print(f"  Evaluations used: {stats.get('evaluations_used', '?')}")

for i, p in enumerate(retained):
    pathway_id = p if isinstance(p, str) else p.get("pathway_id", f"path_{i}")
    print(f"\n  Pathway {pathway_id}:")
    if isinstance(p, dict):
        states = p.get("states", [])
        sequence = [s.get("species_label", "?") for s in states] if states else []
        print(f"    Sequence: {' → '.join(sequence)}")
        print(f"    Max energy: {p.get('max_step_energy_eV', '?')} eV")

# Step 4: Export shortlist for NEB
print("\nStep 4: Exporting pathway shortlist...")
shortlist = tools.extract_pathway_shortlist_tool(
    search_result=result,
    output_base_dir=OUTPUT_DIR,
    run_id="mech_search",
)

# ── Auto-Visualization ────────────────────────────────────────
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter

# Plot best pathway energy diagram
if retained:
    best = retained[0] if isinstance(retained[0], dict) else {}
    states = best.get("states", [])
    if states:
        energies = [s.get("free_energy_eV", 0.0) for s in states]
        labels = [s.get("species_label", f"S{i}") for i, s in enumerate(states)]

        plotter = EnergyDiagramPlotter()
        plotter.plot(energies=energies, labels=labels)
        diagram_path = os.path.join(OUTPUT_DIR, "best_pathway_energy.png")
        plotter.save(diagram_path)
        print(f"Saved: {diagram_path}")

print(f"\nAll outputs in: {OUTPUT_DIR}")
