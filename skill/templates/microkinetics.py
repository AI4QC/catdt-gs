"""CatDT Template: Microkinetic Modeling (CatMAP)

Computes turnover frequency (TOF), selectivity, and coverage from
elementary reaction barriers. Supports both thermal and electrochemical.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

TEMPERATURE_K = 500.0
PRESSURES = {"CO_g": 1.0, "H2_g": 1.0}  # atm
OUTPUT_DIR = os.path.join("output", "kmc", datetime.now().strftime("%Y%m%d_%H%M%S"))

# Define reactions and energetics (replace with your data)
REACTIONS = [
    {"expression": "*_s + CO_g -> CO*", "activation_energy": 0.0},
    {"expression": "CO* + H* <-> CO-H* -> CHO*", "activation_energy": 0.85},
    {"expression": "CHO* + H* <-> CHO-H* -> CH2O*", "activation_energy": 0.72},
    {"expression": "CH2O* -> CH2O_g + *_s", "activation_energy": 0.45},
]

GAS_SPECIES = {
    "CO_g": {"formation_energy": 0.0, "frequencies": [2170]},
    "H2_g": {"formation_energy": 0.0, "frequencies": [4401]},
    "CH2O_g": {"formation_energy": -1.5, "frequencies": [2843, 1746, 1500, 1249, 1167]},
}

ADSORBATE_ENERGIES = {
    "CO*": -0.8,
    "H*": -0.3,
    "CHO*": -0.5,
    "CH2O*": -0.2,
}

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from core.kmc.catmap_predictor import CatMAPPredictor

kmc = CatMAPPredictor(
    catmap_root=os.path.join(CATDT_ROOT, "deps/catmap"),
    work_dir=OUTPUT_DIR,
)

# Add gas species
for name, data in GAS_SPECIES.items():
    kmc.add_gas(name, data["formation_energy"], data.get("frequencies", []))

# Add adsorbates
for name, energy in ADSORBATE_ENERGIES.items():
    kmc.add_adsorbate(name, {"default": energy})

# Add reactions
for rxn in REACTIONS:
    kmc.add_reaction(rxn["expression"])

# Set conditions
kmc.set_conditions(temperature=TEMPERATURE_K, pressures=PRESSURES)

# Run
result = kmc.run()

# Print results
print(f"Temperature: {TEMPERATURE_K} K")
print(f"Pressures: {PRESSURES}")
print(f"TOF: {result.tof:.4e} s^-1 site^-1")
print(f"\nCoverages:")
for species, cov in result.coverages.items():
    print(f"  {species}: {cov:.4f}")
if hasattr(result, "rate_controlling_step"):
    print(f"\nRate-controlling step: {result.rate_controlling_step}")

# ── Auto-Visualization ────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Coverage bar chart
if result.coverages:
    fig, ax = plt.subplots(figsize=(8, 5))
    species = list(result.coverages.keys())
    coverages = list(result.coverages.values())
    ax.bar(species, coverages, color="steelblue")
    ax.set_ylabel("Coverage")
    ax.set_title(f"Surface Coverage at {TEMPERATURE_K} K")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "coverage.png"), dpi=200)
    plt.close(fig)
    print(f"Saved: coverage.png")

print(f"\nAll outputs in: {OUTPUT_DIR}")
