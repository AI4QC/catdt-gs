"""CatDT Template: Structure Relaxation / Single-Point Energy / MD

Uses FairChem UMA model for:
  - Single-point energy calculation
  - Structure optimization (relaxation) with trajectory visualization
  - Molecular dynamics simulation with trajectory GIF

Auto-generates before/after PNGs, trajectory GIF, and energy evolution plot.
"""

import os, sys
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

STRUCTURE_PATH = "structure.vasp"
MODEL = "uma-s-1p1"                # "uma-s-1p1" or "uma-m-1p1"
OUTPUT_DIR = os.path.join("output", "relax", datetime.now().strftime("%Y%m%d_%H%M%S"))

# Mode: "single_point", "relax", or "md"
MODE = "relax"

# Relaxation parameters
RELAX_FMAX = 0.05                  # eV/A
RELAX_MAX_STEPS = 300

# MD parameters (only used if MODE="md")
MD_TEMPERATURE_K = 500.0
MD_TIMESTEP_FS = 1.0
MD_STEPS = 1000
MD_SAVE_INTERVAL = 10

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from ase.io import read, write, Trajectory
from core.pathway.fairchem_predictor import FairchemPredictor

fairchem = FairchemPredictor(
    fairchem_root=os.path.join(CATDT_ROOT, "deps/fairchem"),
    model_name=MODEL,
    use_gpu=True,
)

structure = read(STRUCTURE_PATH)

if MODE == "single_point":
    # ── Single-point energy ──
    result = fairchem.predict_energy(structure, relax=False)
    print(f"Energy: {result.energy:.4f} eV")
    print(f"Max force: {result.forces_max:.4f} eV/A")

elif MODE == "relax":
    # ── Structure optimization ──
    from ase.optimize import LBFGS

    initial = structure.copy()
    calc = fairchem._get_calculator()
    structure.calc = calc

    # Save trajectory
    traj_path = os.path.join(OUTPUT_DIR, "relax.traj")
    traj = Trajectory(traj_path, "w")

    opt = LBFGS(structure, trajectory=traj)
    opt.run(fmax=RELAX_FMAX, steps=RELAX_MAX_STEPS)
    traj.close()

    print(f"Initial energy: {initial.get_potential_energy():.4f} eV" if initial.calc else "")
    print(f"Final energy: {structure.get_potential_energy():.4f} eV")
    print(f"Steps: {opt.nsteps}")

    # Save relaxed structure
    write(os.path.join(OUTPUT_DIR, "relaxed.vasp"), structure, format="vasp")

    # ── Auto-Visualization ──
    from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    viz = CatalystSurfaceVisualizer(quality="high")

    # Before/after PNGs
    viz.visualize_structure(initial, os.path.join(OUTPUT_DIR, "before.png"), title="Before Relaxation")
    viz.visualize_structure(structure, os.path.join(OUTPUT_DIR, "after.png"), title="After Relaxation")

    # Trajectory GIF
    frames = read(traj_path, index=":")
    if len(frames) > 1:
        # Sample for manageable GIF
        step = max(1, len(frames) // 20)
        sampled = frames[::step] + [frames[-1]]
        viz.visualize_trajectory(
            sampled,
            output_file=os.path.join(OUTPUT_DIR, "relaxation.gif"),
            fps=5,
        )
        print(f"Saved: relaxation.gif ({len(sampled)} frames)")

    # Energy evolution plot
    energies = [f.get_potential_energy() for f in frames if f.calc is not None]
    if energies:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(energies, "b-o", markersize=3)
        ax.set_xlabel("Optimization Step")
        ax.set_ylabel("Energy (eV)")
        ax.set_title("Relaxation Energy Convergence")
        fig.tight_layout()
        fig.savefig(os.path.join(OUTPUT_DIR, "energy_convergence.png"), dpi=200)
        plt.close(fig)
        print(f"Saved: energy_convergence.png")

elif MODE == "md":
    # ── Molecular Dynamics ──
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase import units

    calc = fairchem._get_calculator()
    structure.calc = calc

    MaxwellBoltzmannDistribution(structure, temperature_K=MD_TEMPERATURE_K)

    traj_path = os.path.join(OUTPUT_DIR, "md.traj")
    traj = Trajectory(traj_path, "w")

    dyn = Langevin(structure, MD_TIMESTEP_FS * units.fs, temperature_K=MD_TEMPERATURE_K, friction=0.01)
    dyn.attach(traj.write, interval=MD_SAVE_INTERVAL)

    # Collect energies
    energies, temperatures = [], []
    def record():
        energies.append(structure.get_potential_energy())
        temperatures.append(structure.get_temperature())
    dyn.attach(record, interval=MD_SAVE_INTERVAL)

    print(f"Running MD: {MD_STEPS} steps at {MD_TEMPERATURE_K} K...")
    dyn.run(MD_STEPS)
    traj.close()

    # ── Auto-Visualization ──
    from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    viz = CatalystSurfaceVisualizer(quality="high")

    # Trajectory GIF
    frames = read(traj_path, index=":")
    step = max(1, len(frames) // 30)
    sampled = frames[::step]
    viz.visualize_trajectory(
        sampled,
        output_file=os.path.join(OUTPUT_DIR, "md_trajectory.gif"),
        fps=5,
    )
    print(f"Saved: md_trajectory.gif ({len(sampled)} frames)")

    # T and E plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    time_ps = [i * MD_SAVE_INTERVAL * MD_TIMESTEP_FS / 1000 for i in range(len(temperatures))]
    ax1.plot(time_ps, temperatures, "r-")
    ax1.set_ylabel("Temperature (K)")
    ax1.axhline(MD_TEMPERATURE_K, color="k", linestyle="--", alpha=0.5)
    ax2.plot(time_ps, energies, "b-")
    ax2.set_ylabel("Potential Energy (eV)")
    ax2.set_xlabel("Time (ps)")
    fig.suptitle(f"MD at {MD_TEMPERATURE_K} K")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "md_thermo.png"), dpi=200)
    plt.close(fig)
    print(f"Saved: md_thermo.png")

print(f"\nAll outputs in: {OUTPUT_DIR}")
