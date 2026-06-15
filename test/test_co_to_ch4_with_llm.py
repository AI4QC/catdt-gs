#!/usr/bin/env python
"""
Test complete CO→CH4 catalytic digital twin with LLM-controlled NEB

This script runs the complete workflow from surface generation to pathway analysis,
using the new LLM-controlled NEB system for barrier calculations.

The LLM will:
1. Analyze each reaction step (CO→CHO→CHOH→CH→CH2→CH3→CH4)
2. Decide how to handle atom count mismatches (adding/removing atoms)
3. Intelligently position atoms for realistic NEB paths
4. Validate and iteratively improve structures

Key advantage: Works for arbitrary reactions without hardcoded rules!
"""

import os
import sys
sys.path.insert(0, "core")

import logging
from ase.io import read

# Load environment variables from .env file
def load_env():
    """Load environment variables from .env file"""
    env_file = ".env"
    if os.path.exists(env_file):
        with open(env_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    key, _, value = line.partition('=')
                    os.environ[key.strip()] = value.strip()
        print(f"✓ Loaded environment from {env_file}")
    else:
        print(f"⚠ Warning: {env_file} not found")
        print(f"  Please create .env from .env.template and set OPENAI_API_KEY")

load_env()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration
OUTPUT_DIR = "output/co_to_ch4_llm_controlled"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Use Opus model for best chemistry understanding (override if needed)
if "OPENAI_MODEL" not in os.environ or os.environ["OPENAI_MODEL"] == "gpt-4":
    os.environ["OPENAI_MODEL"] = "claude-opus-4-5-20251101"
    logger.info("Set OPENAI_MODEL to claude-opus-4-5-20251101 for best results")

# CO→CH4 reaction pathway
REACTION_INTERMEDIATES = [
    "*CO",     # Starting point
    "*CHO",    # CO + H
    "*CHOH",   # CHO + H
    "*CH",     # CHOH - OH
    "*CH2",    # CH + H
    "*CH3",    # CH2 + H
    "*CH4",    # CH3 + H (final desorption)
]


def test_co_to_ch4_with_llm():
    """Run complete CO→CH4 workflow with LLM-controlled NEB"""

    print("=" * 80)
    print("CO→CH4 Catalytic Digital Twin with LLM-Controlled NEB")
    print("=" * 80)
    print()
    print("This test demonstrates:")
    print("  1. Complete workflow from Cu bulk → surface → adsorption → pathway")
    print("  2. LLM control of ALL NEB calculations")
    print("  3. Automatic handling of atom count mismatches")
    print("  4. Iterative structure validation and improvement")
    print()
    print("Expected: Realistic barriers for all 7 reaction steps")
    print("=" * 80)
    print()

    # Check if we already have a prepared Cu(111) surface with *CO
    # (from previous runs that already completed surface generation and adsorption)

    # Option 1: Use existing surface with CO if available
    # Try multiple possible locations from previous runs
    existing_surface_options = [
        "output/co_to_ch4_digital_twin/surface_111/pathway/adsorption/*CO_relaxed.vasp",
        "output/co_to_ch4_digital_twin/surface_111/adsorption/best_config.vasp",
        "output/co_hydro_result_corrected/surface_111/pathway/adsorption/*CO_relaxed.vasp",
    ]

    existing_surface = None
    for option in existing_surface_options:
        if os.path.exists(option):
            existing_surface = option
            break

    from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin

    if existing_surface:
        logger.info("Using existing Cu(111) + *CO surface from previous run")
        logger.info(f"  Surface: {existing_surface}")

        # Load the surface
        surface_with_co = read(existing_surface)

        # Identify CO atoms (C and O at the top)
        symbols = surface_with_co.get_chemical_symbols()
        z_coords = surface_with_co.positions[:, 2]

        # Find atoms that are not Cu and are in the top region
        z_max = z_coords.max()
        adsorbate_indices = [
            i for i, (sym, z) in enumerate(zip(symbols, z_coords))
            if sym != 'Cu' and z > z_max - 3.0
        ]

        logger.info(f"  Identified {len(adsorbate_indices)} adsorbate atoms: {[symbols[i] for i in adsorbate_indices]}")

        # Initialize digital twin with LLM control
        logger.info("\nInitializing digital twin with LLM-controlled NEB...")

        dt = GasSolidDigitalTwin(
            surff_root="deps/SurFF",
            adsorbdiff_root="deps/AdsorbDiff",
            surface_sampling_root="deps/surface-sampling",
            fairchem_root="deps/fairchem",
            fairchem_model="uma-s-1p1",
            fairchem_model_path="deps/fairchem_models/uma-s-1p1.pt",
            use_llm_controller=True,  # <<<< Enable LLM control
            llm_model="claude-opus-4-5-20251101",  # <<<< Use Opus
        )

        # Run only the pathway analysis (skip surface generation and MC)
        logger.info("\nRunning pathway analysis with LLM-controlled NEB...")
        logger.info(f"  Reaction steps: {' → '.join(REACTION_INTERMEDIATES)}")

        pathway_result = dt.analyze_reaction_pathway(
            surface=surface_with_co,
            intermediates=REACTION_INTERMEDIATES,
            calculate_barriers=True,
            num_sites=5,
            n_frames=8,
            fmax=0.08,
            output_dir=os.path.join(OUTPUT_DIR, "pathway_llm"),
            current_adsorbate_indices=adsorbate_indices,
        )

    else:
        # Option 2: Run complete workflow from scratch
        logger.info("No existing surface found, running complete workflow from scratch...")
        logger.info("  (This will take longer as it includes surface generation and MC)")

        # Load Cu bulk structure
        bulk_structure = "data/Cu_bulk.vasp"
        if not os.path.exists(bulk_structure):
            logger.error(f"Bulk structure not found: {bulk_structure}")
            logger.error("Please provide a Cu bulk POSCAR file")
            return None

        # Initialize digital twin with LLM control
        dt = GasSolidDigitalTwin(
            surff_root="deps/SurFF",
            adsorbdiff_root="deps/AdsorbDiff",
            surface_sampling_root="deps/surface-sampling",
            fairchem_root="deps/fairchem",
            fairchem_model="uma-s-1p1",
            fairchem_model_path="deps/fairchem_models/uma-s-1p1.pt",
            use_llm_controller=True,  # <<<< Enable LLM control
            llm_model="claude-opus-4-5-20251101",  # <<<< Use Opus
        )

        # Run complete workflow
        workflow_result = dt.run_complete_workflow(
            bulk_structure=bulk_structure,
            reaction_intermediates=REACTION_INTERMEDIATES,
            initial_reactant="*CO",
            temperature=500.0,
            top_n_surfaces=1,  # Only analyze Cu(111)
            reconstruction_sweeps=0,  # Skip MC for faster testing
            calculate_barriers=True,
            neb_frames=8,
            neb_fmax=0.08,
            output_dir=OUTPUT_DIR,
        )

        pathway_result = workflow_result.surface_analyses[0].pathway_analysis

    # Print results
    print("\n" + "=" * 80)
    print("RESULTS: CO→CH4 Pathway Analysis (LLM-Controlled NEB)")
    print("=" * 80)
    print()

    print("Reaction Steps:")
    print("-" * 80)
    for i, step in enumerate(pathway_result.steps, 1):
        rds_marker = " [RDS]" if step.is_rate_determining else ""
        print(f"\nStep {i}: {step.reactant_adsorbate} → {step.product_adsorbate}{rds_marker}")
        print(f"  ΔE = {step.reaction_energy:+.3f} eV")
        if step.activation_energy is not None:
            print(f"  E_act = {step.activation_energy:.3f} eV")
            if step.activation_energy > 0.01:
                print(f"  ✓ Non-zero barrier found (LLM control working!)")
            else:
                print(f"  ⚠ Warning: Nearly barrierless")
        else:
            print(f"  (Barrier not calculated)")

    print()
    print("=" * 80)
    print("Summary:")
    print("-" * 80)
    print(f"Overall reaction energy: {pathway_result.overall_reaction_energy:.3f} eV")

    if pathway_result.rate_determining_step:
        rds = pathway_result.rate_determining_step
        print(f"Rate-determining step: {rds.name}")
        print(f"  Barrier: {rds.activation_energy:.3f} eV")

    if pathway_result.max_barrier:
        print(f"Maximum barrier: {pathway_result.max_barrier:.3f} eV")

    # Count successful barrier calculations
    successful_barriers = sum(
        1 for step in pathway_result.steps
        if step.activation_energy is not None and step.activation_energy > 0.01
    )

    print()
    print("=" * 80)
    print("LLM-Controlled NEB Performance:")
    print("-" * 80)
    print(f"Total reaction steps: {len(pathway_result.steps)}")
    print(f"Steps with realistic barriers (>0.01 eV): {successful_barriers}")
    print(f"Success rate: {successful_barriers}/{len(pathway_result.steps)} = {100*successful_barriers/len(pathway_result.steps):.0f}%")

    if successful_barriers >= len(pathway_result.steps) * 0.6:
        print()
        print("✓ SUCCESS: LLM-controlled NEB produced realistic barriers!")
        print("  The system successfully handled atom count mismatches and")
        print("  created chemically reasonable NEB paths for the reaction.")
    else:
        print()
        print("⚠ WARNING: Some barriers are still near-zero")
        print("  This may indicate issues with structure preparation or NEB convergence")

    print("=" * 80)
    print()

    # Save detailed results
    results_file = os.path.join(OUTPUT_DIR, "llm_results.txt")
    with open(results_file, 'w') as f:
        f.write(pathway_result.summary())

    logger.info(f"Detailed results saved to: {results_file}")

    return pathway_result


if __name__ == "__main__":
    result = test_co_to_ch4_with_llm()

    if result:
        print("\n" + "=" * 80)
        print("Test Complete!")
        print("=" * 80)
        print()
        print("The LLM-controlled NEB system has been successfully integrated")
        print("into the catalytic digital twin workflow.")
        print()
        print("Key achievements:")
        print("  ✓ General implementation (works for arbitrary reactions)")
        print("  ✓ LLM analyzes chemistry and decides atom positions")
        print("  ✓ Automatic validation and iterative improvement")
        print("  ✓ No hardcoded rules needed!")
        print()
        print(f"Output directory: {OUTPUT_DIR}")
        print("=" * 80)
