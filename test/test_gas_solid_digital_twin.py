#!/usr/bin/env python
"""
Test script for Gas-Solid Catalysis Digital Twin

This script demonstrates the complete workflow for gas-solid interface
catalysis, from bulk structure to kinetic analysis.

Example: CO oxidation on Pt catalyst
CO + 1/2 O2 -> CO2

Usage:
    python test/test_gas_solid_digital_twin.py
"""

import os
import sys
import logging

# Add core to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))

from ase import Atoms
from ase.build import bulk, surface
from ase.io import write

from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin, run_gas_solid_catalysis


def create_test_bulk():
    """Create a simple Pt bulk structure for testing"""
    pt_bulk = bulk('Pt', 'fcc', a=3.92)
    return pt_bulk


def test_simple_workflow():
    """Test the simplest workflow using convenience function"""
    print("=" * 80)
    print("Test 1: Simple Workflow (CO oxidation on Pt)")
    print("=" * 80)

    # Create test bulk structure
    pt_bulk = create_test_bulk()
    os.makedirs("output/test_gas_solid", exist_ok=True)
    write("output/test_gas_solid/POSCAR_Pt", pt_bulk)

    try:
        # Run simple workflow
        result = run_gas_solid_catalysis(
            bulk_structure="output/test_gas_solid/POSCAR_Pt",
            reaction_intermediates=["*CO", "*O", "*CO2", "*"],
            initial_reactant="*CO",
            temperature=500.0,
            top_n_surfaces=2,  # Only analyze top 2 surfaces for speed
            num_adsorption_sites=5,  # Fewer sites for speed
            reconstruction_sweeps=20,  # Fewer sweeps for speed
            calculate_barriers=False,  # Skip NEB for speed in testing
            output_dir="output/test_gas_solid_simple",
        )

        print("\n" + result.summary())
        print("\n✅ Simple workflow test PASSED")
        return True

    except Exception as e:
        print(f"\n❌ Simple workflow test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_step_by_step():
    """Test step-by-step workflow with full control"""
    print("\n" + "=" * 80)
    print("Test 2: Step-by-Step Workflow")
    print("=" * 80)

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    try:
        # Initialize digital twin
        dt = GasSolidDigitalTwin(
            surff_root="deps/SurFF",
            adsorbdiff_root="deps/AdsorbDiff",
            surface_sampling_root="deps/surface-sampling",
            fairchem_root="deps/fairchem",
            fairchem_model="uma-s-1p1",
            vssr_mc_model="CHGNetNFF",
            use_gpu=False,  # Use CPU for testing
        )

        # Create Pt bulk
        pt_bulk = create_test_bulk()

        # Step 1: Generate surfaces
        print("\n" + "-" * 80)
        print("Step 1: Generating surfaces from bulk...")
        print("-" * 80)
        surff_result, slabs = dt.generate_surfaces_from_bulk(
            bulk_structure=pt_bulk,
            top_n=2,
            output_dir="output/test_gas_solid_steps/surfaces",
        )
        print(f"✅ Generated {len(slabs)} surfaces")

        # Step 2: Predict adsorption sites (on first surface only)
        print("\n" + "-" * 80)
        print("Step 2: Predicting adsorption sites...")
        print("-" * 80)
        adsorption_result = dt.predict_adsorption_sites(
            surface=slabs[0],
            adsorbate="*CO",
            num_sites=5,
            output_dir="output/test_gas_solid_steps/adsorption",
        )
        print(f"✅ Found best site with energy {adsorption_result.best_result.energy:.3f} eV")

        # Step 3: Surface reconstruction
        print("\n" + "-" * 80)
        print("Step 3: Simulating surface reconstruction...")
        print("-" * 80)
        surface_with_adsorbate = adsorption_result.best_result.final_structure
        reconstruction_result = dt.simulate_surface_reconstruction(
            surface=surface_with_adsorbate,
            adsorbates=["C", "O"],  # Elements that can adsorb
            temperature=500.0,
            total_sweeps=20,  # Fewer sweeps for testing
            output_dir="output/test_gas_solid_steps/reconstruction",
        )
        print(f"✅ Sampled {len(reconstruction_result.structures)} structures")
        print(f"   Lowest energy: {reconstruction_result.lowest_energy_structure.energy:.3f} eV")

        # Step 4: Reaction pathway analysis (without barriers for speed)
        print("\n" + "-" * 80)
        print("Step 4: Analyzing reaction pathway...")
        print("-" * 80)
        reconstructed_surface = reconstruction_result.lowest_energy_structure.atoms
        pathway_result = dt.analyze_reaction_pathway(
            surface=reconstructed_surface,
            intermediates=["*CO", "*O", "*CO2", "*"],
            calculate_barriers=False,  # Skip NEB for speed
            num_sites=3,
            output_dir="output/test_gas_solid_steps/pathway",
        )
        print(f"✅ Pathway analysis complete")
        print(f"   Overall ΔE: {pathway_result.overall_reaction_energy:.3f} eV")

        print("\n✅ Step-by-step workflow test PASSED")
        return True

    except Exception as e:
        print(f"\n❌ Step-by-step workflow test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_multiple_surfaces():
    """Test workflow with multiple surface facets"""
    print("\n" + "=" * 80)
    print("Test 3: Multiple Surface Facets Analysis")
    print("=" * 80)

    try:
        pt_bulk = create_test_bulk()
        write("output/test_gas_solid/POSCAR_Pt", pt_bulk)

        result = run_gas_solid_catalysis(
            bulk_structure="output/test_gas_solid/POSCAR_Pt",
            reaction_intermediates=["*O", "*OH", "*"],  # Simpler reaction
            initial_reactant="*O",
            temperature=300.0,
            top_n_surfaces=3,  # Analyze 3 surfaces
            num_adsorption_sites=3,
            reconstruction_sweeps=10,  # Very short for testing
            calculate_barriers=False,
            output_dir="output/test_gas_solid_multi",
        )

        print(f"\n✅ Analyzed {len(result.surface_analyses)} surfaces")
        for i, analysis in enumerate(result.surface_analyses, 1):
            print(f"   Surface {i}: {analysis.miller_index}")
            print(f"      Area fraction: {analysis.area_fraction:.2%}")
            if analysis.pathway_analysis:
                print(f"      Overall ΔE: {analysis.pathway_analysis.overall_reaction_energy:.3f} eV")

        if result.best_surface:
            print(f"\n   Best surface: {result.best_surface.miller_index}")

        print("\n✅ Multiple surfaces test PASSED")
        return True

    except Exception as e:
        print(f"\n❌ Multiple surfaces test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_nh3_synthesis():
    """Test NH3 synthesis reaction on Fe catalyst"""
    print("\n" + "=" * 80)
    print("Test 4: NH3 Synthesis on Fe (N2 + 3H2 -> 2NH3)")
    print("=" * 80)

    try:
        # Create Fe bulk
        fe_bulk = bulk('Fe', 'bcc', a=2.87)
        os.makedirs("output/test_gas_solid", exist_ok=True)
        write("output/test_gas_solid/POSCAR_Fe", fe_bulk)

        result = run_gas_solid_catalysis(
            bulk_structure="output/test_gas_solid/POSCAR_Fe",
            reaction_intermediates=["*N2", "*N", "*NH", "*NH2", "*NH3", "*"],
            initial_reactant="*N2",
            temperature=673.0,  # 400°C, typical for NH3 synthesis
            top_n_surfaces=2,
            num_adsorption_sites=3,
            reconstruction_sweeps=10,
            calculate_barriers=False,
            output_dir="output/test_gas_solid_nh3",
        )

        print("\n" + result.summary(top_n=2))
        print("\n✅ NH3 synthesis test PASSED")
        return True

    except Exception as e:
        print(f"\n❌ NH3 synthesis test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests"""
    print("=" * 80)
    print("Gas-Solid Catalysis Digital Twin - Test Suite")
    print("=" * 80)
    print("\nNote: These tests use reduced parameters for speed.")
    print("For production runs, use larger values for:")
    print("  - top_n_surfaces (5-10)")
    print("  - num_adsorption_sites (10-20)")
    print("  - reconstruction_sweeps (100-200)")
    print("  - calculate_barriers=True (with neb_frames=10-15)")
    print("=" * 80)

    results = []

    # Run tests
    results.append(("Simple Workflow", test_simple_workflow()))
    results.append(("Step-by-Step", test_step_by_step()))
    results.append(("Multiple Surfaces", test_multiple_surfaces()))
    results.append(("NH3 Synthesis", test_nh3_synthesis()))

    # Summary
    print("\n" + "=" * 80)
    print("Test Summary")
    print("=" * 80)
    for name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{name:30s} {status}")

    all_passed = all(r[1] for r in results)
    print("=" * 80)
    if all_passed:
        print("🎉 All tests PASSED!")
    else:
        print("⚠️  Some tests FAILED. Please check the output above.")
    print("=" * 80)

    return all_passed


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
