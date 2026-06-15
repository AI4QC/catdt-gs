"""
Test VSSR-MC Predictor

Test script for the VSSRMCPredictor wrapper class.
Tests automatic offset_data generation and surface reconstruction sampling.
"""

import os
import sys

# Add paths
sys.path.insert(0, "core/reconstruction")
sys.path.insert(0, "deps/surface-sampling")


def test_cu100_surface():
    """Test with a simple Cu(100) surface using CHGNet

    This test demonstrates the simplest usage - just provide a surface
    and the adsorbate elements. offset_data is automatically generated.
    """
    print("=" * 60)
    print("Test 1: Cu(100) surface with auto offset_data")
    print("=" * 60)

    from ase.build import fcc100
    from ase.constraints import FixAtoms

    # Create Cu(100) surface
    slab = fcc100('Cu', size=(2, 2, 4), vacuum=15.0)

    # Set constraints - fix bottom 2 layers
    mask = [atom.position[2] < slab.positions[:, 2].min() + 3.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    print(f"Surface formula: {slab.get_chemical_formula()}")
    print(f"Number of atoms: {len(slab)}")

    from vssr_mc_predictor import VSSRMCPredictor

    predictor = VSSRMCPredictor(
        surface_sampling_root="deps/surface-sampling",
        model_type="CHGNetNFF",
        device="cuda",
        verbose=True,
        keep_files=True,
    )

    # Simple usage - offset_data is auto-generated
    result = predictor.sample(
        surface=slab,
        adsorbates=["Cu"],  # Can add/remove Cu atoms on surface
        canonical=False,  # Semi-grand canonical
        total_sweeps=5,
        sweep_size=5,
        temperature=1.0,
        output_dir="output/vssr/Cu100_auto",
    )

    print("\n" + "=" * 60)
    print("Test 1 Results:")
    print("=" * 60)
    print(result.summary(3))
    print(f"\nAuto-generated offset_data:")
    print(f"  {result.offset_data}")

    return result


def test_srtio3_surface():
    """Test with SrTiO3(001) surface from tutorial data

    Uses the prepared slab from the tutorials.
    """
    print("\n" + "=" * 60)
    print("Test 2: SrTiO3(001) surface with auto offset_data")
    print("=" * 60)

    import pickle

    slab_path = "deps/surface-sampling/tutorials/data/SrTiO3_001/SrTiO3_001_2x2_pristine_slab.pkl"

    if not os.path.exists(slab_path):
        print(f"Tutorial data not found at {slab_path}")
        print("Skipping test 2")
        return None

    with open(slab_path, "rb") as f:
        slab = pickle.load(f)

    print(f"Surface formula: {slab.get_chemical_formula()}")
    print(f"Number of atoms: {len(slab)}")

    from vssr_mc_predictor import VSSRMCPredictor

    predictor = VSSRMCPredictor(
        surface_sampling_root="deps/surface-sampling",
        model_type="CHGNetNFF",
        device="cuda",
        verbose=True,
        keep_files=True,
    )

    # Simple usage - no need to manually specify offset_data
    result = predictor.sample(
        surface=slab,
        adsorbates=["O", "Sr", "Ti"],  # Elements that can be added/removed
        canonical=False,
        total_sweeps=3,
        sweep_size=3,
        temperature=1.0,
        output_dir="output/vssr/SrTiO3_001_auto",
    )

    print("\n" + "=" * 60)
    print("Test 2 Results:")
    print("=" * 60)
    print(result.summary(3))
    print(f"\nAuto-generated offset_data:")
    import json
    print(json.dumps(result.offset_data, indent=2))

    return result


def test_with_manual_offset():
    """Test with manually provided offset_data (for comparison)

    This test shows how to use manually provided offset_data,
    which can be useful when you have pre-computed values.
    """
    print("\n" + "=" * 60)
    print("Test 3: SrTiO3(001) with manual offset_data")
    print("=" * 60)

    import pickle
    import json

    slab_path = "deps/surface-sampling/tutorials/data/SrTiO3_001/SrTiO3_001_2x2_pristine_slab.pkl"
    offset_data_path = "deps/surface-sampling/tutorials/data/SrTiO3_001/nff/offset_data_chgnet_eV.json"

    if not os.path.exists(slab_path) or not os.path.exists(offset_data_path):
        print("Tutorial data not found, skipping test 3")
        return None

    with open(slab_path, "rb") as f:
        slab = pickle.load(f)

    with open(offset_data_path, "r") as f:
        offset_data = json.load(f)

    print(f"Surface formula: {slab.get_chemical_formula()}")
    print(f"Manual offset_data: {offset_data}")

    from vssr_mc_predictor import VSSRMCPredictor

    predictor = VSSRMCPredictor(
        surface_sampling_root="deps/surface-sampling",
        model_type="CHGNetNFF",
        device="cuda",
        verbose=True,
        keep_files=True,
    )

    # Manual offset_data and chem_pots
    result = predictor.sample(
        surface=slab,
        adsorbates=["O", "Sr", "Ti"],
        canonical=False,
        total_sweeps=3,
        sweep_size=3,
        temperature=1.0,
        offset_data=offset_data,  # Manually provided
        chem_pots={"Sr": -2, "Ti": 0, "O": 0},  # Manually provided
        auto_offset=False,  # Disable auto-generation
        output_dir="output/vssr/SrTiO3_001_manual",
    )

    print("\n" + "=" * 60)
    print("Test 3 Results:")
    print("=" * 60)
    print(result.summary(3))

    return result


def test_canonical_sampling():
    """Test canonical (NVT) sampling with fixed number of adsorbates"""
    print("\n" + "=" * 60)
    print("Test 4: Canonical sampling on Cu(100)")
    print("=" * 60)

    from ase.build import fcc100
    from ase.constraints import FixAtoms

    # Create Cu(100) surface
    slab = fcc100('Cu', size=(2, 2, 4), vacuum=15.0)

    # Set constraints
    mask = [atom.position[2] < slab.positions[:, 2].min() + 3.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    print(f"Surface formula: {slab.get_chemical_formula()}")

    from vssr_mc_predictor import VSSRMCPredictor

    predictor = VSSRMCPredictor(
        surface_sampling_root="deps/surface-sampling",
        model_type="CHGNetNFF",
        device="cuda",
        verbose=True,
        keep_files=True,
    )

    # Canonical sampling - fixed number of adsorbates
    result = predictor.sample(
        surface=slab,
        adsorbates=["Cu"],
        canonical=True,  # Fixed composition
        num_adsorbates=4,  # Keep 4 Cu atoms on surface
        total_sweeps=3,
        sweep_size=3,
        temperature=1.0,
        output_dir="output/vssr/Cu100_canonical",
    )

    print("\n" + "=" * 60)
    print("Test 4 Results:")
    print("=" * 60)
    print(result.summary(3))

    return result


def main():
    """Run all tests"""
    print("=" * 60)
    print("VSSR-MC Predictor Test Suite")
    print("=" * 60)

    # Create output directory
    os.makedirs("output/vssr", exist_ok=True)

    results = {}

    # Test 1: Simple Cu surface with auto offset_data
    try:
        results["test1"] = test_cu100_surface()
    except Exception as e:
        print(f"Test 1 failed: {e}")
        import traceback
        traceback.print_exc()
        results["test1"] = None

    # Test 2: SrTiO3 surface with auto offset_data
    try:
        results["test2"] = test_srtio3_surface()
    except Exception as e:
        print(f"Test 2 failed: {e}")
        import traceback
        traceback.print_exc()
        results["test2"] = None

    # Test 3: Manual offset_data
    try:
        results["test3"] = test_with_manual_offset()
    except Exception as e:
        print(f"Test 3 failed: {e}")
        import traceback
        traceback.print_exc()
        results["test3"] = None

    # Test 4: Canonical sampling
    try:
        results["test4"] = test_canonical_sampling()
    except Exception as e:
        print(f"Test 4 failed: {e}")
        import traceback
        traceback.print_exc()
        results["test4"] = None

    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    for test_name, result in results.items():
        status = "PASS" if result is not None else "FAIL"
        print(f"  {test_name}: {status}")

    print(f"\nOutput directory: output/vssr/")


if __name__ == "__main__":
    main()
