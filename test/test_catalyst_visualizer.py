#!/usr/bin/env python
"""
Test script for catalyst_surface_visualizer.py

Tests visualization with structures from data/ and deps/ directories.
Outputs are saved to output/visualization/ directory.
"""

import sys
from pathlib import Path

# Add core to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from ase.io import read

def test_single_structure():
    """Test visualization of a single structure from trajectory."""
    print("\n" + "="*70)
    print("Test 1: Single Structure Visualization")
    print("="*70)

    # Read a single frame from trajectory
    traj_file = "data/ocp/traj/example/101_0.traj"
    atoms = read(traj_file, index=0)

    print(f"Loaded structure from: {traj_file}")
    print(f"Number of atoms: {len(atoms)}")
    print(f"Elements: {', '.join(sorted(set(atoms.get_chemical_symbols())))}")

    # Create output directory
    output_dir = Path("output/visualization")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Visualize with matplotlib (fast, no OVITO required)
    visualizer = CatalystSurfaceVisualizer(
        elevation=30,
        azimuth=45,
        scale=1.0,
        quality='high',
        renderer='matplotlib',
        show_cell=False,
        auto_expand=True
    )

    output_file = output_dir / "single_structure_101_0.png"
    visualizer.visualize_structure(atoms, str(output_file))
    print(f"✓ Saved: {output_file}")


def test_trajectory_animation():
    """Test trajectory animation."""
    print("\n" + "="*70)
    print("Test 2: Trajectory Animation")
    print("="*70)

    # Read trajectory (first 10 frames for speed)
    traj_file = "data/ocp/traj/example/101_0.traj"
    trajectory = read(traj_file, index=':10')

    print(f"Loaded trajectory from: {traj_file}")
    print(f"Number of frames: {len(trajectory)}")
    print(f"Atoms per frame: {len(trajectory[0])}")
    print(f"Elements: {', '.join(sorted(set(trajectory[0].get_chemical_symbols())))}")

    # Create output directory
    output_dir = Path("output/visualization")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create animation
    visualizer = CatalystSurfaceVisualizer(
        elevation=30,
        azimuth=45,
        scale=1.0,
        quality='medium',  # Use medium quality for faster rendering
        renderer='matplotlib',
        show_cell=False,
        auto_expand=True
    )

    output_file = output_dir / "trajectory_101_0.gif"
    visualizer.visualize_trajectory(trajectory, str(output_file), fps=2)
    print(f"✓ Saved: {output_file}")


def test_different_angles():
    """Test visualization from different viewing angles."""
    print("\n" + "="*70)
    print("Test 3: Different Viewing Angles")
    print("="*70)

    # Read a structure
    traj_file = "data/ocp/traj/example/102_0.traj"
    atoms = read(traj_file, index=-1)  # Last frame (relaxed structure)

    print(f"Loaded structure from: {traj_file} (last frame)")
    print(f"Number of atoms: {len(atoms)}")

    # Create output directory
    output_dir = Path("output/visualization")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test different angles
    angles = [
        (30, 45, "default"),
        (0, 0, "top_view"),
        (90, 0, "side_view"),
        (45, 135, "angled")
    ]

    for elevation, azimuth, label in angles:
        visualizer = CatalystSurfaceVisualizer(
            elevation=elevation,
            azimuth=azimuth,
            scale=1.0,
            quality='medium',
            renderer='matplotlib',
            show_cell=False,
            auto_expand=True
        )

        output_file = output_dir / f"angles_102_0_{label}.png"
        visualizer.visualize_structure(atoms, str(output_file))
        print(f"✓ Saved {label} (elev={elevation}°, azim={azimuth}°): {output_file}")


def test_with_cell():
    """Test visualization with cell boundaries."""
    print("\n" + "="*70)
    print("Test 4: Visualization with Cell Boundaries")
    print("="*70)

    # Read a structure
    traj_file = "data/ocp/traj/example/103_0.traj"
    atoms = read(traj_file, index=0)

    print(f"Loaded structure from: {traj_file}")
    print(f"Number of atoms: {len(atoms)}")

    # Create output directory
    output_dir = Path("output/visualization")
    output_dir.mkdir(parents=True, exist_ok=True)

    # With cell
    visualizer = CatalystSurfaceVisualizer(
        elevation=30,
        azimuth=45,
        scale=1.0,
        quality='high',
        renderer='matplotlib',
        show_cell=True,  # Enable cell display (note: matplotlib doesn't show cell)
        auto_expand=True
    )

    output_file = output_dir / "with_cell_103_0.png"
    visualizer.visualize_structure(atoms, str(output_file))
    print(f"✓ Saved: {output_file}")
    print("Note: Cell boundaries require OVITO renderer (tachyon/ospray/opengl/anari)")


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("Catalyst Surface Visualizer Test Suite")
    print("="*70)
    print("\nTesting with structures from data/ocp/traj/example/")
    print("Output directory: output/visualization/")

    try:
        # Test 1: Single structure
        test_single_structure()

        # Test 2: Trajectory animation
        test_trajectory_animation()

        # Test 3: Different angles
        test_different_angles()

        # Test 4: With cell
        test_with_cell()

        print("\n" + "="*70)
        print("All Tests Completed Successfully!")
        print("="*70)
        print("\nGenerated files in output/visualization/:")
        output_dir = Path("output/visualization")
        for f in sorted(output_dir.glob("*")):
            size_mb = f.stat().st_size / 1024 / 1024
            print(f"  - {f.name} ({size_mb:.2f} MB)")

        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
