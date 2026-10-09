#!/usr/bin/env python3
"""
Subprocess script for running VSSR-MC in isolation.

This script runs VSSR-MC in a separate process to avoid GPU/MPI state corruption
from other models (SurFF, AdsorbDiff) that may have been loaded in the main process.
"""

import argparse
import json
import os
import pickle
import sys

# Increase recursion limit at module level
sys.setrecursionlimit(50000)

# Add paths
_this_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.join(_this_dir, "..", "..")
sys.path.insert(0, _project_root)  # project root for `core.reconstruction.*`


def main():
    parser = argparse.ArgumentParser(description="Run VSSR-MC in subprocess")
    parser.add_argument("--surface", required=True, help="Path to pickled surface structure")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    parser.add_argument("--output", required=True, help="Path to output pickle file")
    args = parser.parse_args()

    # Load surface (with adsorbate)
    with open(args.surface, "rb") as f:
        surface = pickle.load(f)

    # Load config
    with open(args.config, "r") as f:
        config = json.load(f)

    # Load clean_slab (without adsorbate) for virtual site generation
    # If clean_slab_path is absent or None, SurfaceSystem will enumerate
    # virtual sites on the seed (with overlayer) — use_seed_for_virtual_sites.
    clean_slab = None
    if config.get("clean_slab_path"):
        with open(config["clean_slab_path"], "rb") as f:
            clean_slab = pickle.load(f)
        print(f"Loaded clean_slab: {clean_slab.get_chemical_formula()} ({len(clean_slab)} atoms)")
    else:
        print("clean_slab_path not provided; virtual sites will be enumerated on the seed structure.")

    # Import VSSR-MC predictor
    from core.reconstruction.vssr_mc_predictor import VSSRMCPredictor

    # Create predictor
    predictor = VSSRMCPredictor(
        surface_sampling_root=config.get("surface_sampling_root", "deps/surface-sampling"),
        model_type=config.get("model_type", "CHGNetNFF"),
        device=config.get("device", "cuda"),
        verbose=config.get("verbose", True),
        keep_files=config.get("keep_files", True),
        potential_she=config.get("potential_she"),
        ph=config.get("ph"),
    )

    try:
        # Run sampling
        sample_kwargs = dict(
            surface=surface,
            clean_slab=clean_slab,
            adsorbates=config.get("adsorbates", []),
            canonical=config.get("canonical", False),
            num_adsorbates=config.get("num_adsorbates", 0),
            adsorbate_counts=config.get("adsorbate_counts"),
            total_sweeps=config.get("total_sweeps", 50),
            sweep_size=config.get("sweep_size", 20),
            temperature=config.get("temperature", 1.0),
            output_dir=config.get("output_dir"),
            surface_indices=config.get("surface_indices"),
            adsorbate_indices=config.get("adsorbate_indices"),
            potential_she=config.get("potential_she"),
            ph=config.get("ph"),
        )
        if config.get("adsorbate_exclusion_radius_A") is not None:
            sample_kwargs["system_settings"] = {
                "adsorbate_exclusion_radius_A": float(
                    config["adsorbate_exclusion_radius_A"]
                )
            }
        if config.get("existing_atom_exclusion_radius_A") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "existing_atom_exclusion_radius_A"
            ] = float(config["existing_atom_exclusion_radius_A"])
        if config.get("min_virtual_site_distance_A") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "min_virtual_site_distance_A"
            ] = float(config["min_virtual_site_distance_A"])
        if config.get("max_virtual_site_distance_to_surface_A") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "max_virtual_site_distance_to_surface_A"
            ] = float(config["max_virtual_site_distance_to_surface_A"])
        if config.get("virtual_site_planar_distance_A") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "planar_distance"
            ] = float(config["virtual_site_planar_distance_A"])
        if config.get("virtual_site_near_reduce") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "near_reduce"
            ] = float(config["virtual_site_near_reduce"])
        if config.get("virtual_site_no_obtuse_hollow") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "no_obtuse_hollow"
            ] = bool(config["virtual_site_no_obtuse_hollow"])
        if config.get("virtual_site_min_count") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "virtual_site_min_count"
            ] = int(config["virtual_site_min_count"])
        if config.get("virtual_site_local_expansion_radius_A") is not None:
            sample_kwargs.setdefault("system_settings", {})[
                "virtual_site_local_expansion_radius_A"
            ] = float(config["virtual_site_local_expansion_radius_A"])
        # Allow passing pre-computed offset_data to ensure consistent
        # energy reference across multiple reconstruction steps
        if config.get("offset_data") is not None:
            sample_kwargs["offset_data"] = config["offset_data"]
            sample_kwargs["auto_offset"] = False
        # Per-element chemical potentials (eV). Non-zero μ_O prevents
        # Ti-only runaway growth on Cu substrates.
        if config.get("chem_pots"):
            sample_kwargs["chem_pots"] = dict(config["chem_pots"])
        result = predictor.sample(**sample_kwargs)

        # Save result - include all VSSRMCResult and SampledStructure fields
        output_data = {
            "surface_name": result.surface_name,
            "model_type": result.model_type,
            "total_sweeps": result.total_sweeps,
            "sweep_size": result.sweep_size,
            "temperature": result.temperature,
            "canonical": result.canonical,
            "structures": [
                {
                    "atoms": s.atoms.copy(),
                    "energy": s.energy,
                    "total_energy": s.total_energy,
                    "sweep_number": s.sweep_number,
                    "num_adsorbates": s.num_adsorbates,
                    "acceptance_rate": s.acceptance_rate,
                }
                for s in result.structures
            ],
            "energy_history": result.energy_history,
            "acceptance_history": result.acceptance_history,
            "adsorption_count_history": result.adsorption_count_history,
            "best_structure": {
                "atoms": result.best_structure.atoms.copy(),
                "energy": result.best_structure.energy,
                "total_energy": result.best_structure.total_energy,
                "sweep_number": result.best_structure.sweep_number,
                "num_adsorbates": result.best_structure.num_adsorbates,
                "acceptance_rate": result.best_structure.acceptance_rate,
            },
            "lowest_energy_structure": {
                "atoms": result.lowest_energy_structure.atoms.copy(),
                "energy": result.lowest_energy_structure.energy,
                "total_energy": result.lowest_energy_structure.total_energy,
                "sweep_number": result.lowest_energy_structure.sweep_number,
                "num_adsorbates": result.lowest_energy_structure.num_adsorbates,
                "acceptance_rate": result.lowest_energy_structure.acceptance_rate,
            },
            "offset_data": result.offset_data,
            "output_dir": result.output_dir,
            "success": True,
        }
        with open(args.output, "wb") as f:
            pickle.dump(output_data, f)

        print(f"VSSR-MC completed successfully. {len(result.structures)} structures sampled.")
        print(f"Lowest energy: {result.lowest_energy_structure.energy:.3f} eV")

    except Exception as e:
        import traceback

        error_data = {
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc(),
        }
        with open(args.output, "wb") as f:
            pickle.dump(error_data, f)
        print(f"VSSR-MC failed: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
