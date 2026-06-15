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
sys.path.insert(0, os.path.join(_this_dir, ".."))  # core/
sys.path.insert(0, _this_dir)  # core/reconstruction/


def patch_chgnet_cutoff():
    """Monkey-patch CHGNet's atom_graph_cutoff to 10.0 to avoid isolated atom issues"""
    import nff.io.chgnet as chgnet_io
    from chgnet.graph import CrystalGraphConverter
    from chgnet.data.dataset import StructureData

    # Store original function
    original_convert_data_batch = chgnet_io.convert_data_batch

    def patched_convert_data_batch(data_batch, cutoff=5.0, shuffle=True):
        """Patched version that uses higher atom_graph_cutoff"""
        from nff.io import AtomsBatch
        from nff.utils.cuda import batch_detach, detach
        from pymatgen.io.ase import AseAtomsAdaptor
        import torch

        detached_batch = batch_detach(data_batch)
        nxyz = detached_batch["nxyz"]
        atoms_batch = AtomsBatch(
            nxyz[:, 0].long(),
            props=detached_batch,
            positions=nxyz[:, 1:],
            cell=(detached_batch["lattice"][0] if "lattice" in detached_batch else None),
            pbc="lattice" in detached_batch,
            cutoff=cutoff,
            dense_nbrs=False,
        )
        atoms_list = atoms_batch.get_list_atoms()

        pymatgen_structures = [AseAtomsAdaptor.get_structure(ab) for ab in atoms_list]

        energies = torch.atleast_1d(data_batch.get("energy"))
        if energies is not None and len(energies) > 0:
            energies_per_atom = energies
        else:
            energies_per_atom = torch.Tensor([0.0] * len(pymatgen_structures))

        energy_grads = data_batch.get("energy_grad")
        num_atoms = detach(data_batch["num_atoms"]).tolist()
        stresses = data_batch.get("stress")
        magmoms = data_batch.get("magmoms")

        if energy_grads is not None and len(energy_grads) > 0:
            forces = [-x for x in energy_grads] if isinstance(energy_grads, list) else -energy_grads
        else:
            forces = None

        if forces is not None and len(forces) > 0:
            forces = torch.split(torch.atleast_2d(forces), num_atoms)
        else:
            forces = [torch.zeros_like(torch.Tensor(ab.get_positions())) for ab in atoms_list]

        if stresses is not None and len(stresses) > 0:
            stresses = torch.split(torch.atleast_2d(stresses), num_atoms)
        if magmoms is not None and len(magmoms) > 0:
            magmoms = torch.split(torch.atleast_2d(magmoms), num_atoms)

        # Defaults (same as CHGNet 0.3.0). The real fix for isolated atoms
        # (e.g. evaporated adatoms after MD) is to drop them BEFORE VSSR-MC
        # in scripts/pdh_smsi/run_vssr_mc.py — bumping cutoff here slows
        # line_graph_adjacency_list enough to stall a 76-atom slab.
        graph_converter = CrystalGraphConverter(
            atom_graph_cutoff=6.0,
            bond_graph_cutoff=3.0,
        )

        return StructureData(
            structures=pymatgen_structures,
            energies=energies_per_atom,
            forces=forces,
            stresses=stresses,
            magmoms=magmoms,
            shuffle=shuffle,
            graph_converter=graph_converter,  # Pass custom graph_converter
        )

    # Apply patch
    chgnet_io.convert_data_batch = patched_convert_data_batch
    print("Patched CHGNet atom_graph_cutoff to 10.0")


def main():
    parser = argparse.ArgumentParser(description="Run VSSR-MC in subprocess")
    parser.add_argument("--surface", required=True, help="Path to pickled surface structure")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    parser.add_argument("--output", required=True, help="Path to output pickle file")
    args = parser.parse_args()

    # Apply CHGNet graph-cutoff patch before any CHGNet/NFF import triggered
    # by VSSRMCPredictor construction. Without this, larger-lattice slabs
    # (e.g. Cu fcc) can yield isolated-atom graphs at the default 6 Å cutoff,
    # which send StructureData.__getitem__ into a random-retry RecursionError.
    try:
        patch_chgnet_cutoff()
    except Exception as exc:
        print(f"WARN: CHGNet cutoff patch failed ({exc}); proceeding with defaults")

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
    from vssr_mc_predictor import VSSRMCPredictor

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
