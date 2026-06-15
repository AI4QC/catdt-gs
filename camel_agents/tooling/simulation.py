
"""Surface/reconstruction/pathway/KMC simulation tool mixin."""

import json
from typing import Any, Dict, List, Optional
from pathlib import Path

from ase.io import read, write

from core.surface.surff_predictor import PredictionResult as SurFFPredictionResult
from core.reconstruction.adsorbdiff_predictor import PredictionOutput as AdsorbDiffPredictionOutput
from core.reconstruction.vssr_mc_predictor import VSSRMCResult
from core.pathway.pathway_predictor import CompletePathwayResult
from core.kmc.catmap_predictor import CatMAPResult

from .common import DEPS_BASE_PATH, logger

try:
    from camel_agents.visualized_digital_twin import VisualizedGasSolidDigitalTwin
except ImportError:
    VisualizedGasSolidDigitalTwin = type(None)


class SimulationToolsMixin:
    def generate_surfaces(
        self,
        bulk_structure_path: str,
        top_n_surfaces: int = 3,
        run_id: str = "default_run",
        workflow_step: str = "01_surfaces"
    ) -> str:
        """
        Generates the most exposed surfaces from a bulk crystal structure using SurFF.

        Args:
            bulk_structure_path (str): Absolute path to the bulk crystal structure file (e.g., POSCAR).
            top_n_surfaces (int): Number of top most exposed surfaces to generate and analyze.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "01_surfaces").

        Returns:
            str: Path to a pickle file containing the SurFFPredictionResult object,
                 and the output directory where individual slab files are saved.
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Generating surfaces for {bulk_structure_path} into {output_dir}")

        dt = self._get_dt_instance(use_visualization=True) # Use visualized DT to ensure viz manager is ready
        surff_result, slabs = dt.generate_surfaces_from_bulk(
            bulk_structure=bulk_structure_path,
            top_n=top_n_surfaces,
            output_dir=str(output_dir)
        )

        # Visualize generated surfaces
        viz_manager = self._get_viz_manager_instance(output_dir)
        surface_files = {}
        surface_energies = {}
        for analysis_info, slab in zip(surff_result.get_top_n(top_n_surfaces), slabs):
            miller_str = ''.join(map(str, analysis_info.miller_index))
            slab_file = output_dir / f"slab_{miller_str}.vasp"
            # Ensure slab file is written for visualization
            if not slab_file.exists():
                try:
                    slab.write(slab_file)
                except Exception as e:
                    logger.warning(f"Could not write slab {miller_str} to {slab_file}: {e}")
                    continue
            surface_files[miller_str] = str(slab_file)
            surface_energies[miller_str] = analysis_info.surface_energy

        if surface_files:
            # Create a dedicated visualization subdirectory
            viz_subdir = output_dir / "visualizations"
            viz_subdir.mkdir(exist_ok=True)
            temp_viz_manager = self._get_viz_manager_instance(viz_subdir) # Use temporary manager for subdir
            
            viz_paths = temp_viz_manager.visualize_generated_surfaces(
                surface_files,
                surface_energies
            )
            logger.info(f"Visualized {len(viz_paths)} surfaces.")
        
        # Save SurFF result object
        result_path = output_dir / "surff_prediction_result.pkl"
        return self._save_result_to_pickle(surff_result, result_path)

    def predict_adsorption_sites(
        self,
        surface_path: str,
        adsorbate_smi: str,
        num_sites: int = 5,
        run_id: str = "default_run",
        workflow_step: str = "02_adsorption",
        llm_review_context: Optional[str] = None
    ) -> str:
        """
        Predicts optimal adsorption sites for a given adsorbate on a surface.

        Args:
            surface_path (str): Absolute path to the surface structure file (e.g., .vasp).
            adsorbate_smi (str): SMILES string or common name for the adsorbate (e.g., "*CO").
            num_sites (int): Number of adsorption sites to sample.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "02_adsorption").
            llm_review_context (str, optional): Context for LLM if reviewing structures (e.g., "CO oxidation on Pt(111)").

        Returns:
            str: Path to a pickle file containing the AdsorbDiffPredictionOutput object.
                 Also saves the best adsorption configuration to a VASP file in the output directory.
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Predicting adsorption sites for {adsorbate_smi} on {surface_path} into {output_dir}")

        dt = self._get_dt_instance(use_visualization=True)

        # Manually manage initial adsorbate structure, if any, for proper surface expansion
        surface_atoms = read(surface_path)

        # AdsorbML needs many more candidates than AdsorbDiff (heuristic
        # enumeration vs diffusion prior), so use the backend-specific count.
        # For adsorbml we take max(LLM's num_sites, configured adsorbml_num_sites).
        effective_num_sites = num_sites
        if getattr(self, "adsorption_backend", "adsorbml") == "adsorbml":
            configured = int(getattr(self, "adsorbml_num_sites", 20))
            effective_num_sites = max(int(num_sites), configured)
            if effective_num_sites != num_sites:
                logger.info(
                    "AdsorbML: overriding num_sites %d → %d (from adsorbml_num_sites config)",
                    num_sites, effective_num_sites,
                )

        # Resolve adsorbate label to a fairchem Adsorbate (typed container).
        # Fairchem's Adsorbate(adsorbate_smiles_from_db=...) does fuzzy lookup
        # and silently returns a RANDOM molecule if the SMILES isn't in its db
        # (e.g. '*H2COOH*' → C2H4O). Our resolver has a proper alias table,
        # SMILES fallbacks, and RDKit 3D generation, so we pre-resolve here
        # and wrap the result in fairchem's Adsorbate type directly — no
        # atoms.info metadata stashing.
        try:
            from core.pathway.adsorbate_resolver import resolve_adsorbate
            from fairchem.data.oc.core import Adsorbate
            ads_atoms, binding_indices = resolve_adsorbate(adsorbate_smi)
            adsorbate_arg: Any = Adsorbate(
                adsorbate_atoms=ads_atoms,
                adsorbate_binding_indices=list(binding_indices) or [0],
            )
            logger.info(
                "Resolved adsorbate '%s' → %s (binding=%s, n=%d)",
                adsorbate_smi, ads_atoms.get_chemical_formula(), binding_indices, len(ads_atoms),
            )
        except Exception as exc:
            logger.warning(
                "adsorbate_resolver failed for '%s' (%s); falling back to fairchem db lookup",
                adsorbate_smi, exc,
            )
            adsorbate_arg = adsorbate_smi

        adsorb_result = dt.predict_adsorption_sites(
            surface=surface_atoms,
            adsorbate=adsorbate_arg,
            num_sites=effective_num_sites,
            output_dir=str(output_dir)
        )

        # Save best result to VASP
        best_config_path = output_dir / "best_adsorption_config.vasp"
        write(best_config_path, adsorb_result.best_result.final_structure)
        logger.info(f"Saved best adsorption configuration to {best_config_path}")

        # Visualize adsorption sites
        viz_subdir = output_dir / "visualizations"
        viz_subdir.mkdir(exist_ok=True)
        viz_manager = self._get_viz_manager_instance(viz_subdir)
        ads_structures = {
            "best": str(best_config_path)
        }
        # AdsorbDiffOutput contains best_result, so we directly visualize it
        try:
            viz_paths = viz_manager.visualize_adsorption_sites(
                ads_structures,
                { "best": adsorb_result.best_result.energy }, # Use actual energy if available
                Path(surface_path).stem # Use filename as identifier
            )
            logger.info(f"Visualized {len(viz_paths)} adsorption configurations.")
        except Exception as e:
            logger.warning(f"Failed to visualize adsorption sites: {e}")

        # LLM review
        if llm_review_context:
            try:
                # LLM review needs the structure and an image
                image_for_llm_review = viz_subdir / "llm_review_temp_img.png"
                viz_manager.structure_visualizer.visualize_structure(
                    adsorb_result.best_result.final_structure,
                    output_file=str(image_for_llm_review),
                    title=f"LLM review: {adsorbate_smi}",
                    show=False,
                )
                llm_review_result = dt._review_structure_with_llm(
                    adsorb_result.best_result.final_structure,
                    intermediate_name=adsorbate_smi,
                    image_path=str(image_for_llm_review),
                    reaction_context=llm_review_context
                )
                logger.info(f"LLM review result for {adsorbate_smi}: {llm_review_result}")
                # Save LLM review result
                with open(output_dir / "llm_review_adsorption.json", "w") as f:
                    json.dump(llm_review_result, f, indent=2, ensure_ascii=False)
            except Exception as e:
                logger.error(f"Error during LLM review of adsorption site: {e}")


        # Save AdsorbDiff result object
        result_path = output_dir / "adsorbdiff_prediction_output.pkl"
        return self._save_result_to_pickle(adsorb_result, result_path)

    def simulate_surface_reconstruction(
        self,
        surface_with_adsorbate_path: str,
        adsorbates_elements_for_mc: List[str], # E.g., ['Pt', 'O'] from the bulk for VSSR-MC to sample
        temperature_k: float, # Temperature in Kelvin for MC
        total_sweeps: int = 100,
        run_id: str = "default_run",
        workflow_step: str = "03_reconstruction",
        clean_slab_path: Optional[str] = None, # Optional: path to the clean slab if known, for virtual site generation
        surface_indices: Optional[List[int]] = None,  # Slab surface atom indices (free for MC)
        adsorbate_indices: Optional[List[int]] = None,  # Pre-existing adsorbate atom indices (free but not "surface")
        num_adsorbates_override: Optional[int] = None,  # Force canonical MC target count (e.g. N=top-2-layer atoms for alloy reconstruction)
        chem_pots: Optional[Dict[str, float]] = None,  # Per-element chemical potentials (eV). Default 0 for all.
        use_seed_for_virtual_sites: bool = False,  # If True, skip passing clean_slab → sites on seed (with overlayer)
        adsorbate_exclusion_radius_A: Optional[float] = None,
        existing_atom_exclusion_radius_A: Optional[float] = None,
    ) -> str:
        """
        Simulates surface reconstruction using VSSR-MC for a given surface and adsorbate elements.

        Args:
            surface_with_adsorbate_path (str): Absolute path to the surface structure file (e.g., .vasp)
                                            which may already contain adsorbates.
            adsorbates_elements_for_mc (List[str]): List of element symbols (e.g., ['Pt', 'O'])
                                                    that VSSR-MC can add/remove to simulate reconstruction.
                                                    These are typically the elements in the bulk/surface.
            temperature_k (float): Temperature in Kelvin for the MC simulation.
            total_sweeps (int): Number of Monte Carlo sweeps to perform.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "03_reconstruction").
            clean_slab_path (str, optional): Absolute path to a clean slab (without adsorbates).
                                             If provided, virtual sites for VSSR-MC will be generated on this.
                                             This helps prevent placing virtual sites above existing adsorbates.

        Returns:
            str: Path to a pickle file containing the VSSRMCResult object.
                 Also saves the lowest energy reconstructed structure to a VASP file.
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Simulating surface reconstruction for {surface_with_adsorbate_path} into {output_dir}")

        dt = self._get_dt_instance(use_visualization=True)

        surface_atoms = read(surface_with_adsorbate_path)
        clean_slab_atoms = read(clean_slab_path) if clean_slab_path else None
        
        mc_result = dt.simulate_surface_reconstruction(
            surface=surface_atoms,
            adsorbates=adsorbates_elements_for_mc,
            temperature=temperature_k,
            total_sweeps=total_sweeps,
            output_dir=str(output_dir),
            surface_indices=surface_indices,
            adsorbate_indices=adsorbate_indices,
            clean_slab=clean_slab_atoms,
            num_adsorbates_override=num_adsorbates_override,
            chem_pots=chem_pots,
            use_seed_for_virtual_sites=use_seed_for_virtual_sites,
            adsorbate_exclusion_radius_A=adsorbate_exclusion_radius_A,
            existing_atom_exclusion_radius_A=existing_atom_exclusion_radius_A,
        )

        # Save lowest energy reconstructed structure to VASP
        lowest_energy_structure_path = output_dir / "lowest_energy_reconstructed_structure.vasp"
        write(lowest_energy_structure_path, mc_result.lowest_energy_structure.atoms)
        logger.info(f"Saved lowest energy reconstructed structure to {lowest_energy_structure_path}")

        # Visualize MC trajectory (if VisualizedGasSolidDigitalTwin is used)
        if isinstance(dt, VisualizedGasSolidDigitalTwin) and dt.enable_visualization:
            try:
                # The _run_vssr_mc_subprocess creates trajectory files, we need to collect them
                # For simplicity, we just use the raw output of mc_result.structures
                # or point to the trajectory if it's saved.
                # Here, we generate one from the saved structures.
                traj_atoms_list = [s.atoms for s in mc_result.structures]
                
                viz_subdir = output_dir / "visualizations"
                viz_subdir.mkdir(exist_ok=True)
                temp_viz_manager = self._get_viz_manager_instance(viz_subdir)

                # Ensure numpy is imported for linspace
                import numpy as np

                if traj_atoms_list:
                    # Limit frames for reasonable GIF size and generation time
                    max_frames = 50
                    if len(traj_atoms_list) > max_frames:
                        indices = np.linspace(0, len(traj_atoms_list) - 1, max_frames, dtype=int)
                        # Ensure first and last frames are included
                        if 0 not in indices: indices = np.insert(indices, 0, 0)
                        if len(traj_atoms_list) - 1 not in indices: indices = np.append(indices, len(traj_atoms_list) - 1)
                        traj_atoms_list = [traj_atoms_list[i] for i in sorted(list(set(indices)))]
                        
                    mc_gif_path = viz_subdir / f"{Path(surface_with_adsorbate_path).stem}_mc_reconstruction.gif"
                    temp_viz_manager.structure_visualizer.visualize_trajectory(
                        trajectory=traj_atoms_list,
                        output_file=str(mc_gif_path),
                        fps=5,
                        show_progress=False,
                    )
                    logger.info(f"Generated MC trajectory GIF: {mc_gif_path}")
                else:
                    logger.warning("No trajectory atoms to visualize for MC reconstruction.")
            except Exception as e:
                logger.warning(f"Failed to visualize MC trajectory: {e}")
                import traceback
                logger.debug(traceback.format_exc())

        # Save VSSRMCResult object
        result_path = output_dir / "vssr_mc_result.pkl"
        return self._save_result_to_pickle(mc_result, result_path)

    def analyze_reaction_pathway(
        self,
        surface_path: str, # Could be pristine, or with a starting adsorbate
        reaction_intermediates: List[str], # E.g., ["*CO", "*O", "*CO2", "*"]
        initial_adsorbate_present_on_surface_smi: str, # E.g., "*CO" if surface_path has *CO
        calculate_barriers: bool = True,
        num_sites: int = 5,
        run_id: str = "default_run",
        workflow_step: str = "04_pathway_analysis"
    ) -> str:
        """
        Analyzes the reaction pathway by calculating adsorption energies and NEB barriers.

        Args:
            surface_path (str): Absolute path to the reference surface structure file (e.g., lowest energy reconstructed surface).
            reaction_intermediates (List[str]): Ordered list of reaction intermediates (e.g., ["*CO", "*O", "*CO2", "*"]).
            initial_adsorbate_present_on_surface_smi (str): The adsorbate that is already present on the 'surface_path' structure.
                                                          This helps correctly identify surface atoms vs. adsorbate atoms.
                                                          Use "" or None if surface_path is a pristine surface.
            calculate_barriers (bool): Whether to perform NEB calculations for barriers.
            num_sites (int): Number of potential adsorption sites to consider for each intermediate.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "04_pathway_analysis").

        Returns:
            str: Path to a pickle file containing the CompletePathwayResult object.
                 Also saves the reaction pathway visualization (GIF) and energy diagram (PNG).
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Analyzing reaction pathway for intermediates {reaction_intermediates} on {surface_path} into {output_dir}")

        dt = self._get_dt_instance(use_visualization=True)

        # Load initial surface atoms to pass to pathway_predictor to help identify surface atoms
        initial_surface_atoms = read(surface_path)
        
        # Determine current adsorbate indices correctly
        # This logic needs to correctly identify atoms belonging to the initial adsorbate from the surface_path
        current_adsorbate_indices = []
        if initial_adsorbate_present_on_surface_smi and initial_adsorbate_present_on_surface_smi != "*":
            # 🔧 实现吸附物识别：找到最上层的非Pt原子
            import numpy as np
            
            symbols = initial_surface_atoms.get_chemical_symbols()
            positions = initial_surface_atoms.positions
            
            # 找到所有非Pt原子
            non_pt_indices = [i for i, sym in enumerate(symbols) if sym != "Pt"]
            
            if non_pt_indices:
                # 按z坐标排序，取最上层的原子作为吸附物
                z_coords = positions[non_pt_indices, 2]
                z_threshold = z_coords.max() - 2.0  # 最上层2Å内的原子
                
                adsorbate_mask = z_coords > z_threshold
                current_adsorbate_indices = [non_pt_indices[i] for i, is_ads in enumerate(adsorbate_mask) if is_ads]
                
                logger.info(f"Identified {len(current_adsorbate_indices)} adsorbate atoms from surface: {initial_adsorbate_present_on_surface_smi}")
                logger.info(f"  Adsorbate atom indices: {current_adsorbate_indices}")
            
        pathway_result = dt.analyze_reaction_pathway(
            surface=initial_surface_atoms, # Pass the Atoms object directly
            intermediates=reaction_intermediates,
            calculate_barriers=calculate_barriers,
            num_sites=num_sites,
            output_dir=str(output_dir),
            current_adsorbate_indices=current_adsorbate_indices
        )

        # Visualize reaction pathway (GIF) and energy diagram (PNG)
        if isinstance(dt, VisualizedGasSolidDigitalTwin) and dt.enable_visualization:
            viz_subdir = output_dir / "visualizations"
            viz_subdir.mkdir(exist_ok=True)
            viz_manager = self._get_viz_manager_instance(viz_subdir)
            
            # --- Pathway GIF ---
            pathway_structures_for_viz = {}
            pathway_adsorption_dir = output_dir / "adsorption" # Pathway predictor saves individual steps here
            
            # Collect structures in the correct order for GIF
            for inter in reaction_intermediates:
                inter_file = pathway_adsorption_dir / f"{inter}_relaxed.vasp"
                if inter_file.exists():
                    pathway_structures_for_viz[inter] = read(str(inter_file)) # Load Atoms object
            
            if len(pathway_structures_for_viz) >= 2:
                try:
                    pathway_gif = viz_manager.visualize_reaction_pathway(
                        pathway_structures_for_viz,
                        sequence=reaction_intermediates,
                        surface_name=Path(surface_path).stem,
                        fps=1,
                        hold_frames=5
                    )
                    logger.info(f"Generated reaction pathway GIF: {pathway_gif}")
                except Exception as e:
                    logger.warning(f"Failed to visualize reaction pathway GIF: {e}")

            # --- Energy Diagram PNG ---
            if calculate_barriers and pathway_result:
                try:
                    # Extract energy profile as done in visualized_digital_twin.py
                    energies = []
                    labels = []
                    is_ts_list = []

                    # The pathway_result's energy_profile attribute directly provides what's needed for visualization
                    # It's a list of tuples: (label, energy, is_transition_state)
                    for label, energy, is_ts in pathway_result.energy_profile:
                        labels.append(label)
                        energies.append(energy)
                        is_ts_list.append(is_ts)
                    
                    energy_diagram_path = viz_manager.visualize_energy_diagram(
                        energies,
                        labels,
                        is_ts_list,
                        filename="energy_diagram.png",
                        title=f"Reaction Energy Profile on {Path(surface_path).stem}"
                    )
                    logger.info(f"Generated energy diagram: {energy_diagram_path}")
                except Exception as e:
                    logger.warning(f"Failed to visualize energy diagram: {e}")

        # Save CompletePathwayResult object
        result_path = output_dir / "complete_pathway_result.pkl"
        return self._save_result_to_pickle(pathway_result, result_path)

    def run_kmc_simulation(
        self,
        pathway_result_path: str,
        temperature_k: float,
        pressures: Dict[str, float], # E.g., {"CO": 1.0, "O2": 0.2}
        run_id: str = "default_run",
        workflow_step: str = "05_kmc_simulation"
    ) -> str:
        """
        Runs Kinetic Monte Carlo (KMC) simulation based on a previously analyzed reaction pathway.

        Args:
            pathway_result_path (str): Absolute path to a pickle file containing the CompletePathwayResult object.
            temperature_k (float): Temperature in Kelvin for the KMC simulation.
            pressures (Dict[str, float]): Dictionary of gas phase pressures in bar (e.g., {"CO": 1.0, "O2": 0.2}).
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "05_kmc_simulation").

        Returns:
            str: Path to a pickle file containing the CatMAPResult object.
                 Also saves the KMC dynamics visualization (GIF) if available.
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Running KMC simulation for pathway {pathway_result_path} at {temperature_k}K, pressures {pressures} into {output_dir}")

        dt = self._get_dt_instance(use_visualization=True)
        pathway_result: CompletePathwayResult = self._load_result_from_pickle(pathway_result_path)

        kmc_result = dt.run_kmc_simulation(
            pathway_result=pathway_result,
            temperature=temperature_k,
            pressures=pressures,
            output_dir=str(output_dir)
        )

        # Visualize KMC dynamics (GIF)
        if isinstance(dt, VisualizedGasSolidDigitalTwin) and dt.enable_visualization and kmc_result:
            try:
                viz_subdir = output_dir / "visualizations"
                viz_subdir.mkdir(exist_ok=True)
                viz_manager = self._get_viz_manager_instance(viz_subdir)
                
                # We need the structures for each intermediate to create the GIF
                # Fetching the base surface from the pathway_result's initial_atoms if available,
                # otherwise try from the reconstructed path.
                base_surface_atoms = None
                first_step_initial_atoms = None
                if pathway_result.steps:
                    first_step_initial_atoms = getattr(pathway_result.steps[0], "initial_atoms", None)
                candidate_base_surfaces = [
                    Path(pathway_result_path).parent.parent / "03_reconstruction" / "lowest_energy_reconstructed_structure.vasp",
                    Path(getattr(pathway_result, "output_dir", "")) / "04_adsorption_energies" / "surface_relaxed.vasp",
                    Path(pathway_result_path).parent / "04_adsorption_energies" / "surface_relaxed.vasp",
                    Path(pathway_result_path).parent / "adsorption" / "surface_relaxed.vasp",
                ]
                if first_step_initial_atoms is not None:
                    base_surface_atoms = first_step_initial_atoms
                else:
                    for candidate in candidate_base_surfaces:
                        if candidate.exists():
                            base_surface_atoms = read(str(candidate))
                            break

                if base_surface_atoms is None:
                    logger.warning("Could not obtain a base surface structure for KMC visualization.")
                    # Fallback: if no base surface, cannot create dynamic visualization correctly
                    raise ValueError("Base surface structure required for KMC dynamic visualization.")

                # Get all intermediate structures from pathway analysis output
                # This needs to include the base_surface_atoms as the "clean" state if '*' is a state
                intermediate_structures = {}
                adsorption_candidates = [
                    Path(getattr(pathway_result, "output_dir", "")) / "04_adsorption_energies",
                    Path(pathway_result_path).parent / "04_adsorption_energies",
                    Path(pathway_result_path).parent / "adsorption",
                ]
                pathway_adsorption_dir = None
                for candidate in adsorption_candidates:
                    if candidate.exists():
                        pathway_adsorption_dir = candidate
                        break
                for inter in pathway_result.adsorbate_energies.keys():
                    if inter == "*":
                        # For the clean surface state, use the base_surface_atoms
                        intermediate_structures["*"] = base_surface_atoms.copy()
                    else:
                        inter_file = pathway_adsorption_dir / f"{inter}_relaxed.vasp" if pathway_adsorption_dir is not None else None
                        if inter_file is not None and inter_file.exists():
                            intermediate_structures[inter] = read(str(inter_file))
                
                # Ensure numpy is imported for linspace
                import numpy as np

                if hasattr(kmc_result, 'state_history') and hasattr(kmc_result, 'time_history'):
                    state_history = kmc_result.state_history
                    time_history = kmc_result.time_history

                    max_frames = 100
                    total_steps = len(state_history)
                    
                    if total_steps == 0:
                        logger.warning("KMC state history is empty, cannot visualize dynamics.")
                        return self._save_result_to_pickle(kmc_result, output_dir / "catmap_result.pkl")

                    # Select frames, ensuring first and last are always present
                    indices = np.linspace(0, total_steps - 1, min(total_steps, max_frames), dtype=int)
                    
                    trajectory = []
                    titles = []
                    for i, idx in enumerate(indices):
                        state = state_history[idx]
                        time = time_history[idx]
                        
                        if state in intermediate_structures:
                            # Need to copy the Atoms object to avoid modifying it during visualization
                            atoms_to_add = intermediate_structures[state].copy()
                            # Optional: Add title to atoms object for visualization_manager
                            atoms_to_add.info['label'] = f"State: {state}\nTime: {time:.4e} s"
                            trajectory.append(atoms_to_add)
                        else:
                            logger.warning(f"KMC state '{state}' not found in intermediate_structures for visualization.")

                    if len(trajectory) > 0:
                        kmc_gif_path = viz_subdir / "kmc_surface_dynamics.gif"
                        viz_manager.structure_visualizer.visualize_trajectory(
                            trajectory,
                            output_file=str(kmc_gif_path),
                            fps=5,
                            show_progress=False # Disable verbose progress for tool output
                        )
                        logger.info(f"Generated KMC dynamics GIF: {kmc_gif_path}")
                    else:
                        logger.warning(f"No valid trajectory frames for KMC visualization.")
                else:
                    logger.warning("KMC result missing state_history or time_history for visualization.")
            except Exception as e:
                logger.warning(f"Failed to visualize KMC trajectory: {e}")
                import traceback
                logger.debug(traceback.format_exc())

        # Save CatMAPResult object
        result_path = output_dir / "catmap_result.pkl"
        return self._save_result_to_pickle(kmc_result, result_path)

    def compute_adsorption_energies(
        self,
        surface_path: str,
        intermediates: List[str],
        output_dir: str,
        num_sites: int = 5,
        fmax: float = 0.05,
        max_steps: int = 300,
        fixed_site: Optional[List[float]] = None,
        skip_surface_relaxation: bool = False,
    ) -> Dict[str, Any]:
        predictor = self.get_shared_fairchem_predictor(
            cache_key="uma_shared",
            model_name="uma-s-1p1",
            use_gpu=True,
            device="cuda",
            work_subdir="_shared_fairchem_global",
            keep_files=False,
            verbose=False,
        )
        result = predictor.predict_pathway_energies(
            surface=surface_path,
            adsorbates=intermediates,
            num_sites=num_sites,
            fmax=fmax,
            max_steps=max_steps,
            output_dir=output_dir,
            fixed_adsorption_site=fixed_site,
            skip_surface_relaxation=skip_surface_relaxation,
        )
        return {
            "adsorbate_energies": result.adsorbate_energies,
            "adsorption_energies": result.adsorption_energies,
            "best_configurations": result.best_configurations,
            "gas_frequencies": dict(getattr(result, "gas_frequencies", {}) or {}),
        }

    def relax_adsorbate_on_surface(
        self,
        structure_path: str,
        adsorbate_indices: List[int],
        output_path: str,
        fmax: float = 0.05,
        max_steps: int = 200,
        free_surface_indices: Optional[List[int]] = None,
    ) -> Dict[str, Any]:
        """Relax adsorbate atoms on surface. Returns energy, path, and force info.

        By default only ``adsorbate_indices`` are free (all surface atoms frozen).
        Pass ``free_surface_indices`` (e.g. the top 2 layers from
        ``surface.top_layer_indices``) to additionally free those surface atoms —
        this is the Phase-2 default and gives physically reasonable geometries
        for Stage-A intermediate construction and Stage-B staged-atom relaxation.
        """
        from ase.io import read as ase_read, write as ase_write
        from ase.optimize import BFGS
        from ase.constraints import FixAtoms

        predictor = self.get_shared_fairchem_predictor(
            cache_key="uma_shared", model_name="uma-s-1p1",
            use_gpu=True, device="cuda",
            work_subdir="_shared_fairchem_global",
            keep_files=False, verbose=False,
        )
        predictor._load_model()

        atoms = ase_read(structure_path)
        atoms.pbc = [True, True, True]
        atoms.calc = predictor._calculator

        free_set = set(int(i) for i in (adsorbate_indices or []))
        if free_surface_indices:
            free_set.update(int(i) for i in free_surface_indices)
        fixed = [i for i in range(len(atoms)) if i not in free_set]
        atoms.set_constraint(FixAtoms(indices=fixed))

        opt = BFGS(atoms, logfile=output_path.replace('.vasp', '_relax.log'))
        opt.run(fmax=fmax, steps=max_steps)

        energy = float(atoms.get_potential_energy())
        import numpy as np
        fmax_final = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))

        try:
            atoms.set_constraint()
        except Exception:
            pass
        ase_write(output_path, atoms)

        return {
            "energy": energy,
            "fmax": fmax_final,
            "output_path": output_path,
            "n_steps": int(getattr(opt, "nsteps", 0)),
        }
