"""CAMEL tool registry facade for CatDT."""

from typing import Any, Dict, Optional

from ase import Atoms

from camel.toolkits import FunctionTool

from camel_agents.tooling import (
    CATDT_CORE_PATH,
    DEPS_BASE_PATH,
    logger,
    CatDTToolRuntimeBase,
    Agent45MementoToolsMixin,
    SimulationToolsMixin,
    ReportingToolsMixin,
    NEBToolsMixin,
    Agent45WorkflowToolsMixin,
    Agent45GeometryTools,
    MechanismToolsMixin,
)


class CatDTTools(
    CatDTToolRuntimeBase,
    Agent45MementoToolsMixin,
    SimulationToolsMixin,
    ReportingToolsMixin,
    NEBToolsMixin,
    Agent45WorkflowToolsMixin,
    MechanismToolsMixin,
):
    """Unified CatDT tool facade; implementations are split across tooling modules."""

    def __init__(
        self,
        output_base_dir: str = "output/catdt_workflow",
        mc_energy_model: str = "CHGNet",
        potential_she=None,
        ph=None,
        mp_api_key=None,
        adsorption_backend: str = "adsorbml",
        adsorbml_num_sites: int = 20,
        adsorbml_placement_mode: str = "random_site_heuristic_placement",
        adsorbml_interstitial_gap: float = 0.1,
        adsorbml_relax_steps: int = 200,
        adsorbml_relax_fmax: float = 0.05,
    ):
        super().__init__(
            output_base_dir=output_base_dir,
            mc_energy_model=mc_energy_model,
            potential_she=potential_she,
            ph=ph,
            mp_api_key=mp_api_key,
            adsorption_backend=adsorption_backend,
            adsorbml_num_sites=adsorbml_num_sites,
            adsorbml_placement_mode=adsorbml_placement_mode,
            adsorbml_interstitial_gap=adsorbml_interstitial_gap,
            adsorbml_relax_steps=adsorbml_relax_steps,
            adsorbml_relax_fmax=adsorbml_relax_fmax,
        )
        self._workflow_bridge = None
        self._current_surface = None  # Set by CatDTCamelWorkflow.run_on_surface()
        self._init_workflow_bridge()

    def _get_surface(self):
        """Return the CatalystSurface currently owned by the parent workflow, if any."""
        return self._current_surface

    def _init_workflow_bridge(self) -> None:
        """Expose workflow helper methods via tools (for orchestration delegation)."""
        try:
            from camel_agents.tooling.species import WorkflowContextMixin
            from camel_agents.tooling.geometry import WorkflowGeometryMixin

            tools_self = self

            class _WorkflowToolBridge(WorkflowContextMixin, WorkflowGeometryMixin):
                def _get_surface(self_bridge):
                    return tools_self._current_surface

            bridge = _WorkflowToolBridge()
            bridge.tools = self
            bridge._pathway_energy_helper = None
            self._workflow_bridge = bridge
        except Exception as exc:
            logger.warning("Workflow helper bridge init failed: %s", exc)
            self._workflow_bridge = None

    def __getattr__(self, name: str) -> Any:
        bridge = self.__dict__.get("_workflow_bridge")
        if bridge is not None:
            try:
                return getattr(bridge, name)
            except AttributeError:
                pass
        raise AttributeError(f"{self.__class__.__name__} has no attribute {name}")

    def generate_surfaces_tool(
        self,
        bulk_structure_path: str = "",
        top_n_surfaces: int = 3,
        run_id: str = "default_run",
        workflow_step: str = "01_surfaces",
    ) -> Any:
        if not str(bulk_structure_path or "").strip():
            return {
                "status": "FAIL",
                "error": "missing_required_arg: bulk_structure_path",
                "required": ["bulk_structure_path", "top_n_surfaces", "run_id", "workflow_step"],
            }
        return self.generate_surfaces(
            bulk_structure_path=bulk_structure_path,
            top_n_surfaces=top_n_surfaces,
            run_id=run_id,
            workflow_step=workflow_step,
        )

    def predict_adsorption_sites_tool(
        self,
        surface_path: str = "",
        adsorbate_smi: str = "",
        num_sites: int = 5,
        run_id: str = "default_run",
        workflow_step: str = "02_adsorption",
        llm_review_context: Any = None,
    ) -> Any:
        if not str(surface_path or "").strip() or not str(adsorbate_smi or "").strip():
            return {
                "status": "FAIL",
                "error": "missing_required_args: surface_path/adsorbate_smi",
                "required": ["surface_path", "adsorbate_smi", "num_sites", "run_id", "workflow_step"],
            }
        return self.predict_adsorption_sites(
            surface_path=surface_path,
            adsorbate_smi=adsorbate_smi,
            num_sites=num_sites,
            run_id=run_id,
            workflow_step=workflow_step,
            llm_review_context=llm_review_context,
        )

    def simulate_surface_reconstruction_tool(
        self,
        surface_with_adsorbate_path: str = "",
        adsorbates_elements_for_mc: list[str] | None = None,
        temperature_k: float = 500.0,
        total_sweeps: int = 100,
        run_id: str = "default_run",
        workflow_step: str = "03_reconstruction",
        clean_slab_path: str | None = None,
    ) -> Any:
        if not str(surface_with_adsorbate_path or "").strip():
            return {
                "status": "FAIL",
                "error": "missing_required_arg: surface_with_adsorbate_path",
                "required": [
                    "surface_with_adsorbate_path",
                    "adsorbates_elements_for_mc",
                    "temperature_k",
                    "total_sweeps",
                    "run_id",
                    "workflow_step",
                ],
            }
        if not adsorbates_elements_for_mc:
            return {
                "status": "FAIL",
                "error": "missing_required_arg: adsorbates_elements_for_mc",
                "required": [
                    "surface_with_adsorbate_path",
                    "adsorbates_elements_for_mc",
                    "temperature_k",
                    "total_sweeps",
                    "run_id",
                    "workflow_step",
                ],
            }
        return self.simulate_surface_reconstruction(
            surface_with_adsorbate_path=surface_with_adsorbate_path,
            adsorbates_elements_for_mc=adsorbates_elements_for_mc,
            temperature_k=temperature_k,
            total_sweeps=total_sweeps,
            run_id=run_id,
            workflow_step=workflow_step,
            clean_slab_path=clean_slab_path,
        )

    def build_steps_payload_tool(self, step_structures: list[dict] | None = None) -> str:
        """LLM-facing wrapper for step payload generation."""
        if not step_structures:
            return "(no step_structures provided)"
        return self.build_steps_payload(workflow=self._workflow_bridge or self, step_structures=step_structures)

    def programmatic_validation_tool(self, step_structures: list[dict] | None = None) -> dict:
        """LLM-facing wrapper for generic geometric validation."""
        if not step_structures:
            return {
                "status": "FAIL",
                "fatal_issues": ["No step_structures provided for validation."],
                "warning_issues": [],
                "summary": "(missing input)",
                "feedback": "Provide canonical step_structures generated by tooling.",
            }
        return self.programmatic_validation(workflow=self._workflow_bridge or self, step_structures=step_structures)

    def get_step_element_deltas_tool(self) -> list[dict]:
        """Query per-step element deltas: what atoms must be added to which side for NEB.

        Returns a list of steps, each with: step_name, to_add_to_reactant,
        to_add_to_product, total_staged.  Call this FIRST to plan your work.
        """
        return self.get_step_element_deltas(workflow=self._workflow_bridge or self)

    def suggest_staged_positions_tool(
        self,
        requests: list[dict] | None = None,
    ) -> dict:
        """Suggest positions for ALL staged atoms across ALL steps in one call.

        Pass a list of requests, each: {"step_name": "...", "element": "H", "side": "reactant"}.
        Returns suggested [x,y,z] positions near chemically relevant anchor atoms.
        Only call this if you need guidance — you can also propose positions yourself.
        """
        if not requests:
            return {"suggestions": [], "error": "requests list is required."}
        return self.suggest_staged_positions(
            workflow=self._workflow_bridge or self, requests=requests,
        )

    def validate_proposed_positions_tool(
        self,
        steps: list[dict] | None = None,
    ) -> dict:
        """Validate ALL proposed atom positions across ALL steps in one call.

        Pass a list of steps, each: {"step_name": "...", "atoms": [{"species": "H", "side": "reactant", "position": [x,y,z]}]}.
        Returns per-step and per-atom diagnostics: distance checks, overlap checks, and suggestions.
        """
        if not steps:
            return {"ok": False, "per_step": [], "suggestion": "steps list is required."}
        return self.validate_proposed_positions(
            workflow=self._workflow_bridge or self, steps=steps,
        )

    def run_agent45_energy_gate_tool(
        self,
        step_structures: list[dict],
        output_dir: str = "",
        relax_if_needed: bool = True,
        force_warning_threshold: float = 2.5,
        force_fatal_threshold: float = 4.0,
        min_pair_fatal_threshold: float = 0.8,
        max_relax_steps: int = 40,
        relax_fmax: float = 0.25,
    ) -> dict:
        """UMA endpoint audit + constrained staged-atom relaxation for Agent4/5 loop."""
        if step_structures:
            sample = step_structures[0]
            if not isinstance(sample, dict) or not isinstance(sample.get("reactant"), Atoms) or not isinstance(sample.get("product"), Atoms):
                return {
                    "status": "FAIL",
                    "summary": "Invalid payload for run_agent45_energy_gate: expected canonical step_structures with ASE Atoms endpoints.",
                    "fatal_issues": [
                        "run_agent45_energy_gate requires step_structures generated after tool-side endpoint construction (reactant/product must be ASE Atoms).",
                    ],
                    "warning_issues": [],
                    "steps": [],
                    "settings": {
                        "relax_if_needed": bool(relax_if_needed),
                        "force_warning_threshold": float(force_warning_threshold),
                        "force_fatal_threshold": float(force_fatal_threshold),
                        "min_pair_fatal_threshold": float(min_pair_fatal_threshold),
                        "max_relax_steps": int(max_relax_steps),
                        "relax_fmax": float(relax_fmax),
                    },
                }
        return self.run_agent45_energy_gate(
            workflow=self._workflow_bridge or self,
            step_structures=step_structures,
            output_dir=output_dir,
            relax_if_needed=relax_if_needed,
            force_warning_threshold=force_warning_threshold,
            force_fatal_threshold=force_fatal_threshold,
            min_pair_fatal_threshold=min_pair_fatal_threshold,
            max_relax_steps=max_relax_steps,
            relax_fmax=relax_fmax,
        )

    def retrieve_agent45_memento_cases_tool(
        self,
        query_text: str = "",
        top_k: int = 6,
        include_negative: bool = True,
        min_score: float = 0.05,
    ) -> dict:
        if not str(query_text or "").strip():
            return {
                "query_text": "",
                "cases": [],
                "prompt_block": "(no query_text provided)",
                "source": "agent45_memento_casebank",
            }
        return self.retrieve_agent45_memento_cases(
            query_text=query_text,
            top_k=top_k,
            include_negative=include_negative,
            min_score=min_score,
        )

    def record_agent45_memento_case_tool(
        self,
        run_id: str,
        reaction_description: str,
        reaction_type: str,
        intermediates: list[str],
        transition_signature: list[str],
        iteration: int,
        validation_status: str,
        issues: list[str],
        feedback: str,
        design_outline: list[dict],
        energy_gate_report: dict | None = None,
        neb_summary: dict | None = None,
        reward: float | None = None,
        extra_query: str = "",
        step_geometry_snapshot: list[dict] | None = None,
        surface_snapshot: dict | None = None,
    ) -> str:
        return self.record_agent45_memento_case(
            run_id=run_id,
            reaction_description=reaction_description,
            reaction_type=reaction_type,
            intermediates=intermediates,
            transition_signature=transition_signature,
            iteration=iteration,
            validation_status=validation_status,
            issues=issues,
            feedback=feedback,
            design_outline=design_outline,
            energy_gate_report=energy_gate_report,
            neb_summary=neb_summary,
            reward=reward,
            extra_query=extra_query,
            step_geometry_snapshot=step_geometry_snapshot,
            surface_snapshot=surface_snapshot,
        )

    def retrieve_agent45_knowledge_items_tool(
        self,
        query_text: str = "",
        top_k: int = 4,
        min_score: float = 0.03,
    ) -> dict:
        if not str(query_text or "").strip():
            return {
                "query_text": "",
                "items": [],
                "prompt_block": "(no query_text provided)",
                "source": "agent45_knowledge_bank",
            }
        return self.retrieve_agent45_knowledge_items(
            query_text=query_text,
            top_k=top_k,
            min_score=min_score,
        )

    def retrieve_agent45_skill_items_tool(
        self,
        query_text: str = "",
        top_k: int = 4,
        min_score: float = 0.03,
    ) -> dict:
        if not str(query_text or "").strip():
            return {
                "query_text": "",
                "items": [],
                "prompt_block": "(no query_text provided)",
                "source": "agent45_skill_bank",
            }
        return self.retrieve_agent45_skill_items(
            query_text=query_text,
            top_k=top_k,
            min_score=min_score,
        )

    def train_agent45_parametric_retriever_tool(
        self,
        min_cases: int = 30,
        max_cases: int = 1500,
        max_pairs: int = 20000,
        max_pos_per_query: int = 4,
        max_neg_per_query: int = 8,
        epochs: int = 1,
        batch_size: int = 16,
        learning_rate: float = 2e-5,
        max_len: int = 256,
        val_ratio: float = 0.15,
        seed: int = 42,
        model_name: str = "",
        device: str = "",
        resume_from_checkpoint: bool = False,
    ) -> dict:
        return self.train_agent45_parametric_retriever(
            min_cases=min_cases,
            max_cases=max_cases,
            max_pairs=max_pairs,
            max_pos_per_query=max_pos_per_query,
            max_neg_per_query=max_neg_per_query,
            epochs=epochs,
            batch_size=batch_size,
            learning_rate=learning_rate,
            max_len=max_len,
            val_ratio=val_ratio,
            seed=seed,
            model_name=model_name,
            device=device,
            resume_from_checkpoint=resume_from_checkpoint,
        )

    def get_agent45_parametric_retriever_status_tool(self) -> dict:
        return self.get_agent45_parametric_retriever_status()

    # ---- Mechanism search tool wrappers ----

    def initialize_mechanism_context_tool(
        self,
        initial_state: Optional[str] = "",
        target_state: Optional[str] = "",
        surface_id: Optional[str] = "",
        surface_facet: Optional[str] = "",
        bulk_formula: Optional[str] = "",
        surface_structure_path: Optional[str] = "",
        clean_slab_path: Optional[str] = "",
        environment: dict | None = None,
        known_intermediates: list[str] | None = None,
        constraints: Optional[str] = "",
        fixed_adsorption_site: list[float] | None = None,
    ) -> dict:
        """Initialize mechanism search context from reaction specification.

        Accepts None for any string argument (the LLM frequently emits JSON
        null for optional fields); they are coerced to "" for downstream use.

        Args:
            fixed_adsorption_site: Optional [x,y,z] site from Agent2/3.
                If set, all intermediates are evaluated at this fixed site,
                skipping the expensive adsorption site search.
        """
        # Normalize all string-typed params: the OpenAI function schema marks them
        # as string-or-null, so the LLM may legitimately send null.
        initial_state = initial_state or ""
        target_state = target_state or ""
        surface_id = surface_id or ""
        surface_facet = surface_facet or ""
        bulk_formula = bulk_formula or ""
        surface_structure_path = surface_structure_path or ""
        clean_slab_path = clean_slab_path or ""
        constraints = constraints or ""

        if not initial_state or not target_state:
            return {"error": "initial_state and target_state are required"}
        return self.initialize_mechanism_context(
            initial_state=initial_state,
            target_state=target_state,
            surface_id=surface_id,
            surface_facet=surface_facet,
            bulk_formula=bulk_formula,
            surface_structure_path=surface_structure_path,
            clean_slab_path=clean_slab_path,
            environment=environment,
            known_intermediates=known_intermediates,
            constraints=constraints,
            fixed_adsorption_site=fixed_adsorption_site,
        )

    def generate_candidate_steps_tool(
        self,
        use_care: bool = True,
        use_generic_ops: bool = True,
        electro: bool = False,
    ) -> dict:
        """Generate candidate elementary steps from enabled backends."""
        return self.generate_candidate_steps(
            use_care=use_care,
            use_generic_ops=use_generic_ops,
            electro=electro,
        )

    def run_mechanism_search_tool(
        self,
        max_depth: int = 8,
        beam_width: int = 8,
        delta_keep: float = 0.10,
        delta_prune: float = 0.30,
        max_state_evaluations: int = 40,
        energy_backend: str = "thermal_uma",
        exploration_mode: str = "agent_guided",
        systematic_strategy: str = "pruning",
        mcts_iterations: int = 100,
        mcts_c: float = 1.4,
        mcts_rollout_depth: int = 4,
        mcts_rollout_bias_beta: float = 1.0,
        mcts_expansion_bias_beta: float = 0.5,
        mcts_distance_penalty_eV: float = 0.3,
    ) -> dict:
        """Run mechanism search using selected exploration module(s).

        Args:
            exploration_mode: 'agent_guided', 'systematic', 'both_sequential', 'both_parallel'
            systematic_strategy: Phase-3 algorithm. One of:
                'quick'   — score Stage-1 paths by max-intermediate UMA
                            energy (fast, no graph re-traversal).
                'pruning' — layer-by-layer beam BFS on CRN + sibling
                            pruning + distance-bucketed round-robin beam.
                'mcts'    — UCB1 tree search; shares intermediate cache.
        """
        return self.run_mechanism_search(
            max_depth=max_depth,
            beam_width=beam_width,
            delta_keep=delta_keep,
            delta_prune=delta_prune,
            max_state_evaluations=max_state_evaluations,
            energy_backend=energy_backend,
            exploration_mode=exploration_mode,
            systematic_strategy=systematic_strategy,
            mcts_iterations=mcts_iterations,
            mcts_c=mcts_c,
            mcts_rollout_depth=mcts_rollout_depth,
            mcts_rollout_bias_beta=mcts_rollout_bias_beta,
            mcts_expansion_bias_beta=mcts_expansion_bias_beta,
            mcts_distance_penalty_eV=mcts_distance_penalty_eV,
        )

    def extract_pathway_shortlist_tool(
        self,
        search_result: dict | None = None,
        output_base_dir: str = "",
        run_id: str = "",
    ) -> dict:
        """Export retained pathways and prepare Agent4/5 handoff payload."""
        if search_result is None:
            return {"error": "search_result is required"}
        return self.extract_pathway_shortlist(
            search_result=search_result,
            output_base_dir=output_base_dir,
            run_id=run_id,
        )

    def relax_adsorbate_tool(
        self,
        structure_path: str = "",
        adsorbate_indices: list[int] | None = None,
        output_path: str = "",
        fmax: float = 0.05,
        max_steps: int = 200,
    ) -> dict:
        """Relax adsorbate atoms on a surface structure (surface atoms fixed).

        Use after adding/removing atoms to get a physically reasonable geometry.
        Returns energy, fmax, and path to relaxed structure.
        """
        if not structure_path or not output_path:
            return {"error": "structure_path and output_path are required"}
        if not adsorbate_indices:
            return {"error": "adsorbate_indices list is required"}
        return self.relax_adsorbate_on_surface(
            structure_path=structure_path,
            adsorbate_indices=adsorbate_indices,
            output_path=output_path,
            fmax=fmax,
            max_steps=max_steps,
        )

    def to_camel_tools(self) -> Dict[str, FunctionTool]:
        """Expose CatDT tool methods as CAMEL `FunctionTool` registry."""
        return {
            "generate_surfaces": FunctionTool(self.generate_surfaces_tool),
            "predict_adsorption_sites": FunctionTool(self.predict_adsorption_sites_tool),
            "simulate_surface_reconstruction": FunctionTool(self.simulate_surface_reconstruction_tool),
            "run_neb_for_steps": FunctionTool(self.run_neb_for_steps),
            "compute_adsorption_energies": FunctionTool(self.compute_adsorption_energies),
            "run_kmc_simulation": FunctionTool(self.run_kmc_simulation),
            "generate_final_report": FunctionTool(self.generate_final_report),
            "validate_neb_endpoints": FunctionTool(self.validate_neb_endpoints),
            "run_neb_with_validation": FunctionTool(self.run_neb_with_validation),
            "check_checkpoint_exists": FunctionTool(self.check_checkpoint_exists),
            "list_available_checkpoints": FunctionTool(self.list_available_checkpoints),
            "export_step_structures": FunctionTool(self.export_step_structures),
            "build_steps_payload": FunctionTool(self.build_steps_payload_tool),
            "programmatic_validation": FunctionTool(self.programmatic_validation_tool),
            "get_step_element_deltas": FunctionTool(self.get_step_element_deltas_tool),
            "suggest_staged_positions": FunctionTool(self.suggest_staged_positions_tool),
            "validate_proposed_positions": FunctionTool(self.validate_proposed_positions_tool),
            "relax_adsorbate": FunctionTool(self.relax_adsorbate_tool),
            "run_agent45_energy_gate": FunctionTool(self.run_agent45_energy_gate_tool),
            "retrieve_agent45_memento_cases": FunctionTool(self.retrieve_agent45_memento_cases_tool),
            "retrieve_agent45_knowledge_items": FunctionTool(self.retrieve_agent45_knowledge_items_tool),
            "retrieve_agent45_skill_items": FunctionTool(self.retrieve_agent45_skill_items_tool),
            "record_agent45_memento_case": FunctionTool(self.record_agent45_memento_case_tool),
            "train_agent45_parametric_retriever": FunctionTool(self.train_agent45_parametric_retriever_tool),
            "get_agent45_parametric_retriever_status": FunctionTool(self.get_agent45_parametric_retriever_status_tool),
            # Mechanism search tools
            "initialize_mechanism_context": FunctionTool(self.initialize_mechanism_context_tool),
            "generate_candidate_steps": FunctionTool(self.generate_candidate_steps_tool),
            "run_mechanism_search": FunctionTool(self.run_mechanism_search_tool),
            "extract_pathway_shortlist": FunctionTool(self.extract_pathway_shortlist_tool),
        }


__all__ = [
    "CATDT_CORE_PATH",
    "DEPS_BASE_PATH",
    "logger",
    "Agent45GeometryTools",
    "CatDTTools",
]
