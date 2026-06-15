"""Mechanism search tools mixin for CatDT.

Exposes the core mechanism search capabilities as tool methods that can be
registered as CAMEL FunctionTools.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class MechanismToolsMixin:
    """Mixin providing mechanism search tool methods for CatDTTools."""

    _mechanism_context: Optional[Dict[str, Any]] = None
    _mechanism_care_bridge: Any = None
    _mechanism_care_network: Any = None

    def initialize_mechanism_context(
        self,
        initial_state: str,
        target_state: str,
        surface_id: str = "",
        surface_facet: str = "",
        bulk_formula: str = "",
        surface_structure_path: str = "",
        clean_slab_path: str = "",
        environment: Optional[Dict[str, Any]] = None,
        known_intermediates: Optional[List[str]] = None,
        constraints: str = "",
        fixed_adsorption_site: Optional[List[float]] = None,
    ) -> Dict[str, Any]:
        """Initialize the mechanism search context.

        Returns a structured context dict and a CARE domain report.
        """
        from camel_agents.mechanism_schemas import MechanismContext
        from camel_agents.tooling.care_bridge import CAREBridge

        ctx = MechanismContext(
            initial_state=initial_state,
            target_state=target_state,
            surface_id=surface_id,
            surface_facet=surface_facet,
            bulk_formula=bulk_formula,
            surface_structure_path=surface_structure_path or "",
            clean_slab_path=clean_slab_path or "",
            environment=environment or {},
            known_intermediates=known_intermediates or [],
            constraints=constraints,
            fixed_adsorption_site=fixed_adsorption_site,
        )
        self._mechanism_context = ctx.model_dump()

        # Check CARE domain
        bridge = CAREBridge(num_cpu=1)
        all_species = [initial_state, target_state] + (known_intermediates or [])
        domain_report = bridge.check_domain(species_labels=all_species)

        self._mechanism_care_bridge = bridge

        return {
            "status": "OK",
            "context": self._mechanism_context,
            "care_domain_report": domain_report.model_dump(),
        }

    def generate_candidate_steps(
        self,
        use_care: bool = True,
        use_generic_ops: bool = True,
        care_seed_species: Optional[List[str]] = None,
        electro: bool = False,
    ) -> Dict[str, Any]:
        """Generate candidate elementary steps from all enabled backends.

        Must call ``initialize_mechanism_context`` first.
        """
        if self._mechanism_context is None:
            return {"error": "Call initialize_mechanism_context first", "candidates": []}

        from core.pathway.candidate_generators import (
            CandidateGeneratorRouter,
            parse_species_elements,
        )
        from camel_agents.tooling.care_bridge import CAREBridge

        ctx = self._mechanism_context
        initial = ctx["initial_state"]
        elements = parse_species_elements(initial)

        # Optionally generate CARE network
        care_network = None
        if use_care and self._mechanism_care_bridge is not None:
            bridge = self._mechanism_care_bridge
            domain = bridge.check_domain(
                species_labels=[initial, ctx["target_state"]] + ctx.get("known_intermediates", []),
            )
            if domain.in_domain:
                seeds = care_seed_species or self._infer_care_seeds(ctx)
                bp_result = bridge.generate_blueprint(
                    seed_species=seeds,
                    electro=electro,
                )
                care_network = bp_result.get("network")
                self._mechanism_care_network = care_network
                if bp_result.get("error"):
                    logger.warning("CARE blueprint: %s", bp_result["error"])

        router = CandidateGeneratorRouter(
            enable_care=use_care and care_network is not None,
            enable_generic=use_generic_ops,
            care_bridge=self._mechanism_care_bridge,
        )

        candidates = router.generate_candidates(
            species_label=initial,
            elements=elements,
            care_network=care_network,
        )

        return {
            "status": "OK",
            "n_candidates": len(candidates),
            "candidates": candidates,
            "backends_used": {
                "care": use_care and care_network is not None,
                "generic_ops": use_generic_ops,
            },
        }

    def estimate_state_free_energy(
        self,
        species_label: str,
        backend: str = "thermal_uma",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """Estimate free energy for a single catalytic state."""
        from core.pathway.free_energy_router import FreeEnergyRouter

        surface_path = None
        if self._mechanism_context:
            surface_path = self._mechanism_context.get("surface_structure_path")

        from core.pathway.free_energy_router import ThermalStateEnergyBackend
        thermal_backend = ThermalStateEnergyBackend(tools=self)
        router = FreeEnergyRouter(thermal_backend=thermal_backend, default_backend=backend)
        return router.estimate_state_free_energy(
            species_label=species_label,
            surface_path=surface_path,
            **kwargs,
        )

    def run_mechanism_search(
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
        enable_llm_filter: bool = False,
    ) -> Dict[str, Any]:
        """Run mechanism search using the selected exploration module(s).

        Args:
            exploration_mode: ``agent_guided`` | ``systematic`` |
                ``both_sequential`` | ``both_parallel``.
            systematic_strategy: within systematic mode, which Phase-3
                algorithm to use:
                * ``quick``  — score Stage-1 paths by max-intermediate
                  UMA energy with shared cache (fast, no graph traversal).
                * ``pruning`` — layer-by-layer beam BFS over the CRN
                  with per-parent sibling pruning + distance-bucketed
                  round-robin beam (default).
                * ``mcts``    — Monte Carlo Tree Search with UCB1.
            mcts_iterations / mcts_c / mcts_rollout_depth: MCTS tuning.
        """
        if self._mechanism_context is None:
            return {"error": "Call initialize_mechanism_context first"}

        from core.pathway.energy_cache import EnergyCache

        ctx = self._mechanism_context
        all_pathways: List[Dict[str, Any]] = []
        all_stats: Dict[str, Any] = {
            "exploration_mode": exploration_mode,
            "systematic_strategy": systematic_strategy,
        }

        # Shared cache across all strategies in this run → any intermediate
        # evaluated once (by UMA) is reused by every subsequent pathway.
        surface_key = (
            f"{ctx.get('surface_id', '')}/{ctx.get('surface_facet', '')}"
            or (ctx.get("surface_structure_path") or "default")
        )
        energy_cache = EnergyCache(surface_key=surface_key)

        # --- Module 1: Agent-guided pathway exploration ---
        if exploration_mode in ("agent_guided", "both_sequential", "both_parallel"):
            guided_result = self._run_agent_guided_search(
                ctx, energy_backend, energy_cache=energy_cache,
            )
            if guided_result.get("pathways"):
                for pw in guided_result["pathways"]:
                    pw["pathway_id"] = f"guided_{pw.get('pathway_id', '')}"
                all_pathways.extend(guided_result["pathways"])
            all_stats["agent_guided"] = guided_result.get("stats", {})

        # --- Module 2: Systematic CRN exploration (pruning | mcts) ---
        if exploration_mode in ("systematic", "both_sequential", "both_parallel"):
            systematic_result = self._run_systematic_search(
                ctx,
                max_depth=max_depth,
                beam_width=beam_width,
                delta_keep=delta_keep,
                delta_prune=delta_prune,
                max_state_evaluations=max_state_evaluations,
                energy_backend=energy_backend,
                llm_filter_top_n=getattr(self, "_systematic_llm_filter_top_n", 10),
                enable_llm_filter=enable_llm_filter,
                strategy=systematic_strategy,
                mcts_iterations=mcts_iterations,
                mcts_c=mcts_c,
                mcts_rollout_depth=mcts_rollout_depth,
                mcts_rollout_bias_beta=mcts_rollout_bias_beta,
                mcts_expansion_bias_beta=mcts_expansion_bias_beta,
                mcts_distance_penalty_eV=mcts_distance_penalty_eV,
                energy_cache=energy_cache,
            )
            if systematic_result.get("pathways"):
                for pw in systematic_result["pathways"]:
                    pw["pathway_id"] = f"systematic_{pw.get('pathway_id', '')}"
                all_pathways.extend(systematic_result["pathways"])
            all_stats["systematic"] = systematic_result.get("stats", {})

        all_stats["energy_cache"] = energy_cache.stats()

        result = {
            "status": "OK",
            "n_pathways": len(all_pathways),
            "pathways": all_pathways,
            "stats": all_stats,
        }
        # Store for workflow retrieval after agent tool call
        self._mechanism_last_search_result = result
        return result

    # ---- Module 1: Agent-guided search ----

    def _run_agent_guided_search(
        self, ctx: Dict[str, Any], energy_backend: str,
        energy_cache: Any = None,
    ) -> Dict[str, Any]:
        """Agent-guided pathway exploration.

        - If user provided ``known_intermediates``: use them directly as a single path.
        - If only initial + target given: call LLM to recommend multiple competing paths,
          evaluate each with UMA, prune by energy.
        """
        initial = ctx["initial_state"]
        target = ctx["target_state"]
        known = ctx.get("known_intermediates", [])
        surface_path = ctx.get("surface_structure_path")
        fixed_site = ctx.get("fixed_adsorption_site")
        surface_id = ctx.get("surface_id", "")
        surface_facet = ctx.get("surface_facet", "")
        constraints = ctx.get("constraints", "")

        if known:
            # User gave intermediates → single path, no LLM needed
            sequences = [[initial] + [k for k in known if k != initial and k != target] + [target]]
            logger.info("Agent-guided: using %d user-supplied intermediates", len(known))
        else:
            # No intermediates → call LLM to propose multiple competing paths
            logger.info("Agent-guided: calling LLM to recommend pathways for %s → %s", initial, target)
            sequences = self._llm_recommend_pathways(
                initial=initial,
                target=target,
                surface_id=surface_id,
                surface_facet=surface_facet,
                constraints=constraints,
            )
            if not sequences:
                return {"pathways": [], "stats": {"error": "LLM returned no pathways"}}
            logger.info("Agent-guided: LLM recommended %d pathway(s)", len(sequences))

        # Evaluate each pathway with UMA
        from core.pathway.candidate_generators import AgentGuidedPathwayGenerator
        from core.pathway.free_energy_router import FreeEnergyRouter, ThermalStateEnergyBackend

        thermal_backend = ThermalStateEnergyBackend(tools=self)
        router = FreeEnergyRouter(thermal_backend=thermal_backend, default_backend=energy_backend)
        generator = AgentGuidedPathwayGenerator()

        all_pathways: List[Dict[str, Any]] = []
        total_evaluated = 0
        # Shared cache threaded from run_mechanism_search; fall back to local.
        if energy_cache is None:
            from core.pathway.energy_cache import EnergyCache
            energy_cache = EnergyCache(surface_key=surface_path or "default")

        for pw_idx, sequence in enumerate(sequences):
            steps = generator.generate_pathway_steps(sequence)
            if not steps:
                continue

            states_data = []
            for i, label in enumerate(sequence):
                hit, cached = energy_cache.get(label)
                if hit:
                    energy = cached
                elif "(g)" in label or "+" in label:
                    energy = None
                    energy_cache.put(label, energy)
                else:
                    est = router.estimate_state_free_energy(
                        label, surface_path=surface_path, fixed_site=fixed_site,
                    )
                    energy = est.get("free_energy_eV")
                    energy_cache.put(label, energy)
                    total_evaluated += 1
                    try:
                        import torch
                        if torch.cuda.is_available():
                            torch.cuda.empty_cache()
                    except Exception:
                        pass

                states_data.append({
                    "species_label": label,
                    "free_energy_eV": energy,
                    "state_id": f"guided_p{pw_idx}_s{i}",
                })

            steps_data = []
            for i, step in enumerate(steps):
                steps_data.append({
                    "step_name": f"{sequence[i]}_to_{sequence[i + 1]}",
                    "operation_type": step.get("operation_type", ""),
                    "bond_changes": step.get("bond_changes", ""),
                })

            # Compute pathway total energy change for ranking
            energies = [s["free_energy_eV"] for s in states_data if s["free_energy_eV"] is not None]
            max_energy = max(energies) if energies else None

            all_pathways.append({
                "pathway_id": f"path_{pw_idx:03d}",
                "is_retained": True,
                "states": states_data,
                "steps": steps_data,
                "max_energy_eV": max_energy,
            })

        # Rank by max intermediate energy (lower = better)
        all_pathways.sort(key=lambda p: p.get("max_energy_eV") or float("inf"))

        return {
            "pathways": all_pathways,
            "stats": {
                "module": "agent_guided",
                "n_pathways_proposed": len(sequences),
                "n_pathways_evaluated": len(all_pathways),
                "n_states_evaluated": total_evaluated,
                "used_llm": len(known) == 0,
            },
        }

    def _llm_recommend_pathways(
        self,
        initial: str,
        target: str,
        surface_id: str = "",
        surface_facet: str = "",
        constraints: str = "",
        n_pathways: int = 3,
    ) -> List[List[str]]:
        """Call LLM (as a proper CAMEL agent) to recommend multiple competing pathways.

        Uses the same CamelWorkflowAgent infrastructure as Agent1-7.
        The agent can also call mechanism tools if needed.
        """
        import json as _json

        try:
            from camel_agents.prompts import TaskPrompts, MECHANISM_CONTEXT_ROUTING_AGENT
            from camel_agents.camel_model_backend import get_camel_model_backend
            from camel_agents.runtime import CamelWorkflowAgent
            try:
                from dotenv import load_dotenv
                load_dotenv()
            except Exception:
                pass

            surface_desc = (
                f"{surface_id}({surface_facet})" if surface_facet
                else (surface_id or "metal surface")
            )
            prompt = TaskPrompts.mechanism_recommend_pathways(
                initial_state=initial,
                target_state=target,
                surface_desc=surface_desc,
                constraints=constraints,
                n_pathways=n_pathways,
            )

            llm = get_camel_model_backend(temperature=0.7, max_tokens=2000)
            agent = CamelWorkflowAgent(
                role_name=MECHANISM_CONTEXT_ROUTING_AGENT.role,
                goal="Recommend multiple competing reaction pathways",
                backstory=MECHANISM_CONTEXT_ROUTING_AGENT.backstory,
                model_backend=llm,
                tools=[],  # pure reasoning, no tool calls
                step_timeout=120,
                message_window_size=4,
                token_limit=16000,
            )
            output = agent.run(prompt, allow_tool_calls=False)
            content = str(getattr(output, "raw", "") or "")

            pathways = self._parse_pathway_json(content, initial, target)
            if pathways:
                return pathways

            logger.warning("LLM pathway recommendation: could not parse response from: %s", content[:200])
        except Exception as exc:
            logger.warning("LLM pathway recommendation failed: %s", exc)

        return self._fallback_generic_pathway(initial, target)

    @staticmethod
    def _parse_pathway_json(
        content: str, initial: str, target: str,
    ) -> List[List[str]]:
        """Extract pathway JSON from LLM response text."""
        import json as _json
        import re

        # Try to find JSON array in the response
        # Look for [[...], [...], ...] pattern
        for match in re.finditer(r'\[[\s\S]*\]', content):
            try:
                data = _json.loads(match.group())
                if isinstance(data, list) and all(isinstance(p, list) for p in data):
                    # Validate each pathway has at least 3 steps and starts/ends correctly
                    valid = []
                    for p in data:
                        if len(p) >= 3 and all(isinstance(s, str) for s in p):
                            valid.append(p)
                    if valid:
                        return valid
            except _json.JSONDecodeError:
                continue
        return []

    @staticmethod
    def _fallback_generic_pathway(initial: str, target: str) -> List[List[str]]:
        """When LLM is unavailable, return a minimal generic path."""
        from core.pathway.candidate_generators import parse_species_elements

        r_elem = parse_species_elements(initial)
        t_elem = parse_species_elements(target)

        # Simple heuristic: if target has more H than initial, it's a reduction
        h_diff = t_elem.get("H", 0) - r_elem.get("H", 0)
        o_diff = t_elem.get("O", 0) - r_elem.get("O", 0)

        intermediates = [initial]

        # Build a minimal path by stepwise H-addition / O-removal
        current = dict(r_elem)
        for _ in range(20):  # safety limit
            if current == t_elem:
                break
            # Remove O first (if target has fewer O)
            if current.get("O", 0) > t_elem.get("O", 0):
                current["O"] = current.get("O", 0) - 1
                if current["O"] == 0:
                    del current["O"]
                label = _build_fallback_label(current)
                intermediates.append(label)
                continue
            # Then add H (if target has more H)
            if current.get("H", 0) < t_elem.get("H", 0):
                current["H"] = current.get("H", 0) + 1
                label = _build_fallback_label(current)
                intermediates.append(label)
                continue
            break

        if intermediates[-1] != target:
            intermediates.append(target)

        return [intermediates]

    # ---- Module 2: Systematic CRN search ----

    def _run_systematic_search(
        self,
        ctx: Dict[str, Any],
        max_depth: int = 8,
        beam_width: int = 8,
        delta_keep: float = 0.10,
        delta_prune: float = 0.30,
        max_state_evaluations: int = 40,
        energy_backend: str = "thermal_uma",
        llm_filter_top_n: int = 10,
        enable_llm_filter: bool = False,
        strategy: str = "pruning",
        mcts_iterations: int = 100,
        mcts_c: float = 1.4,
        mcts_rollout_depth: int = 4,
        mcts_rollout_bias_beta: float = 1.0,
        mcts_expansion_bias_beta: float = 0.5,
        mcts_distance_penalty_eV: float = 0.3,
        energy_cache: Any = None,
    ) -> Dict[str, Any]:
        """Systematic mode: build CRN once, then dispatch to one of three
        Phase-3 strategies on the same Stage-1 CRN graph.

          Stage 1  — ``enumerate_all_paths_bfs`` / ``build_crn_graph``:
                     pure graph construction, no energy evaluation.
          Stage 1.5 — (optional, ``enable_llm_filter=False`` by default)
                     LLM narrows the CRN to chemically plausible species
                     by filtering the raw DFS-extracted paths.
          Stage 2 / Phase 3 — three strategies:
            * ``quick``   : score the Stage-1 paths by max-intermediate
                            UMA energy with a shared cache and return the
                            sorted list. Fast, no graph re-traversal.
            * ``pruning`` : layer-by-layer beam BFS over the CRN graph
                            with per-parent sibling comparison
                            (delta_keep / delta_prune) and a
                            distance-bucketed round-robin beam cap that
                            guarantees target-parents survive.
            * ``mcts``    : UCB1 tree search on the same CRN graph.

        All three strategies share the same Stage-1 CRN so the
        comparison is apples-to-apples.
        """
        # ---------- Shared Stage-1 setup (all 3 strategies) ----------
        from core.pathway.search_engine import MechanismSearchEngine
        from core.pathway.pruning_policy import PruningPolicy
        from core.pathway.candidate_generators import (
            CandidateGeneratorRouter, AgentGuidedPathwayGenerator,
            parse_species_elements,
        )
        from core.pathway.free_energy_router import FreeEnergyRouter, ThermalStateEnergyBackend

        initial = ctx["initial_state"]
        target = ctx["target_state"]
        surface_path = ctx.get("surface_structure_path")
        fixed_site = ctx.get("fixed_adsorption_site")

        logger.info("Systematic Stage 1: build CRN graph (no UMA)")

        router = CandidateGeneratorRouter(
            enable_care=self._mechanism_care_network is not None,
            enable_generic=True,
            care_bridge=self._mechanism_care_bridge,
        )
        initial_elements = parse_species_elements(initial)
        target_elements = parse_species_elements(target)
        # Auto-size the CRN search window to the reactant + target
        # formulas: narrows gas-partner elements AND rescales
        # max_carbon / max_heavy_atoms / max_partner_heavy_atoms so the
        # same code handles any reaction (FT, OER, NRR, …) without
        # manual tuning.
        router.configure_task(
            reactant_elements=initial_elements,
            target_elements=target_elements,
        )
        engine = MechanismSearchEngine(
            generator=router,
            energy_router=FreeEnergyRouter(default_backend="heuristic"),
            policy=PruningPolicy(max_depth=max_depth),
        )
        engine.initialize_root(initial, initial_elements)

        # Build CRN ONCE (filtered to target-reachable subgraph). All
        # three Phase-3 strategies (quick / pruning / mcts) take this
        # exact dict as input, so the comparison is strictly apples to
        # apples.
        crn_graph = engine.build_crn_graph(
            max_depth=max_depth,
            target_label=target,
            target_elements=target_elements,
        )
        logger.info(
            "Stage 1 CRN: %d nodes (target-reachable), %d edge lists",
            len(set(crn_graph.keys()) | {e["label"] for adj in crn_graph.values() for e in adj}),
            len(crn_graph),
        )

        # Derive the path list from the same CRN dict (DFS bounded by
        # max_paths=50, max_depth, max_dfs_steps).
        raw_path_dicts = engine.extract_paths_from_graph(
            crn_graph,
            root_label=initial,
            root_elements=initial_elements,
            target_label=target,
            target_elements=target_elements,
            max_depth=max_depth,
            max_paths=50,
        )
        n_raw = len(raw_path_dicts)
        logger.info("Stage 1 path extraction: %d raw pathways to target", n_raw)

        if not raw_path_dicts:
            return {
                "pathways": [],
                "stats": {"module": "systematic", "phase1_paths": 0},
            }

        raw_sequences: List[List[str]] = []
        for path in raw_path_dicts:
            seq = [entry["species_label"] for entry in path]
            if len(seq) >= 2 and seq not in raw_sequences:
                raw_sequences.append(seq)

        # MCTS dispatch — pass the same crn_graph so all three see
        # identical input.
        if strategy == "mcts":
            return self._run_systematic_mcts(
                ctx,
                max_depth=max_depth,
                max_state_evaluations=max_state_evaluations,
                energy_backend=energy_backend,
                mcts_iterations=mcts_iterations,
                mcts_c=mcts_c,
                mcts_rollout_depth=mcts_rollout_depth,
                mcts_rollout_bias_beta=mcts_rollout_bias_beta,
                mcts_expansion_bias_beta=mcts_expansion_bias_beta,
                mcts_distance_penalty_eV=mcts_distance_penalty_eV,
                energy_cache=energy_cache,
                shared_crn_graph=crn_graph,
                shared_router=router,
            )

        # ================================================================
        # Phase 2: LLM chemical plausibility filter (optional, off by default)
        # ================================================================
        if not enable_llm_filter:
            logger.info(
                "Systematic Phase 2: LLM filter disabled (enable_llm_filter=False); "
                "passing all %d raw pathways to Stage 2.",
                len(raw_sequences),
            )
            filtered_sequences = raw_sequences
        elif len(raw_sequences) > llm_filter_top_n:
            logger.info(
                "Systematic Phase 2: LLM filter %d → top %d pathways",
                len(raw_sequences), llm_filter_top_n,
            )
            filtered_sequences = self._llm_filter_pathways(
                raw_sequences, initial, target, ctx, top_n=llm_filter_top_n,
            )
        else:
            logger.info(
                "Systematic Phase 2: %d pathways ≤ top_n=%d, skip LLM filter",
                len(raw_sequences), llm_filter_top_n,
            )
            filtered_sequences = raw_sequences

        n_filtered = len(filtered_sequences)
        logger.info("Systematic Phase 2 done: %d pathways retained", n_filtered)

        # Shared Phase 3 setup: UMA backend + energy cache + helpers.
        thermal_backend = ThermalStateEnergyBackend(tools=self)
        uma_router = FreeEnergyRouter(
            thermal_backend=thermal_backend, default_backend=energy_backend,
        )
        generator = AgentGuidedPathwayGenerator()

        # Use the caller-supplied EnergyCache when available so that
        # pruning / quick share the same cache as mcts (and persist
        # across both_sequential runs). Fall back to a fresh one when
        # called outside the dispatch.
        from core.pathway.energy_cache import EnergyCache
        if energy_cache is None:
            energy_cache = EnergyCache(surface_key=surface_path or "default")
        total_evaluated = 0

        # Gas-phase free-energy DB: enables stoichiometry-corrected ΔG
        # for every CRN step. Without it, scoring degrades to absolute
        # surface G (size-biased, comparison across paths is unfair).
        from core.pathway.catdt_molecule_db import load_default as _load_db
        try:
            gas_db = _load_db()
            if not gas_db.molecules:
                gas_db = None
                logger.warning(
                    "CatDT molecule DB is empty; ΔG corrections unavailable. "
                    "Build with `python scripts/build_catdt_molecule_db.py`."
                )
        except Exception as exc:
            logger.warning("Failed to load CatDT molecule DB (%s); using naïve ΔG.", exc)
            gas_db = None

        env = ctx.get("environment", {}) or {}
        T_env = float(env.get("T", 298.0))
        P_env = float(env.get("P", 1e5))

        def _surface_thermal_correction(label: str) -> float:
            """ZPE + thermal H − TS from the matching gas molecule;
            zero if no formula match in DB."""
            if gas_db is None:
                return 0.0
            elem = parse_species_elements(label)
            sp = gas_db.species_for_formula(elem)
            if sp is None:
                return 0.0
            try:
                return (
                    gas_db.G(sp, T=T_env, P=P_env)
                    - gas_db.molecules[sp]["E_elec_eV"]
                )
            except Exception:
                return 0.0

        # UMA's compute_adsorption_energies uses hard-coded per-atom
        # reference energies internally:
        #   E_ads_atomic(*X) = E(slab+X) - E(slab) - Σ atomic_ref[e]·n_e
        # We convert to a molecular reference (UMA-computed gas molecule
        # E_elec from the CatDT DB) so different species are comparable
        # on a single thermodynamic scale:
        #   E_ads_mol(*X) = E_ads_atomic(*X) + (Σ atomic_ref) − E_gas(X)
        # When no DB match exists (radical fragment), keep atomic ref.
        SANITY_ABSG_EV = 30.0
        TYPICAL_ADS_EV = -1.5
        # OC20 atomic-reference table (mirrors fairchem_predictor.py).
        ATOMIC_REF_OC20 = {
            "H": -3.477, "C": -7.282, "N": -8.083,
            "O": -7.204, "F": -4.891, "S": -4.659,
        }

        def _atomic_ref_sum(elem: Dict[str, int]) -> float:
            return sum(n * ATOMIC_REF_OC20.get(e, 0.0) for e, n in elem.items())

        def _molecular_ref_E_elec(label: str) -> Optional[float]:
            """UMA E_elec of the matching gas-phase molecule from the DB."""
            if gas_db is None:
                return None
            elem = parse_species_elements(label)
            sp = gas_db.species_for_formula(elem)
            if sp is None:
                return None
            return gas_db.molecules[sp].get("E_elec_eV")

        def _gas_reference_G(label: str) -> Optional[float]:
            """μ_gas(species formula) + typical surface stabilization."""
            if gas_db is None:
                return None
            try:
                elem = parse_species_elements(label)
                return gas_db.mu_for_delta(elem, T=T_env, P=P_env) + TYPICAL_ADS_EV
            except Exception:
                return None

        def _eval_energy(label: str) -> Optional[float]:
            nonlocal total_evaluated
            hit, cached = energy_cache.get(label)
            if hit:
                return cached
            if "(g)" in label or "+" in label:
                if gas_db is not None:
                    try:
                        elem = parse_species_elements(label)
                        energy = gas_db.mu_for_delta(elem, T=T_env, P=P_env)
                    except Exception:
                        energy = None
                else:
                    energy = None
                energy_cache.put(label, energy)
                return energy
            est = uma_router.estimate_state_free_energy(
                label, surface_path=surface_path, fixed_site=fixed_site,
            )
            energy = est.get("free_energy_eV")
            if energy is None:
                fb = _gas_reference_G(label)
                if fb is not None:
                    energy = fb
            else:
                # G_full(slab+X) = E_total(slab+X) + thermal_correction(X)
                # E_total is the raw UMA energy (typically -200 to
                # -500 eV for slab+ads). E(slab) is constant and cancels
                # in ΔG between two surface species, so it doesn't enter
                # the calculation.
                energy = energy + _surface_thermal_correction(label)
            energy_cache.put(label, energy)
            total_evaluated += 1
            try:
                import torch
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            except Exception:
                pass
            return energy

        def _step_dG(parent_label, child_label, parent_elements, child_elements):
            """ΔG for the elementary step parent→child (eV).

            Uses stoichiometry-corrected formula when ``gas_db`` is
            available; otherwise falls back to naïve G(child)−G(parent).
            """
            pG = _eval_energy(parent_label)
            cG = _eval_energy(child_label)
            if pG is None or cG is None:
                return None
            if gas_db is None:
                return cG - pG
            dG = uma_router.step_dG(
                parent_G=pG, child_G=cG,
                parent_elements=parent_elements,
                child_elements=child_elements,
                gas_db=gas_db, T=T_env, P=P_env,
            )
            return dG if dG is not None else (cG - pG)

        target_norm = target.strip().lower().replace("*", "").replace(" ", "")
        for suffix in ("(g)", "(s)", "(l)"):
            target_norm = target_norm.replace(suffix, "")

        def _is_goal(label: str, elements: Dict[str, int]) -> bool:
            norm = label.strip().lower().replace("*", "").replace(" ", "")
            for suffix in ("(g)", "(s)", "(l)"):
                norm = norm.replace(suffix, "")
            if norm == target_norm:
                return True
            if target_elements and elements == target_elements:
                return True
            return False

        # ================================================================
        # Phase 3 — strategy dispatch
        # ================================================================
        if strategy == "quick":
            # Quick path: just score each Stage-1 path by max-intermediate
            # energy with shared cache, sort, and return.
            logger.info(
                "Systematic Phase 3 [quick]: UMA evaluation along %d paths",
                len(filtered_sequences),
            )
            pathways_data: List[Dict[str, Any]] = []
            for seq_id, sequence in enumerate(filtered_sequences):
                if total_evaluated >= max_state_evaluations:
                    logger.info(
                        "Phase 3 halted at path %d/%d: eval budget %d reached.",
                        seq_id, len(filtered_sequences), max_state_evaluations,
                    )
                    break
                steps = generator.generate_pathway_steps(sequence)
                states_data = []
                seq_elements = []
                for i, label in enumerate(sequence):
                    energy = _eval_energy(label)
                    elem = parse_species_elements(label)
                    seq_elements.append(elem)
                    states_data.append({
                        "species_label": label,
                        "free_energy_eV": energy,
                        "state_id": f"sys_p{seq_id}_s{i}",
                    })
                # Stoichiometry-corrected ΔG per step
                step_dGs: List[Optional[float]] = []
                for i in range(len(sequence) - 1):
                    dG = _step_dG(
                        sequence[i], sequence[i + 1],
                        seq_elements[i], seq_elements[i + 1],
                    )
                    step_dGs.append(dG)
                steps_data = []
                for i, step in enumerate(steps):
                    steps_data.append({
                        "step_name": f"{sequence[i]}_to_{sequence[i + 1]}",
                        "operation_type": step.get("operation_type", ""),
                        "bond_changes": step.get("bond_changes", ""),
                        "dG_step_eV": step_dGs[i] if i < len(step_dGs) else None,
                    })
                valid_dGs = [d for d in step_dGs if d is not None]
                max_step_dG = max(valid_dGs) if valid_dGs else None
                energies = [
                    s["free_energy_eV"] for s in states_data
                    if s["free_energy_eV"] is not None
                ]
                max_energy = max(energies) if energies else None
                pathways_data.append({
                    "pathway_id": f"path_{seq_id:03d}",
                    "is_retained": True,
                    "states": states_data,
                    "steps": steps_data,
                    "max_step_dG_eV": max_step_dG,   # primary ranking metric
                    "max_energy_eV": max_energy,     # legacy / debug
                })

            # Sort by max step ΔG (primary); fall back to max abs G when
            # ΔG unavailable.
            pathways_data.sort(
                key=lambda p: (
                    p.get("max_step_dG_eV") if p.get("max_step_dG_eV") is not None
                    else float("inf")
                )
            )

            return {
                "pathways": pathways_data,
                "stats": {
                    "module": "systematic",
                    "strategy": "quick",
                    "phase1_raw_paths": n_raw,
                    "phase2_filtered": n_filtered,
                    "phase3_uma_evaluated": total_evaluated,
                    "phase3_goal_paths": len(pathways_data),
                    "llm_filter_enabled": enable_llm_filter,
                },
            }

        # ================================================================
        # Phase 3 — pruning strategy: layer-by-layer beam BFS over CRN.
        # Build the full CRN adjacency dict (filtered to target-reachable
        # nodes), pre-compute dist_to_target via reverse BFS, then walk
        # the graph level by level applying per-parent sibling pruning
        # plus a distance-bucketed round-robin beam cap that guarantees
        # every distance class — including target-parents — keeps at
        # least one representative.
        # ================================================================
        logger.info(
            "Systematic Phase 3 [pruning]: beam BFS on shared CRN + sibling pruning "
            "(delta_keep=%.2f, delta_prune=%.2f, %d edge lists)",
            delta_keep, delta_prune, len(crn_graph),
        )

        # Stage 1.5 — optional LLM filter narrows allowed species.
        allowed_species: Optional[Set[str]] = None
        if enable_llm_filter:
            allowed_species = set()
            for seq in filtered_sequences:
                allowed_species.update(seq)
            logger.info(
                "Stage 1.5 LLM filter active: %d species allowed",
                len(allowed_species),
            )

        from core.pathway.pruning_policy import PruningPolicy
        policy = PruningPolicy(
            delta_keep=delta_keep, delta_prune=delta_prune,
            beam_width=beam_width, max_depth=max_depth,
            max_state_evaluations=max_state_evaluations,
        )
        prune_log: List[Dict[str, Any]] = []

        # Pre-compute dist_to_target on the filtered CRN (reverse BFS).
        # Used by min-dist sibling exemption + distance-bucket beam to
        # guarantee target-progress nodes survive energy pruning.
        dist_to_target = engine.compute_dist_to_target(
            crn_graph, target_label=target, target_elements=target_elements,
        )
        logger.info(
            "Pruning beam pre-computed dist_to_target for %d nodes; "
            "%d target-parents (dist=1).",
            len(dist_to_target),
            sum(1 for d in dist_to_target.values() if d == 1),
        )

        root_state = {"label": initial, "elements": initial_elements, "path": [initial]}
        active: List[Dict[str, Any]] = [root_state]
        goal_paths: List[List[str]] = []

        for depth_level in range(1, max_depth + 1):
            if total_evaluated >= max_state_evaluations:
                logger.info(
                    "Pruning halted at depth %d: eval budget %d reached.",
                    depth_level, max_state_evaluations,
                )
                break

            children_by_parent: Dict[int, List[Dict[str, Any]]] = {}
            for parent_id, parent_state in enumerate(active):
                parent_label = parent_state["label"]
                parent_path = parent_state["path"]
                edges = crn_graph.get(parent_label, [])
                for edge in edges:
                    child_label = edge["label"]
                    child_elem = edge["elements"]
                    if child_label in parent_path:
                        continue
                    if allowed_species is not None and child_label not in allowed_species:
                        continue
                    children_by_parent.setdefault(parent_id, []).append({
                        "label": child_label,
                        "elements": child_elem,
                        "parent_id": parent_id,
                        "path": parent_path + [child_label],
                    })

            if not children_by_parent:
                break

            new_active: List[Dict[str, Any]] = []
            for parent_id, child_list in children_by_parent.items():
                # Compare siblings on stoichiometry-corrected ΔG_step
                # (parent → child) rather than child's absolute G —
                # absolute G is biased by molecule size and unfair when
                # different element compositions appear at the same depth.
                parent_state = active[parent_id]
                parent_label = parent_state["label"]
                parent_elem = parent_state["elements"]
                sibling_estimates = []
                for cidx, child in enumerate(child_list):
                    dG = _step_dG(
                        parent_label, child["label"],
                        parent_elem, child["elements"],
                    )
                    sibling_estimates.append({
                        "state_id": f"{depth_level}:{parent_id}:{cidx}",
                        "free_energy_eV": dG,   # ΔG_step (kept key name for policy compat)
                    })
                    if total_evaluated >= max_state_evaluations:
                        break

                none_sibs = [
                    s for s in sibling_estimates
                    if s.get("free_energy_eV") is None
                ]
                valid_sibs = [
                    s for s in sibling_estimates
                    if s.get("free_energy_eV") is not None
                ]
                decision = policy.compare_siblings(valid_sibs) if valid_sibs else {
                    "kept": [], "pruned": [], "marginal": [], "decisions": {},
                }
                kept = set(decision["kept"]) | set(decision["marginal"])
                kept |= {s["state_id"] for s in none_sibs}
                pruned = set(decision["pruned"])

                # Target-progress exemption: among siblings of a parent,
                # find the minimum dist_to_target; every child at that
                # minimum distance is unconditionally kept (regardless of
                # energy). Siblings at higher dist still get energy-
                # pruned. Rationale: "if you're making the most progress
                # toward target, energy is secondary; only siblings that
                # stagnate or regress get cut on energy."
                sib_dists = [
                    dist_to_target.get(c["label"], 10**9)
                    for c in child_list
                ]
                if sib_dists:
                    min_d = min(sib_dists)
                    for cidx, child in enumerate(child_list):
                        if dist_to_target.get(child["label"], 10**9) == min_d:
                            sid = f"{depth_level}:{parent_id}:{cidx}"
                            kept.add(sid)
                            pruned.discard(sid)

                for cidx, child in enumerate(child_list):
                    sid = f"{depth_level}:{parent_id}:{cidx}"
                    if sid in pruned:
                        prune_log.append({
                            "depth": depth_level,
                            "label": child["label"],
                            "reason": decision["decisions"].get(sid, "pruned"),
                        })
                        continue
                    if sid not in kept:
                        continue
                    if _is_goal(child["label"], child["elements"]):
                        goal_paths.append(child["path"])
                    else:
                        new_active.append(child)

            # Reachability filter: drop nodes whose dist_to_target exceeds
            # remaining hops.
            remaining_depth = max_depth - depth_level
            reachable_active: List[Dict[str, Any]] = []
            for item in new_active:
                d = dist_to_target.get(item["label"])
                if d is None or d <= remaining_depth:
                    reachable_active.append(item)
                else:
                    prune_log.append({
                        "depth": depth_level,
                        "label": item["label"],
                        "reason": f"dist_to_target={d} > remaining_depth={remaining_depth}",
                    })
            new_active = reachable_active

            # Distance-bucketed round-robin beam cap.
            if len(new_active) > beam_width:
                by_dist: Dict[int, List[Dict[str, Any]]] = {}
                for item in new_active:
                    d = dist_to_target.get(item["label"], 10**9)
                    by_dist.setdefault(d, []).append(item)
                # Within each distance bucket, rank by ΔG_step from
                # the item's own parent (consistent with sibling pruning).
                def _item_score(it):
                    parent = active[it["parent_id"]]
                    dG = _step_dG(
                        parent["label"], it["label"],
                        parent["elements"], it["elements"],
                    )
                    return float("inf") if dG is None else dG
                for d, lst in by_dist.items():
                    lst.sort(key=_item_score)
                kept_items: List[Dict[str, Any]] = []
                bucket_keys = sorted(by_dist.keys())
                idx = {d: 0 for d in bucket_keys}
                while len(kept_items) < beam_width and any(
                    idx[d] < len(by_dist[d]) for d in bucket_keys
                ):
                    for d in bucket_keys:
                        if idx[d] < len(by_dist[d]):
                            kept_items.append(by_dist[d][idx[d]])
                            idx[d] += 1
                            if len(kept_items) >= beam_width:
                                break
                kept_set = {id(it) for it in kept_items}
                for item in new_active:
                    if id(item) not in kept_set:
                        prune_log.append({
                            "depth": depth_level,
                            "label": item["label"],
                            "reason": (
                                f"beam_cap dist={dist_to_target.get(item['label'], '?')}"
                            ),
                        })
                active = kept_items
            else:
                active = new_active
            if not active:
                break

        # Deduplicate goal paths (label sequence) and emit pathways.
        seen_seq: Set[Tuple[str, ...]] = set()
        unique_goal_paths: List[List[str]] = []
        for p in goal_paths:
            sig = tuple(p)
            if sig in seen_seq:
                continue
            seen_seq.add(sig)
            unique_goal_paths.append(p)

        pathways_data = []
        for seq_id, sequence in enumerate(unique_goal_paths):
            steps = generator.generate_pathway_steps(sequence)
            seq_elements = [parse_species_elements(lbl) for lbl in sequence]
            states_data = [
                {
                    "species_label": label,
                    "free_energy_eV": energy_cache.get(label)[1],
                    "state_id": f"sys_p{seq_id}_s{i}",
                }
                for i, label in enumerate(sequence)
            ]
            step_dGs: List[Optional[float]] = []
            for i in range(len(sequence) - 1):
                step_dGs.append(_step_dG(
                    sequence[i], sequence[i + 1],
                    seq_elements[i], seq_elements[i + 1],
                ))
            steps_data = []
            for i, step in enumerate(steps):
                steps_data.append({
                    "step_name": f"{sequence[i]}_to_{sequence[i + 1]}",
                    "operation_type": step.get("operation_type", ""),
                    "bond_changes": step.get("bond_changes", ""),
                    "dG_step_eV": step_dGs[i] if i < len(step_dGs) else None,
                })
            valid_dGs = [d for d in step_dGs if d is not None]
            max_step_dG = max(valid_dGs) if valid_dGs else None
            energies = [
                s["free_energy_eV"] for s in states_data
                if s["free_energy_eV"] is not None
            ]
            max_energy = max(energies) if energies else None
            pathways_data.append({
                "pathway_id": f"path_{seq_id:03d}",
                "is_retained": True,
                "states": states_data,
                "steps": steps_data,
                "max_step_dG_eV": max_step_dG,
                "max_energy_eV": max_energy,
            })

        pathways_data.sort(
            key=lambda p: (
                p.get("max_step_dG_eV") if p.get("max_step_dG_eV") is not None
                else float("inf")
            )
        )

        return {
            "pathways": pathways_data,
            "stats": {
                "module": "systematic",
                "strategy": "pruning",
                "phase1_raw_paths": n_raw,
                "phase2_filtered": n_filtered if enable_llm_filter else n_raw,
                "phase3_uma_evaluated": total_evaluated,
                "phase3_goal_paths": len(unique_goal_paths),
                "phase3_pruned_nodes": len(prune_log),
                "delta_keep_eV": delta_keep,
                "delta_prune_eV": delta_prune,
                "llm_filter_enabled": enable_llm_filter,
            },
        }

    # ---- Module 2b: Systematic MCTS search ----

    def _run_systematic_mcts(
        self,
        ctx: Dict[str, Any],
        max_depth: int = 8,
        max_state_evaluations: int = 40,
        energy_backend: str = "thermal_uma",
        mcts_iterations: int = 100,
        mcts_c: float = 1.4,
        mcts_rollout_depth: int = 4,
        mcts_rollout_bias_beta: float = 1.0,
        mcts_expansion_bias_beta: float = 0.5,
        mcts_distance_penalty_eV: float = 0.3,
        energy_cache: Any = None,
        shared_crn_graph: Optional[Dict[str, List[Dict[str, Any]]]] = None,
        shared_router: Any = None,
    ) -> Dict[str, Any]:
        """MCTS over the catalytic state CRN.

        Uses UCB1 to balance low-energy exploitation against exploration of
        branches that look bad early but may lead to low barriers later
        ("hard up front, easy after"). State energies come from UMA and
        are memoised in ``energy_cache`` so the same intermediate is never
        evaluated twice across rollouts or pathways.
        """
        from core.pathway.mcts_search import MCTSSearchEngine
        from core.pathway.candidate_generators import (
            CandidateGeneratorRouter, AgentGuidedPathwayGenerator,
            parse_species_elements,
        )
        from core.pathway.free_energy_router import (
            FreeEnergyRouter, ThermalStateEnergyBackend,
        )
        from core.pathway.energy_cache import EnergyCache

        initial = ctx["initial_state"]
        target = ctx["target_state"]
        surface_path = ctx.get("surface_structure_path")
        fixed_site = ctx.get("fixed_adsorption_site")

        if energy_cache is None:
            energy_cache = EnergyCache(surface_key=surface_path or "default")

        # MCTS is always invoked via ``_run_systematic_search`` which
        # builds the CRN once and passes it in. Direct invocation
        # without ``shared_crn_graph`` is unsupported.
        if shared_crn_graph is None or shared_router is None:
            raise ValueError(
                "_run_systematic_mcts requires shared_crn_graph + shared_router; "
                "invoke via _run_systematic_search instead."
            )
        router = shared_router
        crn_graph = shared_crn_graph
        logger.info(
            "MCTS using shared CRN: %d nodes, %d edges",
            len(crn_graph), sum(len(v) for v in crn_graph.values()),
        )

        thermal_backend = ThermalStateEnergyBackend(tools=self)
        uma_router = FreeEnergyRouter(
            thermal_backend=thermal_backend, default_backend=energy_backend,
        )

        # Gas-phase DB (loaded ONCE up-front so the energy + step-ΔG
        # functions below can capture it via closure).
        from core.pathway.catdt_molecule_db import load_default as _load_db
        try:
            _gas_db = _load_db()
            if not _gas_db.molecules:
                _gas_db = None
        except Exception:
            _gas_db = None
        env = ctx.get("environment", {}) or {}
        T_env = float(env.get("T", 298.0))
        P_env = float(env.get("P", 1e5))

        def _surface_thermal_correction_mcts(label: str) -> float:
            if _gas_db is None:
                return 0.0
            elem = parse_species_elements(label)
            sp = _gas_db.species_for_formula(elem)
            if sp is None:
                return 0.0
            try:
                return (
                    _gas_db.G(sp, T=T_env, P=P_env)
                    - _gas_db.molecules[sp]["E_elec_eV"]
                )
            except Exception:
                return 0.0

        SANITY_ABSG_EV = 30.0
        TYPICAL_ADS_EV = -1.5
        ATOMIC_REF_OC20 = {
            "H": -3.477, "C": -7.282, "N": -8.083,
            "O": -7.204, "F": -4.891, "S": -4.659,
        }

        def _atomic_ref_sum_mcts(elem):
            return sum(n * ATOMIC_REF_OC20.get(e, 0.0) for e, n in elem.items())

        def _molecular_ref_E_elec_mcts(label):
            if _gas_db is None:
                return None
            elem = parse_species_elements(label)
            sp = _gas_db.species_for_formula(elem)
            if sp is None:
                return None
            return _gas_db.molecules[sp].get("E_elec_eV")

        def _gas_reference_G_mcts(label: str) -> Optional[float]:
            if _gas_db is None:
                return None
            try:
                elem = parse_species_elements(label)
                return _gas_db.mu_for_delta(elem, T=T_env, P=P_env) + TYPICAL_ADS_EV
            except Exception:
                return None

        def _energy_fn(label: str) -> Optional[float]:
            try:
                est = uma_router.estimate_state_free_energy(
                    label, surface_path=surface_path, fixed_site=fixed_site,
                )
                val = est.get("free_energy_eV")
            except Exception as exc:
                logger.warning("MCTS energy eval failed for %s: %s", label, exc)
                val = None
            if val is None:
                fb = _gas_reference_G_mcts(label)
                if fb is not None:
                    val = fb
            else:
                # G_full = E_total(slab+X) + thermal_correction(X). E(slab)
                # cancels in ΔG between two surface species.
                val = val + _surface_thermal_correction_mcts(label)
            try:
                import torch
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            except Exception:
                pass
            return val

        def _step_dg_fn(p_label, c_label, p_elem, c_elem):
            """ΔG_step bottleneck input for MCTS reward."""
            hit_p, pG = energy_cache.get(p_label)
            if not hit_p:
                pG = _energy_fn(p_label); energy_cache.put(p_label, pG)
            hit_c, cG = energy_cache.get(c_label)
            if not hit_c:
                cG = _energy_fn(c_label); energy_cache.put(c_label, cG)
            if pG is None or cG is None:
                return None
            if _gas_db is None:
                return cG - pG
            return uma_router.step_dG(
                parent_G=pG, child_G=cG,
                parent_elements=p_elem, child_elements=c_elem,
                gas_db=_gas_db, T=T_env, P=P_env,
            ) or (cG - pG)

        engine = MCTSSearchEngine(
            generator=router,
            crn_graph=crn_graph,
            energy_cache=energy_cache,
            energy_fn=_energy_fn,
            step_dg_fn=_step_dg_fn,
            c_uct=mcts_c,
            rollout_depth=mcts_rollout_depth,
            max_depth=max_depth,
            max_iterations=mcts_iterations,
            max_state_evaluations=max_state_evaluations,
            rollout_bias_beta=mcts_rollout_bias_beta,
            expansion_bias_beta=mcts_expansion_bias_beta,
            distance_penalty_eV=mcts_distance_penalty_eV,
        )

        logger.info(
            "Systematic MCTS: iters=%d c=%.2f rollout=%d max_depth=%d "
            "budget=%d rollout_beta=%.2f expansion_beta=%.2f dist_pen=%.3f",
            mcts_iterations, mcts_c, mcts_rollout_depth, max_depth,
            max_state_evaluations,
            mcts_rollout_bias_beta, mcts_expansion_bias_beta,
            mcts_distance_penalty_eV,
        )

        mcts_result = engine.run(
            initial_label=initial,
            target_label=target,
            initial_elements=parse_species_elements(initial),
            target_elements=parse_species_elements(target),
        )
        mcts_paths = mcts_result.get("paths", [])
        logger.info("Systematic MCTS done: %d paths to target", len(mcts_paths))

        # Re-express MCTS paths in the shared pathway schema and confirm
        # energies for any intermediates touched only during rollout.
        generator = AgentGuidedPathwayGenerator()
        pathways_data: List[Dict[str, Any]] = []
        total_evaluated_before = energy_cache.stats()["misses"]

        for pw_idx, chain in enumerate(mcts_paths):
            sequence = [n["species_label"] for n in chain]
            if len(sequence) < 2:
                continue
            steps = generator.generate_pathway_steps(sequence)
            if not steps:
                continue

            states_data = []
            for i, label in enumerate(sequence):
                hit, cached = energy_cache.get(label)
                if hit:
                    energy = cached
                elif "(g)" in label or "+" in label:
                    energy = None
                    energy_cache.put(label, energy)
                else:
                    energy = _energy_fn(label)
                    energy_cache.put(label, energy)

                states_data.append({
                    "species_label": label,
                    "free_energy_eV": energy,
                    "state_id": f"mcts_p{pw_idx}_s{i}",
                })

            seq_elements = [parse_species_elements(lbl) for lbl in sequence]
            step_dGs: List[Optional[float]] = []
            for i in range(len(sequence) - 1):
                step_dGs.append(_step_dg_fn(
                    sequence[i], sequence[i + 1],
                    seq_elements[i], seq_elements[i + 1],
                ))
            steps_data = []
            for i, step in enumerate(steps):
                steps_data.append({
                    "step_name": f"{sequence[i]}_to_{sequence[i + 1]}",
                    "operation_type": step.get("operation_type", ""),
                    "bond_changes": step.get("bond_changes", ""),
                    "dG_step_eV": step_dGs[i] if i < len(step_dGs) else None,
                })
            valid_dGs = [d for d in step_dGs if d is not None]
            max_step_dG = max(valid_dGs) if valid_dGs else None
            energies = [s["free_energy_eV"] for s in states_data if s["free_energy_eV"] is not None]
            max_energy = max(energies) if energies else None

            pathways_data.append({
                "pathway_id": f"path_{pw_idx:03d}",
                "is_retained": True,
                "states": states_data,
                "steps": steps_data,
                "max_step_dG_eV": max_step_dG,
                "max_energy_eV": max_energy,
            })

        pathways_data.sort(
            key=lambda p: (
                p.get("max_step_dG_eV") if p.get("max_step_dG_eV") is not None
                else float("inf")
            )
        )

        stats = {
            "module": "systematic",
            "strategy": "mcts",
            "n_paths_found": len(mcts_paths),
            "n_pathways_emitted": len(pathways_data),
            "mcts_stats": mcts_result.get("stats", {}),
            "finalize_evaluations": energy_cache.stats()["misses"] - total_evaluated_before,
        }
        return {"pathways": pathways_data, "stats": stats}

    def _llm_filter_pathways(
        self,
        sequences: List[List[str]],
        initial: str,
        target: str,
        ctx: Dict[str, Any],
        top_n: int = 10,
    ) -> List[List[str]]:
        """Use LLM to rank pathway sequences by chemical plausibility.

        Returns the top-N most promising sequences.
        """
        import json as _json

        paths_text = ""
        for i, seq in enumerate(sequences):
            paths_text += f"  Path {i+1}: {' → '.join(seq)}\n"

        surface_id = ctx.get("surface_id", "")
        surface_facet = ctx.get("surface_facet", "")
        surface_desc = f"{surface_id}({surface_facet})" if surface_facet else (surface_id or "metal surface")

        prompt = f"""You are a heterogeneous catalysis expert. Evaluate these candidate reaction pathways and select the {top_n} most chemically plausible ones.

Reaction: {initial} → {target} on {surface_desc}

Candidate pathways (from systematic bond-operation enumeration):
{paths_text}

Ranking criteria:
1. Thermodynamic feasibility: avoid very high-energy intermediates (e.g. bare *C on most surfaces)
2. Kinetic accessibility: each step should be a single elementary reaction (one bond break/form)
3. Literature precedent: prefer pathways documented for this reaction/surface family
4. Step count: shorter pathways preferred (fewer barriers to cross)
5. Avoid unrealistic intermediates (e.g. *CH4O is unlikely on most surfaces)

Output ONLY a JSON array of the selected pathway indices (1-based), sorted best-first.
Example: [3, 1, 7, 5]

Select up to {top_n} pathways:"""

        try:
            from camel_agents.camel_model_backend import get_camel_model_backend
            from camel_agents.runtime import CamelWorkflowAgent
            from camel_agents.prompts import MECHANISM_CONTEXT_ROUTING_AGENT
            try:
                from dotenv import load_dotenv; load_dotenv()
            except Exception:
                pass

            llm = get_camel_model_backend(temperature=0.3, max_tokens=500)
            agent = CamelWorkflowAgent(
                role_name=MECHANISM_CONTEXT_ROUTING_AGENT.role,
                goal="Filter pathways by chemical plausibility",
                backstory="You are a catalysis expert ranking reaction pathways.",
                model_backend=llm,
                tools=[], step_timeout=60,
                message_window_size=4, token_limit=16000,
            )
            output = agent.run(prompt, allow_tool_calls=False)
            content = str(getattr(output, "raw", "") or "")

            # Parse indices from response
            import re
            for match in re.finditer(r'\[[\d,\s]+\]', content):
                try:
                    indices = _json.loads(match.group())
                    if isinstance(indices, list) and all(isinstance(x, int) for x in indices):
                        # Convert 1-based to 0-based, filter valid
                        selected = []
                        for idx in indices:
                            i = idx - 1
                            if 0 <= i < len(sequences) and sequences[i] not in selected:
                                selected.append(sequences[i])
                            if len(selected) >= top_n:
                                break
                        if selected:
                            logger.info("LLM filter selected %d/%d pathways", len(selected), len(sequences))
                            return selected
                except _json.JSONDecodeError:
                    continue

            logger.warning("LLM filter: could not parse response, returning first %d", top_n)
        except Exception as exc:
            logger.warning("LLM filter failed: %s, returning first %d", exc, top_n)

        return sequences[:top_n]

    def extract_pathway_shortlist(
        self,
        search_result: Optional[Dict[str, Any]] = None,
        output_base_dir: str = "",
        run_id: str = "",
    ) -> Dict[str, Any]:
        """Export retained pathways to disk and prepare Agent4/5 handoff.

        Args:
            search_result: output from ``run_mechanism_search``
            output_base_dir: base output directory
            run_id: workflow run ID
        """
        if search_result is None:
            return {"error": "No search result provided"}

        from core.pathway.pathway_bundle import PathwayBundleBuilder

        if not output_base_dir and self._mechanism_context:
            output_base_dir = "output/catdt_workflow"
        if not run_id:
            run_id = "mechanism_search"

        builder = PathwayBundleBuilder(
            output_base_dir=output_base_dir,
            run_id=run_id,
        )

        retained = [p for p in search_result.get("pathways", []) if p.get("is_retained", True)]

        manifest = builder.export_multi_pathway_bundle(retained)

        # Build Agent4/5 handoff: for each retained pathway, create step-level inputs
        agent45_payloads: List[Dict[str, Any]] = []
        for pw in retained:
            steps = pw.get("steps", [])
            states = pw.get("states", [])
            if len(states) >= 2:
                payload = {
                    "pathway_id": pw["pathway_id"],
                    "intermediates": [s["species_label"] for s in states],
                    "steps": [
                        {
                            "step_name": step.get("step_name", ""),
                            "reactant_formula": states[i]["species_label"] if i < len(states) else "",
                            "product_formula": states[i + 1]["species_label"] if i + 1 < len(states) else "",
                        }
                        for i, step in enumerate(steps)
                    ],
                }
                agent45_payloads.append(payload)

        shortlist_result = {
            "status": "OK",
            "manifest": manifest,
            "agent45_payloads": agent45_payloads,
            "n_retained": len(retained),
        }
        # Store for workflow retrieval after agent tool call
        self._mechanism_last_shortlist = shortlist_result
        return shortlist_result

    # ---- internal helpers ----

    @staticmethod
    def _infer_care_seeds(ctx: Dict[str, Any]) -> List[str]:
        """Infer CARE seed species from mechanism context."""
        seeds = set()
        for label in [ctx.get("initial_state", ""), ctx.get("target_state", "")] + ctx.get("known_intermediates", []):
            cleaned = label.replace("*", "").replace("(g)", "").replace("(s)", "").strip()
            if cleaned:
                seeds.add(cleaned)
        return sorted(seeds)


def _build_fallback_label(elements: Dict[str, int]) -> str:
    """Build *CHO style label from element dict (fallback helper)."""
    parts = []
    for elem in ["C", "H"]:
        count = elements.get(elem, 0)
        if count > 0:
            parts.append(elem if count == 1 else f"{elem}{count}")
    for elem in sorted(elements.keys()):
        if elem in ("C", "H"):
            continue
        count = elements[elem]
        if count > 0:
            parts.append(elem if count == 1 else f"{elem}{count}")
    return f"*{''.join(parts)}"


__all__ = ["MechanismToolsMixin"]
