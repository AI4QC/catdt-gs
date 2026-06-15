"""Mechanism search stage runner.

Orchestrates Agent M1 and M2 (or runs tool-only mode) between Agent3 and
Agent4/5 in the CatDT workflow.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class MechanismStageRunner:
    """Runs the mechanism search stage after Agent3 completes.

    Can operate in two modes:
    1. **Agent mode**: Agent M1 and M2 coordinate via LLM + tool calls.
    2. **Tool-only mode**: Deterministic pipeline, no LLM involved.
       (Default for Phase 1 — simpler and more predictable.)
    """

    def __init__(self, tools: Any = None, use_agents: bool = False):
        self.tools = tools
        self.use_agents = use_agents

    def run_after_agent3(
        self,
        state: Any,
        config: Any,
        surface: Any = None,
    ) -> Dict[str, Any]:
        """Execute the mechanism search stage.

        Args:
            state: WorkflowState after Agent3 completion (LLM state only)
            config: CatDTConfig with mechanism_* fields
            surface: CatalystSurface holding physical state / file paths

        Returns:
            MultiPathwaySearchResult-compatible dict
        """
        if self.tools is None:
            return {"error": "tools not initialized"}

        # Extract context from workflow state
        initial_state = ""
        target_state = ""
        known_intermediates: List[str] = []

        if state.reaction_context:
            intermediates = list(state.reaction_context.intermediates or [])
            if len(intermediates) >= 2:
                initial_state = intermediates[0]
                target_state = intermediates[-1]
                known_intermediates = intermediates[1:-1]

        if not initial_state or not target_state:
            return {"error": "Cannot determine initial/target states from reaction context"}

        surface_id = getattr(surface, "surface_id", "") if surface is not None else ""
        surface_facet = ""
        surface_path_for_facet = getattr(surface, "surface_path", None) if surface is not None else None
        if surface_path_for_facet:
            import re
            m = re.search(r"(\d{3,4})", str(surface_path_for_facet))
            if m:
                surface_facet = m.group(1)

        recon_path = getattr(surface, "reconstructed_surface_path", None) if surface is not None else None
        clean_path = getattr(surface, "clean_slab_path", None) if surface is not None else None

        # Step 1: Initialize context
        logger.info(
            "Mechanism search: %s → %s (intermediates: %s)",
            initial_state, target_state, known_intermediates,
        )

        ctx_result = self.tools.initialize_mechanism_context(
            initial_state=initial_state,
            target_state=target_state,
            surface_id=surface_id,
            surface_facet=surface_facet,
            surface_structure_path=recon_path or surface_path_for_facet or "",
            clean_slab_path=clean_path or "",
            known_intermediates=known_intermediates,
        )

        if ctx_result.get("error"):
            logger.warning("Mechanism context init failed: %s", ctx_result["error"])
            return ctx_result

        # Step 2: Generate candidate steps
        candidates_result = self.tools.generate_candidate_steps(
            use_care=config.enable_mechanism_search,
            use_generic_ops=True,
        )
        logger.info(
            "Mechanism search: %d candidates generated",
            candidates_result.get("n_candidates", 0),
        )

        # Step 3: Run search
        exploration_mode = getattr(config, "mechanism_exploration_mode", "agent_guided")
        search_result = self.tools.run_mechanism_search(
            max_depth=config.mechanism_max_depth,
            beam_width=config.mechanism_beam_width,
            delta_keep=config.mechanism_delta_keep,
            delta_prune=config.mechanism_delta_prune,
            max_state_evaluations=config.mechanism_max_state_evaluations,
            exploration_mode=exploration_mode,
            systematic_strategy=getattr(
                config, "mechanism_systematic_strategy", "pruning",
            ),
            mcts_iterations=getattr(config, "mechanism_mcts_iterations", 100),
            mcts_c=getattr(config, "mechanism_mcts_c", 1.4),
            mcts_rollout_depth=getattr(config, "mechanism_mcts_rollout_depth", 4),
            mcts_rollout_bias_beta=getattr(
                config, "mechanism_mcts_rollout_bias_beta", 1.0,
            ),
            mcts_expansion_bias_beta=getattr(
                config, "mechanism_mcts_expansion_bias_beta", 0.5,
            ),
            mcts_distance_penalty_eV=getattr(
                config, "mechanism_mcts_distance_penalty_eV", 0.3,
            ),
        )
        logger.info(
            "Mechanism search complete: %d pathways found, stats=%s",
            search_result.get("n_pathways", 0),
            search_result.get("stats", {}),
        )

        # Step 4: Export and prepare handoff
        shortlist = self.tools.extract_pathway_shortlist(
            search_result=search_result,
            output_base_dir=state.output_base_dir,
            run_id=state.run_id,
        )

        return {
            "status": "OK",
            "context": ctx_result.get("context"),
            "search_result": search_result,
            "shortlist": shortlist,
        }

    def handoff_to_agent45_payload(
        self,
        search_output: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Convert mechanism search output to Agent4/5-ready payload.

        Returns:
            Dict with ``pathways`` list, each containing ``intermediates``
            and ``steps`` that Agent4/5 can consume.
        """
        shortlist = search_output.get("shortlist", {})
        payloads = shortlist.get("agent45_payloads", [])

        if not payloads:
            return {"pathways": [], "n_pathways": 0}

        return {
            "pathways": payloads,
            "n_pathways": len(payloads),
            "manifest": shortlist.get("manifest", {}),
        }


__all__ = ["MechanismStageRunner"]
