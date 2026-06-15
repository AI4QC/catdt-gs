"""Pruning policy for mechanism search beam management.

Implements sibling comparison, beam selection, and budget control.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class PruningPolicy:
    """Configurable policy for pathway branch pruning.

    The core logic:
    - If sibling energy gap <= ``delta_keep`` → keep both
    - If sibling energy gap >= ``delta_prune`` → prune the higher one
    - In between → keep but mark as marginal
    - After sibling comparison, apply beam width to limit frontier size
    """

    delta_keep: float = 0.10   # eV
    delta_prune: float = 0.30  # eV
    beam_width: int = 8
    max_depth: int = 8
    max_state_evaluations: int = 40

    def compare_siblings(
        self,
        sibling_estimates: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Compare sibling candidates and decide keep/prune.

        Args:
            sibling_estimates: list of dicts, each with at least
                ``state_id`` and ``free_energy_eV``.

        Returns:
            Dict with ``kept``, ``pruned``, ``marginal`` lists of state_ids,
            plus ``decisions`` mapping state_id → reason.
        """
        if not sibling_estimates:
            return {"kept": [], "pruned": [], "marginal": [], "decisions": {}}

        # Sort by energy (None → infinity)
        def _energy(e: Dict) -> float:
            v = e.get("free_energy_eV")
            return float(v) if v is not None else float("inf")

        sorted_siblings = sorted(sibling_estimates, key=_energy)
        best_energy = _energy(sorted_siblings[0])

        kept: List[str] = []
        pruned: List[str] = []
        marginal: List[str] = []
        decisions: Dict[str, str] = {}

        for est in sorted_siblings:
            sid = est.get("state_id", "?")
            energy = _energy(est)
            gap = energy - best_energy

            if energy == float("inf"):
                pruned.append(sid)
                decisions[sid] = "pruned: no energy estimate"
            elif gap <= self.delta_keep:
                kept.append(sid)
                decisions[sid] = f"kept: gap={gap:.3f} eV <= delta_keep={self.delta_keep}"
            elif gap >= self.delta_prune:
                pruned.append(sid)
                decisions[sid] = f"pruned: gap={gap:.3f} eV >= delta_prune={self.delta_prune}"
            else:
                marginal.append(sid)
                decisions[sid] = f"marginal: delta_keep < gap={gap:.3f} eV < delta_prune"

        return {
            "kept": kept,
            "pruned": pruned,
            "marginal": marginal,
            "decisions": decisions,
        }

    def select_frontier(
        self,
        nodes: List[Dict[str, Any]],
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Select the top ``beam_width`` nodes from a frontier.

        Args:
            nodes: list of dicts with at least ``state_id`` and ``free_energy_eV``.

        Returns:
            (selected, dropped) tuple.
        """
        def _energy(n: Dict) -> float:
            v = n.get("free_energy_eV")
            return float(v) if v is not None else float("inf")

        sorted_nodes = sorted(nodes, key=_energy)
        selected = sorted_nodes[: self.beam_width]
        dropped = sorted_nodes[self.beam_width :]
        return selected, dropped

    def should_stop(
        self,
        depth: int,
        evaluations_used: int,
    ) -> Tuple[bool, str]:
        """Check whether the search should terminate."""
        if depth >= self.max_depth:
            return True, f"max_depth={self.max_depth} reached"
        if evaluations_used >= self.max_state_evaluations:
            return True, f"max_state_evaluations={self.max_state_evaluations} reached"
        return False, ""


__all__ = ["PruningPolicy"]
