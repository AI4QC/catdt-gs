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

        All numeric ``free_energy_eV`` values are assumed to share ONE
        convention (corrected adsorption free energy ΔG_ads, or step ΔG
        derived from it — enforced upstream by the free-energy router).
        Estimates with ``free_energy_eV=None`` (or whose ``convention``
        field disagrees with the group's) are INCOMPARABLE: they are never
        used as the prune reference and are never pruned against numeric
        values — they are kept with prune_reason ``energy_unavailable_kept``.

        Args:
            sibling_estimates: list of dicts, each with at least
                ``state_id`` and ``free_energy_eV``.

        Returns:
            Dict with ``kept``, ``pruned``, ``marginal`` lists of state_ids,
            plus ``decisions`` mapping state_id → reason.
        """
        if not sibling_estimates:
            return {"kept": [], "pruned": [], "marginal": [], "decisions": {}}

        kept: List[str] = []
        pruned: List[str] = []
        marginal: List[str] = []
        decisions: Dict[str, str] = {}

        # Split incomparable (None-energy or convention-mismatched) siblings
        # from numeric ones. Incomparable siblings are always kept.
        conventions = {
            e.get("convention")
            for e in sibling_estimates
            if e.get("free_energy_eV") is not None and e.get("convention") is not None
        }
        ref_convention = next(iter(conventions)) if len(conventions) == 1 else None

        numeric: List[Dict[str, Any]] = []
        for est in sibling_estimates:
            sid = est.get("state_id", "?")
            v = est.get("free_energy_eV")
            conv = est.get("convention")
            if v is None:
                kept.append(sid)
                decisions[sid] = "energy_unavailable_kept"
                logger.info(
                    "compare_siblings: %s has no energy estimate — "
                    "incomparable, kept.", sid,
                )
                continue
            if len(conventions) > 1 and conv != ref_convention:
                kept.append(sid)
                decisions[sid] = "energy_unavailable_kept"
                logger.warning(
                    "compare_siblings: %s carries convention %r != group %r — "
                    "incomparable, kept.", sid, conv, ref_convention,
                )
                continue
            numeric.append(est)

        if not numeric:
            return {
                "kept": kept,
                "pruned": pruned,
                "marginal": marginal,
                "decisions": decisions,
            }

        def _energy(e: Dict) -> float:
            return float(e["free_energy_eV"])

        sorted_siblings = sorted(numeric, key=_energy)
        best_energy = _energy(sorted_siblings[0])

        for est in sorted_siblings:
            sid = est.get("state_id", "?")
            energy = _energy(est)
            gap = energy - best_energy

            if gap <= self.delta_keep:
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

        This is a capacity cap, not an energy comparison: nodes without an
        energy estimate (``free_energy_eV=None``) sort last and are only
        dropped when the beam overflows.

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
