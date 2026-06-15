from __future__ import annotations

from typing import Sequence

import numpy as np

from .types import CandidateScore


def _minmax(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return values
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if np.isclose(vmax, vmin):
        return np.zeros_like(values)
    return (values - vmin) / (vmax - vmin)


def _pairwise_distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def rank_candidates(
    candidates: Sequence[CandidateScore],
    top_k: int,
    preselect_k: int,
    weights: dict[str, float],
    diversity_lambda: float,
) -> list[CandidateScore]:
    if top_k <= 0:
        return []
    if len(candidates) <= top_k:
        return list(candidates)

    force = np.array([c.uncertainty_force for c in candidates], dtype=float)
    energy = np.array([c.uncertainty_energy for c in candidates], dtype=float)
    potential = np.array([c.uncertainty_potential for c in candidates], dtype=float)

    score_unc = (
        weights.get("force", 0.6) * _minmax(force)
        + weights.get("energy", 0.3) * _minmax(energy)
        + weights.get("potential", 0.1) * _minmax(potential)
    )

    ranked_indices = np.argsort(-score_unc)
    shortlist_indices = ranked_indices[: min(preselect_k, len(ranked_indices))]
    shortlist = [candidates[int(i)] for i in shortlist_indices]

    if any(c.descriptor is None for c in shortlist) or diversity_lambda <= 0.0:
        return shortlist[:top_k]

    selected: list[CandidateScore] = [shortlist[0]]
    while len(selected) < min(top_k, len(shortlist)):
        best_idx = -1
        best_score = -1.0
        for idx, cand in enumerate(shortlist):
            if cand in selected:
                continue
            unc_component = float(score_unc[int(shortlist_indices[idx])])
            dist = min(
                _pairwise_distance(cand.descriptor, chosen.descriptor)  # type: ignore[arg-type]
                for chosen in selected
            )
            combined = (
                1.0 - diversity_lambda
            ) * unc_component + diversity_lambda * dist
            if combined > best_score:
                best_score = combined
                best_idx = idx
        if best_idx < 0:
            break
        selected.append(shortlist[best_idx])

    return selected
