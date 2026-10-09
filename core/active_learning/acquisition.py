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

    # Track selections by shortlist index: `cand in selected` would invoke
    # the dataclass __eq__ on ndarray fields and raise ValueError.
    selected_indices: list[int] = [0]
    while len(selected_indices) < min(top_k, len(shortlist)):
        selected_set = set(selected_indices)
        remaining = [i for i in range(len(shortlist)) if i not in selected_set]
        if not remaining:
            break
        dists = np.array(
            [
                min(
                    _pairwise_distance(
                        shortlist[i].descriptor,  # type: ignore[arg-type]
                        shortlist[j].descriptor,  # type: ignore[arg-type]
                    )
                    for j in selected_indices
                )
                for i in remaining
            ],
            dtype=float,
        )
        # Normalize distances over the current shortlist so they are on the
        # same [0, 1] scale as the min-max-normalized uncertainty score.
        dists_norm = _minmax(dists)
        best_idx = -1
        best_score = -1.0
        for pos, idx in enumerate(remaining):
            unc_component = float(score_unc[int(shortlist_indices[idx])])
            combined = (
                1.0 - diversity_lambda
            ) * unc_component + diversity_lambda * float(dists_norm[pos])
            if combined > best_score:
                best_score = combined
                best_idx = idx
        if best_idx < 0:
            break
        selected_indices.append(best_idx)

    return [shortlist[i] for i in selected_indices]
