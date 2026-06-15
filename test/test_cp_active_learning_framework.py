from __future__ import annotations

import pathlib
import sys

import numpy as np

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from core.active_learning.acquisition import rank_candidates
from core.active_learning.types import CandidateScore


def test_rank_candidates_balances_uncertainty_and_diversity() -> None:
    candidates = [
        CandidateScore(
            structure_id="a",
            uncertainty_force=0.9,
            uncertainty_energy=0.7,
            uncertainty_potential=0.1,
            descriptor=np.array([0.0, 0.0]),
        ),
        CandidateScore(
            structure_id="b",
            uncertainty_force=0.8,
            uncertainty_energy=0.8,
            uncertainty_potential=0.2,
            descriptor=np.array([1.0, 1.0]),
        ),
        CandidateScore(
            structure_id="c",
            uncertainty_force=0.85,
            uncertainty_energy=0.85,
            uncertainty_potential=0.3,
            descriptor=np.array([0.95, 0.95]),
        ),
        CandidateScore(
            structure_id="d",
            uncertainty_force=0.1,
            uncertainty_energy=0.1,
            uncertainty_potential=0.1,
            descriptor=np.array([3.0, 3.0]),
        ),
    ]

    selected = rank_candidates(
        candidates,
        top_k=2,
        preselect_k=4,
        weights={"force": 0.6, "energy": 0.3, "potential": 0.1},
        diversity_lambda=0.3,
    )

    ids = [x.structure_id for x in selected]
    assert len(ids) == 2
    assert "c" in ids
    assert "d" not in ids


def test_rank_candidates_works_without_descriptors() -> None:
    candidates = [
        CandidateScore(
            structure_id="x",
            uncertainty_force=0.3,
            uncertainty_energy=0.2,
            uncertainty_potential=0.1,
            descriptor=None,
        ),
        CandidateScore(
            structure_id="y",
            uncertainty_force=0.9,
            uncertainty_energy=0.9,
            uncertainty_potential=0.1,
            descriptor=None,
        ),
    ]
    selected = rank_candidates(
        candidates,
        top_k=1,
        preselect_k=2,
        weights={"force": 0.5, "energy": 0.4, "potential": 0.1},
        diversity_lambda=0.5,
    )
    assert [x.structure_id for x in selected] == ["y"]
