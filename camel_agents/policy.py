"""Evolvability policy for CAMEL workflow strategy adaptation."""

from __future__ import annotations

import base64
import json
import logging
import pickle
import random
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class EvolvablePolicy:
    """Simple bandit-style evolvability policy updated from workflow rewards."""

    def __init__(
        self,
        policy_path: Path,
        available_strategies: Optional[List[str]] = None,
        epsilon: float = 0.2,
        learning_rate: float = 0.3,
        seed: int = 42,
    ):
        self.policy_path = Path(policy_path)
        self.available_strategies = available_strategies or ["balanced", "conservative", "exploratory"]
        self.epsilon = float(epsilon)
        self.learning_rate = float(learning_rate)
        self.seed = int(seed)
        self.q_values: Dict[str, float] = {name: 0.0 for name in self.available_strategies}
        self.visits: Dict[str, int] = {name: 0 for name in self.available_strategies}
        self._rng = random.Random(self.seed)
        self._load()

    def _load(self) -> None:
        if not self.policy_path.exists():
            return
        try:
            payload = json.loads(self.policy_path.read_text(encoding="utf-8"))
            q_values = payload.get("q_values", {})
            visits = payload.get("visits", {})
            for name in self.available_strategies:
                if name in q_values:
                    self.q_values[name] = float(q_values[name])
                if name in visits:
                    self.visits[name] = int(visits[name])
            rng_state = payload.get("rng_state")
            if isinstance(rng_state, list):
                # New JSON-native format: nested lists of ints (see _save).
                self._rng.setstate(self._rng_state_from_json(rng_state))
            elif isinstance(rng_state, str) and rng_state.strip():
                # Legacy base64-pickle format (older policy files).
                self._rng.setstate(pickle.loads(base64.b64decode(rng_state.strip().encode("ascii"))))
        except Exception as exc:
            logger.warning("Failed to load evolvable policy: %s", exc)

    @staticmethod
    def _rng_state_to_json(state):
        """Convert random.Random.getstate() (nested tuples of ints) to lists."""
        return [list(item) if isinstance(item, tuple) else item for item in state]

    @staticmethod
    def _rng_state_from_json(state):
        """Inverse of _rng_state_to_json: lists back to the tuples setstate expects."""
        return tuple(tuple(item) if isinstance(item, list) else item for item in state)

    def _save(self) -> None:
        self.policy_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "q_values": self.q_values,
            "visits": self.visits,
            "epsilon": self.epsilon,
            "learning_rate": self.learning_rate,
            "available_strategies": self.available_strategies,
            "seed": self.seed,
            "rng_state": self._rng_state_to_json(self._rng.getstate()),
        }
        self.policy_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    def choose_strategy(self) -> str:
        if self._rng.random() < self.epsilon:
            return self._rng.choice(self.available_strategies)
        return max(self.available_strategies, key=lambda name: self.q_values.get(name, 0.0))

    def update(self, strategy: str, reward: float) -> None:
        if strategy not in self.q_values:
            return
        old = self.q_values[strategy]
        self.q_values[strategy] = old + self.learning_rate * (float(reward) - old)
        self.visits[strategy] = self.visits.get(strategy, 0) + 1
        self._save()
