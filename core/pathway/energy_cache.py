"""Shared intermediate energy cache for mechanism search.

Caches UMA/free-energy evaluations across pathways and search strategies
(agent-guided, systematic-pruning, systematic-mcts). A single cache instance
is threaded through one ``run_mechanism_search`` call so repeated
intermediates are never re-evaluated.

Keys are ``(surface_key, species_label)``; gas-phase and composite labels
(``(g)``, ``+``) are stored with value ``None`` per upstream convention.
"""

from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

logger = logging.getLogger(__name__)


CacheKey = Tuple[str, str]


class EnergyCache:
    """In-memory cache with optional pickle persistence."""

    def __init__(
        self,
        surface_key: str = "",
        persist_path: Optional[str] = None,
    ):
        self.surface_key = surface_key or "default"
        self.persist_path = persist_path
        self._store: Dict[CacheKey, Optional[float]] = {}
        self._hits = 0
        self._misses = 0

        if persist_path:
            self._load()

    # ----- public API -----

    def get(self, species_label: str) -> Tuple[bool, Optional[float]]:
        """Return ``(hit, energy)``. Hit is True even when stored value is None
        (gas-phase / composite labels intentionally cached as None)."""
        key = (self.surface_key, species_label)
        if key in self._store:
            self._hits += 1
            return True, self._store[key]
        self._misses += 1
        return False, None

    def put(self, species_label: str, energy: Optional[float]) -> None:
        key = (self.surface_key, species_label)
        self._store[key] = energy

    def get_or_compute(
        self,
        species_label: str,
        compute_fn: Callable[[str], Optional[float]],
    ) -> Optional[float]:
        """Return cached energy, or call ``compute_fn(label)`` and cache it.

        ``compute_fn`` is invoked only on miss. Callers are responsible for
        deciding which labels are worth computing (gas-phase should return
        None and still be cached).
        """
        hit, val = self.get(species_label)
        if hit:
            return val
        val = compute_fn(species_label)
        self.put(species_label, val)
        return val

    def stats(self) -> Dict[str, int]:
        return {
            "size": len(self._store),
            "hits": self._hits,
            "misses": self._misses,
            "surface_key": self.surface_key,
        }

    def save(self) -> None:
        if not self.persist_path:
            return
        try:
            p = Path(self.persist_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(p, "wb") as f:
                pickle.dump(
                    {"surface_key": self.surface_key, "store": self._store},
                    f,
                )
            logger.info("EnergyCache saved %d entries to %s", len(self._store), p)
        except Exception as exc:
            logger.warning("EnergyCache save failed: %s", exc)

    # ----- internal -----

    def _load(self) -> None:
        p = Path(self.persist_path)
        if not p.exists():
            return
        try:
            with open(p, "rb") as f:
                data = pickle.load(f)
            store = data.get("store", {})
            saved_surface = data.get("surface_key", "")
            if saved_surface and saved_surface != self.surface_key:
                logger.info(
                    "EnergyCache: skip load; saved surface %r != current %r",
                    saved_surface, self.surface_key,
                )
                return
            self._store.update(store)
            logger.info("EnergyCache loaded %d entries from %s", len(store), p)
        except Exception as exc:
            logger.warning("EnergyCache load failed: %s", exc)


__all__ = ["EnergyCache"]
