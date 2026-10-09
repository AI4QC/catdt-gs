"""Shared intermediate energy cache for mechanism search.

Caches ΔG_ads / free-energy evaluations across pathways and search
strategies (agent-guided, systematic-pruning, systematic-mcts). A single
cache instance is threaded through one ``run_mechanism_search`` call so
repeated intermediates are never re-evaluated.

Keys are ``(surface_key, species_label)``; values are floats on the shared
ΔG_ads convention (see ``core.pathway.free_energy_router``) or ``None``
for species that could not be evaluated.

Persistence is JSON (entries are floats / None — no structures are
persisted), tagged with metadata (surface key, model name, temperature,
pressure, convention). A persisted cache whose metadata does not match the
current settings is ignored on load, so energies computed under one
model / set of operating conditions are never silently reused under
another.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

logger = logging.getLogger(__name__)


CacheKey = Tuple[str, str]

_PERSIST_FORMAT = "catdt-energy-cache-json-v2"


class EnergyCache:
    """In-memory cache with optional JSON persistence."""

    def __init__(
        self,
        surface_key: str = "",
        persist_path: Optional[str] = None,
        model_name: Optional[str] = None,
        temperature_K: Optional[float] = None,
        pressure_Pa: Optional[float] = None,
        convention: str = "dG_ads",
    ):
        self.surface_key = surface_key or "default"
        self.persist_path = persist_path
        self.model_name = model_name
        self.temperature_K = temperature_K
        self.pressure_Pa = pressure_Pa
        self.convention = convention
        self._store: Dict[CacheKey, Optional[float]] = {}
        self._hits = 0
        self._misses = 0

        if persist_path:
            self._load()

    # ----- public API -----

    def get(self, species_label: str) -> Tuple[bool, Optional[float]]:
        """Return ``(hit, energy)``. Hit is True even when the stored value
        is None (species evaluated and found unevaluable)."""
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
        deciding which labels are worth computing (unevaluable labels should
        return None and still be cached).
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

    # ----- persistence -----

    def _metadata(self) -> Dict[str, Any]:
        return {
            "surface_key": self.surface_key,
            "model_name": self.model_name,
            "temperature_K": self.temperature_K,
            "pressure_Pa": self.pressure_Pa,
            "convention": self.convention,
        }

    def _metadata_matches(self, saved: Dict[str, Any]) -> bool:
        """A saved cache is reusable only when every metadata field that is
        set on EITHER side matches exactly. ``None == None`` is a match
        (field unknown on both sides), but a value on one side and None on
        the other is a mismatch — we cannot prove the energies were
        computed under the current settings."""
        current = self._metadata()
        for field, cur_val in current.items():
            if saved.get(field) != cur_val:
                return False
        return True

    def save(self) -> None:
        if not self.persist_path:
            return
        try:
            p = Path(self.persist_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            # Persist only entries for this cache's surface_key; values are
            # floats / None (no Atoms or other structures are ever cached).
            store = {
                label: energy
                for (skey, label), energy in self._store.items()
                if skey == self.surface_key
            }
            payload = {
                "format": _PERSIST_FORMAT,
                "metadata": self._metadata(),
                "store": store,
            }
            with open(p, "w") as f:
                json.dump(payload, f, indent=2, sort_keys=True)
            logger.info("EnergyCache saved %d entries to %s", len(store), p)
        except Exception as exc:
            logger.warning("EnergyCache save failed: %s", exc)

    # ----- internal -----

    def _load(self) -> None:
        p = Path(self.persist_path)
        if not p.exists():
            return
        try:
            with open(p) as f:
                data = json.load(f)
            if not isinstance(data, dict) or data.get("format") != _PERSIST_FORMAT:
                logger.info(
                    "EnergyCache: skip load of %s — unsupported / legacy "
                    "format (expected %s).", p, _PERSIST_FORMAT,
                )
                return
            saved_meta = data.get("metadata", {}) or {}
            if not self._metadata_matches(saved_meta):
                logger.info(
                    "EnergyCache: skip load of %s — metadata mismatch "
                    "(saved %r != current %r).", p, saved_meta, self._metadata(),
                )
                return
            store = data.get("store", {}) or {}
            for label, energy in store.items():
                value = float(energy) if energy is not None else None
                self._store[(self.surface_key, str(label))] = value
            logger.info("EnergyCache loaded %d entries from %s", len(store), p)
        except Exception as exc:
            logger.warning("EnergyCache load failed: %s", exc)


__all__ = ["EnergyCache"]
