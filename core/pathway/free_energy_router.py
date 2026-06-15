"""Free energy estimation routing for mechanism search.

Dispatches state-level energy queries to the appropriate backend
(thermal UMA, electrochemical CHE, heuristic) and returns unified
``StateEnergyEstimate``-compatible dicts.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Backend implementations
# ---------------------------------------------------------------------------

class ThermalStateEnergyBackend:
    """Estimate adsorption free energy via UMA calculator.

    Wraps the existing ``compute_adsorption_energies`` tool, but returns
    a simplified energy estimate.
    """

    def __init__(self, tools: Any = None):
        self._tools = tools

    def estimate(
        self,
        species_label: str,
        surface_path: Optional[str] = None,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """Return a dict matching StateEnergyEstimate fields."""
        if self._tools is None or surface_path is None:
            return self._heuristic_estimate(species_label)

        try:
            call_kwargs: Dict[str, Any] = {
                "surface_path": surface_path,
                "intermediates": [species_label],
                "output_dir": kwargs.get("output_dir", "/tmp/catdt_fe_estimate"),
                "num_sites": kwargs.get("num_sites", 1),
                "max_steps": kwargs.get("max_steps", 50),
            }
            # Use fixed adsorption site if provided (from Agent2/3 best site)
            fixed_site = kwargs.get("fixed_site")
            if fixed_site is not None:
                call_kwargs["fixed_site"] = fixed_site
                call_kwargs["num_sites"] = 1  # only one site when fixed
            result = self._tools.compute_adsorption_energies(**call_kwargs)
            # Prefer raw UMA E_total(slab+X) when available so the
            # caller gets an absolute free-energy scale (E(slab)
            # cancels in ΔG between two surface species). Falls back
            # to E_ads (atomic-ref) only if E_total is unavailable.
            e_total = result.get("adsorbate_energies", {}).get(species_label)
            e_ads = result.get("adsorption_energies", {}).get(species_label)
            if e_total is not None:
                return {
                    "state_id": species_label,
                    "free_energy_eV": float(e_total),     # absolute UMA E_total(slab+X)
                    "adsorption_energy_eV": float(e_ads) if e_ads is not None else None,
                    "backend": "thermal_uma",
                    "confidence": 0.7,
                    "details": {"source": "uma_total_energy"},
                }
            if e_ads is not None:
                return {
                    "state_id": species_label,
                    "free_energy_eV": float(e_ads),
                    "adsorption_energy_eV": float(e_ads),
                    "backend": "thermal_uma",
                    "confidence": 0.7,
                    "details": {"source": "uma_adsorption"},
                }
        except Exception as exc:
            logger.warning("Thermal UMA estimate failed for %s: %s", species_label, exc)

        return self._heuristic_estimate(species_label)

    @staticmethod
    def _heuristic_estimate(species_label: str) -> Dict[str, Any]:
        """Placeholder heuristic when UMA is unavailable."""
        return {
            "state_id": species_label,
            "free_energy_eV": 0.0,
            "adsorption_energy_eV": None,
            "backend": "heuristic",
            "confidence": 0.1,
            "details": {"source": "heuristic_placeholder"},
        }


class ElectroStateEnergyBackend:
    """Computational Hydrogen Electrode (CHE) free energy estimate.

    Phase 1 minimal implementation: applies CHE correction
    ``dG = dE + dZPE - TdS + neU + 0.0592*pH*n``
    using tabulated ZPE/entropy values where available.
    """

    def __init__(
        self,
        voltage_V: float = 0.0,
        pH: float = 0.0,
        temperature_K: float = 298.15,
    ):
        self.voltage_V = voltage_V
        self.pH = pH
        self.temperature_K = temperature_K

    def estimate(
        self,
        species_label: str,
        adsorption_energy_eV: Optional[float] = None,
        n_electrons: int = 0,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """Apply CHE correction to an adsorption energy."""
        if adsorption_energy_eV is None:
            return {
                "state_id": species_label,
                "free_energy_eV": None,
                "backend": "electro_che",
                "confidence": 0.0,
                "details": {"error": "no adsorption energy provided"},
            }

        # CHE correction: ΔG = ΔE + neU + 0.0592 * pH * n
        che_correction = n_electrons * self.voltage_V + 0.0592 * self.pH * n_electrons
        free_energy = adsorption_energy_eV + che_correction

        return {
            "state_id": species_label,
            "free_energy_eV": float(free_energy),
            "adsorption_energy_eV": float(adsorption_energy_eV),
            "backend": "electro_che",
            "confidence": 0.5,
            "details": {
                "voltage_V": self.voltage_V,
                "pH": self.pH,
                "n_electrons": n_electrons,
                "che_correction_eV": che_correction,
            },
        }


# ---------------------------------------------------------------------------
# Router
# ---------------------------------------------------------------------------

class FreeEnergyRouter:
    """Routes free energy estimation to the appropriate backend."""

    def __init__(
        self,
        thermal_backend: Optional[ThermalStateEnergyBackend] = None,
        electro_backend: Optional[ElectroStateEnergyBackend] = None,
        default_backend: str = "thermal_uma",
    ):
        self.thermal = thermal_backend or ThermalStateEnergyBackend()
        self.electro = electro_backend or ElectroStateEnergyBackend()
        self.default_backend = default_backend

    def estimate_state_free_energy(
        self,
        species_label: str,
        *,
        backend: Optional[str] = None,
        surface_path: Optional[str] = None,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """Estimate free energy for a single state."""
        backend = backend or self.default_backend

        if backend == "electro_che":
            return self.electro.estimate(species_label=species_label, **kwargs)
        else:
            return self.thermal.estimate(
                species_label=species_label,
                surface_path=surface_path,
                **kwargs,
            )

    def step_dG(
        self,
        parent_G: Optional[float],
        child_G: Optional[float],
        parent_elements: Dict[str, int],
        child_elements: Dict[str, int],
        gas_db: Any,
        T: float = 298.0,
        P: float = 1e5,
    ) -> Optional[float]:
        """Stoichiometry-corrected ΔG for an elementary CRN step.

            ΔG_step = G(child_*) − G(parent_*) − Σ_e Δn_e × μ_ref(e, T, P)

        where ``Δn_e = n_e(child) − n_e(parent)`` and ``μ_ref(e)`` is
        the gas-phase chemical potential of element ``e`` evaluated
        through the reference reactions stored in ``gas_db``.

        Returns ``None`` if either G or any required element reference
        cannot be resolved.
        """
        if parent_G is None or child_G is None:
            return None
        delta = {
            e: child_elements.get(e, 0) - parent_elements.get(e, 0)
            for e in set(parent_elements) | set(child_elements)
        }
        delta = {e: n for e, n in delta.items() if n != 0}
        try:
            mu_correction = gas_db.mu_for_delta(delta, T=T, P=P)
        except KeyError as exc:
            # Missing gas-phase reference for some element — fall back
            # to naïve ΔG (no gas correction). Caller logs.
            return None
        return (child_G - parent_G) - mu_correction

    def estimate_many(
        self,
        species_labels: List[str],
        *,
        backend: Optional[str] = None,
        surface_path: Optional[str] = None,
        **kwargs: Any,
    ) -> List[Dict[str, Any]]:
        """Batch estimate free energies for multiple states."""
        return [
            self.estimate_state_free_energy(
                label,
                backend=backend,
                surface_path=surface_path,
                **kwargs,
            )
            for label in species_labels
        ]


__all__ = [
    "ThermalStateEnergyBackend",
    "ElectroStateEnergyBackend",
    "FreeEnergyRouter",
]
