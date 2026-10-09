"""Free energy estimation routing for mechanism search.

Energy convention (single source of truth)
------------------------------------------
Every species-level value emitted by this module (``free_energy_eV``) is
the **corrected adsorption free energy**

    ΔG_ads(X) = E(slab+X, relaxed) − E(clean slab) − G_gas(X; T, P)

where ``G_gas`` is the gas-phase chemical potential of X's composition at
the operating temperature and pressure, resolved through the CatDT
molecule database (``core.pathway.catdt_molecule_db``) via
``mu_for_delta``: the Gibbs free energy (electronic + ZPE + vibrational /
rotational / translational thermal terms) of the matching gas molecule or
radical when one is stored, with element-reference decomposition (e.g.
μ(H) = ½ G(H₂)) as fallback for compositions without a stored molecule
(monatomic adsorbates such as *H, *O). The per-state reference choice
cancels exactly in step free energies (it is subtracted and re-added by
the same resolver); it only normalises state-level comparisons. All
search algorithms (beam search, MCTS, pruning, ranking) compare values on
this one scale.

Step free energies are differences of these plus gas-side bookkeeping via
the same database (see ``FreeEnergyRouter.step_dG``):

    ΔG_step = ΔG_ads(child) − ΔG_ads(parent)
              + [G_gas(child) − G_gas(parent) − Σ_e Δn_e · μ_e(T, P)]

When an estimate cannot be produced on this scale (no UMA tools, no
surface, no gas reference in the DB), the result carries
``free_energy_eV = None`` with ``confidence = 0.0`` and
``comparable = False`` — callers must treat such states as INCOMPARABLE
(never prune them against numeric values, never use them as a prune
reference). No fabricated placeholder values are ever returned.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


# OC20 per-atom reference energies (eV) used internally by UMA's
# ``compute_adsorption_energies`` (mirrors
# ``fairchem_predictor.FairchemPredictor.ATOMIC_REFERENCE_ENERGIES``):
#   E_ads_atomic(*X) = E(slab+X) − E(slab) − Σ_e n_e · ATOMIC_REF_OC20[e]
# We invert this to recover the slab-referenced total and re-reference it
# against the gas-phase molecule:
#   ΔG_ads(*X) = E_ads_atomic(*X) + Σ_e n_e · ATOMIC_REF_OC20[e] − G_gas(X)
ATOMIC_REF_OC20: Dict[str, float] = {
    "H": -3.477, "C": -7.282, "N": -8.083,
    "O": -7.204, "F": -4.891, "S": -4.659,
}


def _unavailable_estimate(species_label: str, reason: str) -> Dict[str, Any]:
    """Estimate result for a state that cannot be evaluated on the
    ΔG_ads scale. ``free_energy_eV=None`` + ``comparable=False`` mark the
    state as incomparable; downstream pruning/ranking must keep it."""
    return {
        "state_id": species_label,
        "free_energy_eV": None,
        "adsorption_energy_eV": None,
        "backend": "thermal_uma",
        "convention": "dG_ads",
        "comparable": False,
        "confidence": 0.0,
        "details": {"error": reason},
    }


# ---------------------------------------------------------------------------
# Backend implementations
# ---------------------------------------------------------------------------

class ThermalStateEnergyBackend:
    """Estimate the corrected adsorption free energy ΔG_ads via UMA.

    Wraps the existing ``compute_adsorption_energies`` tool. Returns the
    slab-referenced, gas-corrected ΔG_ads defined in the module docstring,
    or an incomparable ``None`` result when any required quantity (UMA
    tools, surface, gas reference) is missing.
    """

    def __init__(
        self,
        tools: Any = None,
        gas_db: Any = None,
        temperature_K: float = 298.0,
        pressure_Pa: float = 1e5,
    ):
        self._tools = tools
        self._gas_db = gas_db
        self.temperature_K = temperature_K
        self.pressure_Pa = pressure_Pa

    def _gas_G_for_elements(self, elements: Dict[str, int]) -> Optional[float]:
        """Gas-phase chemical potential of ``elements`` at (T, P), resolved
        through the molecule DB (whole-molecule G preferred, element-
        reference decomposition fallback — see module docstring)."""
        if self._gas_db is None:
            return None
        formula = {e: n for e, n in elements.items() if n}
        if not formula:
            return None
        try:
            return float(self._gas_db.mu_for_delta(
                formula, T=self.temperature_K, P=self.pressure_Pa,
            ))
        except Exception as exc:
            logger.warning("Gas reference lookup failed for %s: %s", formula, exc)
            return None

    def estimate(
        self,
        species_label: str,
        surface_path: Optional[str] = None,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """Return a dict matching StateEnergyEstimate fields.

        ``free_energy_eV`` is ΔG_ads(X) or None (incomparable).
        """
        if self._tools is None or surface_path is None:
            return _unavailable_estimate(
                species_label, "no UMA tools or surface available",
            )

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
            e_ads = result.get("adsorption_energies", {}).get(species_label)
            if e_ads is None:
                return _unavailable_estimate(
                    species_label, "UMA returned no adsorption energy",
                )

            from core.pathway.candidate_generators import parse_species_elements
            elements = parse_species_elements(species_label)
            missing = [e for e in elements if e not in ATOMIC_REF_OC20]
            if missing:
                return _unavailable_estimate(
                    species_label,
                    f"no OC20 atomic reference for element(s) {missing}",
                )
            ref_sum = sum(
                n * ATOMIC_REF_OC20[e] for e, n in elements.items()
            )
            gas_G = self._gas_G_for_elements(elements)
            if gas_G is None:
                return _unavailable_estimate(
                    species_label,
                    "no gas-phase free-energy reference in CatDT molecule DB",
                )

            # ΔG_ads = E(slab+X) − E(slab) − G_gas(X)
            #        = E_ads_atomic + Σ atomic_ref − G_gas(X)
            dg_ads = float(e_ads) + ref_sum - gas_G
            return {
                "state_id": species_label,
                "free_energy_eV": dg_ads,
                "adsorption_energy_eV": float(e_ads),
                "backend": "thermal_uma",
                "convention": "dG_ads",
                "comparable": True,
                "confidence": 0.7,
                "details": {
                    "source": "uma_dG_ads",
                    "e_ads_atomic_eV": float(e_ads),
                    "atomic_ref_sum_eV": ref_sum,
                    "gas_G_eV": gas_G,
                    "temperature_K": self.temperature_K,
                    "pressure_Pa": self.pressure_Pa,
                },
            }
        except Exception as exc:
            logger.warning("Thermal UMA estimate failed for %s: %s", species_label, exc)
            return _unavailable_estimate(species_label, f"UMA estimate failed: {exc}")


class ElectroStateEnergyBackend:
    """Computational Hydrogen Electrode (CHE) free energy estimate.

    Implements (only) the CHE electrochemical shift on a caller-supplied
    adsorption energy:

        ΔG = ΔE_ads + n·e·U + n · (2.303 · k_B · T / e) · pH

    i.e. the electron-transfer term plus the Nernstian pH shift evaluated
    at ``temperature_K``. ZPE / −TΔS vibrational corrections are NOT
    applied here — they are expected to already be contained in the input
    (e.g. a ΔG_ads value from ``ThermalStateEnergyBackend``).
    """

    # Boltzmann constant in eV/K
    _KB_EV = 8.617333262e-5

    def __init__(
        self,
        voltage_V: float = 0.0,
        pH: float = 0.0,
        temperature_K: float = 298.15,
    ):
        self.voltage_V = voltage_V
        self.pH = pH
        self.temperature_K = temperature_K

    @property
    def nernst_pH_coeff_eV(self) -> float:
        """2.303 · k_B · T / e in eV per pH unit (0.0592 at 298.15 K)."""
        import math
        return math.log(10.0) * self._KB_EV * self.temperature_K

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
                "comparable": False,
                "confidence": 0.0,
                "details": {"error": "no adsorption energy provided"},
            }

        # CHE correction: ΔG = ΔE + neU + (2.303 kB T / e) · pH · n
        pH_coeff = self.nernst_pH_coeff_eV
        che_correction = n_electrons * self.voltage_V + pH_coeff * self.pH * n_electrons
        free_energy = adsorption_energy_eV + che_correction

        return {
            "state_id": species_label,
            "free_energy_eV": float(free_energy),
            "adsorption_energy_eV": float(adsorption_energy_eV),
            "backend": "electro_che",
            "comparable": True,
            "confidence": 0.5,
            "details": {
                "voltage_V": self.voltage_V,
                "pH": self.pH,
                "temperature_K": self.temperature_K,
                "nernst_pH_coeff_eV": pH_coeff,
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
        parent_gas_G: Optional[float] = None,
        child_gas_G: Optional[float] = None,
    ) -> Optional[float]:
        """Step free energy for an elementary CRN step parent → child.

        ``parent_G`` / ``child_G`` are corrected adsorption free energies
        ΔG_ads (module-docstring convention). The exact step free energy is

            ΔG_step = ΔG_ads(child) − ΔG_ads(parent)
                      + [G_gas(child) − G_gas(parent) − Σ_e Δn_e · μ_e(T, P)]

        The bracketed gas-side term re-attaches the gas references removed
        in ΔG_ads and charges/credits the chemical potential of gas species
        exchanged with the reservoir (``Δn_e = n_e(child) − n_e(parent)``,
        μ resolved through ``gas_db.mu_for_delta``, which prefers whole
        gas molecules and falls back to element references). It equals the
        gas-phase reaction free energy "parent(g) + gas sources →
        child(g)" and is exactly 0 for isomerisation / adsorption /
        desorption steps where both sides share one formula.

        ``parent_gas_G`` / ``child_gas_G`` let composite-aware callers
        (e.g. co-adsorbed "*A+*B" states) supply pre-summed gas references;
        by default each side is resolved from its element dict via
        ``gas_db.mu_for_delta`` — the SAME resolver used when the ΔG_ads
        values were built, so the reference cancels exactly.

        Returns ``None`` if either ΔG_ads, either gas reference, or any
        required element chemical potential cannot be resolved — the step
        is then INCOMPARABLE and must not be ranked against numeric values.
        """
        if parent_G is None or child_G is None:
            return None

        p_formula = {e: n for e, n in parent_elements.items() if n}
        c_formula = {e: n for e, n in child_elements.items() if n}

        # Same formula on both sides → gas terms cancel exactly (the
        # resolver is deterministic) and Δn = 0.
        if p_formula == c_formula and parent_gas_G is None and child_gas_G is None:
            return child_G - parent_G

        if gas_db is None:
            return None

        def _resolve_gas_G(formula: Dict[str, int], override: Optional[float]) -> Optional[float]:
            if override is not None:
                return override
            if not formula:
                return 0.0
            try:
                return float(gas_db.mu_for_delta(formula, T=T, P=P))
            except Exception:
                return None

        pgg = _resolve_gas_G(p_formula, parent_gas_G)
        cgg = _resolve_gas_G(c_formula, child_gas_G)
        if pgg is None or cgg is None:
            return None

        delta = {
            e: c_formula.get(e, 0) - p_formula.get(e, 0)
            for e in set(p_formula) | set(c_formula)
        }
        delta = {e: n for e, n in delta.items() if n != 0}
        try:
            mu_correction = gas_db.mu_for_delta(delta, T=T, P=P)
        except KeyError:
            # Missing gas-phase reference for some element — the step
            # cannot be expressed on the shared scale. Caller logs.
            return None
        return (child_G - parent_G) + (cgg - pgg) - mu_correction

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
    "ATOMIC_REF_OC20",
    "ThermalStateEnergyBackend",
    "ElectroStateEnergyBackend",
    "FreeEnergyRouter",
]
