"""
UMA Surface Calculator Adapters for VSSR-MC.

Wraps FairChem's FAIRChemCalculator (UMA model) to work with VSSR-MC's
SurfaceSystem, which requires a `surface_energy` property.

Two modes:
- UMASurfaceCalculator: thermal (gas-solid), same formula as EnsembleNFFSurface
- UMAPourbaixCalculator: electrochemical (liquid-solid), same formula as NFFPourbaix
"""

import logging
from collections import Counter
from typing import Dict, List, Optional

import ase
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

logger = logging.getLogger(__name__)

KB_EV = 8.617333262e-5  # eV/K


class UMASurfaceCalculator(Calculator):
    """Wraps any ASE Calculator (e.g. UMA FAIRChemCalculator) with surface energy.

    Surface energy formula (same as EnsembleNFFSurface, Du et al. 2023):
        E_surf = E_total - N_ref × E_bulk_ref
                 - Σ_i (excess_i × E_bulk_i)
                 - Σ_i (excess_i × μ_i)
    """

    implemented_properties = ("energy", "forces", "surface_energy")

    def __init__(
        self,
        base_calculator: Calculator,
        offset_data: Optional[Dict] = None,
        chem_pots: Optional[Dict[str, float]] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.base = base_calculator
        self.offset_data = offset_data or {}
        self.chem_pots = chem_pots or {}
        self._logger = logging.getLogger(__name__)

    def set(self, **kwargs) -> dict:
        changed = {}
        if "chem_pots" in kwargs:
            self.chem_pots = kwargs.pop("chem_pots")
            changed["chem_pots"] = self.chem_pots
            self._logger.info("chemical potentials: %s are set from parameters", self.chem_pots)
        if "offset_data" in kwargs:
            self.offset_data = kwargs.pop("offset_data")
            changed["offset_data"] = self.offset_data
            self._logger.info("offset data: %s is set from parameters", self.offset_data)
        # Store remaining kwargs in parameters dict for compatibility
        self.parameters.update(kwargs)
        self.parameters.update(changed)
        return changed

    def calculate(
        self,
        atoms: ase.Atoms = None,
        properties=("energy", "forces", "surface_energy"),
        system_changes=all_changes,
    ):
        if atoms is None:
            atoms = self.atoms

        Calculator.calculate(self, atoms, properties, system_changes)

        # Delegate energy/forces to base calculator
        base_props = [p for p in properties if p in ("energy", "forces")]
        if base_props:
            self.base.calculate(atoms, base_props, system_changes)
            self.results.update(self.base.results)

        if "surface_energy" in properties:
            self.results["surface_energy"] = self._compute_surface_energy(atoms)

        if hasattr(atoms, "results"):
            atoms.results.update(self.results)

    def get_potential_energy(self, atoms=None, **kwargs):
        if atoms is None:
            atoms = self.atoms
        # Follow ASE Calculator semantics: recompute when the atoms changed
        # since the last calculation instead of blindly reusing cached results.
        system_changes = self.check_state(atoms)
        if system_changes or "energy" not in self.results:
            self.calculate(atoms, ["energy"], system_changes or all_changes)
        return float(self.results["energy"])

    def _compute_surface_energy(self, atoms: ase.Atoms) -> float:
        """Same formula as EnsembleNFFSurface.get_surface_energy()."""
        e_total = self.get_potential_energy(atoms)

        if not self.offset_data:
            return e_total

        ads_count = Counter(atoms.get_chemical_symbols())
        bulk_energies = self.offset_data.get("bulk_energies", {})
        stoics = self.offset_data.get("stoics", {})
        ref_formula = self.offset_data.get("ref_formula", "")
        ref_element = self.offset_data.get("ref_element", "")

        if not ref_element or ref_element not in ads_count:
            return e_total

        # Subtract bulk reference energy
        bulk_ref_en = ads_count[ref_element] * bulk_energies.get(ref_formula, 0.0)
        for ele in ads_count:
            if ele != ref_element:
                s_ele = stoics.get(ele, 0)
                s_ref = stoics.get(ref_element, 1)
                excess = ads_count[ele] - s_ele / s_ref * ads_count[ref_element]
                bulk_ref_en += excess * bulk_energies.get(ele, 0.0)

        surface_energy = e_total - bulk_ref_en

        # Subtract chemical potential deviation
        for ele in ads_count:
            if ele != ref_element and ele in self.chem_pots:
                s_ele = stoics.get(ele, 0)
                s_ref = stoics.get(ref_element, 1)
                excess = ads_count[ele] - s_ele / s_ref * ads_count[ref_element]
                surface_energy -= excess * self.chem_pots[ele]

        return float(surface_energy)


class UMAPourbaixCalculator(UMASurfaceCalculator):
    """Wraps UMA with Pourbaix grand potential for electrochemical MC.

    Exact same formulas as NFFPourbaix (mcmc/calculators/calculators.py),
    including adsorbate corrections (OH ZPE-TS).

    Accepts PourbaixAtom objects (from mcmc.pourbaix.atoms) or dicts with
    the same keys (atom_std_state_energy, num_e, num_H, delta_G2_std, species_conc).
    """

    implemented_properties = ("energy", "forces", "surface_energy")

    def __init__(
        self,
        base_calculator: Calculator,
        pourbaix_atoms: Optional[Dict] = None,
        potential_she: float = 0.0,
        ph: float = 7.0,
        temperature_k: float = 300.0,
        adsorbate_corrections: Optional[Dict[str, float]] = None,
        **kwargs,
    ):
        super().__init__(base_calculator, **kwargs)
        self.pourbaix_atoms = pourbaix_atoms or {}
        self.phi = potential_she
        self.pH = ph
        self.temp = KB_EV * temperature_k  # kT in eV
        self.adsorbate_corrections = adsorbate_corrections or {}

    def set(self, **kwargs) -> dict:
        changed = super().set(**kwargs)
        if "phi" in kwargs:
            self.phi = kwargs["phi"]
            self._logger.info("potential: %.3f V vs SHE", self.phi)
        if "pH" in kwargs:
            self.pH = kwargs["pH"]
            self._logger.info("pH: %.1f", self.pH)
        if "temperature" in kwargs:
            self.temp = kwargs["temperature"]
            self._logger.info("temperature: %.4f kT", self.temp)
        if "pourbaix_atoms" in kwargs:
            self.pourbaix_atoms = kwargs["pourbaix_atoms"]
            self._logger.info("Pourbaix atoms: %s", list(self.pourbaix_atoms.keys()))
        if "adsorbate_corrections" in kwargs:
            self.adsorbate_corrections = kwargs["adsorbate_corrections"]
        return changed

    @staticmethod
    def _pa_attr(pa, key, default=0.0):
        """Get attribute from PourbaixAtom object or dict."""
        if hasattr(pa, key):
            return getattr(pa, key)
        if isinstance(pa, dict):
            return pa.get(key, default)
        return default

    def _compute_surface_energy(self, atoms: ase.Atoms) -> float:
        """Pourbaix grand potential = -(ΔG₁ + ΔG₂)."""
        return -(self._get_delta_g1(atoms) + self._get_delta_g2(atoms))

    def _get_delta_g1(self, atoms: ase.Atoms) -> float:
        """Dissociation free energy (same as NFFPourbaix.get_delta_G1).

        ΔG₁ = Σ(count × μ_atom_std_state) - E_slab + adsorbate_corrections
        """
        e_slab = self.get_potential_energy(atoms)

        # Sum standard-state chemical potentials
        atoms_count = Counter(atoms.get_chemical_symbols())
        sum_chem_pots = 0.0
        for atom, count in atoms_count.items():
            if atom in self.pourbaix_atoms:
                sum_chem_pots += count * self._pa_attr(
                    self.pourbaix_atoms[atom], "atom_std_state_energy"
                )

        # Adsorbate corrections (OH ZPE-TS), same logic as NFFPourbaix
        try:
            from ase.formula import Formula
        except ImportError as exc:
            # No pymatgen fallback: Composition does not implement the ase
            # Formula API (count/divmod/from_dict) used below.
            raise ImportError(
                "ase.formula.Formula is required for adsorbate corrections "
                "in UMAPourbaixCalculator._get_delta_g1"
            ) from exc
        if self.adsorbate_corrections:
            formula = Formula(atoms.get_chemical_formula())
            for adsorbate, correction in self.adsorbate_corrections.items():
                if "O" in adsorbate and "H" in adsorbate:
                    HO_diff = max(formula["H"] - formula["O"], 0)
                    if HO_diff > 0:
                        from ase.formula import Formula as F
                        formula_dict_to_subtract = (F("H2O") * HO_diff).count()
                        formula_dict = formula.count()
                        formula_dict = {
                            k: formula_dict[k] - formula_dict_to_subtract.get(k, 0)
                            for k in formula_dict
                        }
                        formula = F.from_dict(formula_dict)
                div, _ = divmod(formula, adsorbate)
                e_slab += div * correction

        return sum_chem_pots - e_slab

    def _get_delta_g2(self, atoms: ase.Atoms) -> float:
        """Electrochemical free energy (same as NFFPourbaix.get_delta_G2).

        ΔG₂ = Σ_atoms(ΔG₂_std - n_e·U - 2.3·n_H·kT·pH + kT·ln(a))
        """
        delta_g2 = 0.0
        for symbol in atoms.get_chemical_symbols():
            if symbol in self.pourbaix_atoms:
                pa = self.pourbaix_atoms[symbol]
                n_e = self._pa_attr(pa, "num_e", 0)
                n_h = self._pa_attr(pa, "num_H", 0)
                dg2_std = self._pa_attr(pa, "delta_G2_std", 0.0)
                conc = self._pa_attr(pa, "species_conc", 1e-6)
                delta_g2_non_std = (
                    -n_e * self.phi
                    - np.log(10) * n_h * self.temp * self.pH
                    + self.temp * np.log(max(conc, 1e-30))
                )
                delta_g2 += dg2_std + delta_g2_non_std
        return delta_g2
