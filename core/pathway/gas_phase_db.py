"""Gas-phase free-energy database for UniMech pathway evaluation.

Stores per-molecule electronic + vibrational data computed from a
UMA-relaxed gas-phase structure in a ``box_size`` (default 30 Å) cubic
vacuum box. Free energies are obtained on-the-fly via
``ase.thermochemistry.IdealGasThermo`` at any requested (T, P) — no
discrete grid, fully continuous in both arguments.

Each record stores everything ``IdealGasThermo`` needs:
  - ``E_elec_eV``           : potential energy from UMA relaxation
  - ``vib_energies_eV``     : list of real-mode vibrational quanta (eV)
  - ``geometry``            : 'monatomic' / 'linear' / 'nonlinear'
  - ``symmetry_number``     : rotational σ
  - ``spin``                : number of unpaired electrons
  - ``atom_symbols``        : list[str]
  - ``atom_positions_A``    : list[[x,y,z]] in Å
  - ``formula``             : {element: count}
  - ``smiles``, ``molecular_mass_amu`` : metadata

Element chemical potentials μ_ref(e, T, P) are derived through
``DEFAULT_ELEMENT_REFERENCES`` (overridable via ``element_refs``).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Element-reference mapping
# ---------------------------------------------------------------------------
DEFAULT_ELEMENT_REFERENCES: Dict[str, List[Tuple[str, float]]] = {
    "H":  [("H2", 0.5)],                   # ½ H₂(g)
    "O":  [("H2O", 1.0), ("H2", -1.0)],    # H₂O − H₂  (avoids UMA's O₂ over-binding)
    "C":  [("CO", 1.0)],                   # CO(g) — default for FT/methanation
    "N":  [("N2", 0.5)],                   # ½ N₂(g)
    "S":  [("H2S", 1.0), ("H2", -1.0)],    # H₂S − H₂
    "F":  [("HF", 1.0), ("H2", -0.5)],     # HF − ½ H₂
    "Cl": [("HCl", 1.0), ("H2", -0.5)],    # HCl − ½ H₂
}


# ---------------------------------------------------------------------------
# DB class
# ---------------------------------------------------------------------------
class GasPhaseFreeEnergyDB:
    """Read-side API for the gas-phase free-energy database.

    Construction is cheap (loads JSON into memory). ``G(name, T, P)``
    reconstructs an ``IdealGasThermo`` object from the stored data and
    evaluates Gibbs free energy at the requested point — continuous in
    both T and P. Reconstruction is lazy + memoised per name.
    """

    def __init__(
        self,
        db_path: Path | str,
        element_refs: Optional[Dict[str, List[Tuple[str, float]]]] = None,
    ):
        self.db_path = Path(db_path)
        if self.db_path.exists():
            with open(self.db_path) as f:
                raw = json.load(f)
            self.metadata = raw.get("metadata", {})
            self.molecules: Dict[str, Dict[str, Any]] = raw.get("molecules", {})
        else:
            logger.info("Gas-phase DB empty at %s — initializing.", self.db_path)
            self.metadata = {}
            self.molecules = {}
        self.element_refs = element_refs or DEFAULT_ELEMENT_REFERENCES
        self._thermo_cache: Dict[str, Any] = {}   # name → IdealGasThermo

    # ------------------------------------------------------------------ I/O

    def save(self) -> None:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"metadata": self.metadata, "molecules": self.molecules}
        with open(self.db_path, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
        logger.info(
            "Gas-phase DB saved: %s (%d molecules)",
            self.db_path, len(self.molecules),
        )

    def has(self, name: str) -> bool:
        return name in self.molecules

    def upsert(self, name: str, record: Dict[str, Any]) -> None:
        self.molecules[name] = record
        self._thermo_cache.pop(name, None)

    # ------------------------------------------------------------- Thermo

    def _build_thermo(self, name: str) -> Any:
        """Reconstruct an IdealGasThermo object from the stored record."""
        from ase import Atoms
        from ase.thermochemistry import IdealGasThermo
        m = self.molecules[name]
        symbols = m.get("atom_symbols", [])
        positions = m.get("atom_positions_A", [])
        if not symbols or not positions:
            raise ValueError(
                f"Molecule '{name}' lacks atom_symbols / atom_positions_A — "
                "rebuild via scripts/build_gas_phase_db.py."
            )
        atoms = Atoms(symbols=symbols, positions=positions)
        return IdealGasThermo(
            vib_energies=m.get("vib_energies_eV", []),
            potentialenergy=m["E_elec_eV"],
            atoms=atoms,
            geometry=m.get("geometry", "nonlinear"),
            symmetrynumber=m.get("symmetry_number", 1),
            spin=m.get("spin", 0),
        )

    def _get_thermo(self, name: str) -> Any:
        if name not in self.molecules:
            raise KeyError(f"Gas-phase DB missing '{name}'.")
        if name not in self._thermo_cache:
            self._thermo_cache[name] = self._build_thermo(name)
        return self._thermo_cache[name]

    # -------------------------------------------------------------- Queries

    def G(self, name: str, T: float = 298.0, P: float = 1e5) -> float:
        """Continuous Gibbs free energy of ``name`` (eV) at given T (K), P (Pa)."""
        thermo = self._get_thermo(name)
        return float(thermo.get_gibbs_energy(temperature=T, pressure=P, verbose=False))

    def thermal_correction(self, name: str, T: float = 298.0, P: float = 1e5) -> float:
        """G(T,P) − E_elec — the thermal/entropic correction added to E_elec."""
        return self.G(name, T=T, P=P) - self.molecules[name]["E_elec_eV"]

    def mu_element(
        self,
        element: str,
        T: float = 298.0,
        P: float = 1e5,
    ) -> float:
        """Chemical potential of an element (eV/atom) via reference reactions."""
        if element not in self.element_refs:
            raise KeyError(f"No reference reaction for element '{element}'.")
        return sum(
            coeff * self.G(species, T=T, P=P)
            for species, coeff in self.element_refs[element]
        )

    def mu_for_delta(
        self,
        element_delta: Dict[str, int],
        T: float = 298.0,
        P: float = 1e5,
    ) -> float:
        """Σ Δn_e × μ_ref(e). For ΔG_step bookkeeping."""
        return sum(
            n * self.mu_element(elem, T=T, P=P)
            for elem, n in element_delta.items() if n != 0
        )

    def species_for_formula(self, formula: Dict[str, int]) -> Optional[str]:
        for name, m in self.molecules.items():
            if m.get("formula") == formula:
                return name
        return None


# ---------------------------------------------------------------------------
# Default-path helper
# ---------------------------------------------------------------------------
def default_db_path() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "gas_phase_db.json"


def load_default(
    element_refs: Optional[Dict[str, List[Tuple[str, float]]]] = None,
) -> GasPhaseFreeEnergyDB:
    return GasPhaseFreeEnergyDB(default_db_path(), element_refs=element_refs)


__all__ = [
    "GasPhaseFreeEnergyDB",
    "DEFAULT_ELEMENT_REFERENCES",
    "default_db_path",
    "load_default",
]
