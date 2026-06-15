"""CatDT unified molecule database.

A single source of truth merging:
  1. Fairchem's ``ADSORBATE_PKL`` — 475 adsorbate / surface-radical
     structures (ASE atoms + binding-atom indices) used for placing
     adsorbates on slab surfaces.
  2. The gas-phase free-energy table — stable molecules with electronic
     energy + vibrational frequencies + IdealGasThermo data, used for
     stoichiometry-corrected ΔG bookkeeping in mechanism search.

Each record stores everything ``IdealGasThermo`` (or ``HarmonicThermo``)
needs to evaluate Gibbs free energy at any (T, P), so G(T, P) is a
**continuous function** — no discrete grid.

Schema (per molecule, keyed by canonical name):

    {
      "smiles": "CC",
      "formula": {"C": 2, "H": 6},
      "geometry": "nonlinear",            # 'monatomic' / 'linear' / 'nonlinear'
      "symmetry_number": 6,
      "spin": 0,
      "molecular_mass_amu": 30.07,

      # Structure (UMA-relaxed by default; falls back to fairchem geom if
      # relaxation unavailable).
      "atom_symbols": ["C","C","H","H",...],
      "atom_positions_A": [[x,y,z], ...],
      "binding_indices": [0],             # for adsorbate placement; from fairchem

      # Electronic + vibrational data — only present if UMA computed it.
      "E_elec_eV": -40.218,
      "vib_energies_eV": [...],
      "frequencies_cm-1": [...],
      "ZPE_eV": 1.97,

      # Provenance.
      "source": "fairchem|gas|both"
    }
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
    # Self-consistent elemental chemical potentials. Each μ(e) is
    # G(reference_compound) minus all other elements' μ contributions, so
    # Σ_e n_e × μ(e) = G(compound) for every reference compound.
    "H":  [("H2", 0.5)],                                # ½ H₂(g)
    "N":  [("N2", 0.5)],                                # ½ N₂(g)
    "O":  [("H2O", 1.0), ("H2", -1.0)],                 # H₂O − H₂  (avoids UMA's O₂ over-binding)
    "C":  [("CO", 1.0), ("H2O", -1.0), ("H2", 1.0)],    # CO − μ(O) = CO − H₂O + H₂  (FT default: C source = CO)
    "S":  [("H2S", 1.0), ("H2", -1.0)],                 # H₂S − H₂
    "F":  [("HF", 1.0), ("H2", -0.5)],                  # HF − ½ H₂
    "Cl": [("HCl", 1.0), ("H2", -0.5)],                 # HCl − ½ H₂
}


# ---------------------------------------------------------------------------
# DB class
# ---------------------------------------------------------------------------
class CatDTMoleculeDB:
    """Unified read-side API for the CatDT molecule database.

    Construction is cheap (loads JSON into memory). ``G(name, T, P)``
    reconstructs an ``IdealGasThermo`` object from stored data and
    evaluates Gibbs free energy at the requested (T, P) — continuous in
    both arguments, lazy + memoised per name.

    For records with only structure (fairchem adsorbates without a
    UMA-computed thermo entry), ``G(...)`` raises ``ValueError``;
    ``has_thermo(name)`` checks availability first.
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
            logger.info("CatDT molecule DB empty at %s — initializing.", self.db_path)
            self.metadata = {}
            self.molecules = {}
        self.element_refs = element_refs or DEFAULT_ELEMENT_REFERENCES
        self._thermo_cache: Dict[str, Any] = {}

    # ------------------------------------------------------------------ I/O

    def save(self) -> None:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"metadata": self.metadata, "molecules": self.molecules}
        with open(self.db_path, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
        logger.info(
            "CatDT molecule DB saved: %s (%d molecules)",
            self.db_path, len(self.molecules),
        )

    @staticmethod
    def _candidate_keys(name: str) -> List[str]:
        """Generate lookup variants so callers don't need to track the
        '*' prefix convention. Order matters: prefer exact match, then
        starred fairchem form, then bare-stripped form."""
        out = [name]
        if name.startswith("*"):
            out.append(name[1:])
        else:
            out.append("*" + name)
        return out

    def _resolve(self, name: str) -> Optional[str]:
        for k in self._candidate_keys(name):
            if k in self.molecules:
                return k
        return None

    def has(self, name: str) -> bool:
        return self._resolve(name) is not None

    def has_thermo(self, name: str) -> bool:
        k = self._resolve(name)
        if k is None:
            return False
        m = self.molecules[k]
        return (
            "E_elec_eV" in m
            and "vib_energies_eV" in m
            and m.get("atom_symbols")
            and m.get("atom_positions_A")
        )

    def upsert(self, name: str, record: Dict[str, Any]) -> None:
        existing = self.molecules.get(name, {})
        existing.update(record)
        self.molecules[name] = existing
        self._thermo_cache.pop(name, None)

    # -------------------------------------------------------- Structure API

    def get_atoms(self, name: str):
        """Return an ASE Atoms object (positions + symbols) for adsorbate
        placement / NEB seed structures."""
        from ase import Atoms
        k = self._resolve(name)
        if k is None:
            raise KeyError(f"DB missing '{name}'.")
        m = self.molecules[k]
        symbols = m.get("atom_symbols")
        positions = m.get("atom_positions_A")
        if not symbols or not positions:
            raise ValueError(f"Molecule '{name}' has no structural data.")
        return Atoms(symbols=symbols, positions=positions)

    def get_binding_indices(self, name: str) -> List[int]:
        k = self._resolve(name)
        return list(self.molecules.get(k, {}).get("binding_indices", []))

    # ---------------------------------------------------------- Thermo API

    def _build_thermo(self, key: str) -> Any:
        from ase import Atoms
        from ase.thermochemistry import IdealGasThermo
        m = self.molecules[key]
        atoms = Atoms(
            symbols=m["atom_symbols"],
            positions=m["atom_positions_A"],
        )
        return IdealGasThermo(
            vib_energies=m["vib_energies_eV"],
            potentialenergy=m["E_elec_eV"],
            atoms=atoms,
            geometry=m.get("geometry", "nonlinear"),
            symmetrynumber=m.get("symmetry_number", 1),
            spin=m.get("spin", 0),
        )

    def _get_thermo(self, name: str) -> Any:
        key = self._resolve(name)
        if key is None:
            raise KeyError(f"DB missing '{name}'.")
        if not self.has_thermo(key):
            raise ValueError(f"'{name}' has no thermo data — rebuild required.")
        if key not in self._thermo_cache:
            self._thermo_cache[key] = self._build_thermo(key)
        return self._thermo_cache[key]

    def G(self, name: str, T: float = 298.0, P: float = 1e5) -> float:
        """Continuous Gibbs free energy (eV) at any (T, P)."""
        thermo = self._get_thermo(name)
        return float(thermo.get_gibbs_energy(
            temperature=T, pressure=P, verbose=False,
        ))

    def thermal_correction(self, name: str, T: float = 298.0, P: float = 1e5) -> float:
        key = self._resolve(name)
        if key is None:
            raise KeyError(f"DB missing '{name}'.")
        return self.G(name, T=T, P=P) - self.molecules[key]["E_elec_eV"]

    def mu_element(
        self, element: str, T: float = 298.0, P: float = 1e5,
    ) -> float:
        if element not in self.element_refs:
            raise KeyError(f"No reference reaction for element '{element}'.")
        return sum(
            coeff * self.G(species, T=T, P=P)
            for species, coeff in self.element_refs[element]
        )

    def _atomic_mu_for_delta(
        self,
        element_delta: Dict[str, int],
        T: float = 298.0,
        P: float = 1e5,
    ) -> float:
        """Element-decomposition fallback (legacy behavior)."""
        return sum(
            n * self.mu_element(elem, T=T, P=P)
            for elem, n in element_delta.items() if n != 0
        )

    def _match_single_molecule(
        self, formula: Dict[str, int],
    ) -> Optional[Tuple[str, float]]:
        """If ``formula = c × stored_formula`` for some molecule in the DB
        (any source, c rational > 0), return ``(name, c)``. Else None.

        The whole 569-molecule DB is searched — including fairchem
        radicals like *CH, *CHO, *OH, *CH₂. This lets transfer-deltas
        with no closed-shell stable gas analog (e.g. Δn={C:1,H:1}) still
        get a proper free-energy chemical potential from the matching
        radical's UMA-relaxed gas G (with ZPE + vib entropy), rather
        than collapsing to per-element decomposition.

        Preference order: closed-shell over open-shell, stable
        ('gas'/'both' source) over fairchem radical, integer scale over
        fractional, smaller scale, and finally larger molecule (= more
        elementary single-molecule transfer).
        """
        formula = {e: n for e, n in formula.items() if n != 0}
        if not formula:
            return None
        best: Optional[Tuple[str, float]] = None
        best_score: Optional[Tuple[int, int, int, float, int]] = None
        for name, m in self.molecules.items():
            if not self.has_thermo(name):
                continue
            F = {e: n for e, n in m.get("formula", {}).items() if n != 0}
            if set(F.keys()) != set(formula.keys()):
                continue
            ratios = [formula[e] / F[e] for e in formula]
            ref = ratios[0]
            if ref <= 0 or any(abs(r - ref) > 1e-6 for r in ratios):
                continue
            atom_count = sum(F.values())
            spin = int(m.get("spin", 0))
            source = m.get("source", "fairchem")
            # Priority (lower wins):
            #   1) integer scale  > fractional (transferring whole molecule
            #      is more physical than fraction of one).
            #   2) smaller scale  (1×CH > 6×... > etc.)
            #   3) closed-shell   over open-shell.
            #   4) stable source  ('gas'/'both') over fairchem-only radical.
            #   5) larger molecule per scale unit (more elementary per atom).
            int_score = 0 if abs(ref - round(ref)) < 1e-6 else (
                1 if abs(2 * ref - round(2 * ref)) < 1e-6 else 2
            )
            source_score = 0 if source in ("gas", "both") else 1
            score = (int_score, ref, spin, source_score, -atom_count)
            if best is None or score < best_score:
                best = (name, float(ref))
                best_score = score
        return best

    def mu_for_delta(
        self,
        element_delta: Dict[str, int],
        T: float = 298.0,
        P: float = 1e5,
    ) -> float:
        """Stoichiometry-aware μ correction for an elementary step.

        Fully general: tries to identify the actual gas-phase molecule
        being consumed / produced from the step's element delta, falling
        back to per-element atom-reference decomposition only when the
        delta has no single-molecule analog (mixed-sign cases or
        radical residues).

        Recognized cases:
          Δn={H:2}    → consumes 1 H₂(g)
          Δn={H:1}    → consumes 0.5 H₂(g)
          Δn={C:1,O:1}→ consumes 1 CO(g)
          Δn={C:1,H:4}→ consumes 1 CH₄(g)
          Δn={H:2,O:1}→ consumes 1 H₂O(g)
          Δn={H:-2,O:-1}→ produces 1 H₂O(g)
          (any other formula matching a stored molecule, with rational
          scale)

        Otherwise falls back to ``_atomic_mu_for_delta``.
        """
        delta = {e: n for e, n in element_delta.items() if n != 0}
        if not delta:
            return 0.0

        pos = {e: n for e, n in delta.items() if n > 0}
        neg = {e: -n for e, n in delta.items() if n < 0}

        # Pure-sign delta: try a single-molecule match first.
        if pos and not neg:
            match = self._match_single_molecule(pos)
            if match is not None:
                name, scale = match
                return scale * self.G(name, T=T, P=P)
        if neg and not pos:
            match = self._match_single_molecule(neg)
            if match is not None:
                name, scale = match
                return -scale * self.G(name, T=T, P=P)

        # Mixed-sign delta: try to identify pos and neg halves separately;
        # any residue falls back to atom decomposition.
        mu = 0.0
        residue: Dict[str, int] = {}
        if pos:
            match = self._match_single_molecule(pos)
            if match is not None:
                name, scale = match
                mu += scale * self.G(name, T=T, P=P)
            else:
                for e, n in pos.items():
                    residue[e] = residue.get(e, 0) + n
        if neg:
            match = self._match_single_molecule(neg)
            if match is not None:
                name, scale = match
                mu -= scale * self.G(name, T=T, P=P)
            else:
                for e, n in neg.items():
                    residue[e] = residue.get(e, 0) - n

        if residue:
            mu += self._atomic_mu_for_delta(residue, T=T, P=P)
        return mu

    def species_for_formula(
        self,
        formula: Dict[str, int],
        require_thermo: bool = True,
    ) -> Optional[str]:
        """Find a stored molecule whose formula matches.
        When ``require_thermo`` is True, only entries with thermal data
        are eligible (so callers don't trip on radicals lacking it).
        Prefers stable closed-shell entries (spin=0) when multiple match."""
        candidates: List[str] = []
        for name, m in self.molecules.items():
            if m.get("formula") != formula:
                continue
            if require_thermo and not self.has_thermo(name):
                continue
            candidates.append(name)
        if not candidates:
            return None
        # Prefer closed-shell + 'gas' / 'both' source over fairchem radicals.
        candidates.sort(key=lambda n: (
            self.molecules[n].get("spin", 0),
            0 if self.molecules[n].get("source") in ("gas", "both") else 1,
        ))
        return candidates[0]


# ---------------------------------------------------------------------------
# Default-path helper
# ---------------------------------------------------------------------------
def default_db_path() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "catdt_molecule_db.json"


def load_default(
    element_refs: Optional[Dict[str, List[Tuple[str, float]]]] = None,
) -> CatDTMoleculeDB:
    return CatDTMoleculeDB(default_db_path(), element_refs=element_refs)


__all__ = [
    "CatDTMoleculeDB",
    "DEFAULT_ELEMENT_REFERENCES",
    "default_db_path",
    "load_default",
]
