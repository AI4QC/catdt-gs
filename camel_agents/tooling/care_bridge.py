"""CARE bridge layer for CatDT.

Wraps CARE's blueprint expansion and template chains, converting CARE
objects into CatDT's own ``CatalyticStateRecord`` / ``ElementaryStepCandidate``
representations.  CARE is treated as a pluggable backend — if it is
unavailable or the query falls outside its domain, the bridge returns an
explicit rejection rather than crashing.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy CARE imports — fail gracefully
# ---------------------------------------------------------------------------

_CARE_AVAILABLE: Optional[bool] = None


def _check_care_available() -> bool:
    global _CARE_AVAILABLE
    if _CARE_AVAILABLE is not None:
        return _CARE_AVAILABLE
    try:
        import care.crn.utils.blueprint  # noqa: F401
        _CARE_AVAILABLE = True
    except Exception:
        _CARE_AVAILABLE = False
    return _CARE_AVAILABLE


# CARE's hard-coded element space
_CARE_SUPPORTED_ELEMENTS = frozenset({"C", "H", "O", "N"})


# ---------------------------------------------------------------------------
# Domain report
# ---------------------------------------------------------------------------

class CAREDomainReport(BaseModel):
    """Whether the query falls inside CARE's supported domain."""
    in_domain: bool = False
    unsupported_elements: List[str] = Field(default_factory=list)
    reason: str = ""
    care_available: bool = False


# ---------------------------------------------------------------------------
# Bridge class
# ---------------------------------------------------------------------------

class CAREBridge:
    """Adapter between CARE reaction-network generation and CatDT schemas."""

    def __init__(self, num_cpu: int = 1, show_progress: bool = False):
        self.num_cpu = num_cpu
        self.show_progress = show_progress

    # ----- domain check -----

    def check_domain(
        self,
        species_elements: Optional[set] = None,
        species_labels: Optional[List[str]] = None,
    ) -> CAREDomainReport:
        """Check whether the given species fall inside CARE's supported domain."""
        care_available = _check_care_available()

        if species_elements is None:
            species_elements = set()
        if species_labels:
            for label in species_labels:
                # Extract elements from labels like "*CO", "H2(g)", "CH3OH"
                cleaned = label.replace("*", "").replace("(g)", "").replace("(s)", "")
                import re
                found = re.findall(r"[A-Z][a-z]?", cleaned)
                species_elements.update(found)

        # Remove metal surface elements (not part of CARE's species space)
        adsorbate_elements = species_elements - {
            "Cu", "Ag", "Au", "Pt", "Pd", "Ni", "Fe", "Co", "Rh", "Ir", "Ru",
            "Ti", "Zr", "Sr", "Ba", "Mn", "Zn", "Sn", "Al", "Ga", "In",
        }

        unsupported = sorted(adsorbate_elements - _CARE_SUPPORTED_ELEMENTS)

        if not care_available:
            return CAREDomainReport(
                in_domain=False,
                unsupported_elements=unsupported,
                reason="CARE library is not available in this environment",
                care_available=False,
            )

        if unsupported:
            return CAREDomainReport(
                in_domain=False,
                unsupported_elements=unsupported,
                reason=f"Elements {unsupported} are outside CARE's C/H/O/N domain",
                care_available=True,
            )

        return CAREDomainReport(in_domain=True, care_available=True)

    # ----- blueprint generation -----

    def generate_blueprint(
        self,
        seed_species: List[str],
        electro: bool = False,
        additional_rxns: bool = True,
        ncc: Optional[int] = None,
        noc: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Generate a CARE reaction network blueprint from seed species.

        Returns a dict with keys: ``network``, ``intermediates``, ``reactions``,
        ``stats``, or ``error`` on failure.
        """
        if not _check_care_available():
            return {"error": "CARE not available", "network": None}

        try:
            from care.crn.utils.blueprint import gen_blueprint

            network = gen_blueprint(
                cs=seed_species,
                ncc=ncc,
                noc=noc,
                electro=electro,
                additional_rxns=additional_rxns,
                num_cpu=self.num_cpu,
                show_progress=self.show_progress,
            )

            return {
                "network": network,
                "intermediates": self._extract_intermediates(network),
                "reactions": self._extract_reactions(network),
                "stats": {
                    "n_intermediates": sum(1 for n in network.nodes if self._is_intermediate(n)),
                    "n_reactions": sum(1 for n in network.nodes if not self._is_intermediate(n)),
                },
            }
        except Exception as exc:
            logger.warning("CARE blueprint generation failed: %s", exc)
            return {"error": str(exc), "network": None}

    # ----- individual template expansion -----

    def expand_dissociation(self, chemical_space: List[str]) -> Dict[str, Any]:
        """Run CARE dissociation template on a list of species SMILES/formulas."""
        if not _check_care_available():
            return {"error": "CARE not available", "intermediates": {}, "reactions": []}

        try:
            from care.crn.templates.dissociation import dissociate

            intermediates, reactions = dissociate(
                chemical_space=chemical_space,
                ncpus=self.num_cpu,
                show_progress=self.show_progress,
            )
            return {
                "intermediates": {k: self._intermediate_to_dict(v) for k, v in intermediates.items()},
                "reactions": [self._reaction_to_dict(r) for r in reactions],
            }
        except Exception as exc:
            logger.warning("CARE dissociation expansion failed: %s", exc)
            return {"error": str(exc), "intermediates": {}, "reactions": []}

    def expand_adsorption(self, intermediates_dict: Any) -> List[Dict[str, Any]]:
        """Run CARE adsorption template on intermediates."""
        if not _check_care_available():
            return []

        try:
            from care.crn.templates.adsorption import gen_adsorption_reactions

            reactions = gen_adsorption_reactions(
                intermediates=intermediates_dict,
                num_cpu=self.num_cpu,
                show_progress=self.show_progress,
            )
            return [self._reaction_to_dict(r) for r in reactions]
        except Exception as exc:
            logger.warning("CARE adsorption expansion failed: %s", exc)
            return []

    def expand_rearrangement(self, intermediates_dict: Any) -> List[Dict[str, Any]]:
        """Run CARE rearrangement template on intermediates."""
        if not _check_care_available():
            return []

        try:
            from care.crn.templates.rearrengement import gen_rearrangement_reactions

            reactions = gen_rearrangement_reactions(
                intermediates=intermediates_dict,
                num_cpu=self.num_cpu,
                show_progress=self.show_progress,
            )
            return [self._reaction_to_dict(r) for r in reactions]
        except Exception as exc:
            logger.warning("CARE rearrangement expansion failed: %s", exc)
            return []

    def expand_pcet(self, intermediates_dict: Any, reactions_list: Any) -> List[Dict[str, Any]]:
        """Run CARE PCET template on intermediates + existing reactions."""
        if not _check_care_available():
            return []

        try:
            from care.crn.templates.pcet import gen_pcet_reactions

            pcet_reactions = gen_pcet_reactions(
                intermediates=intermediates_dict,
                reactions=reactions_list,
                show_progress=self.show_progress,
            )
            return [self._reaction_to_dict(r) for r in pcet_reactions]
        except Exception as exc:
            logger.warning("CARE PCET expansion failed: %s", exc)
            return []

    # ----- conversion to CatDT schemas -----

    def reaction_network_to_candidates(
        self, network: Any,
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Convert a CARE ReactionNetwork to CatDT state records + step candidates.

        Returns:
            (state_records, step_candidates) — both as list of dicts
            matching CatalyticStateRecord / ElementaryStepCandidate field names.
        """
        state_records = []
        step_candidates = []

        intermediates = self._extract_intermediates(network)
        for key, info in intermediates.items():
            state_records.append({
                "state_id": f"care_{key}",
                "species_label": info.get("formula", key),
                "phase": info.get("phase", "ads"),
                "elements": info.get("elements", {}),
                "extra": {"care_code": key, "care_smiles": info.get("smiles", "")},
            })

        reactions = self._extract_reactions(network)
        for idx, rxn in enumerate(reactions):
            step_candidates.append({
                "step_id": f"care_step_{idx:04d}",
                "reactant_state_id": f"care_{rxn.get('reactant_key', '')}",
                "product_state_id": f"care_{rxn.get('product_key', '')}",
                "operation_type": rxn.get("reaction_type", "other"),
                "bond_changes": rxn.get("bond_changes", ""),
                "source_backend": "care",
                "extra": rxn,
            })

        return state_records, step_candidates

    # ----- internal helpers -----

    @staticmethod
    def _is_intermediate(node: Any) -> bool:
        try:
            from care.crn.intermediate import Intermediate
            return isinstance(node, Intermediate)
        except ImportError:
            return hasattr(node, "formula") and hasattr(node, "phase")

    def _extract_intermediates(self, network: Any) -> Dict[str, Dict[str, Any]]:
        result = {}
        for node in network.nodes:
            if self._is_intermediate(node):
                key = getattr(node, "code", None) or str(id(node))
                result[key] = self._intermediate_to_dict(node)
        return result

    def _extract_reactions(self, network: Any) -> List[Dict[str, Any]]:
        result = []
        for node in network.nodes:
            if not self._is_intermediate(node):
                rxn_dict = self._reaction_to_dict(node)
                result.append(rxn_dict)
        return result

    @staticmethod
    def _intermediate_to_dict(inter: Any) -> Dict[str, Any]:
        formula = getattr(inter, "formula", "") or ""
        phase = getattr(inter, "phase", "ads") or "ads"
        smiles = getattr(inter, "smiles", "") or ""
        code = getattr(inter, "code", "") or ""

        # Parse element counts from formula
        import re
        elements: Dict[str, int] = {}
        for match in re.finditer(r"([A-Z][a-z]?)(\d*)", formula):
            elem = match.group(1)
            count = int(match.group(2) or 1)
            if elem:
                elements[elem] = elements.get(elem, 0) + count

        return {
            "formula": formula,
            "phase": phase,
            "smiles": smiles,
            "code": code,
            "elements": elements,
        }

    @staticmethod
    def _reaction_to_dict(rxn: Any) -> Dict[str, Any]:
        rxn_type = getattr(rxn, "r_type", "") or getattr(rxn, "reaction_type", "") or "other"

        reactants = []
        products = []
        reactant_key = ""
        product_key = ""
        bond_changes = ""

        if hasattr(rxn, "reactants") and hasattr(rxn, "products"):
            for r in (rxn.reactants or []):
                label = getattr(r, "formula", str(r))
                reactants.append(label)
                if not reactant_key:
                    reactant_key = getattr(r, "code", label)
            for p in (rxn.products or []):
                label = getattr(p, "formula", str(p))
                products.append(label)
                if not product_key:
                    product_key = getattr(p, "code", label)

        if hasattr(rxn, "bb"):
            bb = rxn.bb
            bond_changes = f"break {getattr(bb, 'atom1', '?')}-{getattr(bb, 'atom2', '?')}"
        elif hasattr(rxn, "bond_changes"):
            bond_changes = str(rxn.bond_changes)

        return {
            "reaction_type": str(rxn_type),
            "reactants": reactants,
            "products": products,
            "reactant_key": reactant_key,
            "product_key": product_key,
            "bond_changes": bond_changes,
        }


# ---------------------------------------------------------------------------
# Module-level convenience functions
# ---------------------------------------------------------------------------

def care_check_domain(species_labels: List[str]) -> CAREDomainReport:
    """Quick domain check without instantiating the full bridge."""
    return CAREBridge().check_domain(species_labels=species_labels)


def care_generate_crn_blueprint(
    seed_species: List[str],
    electro: bool = False,
    num_cpu: int = 1,
) -> Dict[str, Any]:
    """One-shot CRN blueprint generation."""
    return CAREBridge(num_cpu=num_cpu).generate_blueprint(
        seed_species=seed_species,
        electro=electro,
    )


__all__ = [
    "CAREBridge",
    "CAREDomainReport",
    "care_check_domain",
    "care_generate_crn_blueprint",
]
