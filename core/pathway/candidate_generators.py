"""Candidate generation: two independent exploration modules.

Module 1 — **AgentGuidedPathwayGenerator**
    LLM (or user) supplies a recommended pathway as an ordered list of
    intermediate labels.  This module validates the sequence, fills in
    missing steps, and produces candidate steps to evaluate.
    Fast: only evaluates species on the suggested path.
    Works for ANY reaction on ANY surface — the LLM proposes, UMA evaluates.

Module 2 — **SystematicCRNExplorer**
    Builds a reaction network from scratch using RDKit recursive bond
    operations (same principle as CARE, but element-agnostic).
    Thorough: enumerates all reachable states up to search budget.
    Discovers pathways the LLM might miss.

Both modules produce the same output format: lists of candidate dicts
with ``product_label``, ``product_elements``, ``operation_type``, etc.
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# RDKit availability
# ---------------------------------------------------------------------------

_RDKIT_AVAILABLE: Optional[bool] = None


def _check_rdkit() -> bool:
    global _RDKIT_AVAILABLE
    if _RDKIT_AVAILABLE is not None:
        return _RDKIT_AVAILABLE
    try:
        from rdkit import Chem  # noqa: F401
        _RDKIT_AVAILABLE = True
    except ImportError:
        _RDKIT_AVAILABLE = False
    return _RDKIT_AVAILABLE


# ═══════════════════════════════════════════════════════════════════════════
# Module 1: Agent-Guided Pathway Generator
# ═══════════════════════════════════════════════════════════════════════════

class AgentGuidedPathwayGenerator:
    """Generate candidate pathways from LLM / user recommended sequences.

    Input: an ordered list of intermediate labels, e.g.
        ["*CO2", "*COOH", "*CO", "*CHO", "*CH2O", "*CH3O", "*CH3OH", "*CH3", "*CH4"]

    Output: a list of (reactant_label, product_label, operation_type) steps
    that the search engine can directly evaluate with UMA.

    This works for ANY catalytic reaction — the chemical knowledge comes
    from the LLM or user, not from hard-coded rules.
    """

    def generate_pathway_steps(
        self,
        intermediates: List[str],
    ) -> List[Dict[str, Any]]:
        """Convert an ordered intermediate sequence into candidate steps."""
        if len(intermediates) < 2:
            return []

        steps: List[Dict[str, Any]] = []
        for i in range(len(intermediates) - 1):
            reactant = intermediates[i]
            product = intermediates[i + 1]
            r_elements = parse_species_elements(reactant)
            p_elements = parse_species_elements(product)
            op_type = _infer_operation_type(r_elements, p_elements, reactant, product)

            steps.append({
                "reactant_label": reactant,
                "product_label": product,
                "product_elements": p_elements,
                "operation_type": op_type,
                "bond_changes": f"{reactant} → {product}",
                "source_backend": "agent_guided",
            })

        return steps

    def generate_branching_candidates(
        self,
        recommended_pathways: List[List[str]],
    ) -> List[List[Dict[str, Any]]]:
        """Generate steps for multiple recommended pathways (competing routes)."""
        return [self.generate_pathway_steps(pw) for pw in recommended_pathways]


def _infer_operation_type(
    r_elem: Dict[str, int], p_elem: Dict[str, int],
    r_label: str, p_label: str,
) -> str:
    """Infer elementary step type from element changes."""
    r_total = sum(r_elem.values())
    p_total = sum(p_elem.values())
    h_diff = p_elem.get("H", 0) - r_elem.get("H", 0)
    o_diff = p_elem.get("O", 0) - r_elem.get("O", 0)
    c_diff = p_elem.get("C", 0) - r_elem.get("C", 0)

    if "(g)" in p_label and "*" in r_label:
        return "desorption"
    if "*" in p_label and "(g)" in r_label:
        return "adsorption"
    if h_diff > 0 and o_diff == 0 and c_diff == 0:
        return "hydrogenation"
    if h_diff < 0 and o_diff == 0 and c_diff == 0:
        return "dehydrogenation"
    if p_total < r_total:
        return "dissociation"
    if h_diff > 0 and o_diff != 0:
        return "pcet"
    return "other"


# ═══════════════════════════════════════════════════════════════════════════
# Module 2: Systematic CRN Explorer
# ═══════════════════════════════════════════════════════════════════════════

class SystematicCRNExplorer:
    """Build reaction network from scratch via recursive bond operations.

    Same principle as CARE's dissociate() pipeline, but:
    - Element-agnostic (works for any elements RDKit can represent)
    - Generates bond-breaking, hydrogenation, AND C-C coupling candidates
    - Supports C2–C6 products via coupling of surface intermediates
    - Integrates with CARE bridge when in-domain
    - Produces candidates one layer at a time for the search engine
    """

    # Default gas-phase reactants always implicitly available in a
    # catalytic environment. We dissociate each one at search
    # initialization to seed atomic / small-fragment partners (*H, *O,
    # *OH, *CO, *N, ...) into the species registry. This makes all
    # "+H / +O / +OH / +CO / Eley-Rideal / PCET" operations emerge from
    # the same generic bond-formation primitive instead of hardcoded
    # special cases. Users do not need to provide this list.
    DEFAULT_GAS_PARTNERS: Tuple[str, ...] = (
        "[H][H]",    # H2
        "O=O",       # O2
        "O",         # H2O
        "[C-]#[O+]", # CO
        "O=C=O",     # CO2
        "N#N",       # N2
    )

    def __init__(
        self,
        care_bridge: Any = None,
        max_carbon: int = 6,
        max_heavy_atoms: int = 12,
        max_partner_heavy_atoms: int = 3,
        gas_partners: Optional[Tuple[str, ...]] = None,
    ):
        self.care_bridge = care_bridge
        self.max_carbon = max_carbon
        self.max_heavy_atoms = max_heavy_atoms  # max non-H atoms in association products
        # For association, only try partners with at most this many heavy
        # atoms. Large-core × large-partner coupling is skipped because
        # (a) BFS from the same small-partner pool reaches those products
        # via multiple small steps anyway, (b) iterating over a runaway
        # registry makes Stage-1 CRN construction quadratic in its own
        # output. Covers common small fragments (*H, *O, *OH, *CO, *CO2,
        # *CH, *CH2, *CH3, *NH, *NH2, *C2H5, *C3H7, …).
        self.max_partner_heavy_atoms = max_partner_heavy_atoms
        self._smiles_cache: Dict[str, str] = {}
        # Registry of ALL discovered surface species for association reactions
        self._species_registry: Dict[str, Dict[str, Any]] = {}  # label → {elements, smiles}
        # Registry of known formulas for isomerization detection
        self._formula_registry: Dict[str, List[str]] = {}  # formula → [label1, label2, ...]
        # Gas partners are seeded lazily on the first expansion so the
        # constructor stays cheap.
        self._gas_partner_smiles: Tuple[str, ...] = (
            tuple(gas_partners) if gas_partners is not None
            else self.DEFAULT_GAS_PARTNERS
        )
        self._gas_partners_seeded: bool = False
        # Elements declared relevant for this search (union of
        # reactant/target element sets). Set via configure_task() before
        # the first expansion; restricts default gas seeding to partners
        # whose atoms lie in this set.
        self._relevant_elements: Optional[Set[str]] = None

    def configure_task(
        self,
        relevant_elements: Optional[Set[str]] = None,
        reactant_elements: Optional[Dict[str, int]] = None,
        target_elements: Optional[Dict[str, int]] = None,
        margin: int = 2,
    ) -> None:
        """Declare the reaction task, optionally auto-sizing size limits.

        Two ways to call:
          1. ``configure_task(relevant_elements={"C","H","O"})`` — only
             filter default gas partners to those elements; keep the
             default size limits.
          2. ``configure_task(reactant_elements=..., target_elements=...)``
             — additionally auto-infer ``max_carbon`` / ``max_heavy_atoms``
             / ``max_partner_heavy_atoms`` from the reactant and target
             formulas plus a small margin. Works for any reaction (FT,
             OER, NRR, …) without tuning: the search window stretches to
             cover the target's size without exploding the CRN.

        All logic is element-agnostic (set intersection + atom counting)
        and adds no reaction-specific chemistry.
        """
        if relevant_elements is not None:
            elements = set(relevant_elements)
        elif reactant_elements is not None or target_elements is not None:
            elements = set(reactant_elements or {}) | set(target_elements or {})
        else:
            elements = set()
        self._relevant_elements = elements or None

        if reactant_elements is not None and target_elements is not None:
            def _heavy(ec: Dict[str, int]) -> int:
                return sum(c for e, c in ec.items() if e != "H")

            max_c = max(
                reactant_elements.get("C", 0),
                target_elements.get("C", 0),
            )
            max_heavy = max(_heavy(reactant_elements), _heavy(target_elements))
            self.max_carbon = max_c + margin
            self.max_heavy_atoms = max_heavy + margin
            # Partner ≤ ceil(target_heavy / 2): enough for A+A dimerisation
            # and half-and-half coupling, but bounds the per-expansion work.
            self.max_partner_heavy_atoms = max(3, (max_heavy + 1) // 2)

    @staticmethod
    def _fragment_gas_molecule(gas_smiles: str) -> List[Tuple[str, Dict[str, int], str]]:
        """Enumerate every fragment(+ H-containing fragment) that can come
        from breaking any single bond in ``gas_smiles``. Returns
        ``[(label, element_counts, smiles_or_empty), ...]`` with fragments
        recursively also dissociated, so e.g. H2O gives {*OH, *H, *O}.
        """
        if not _check_rdkit():
            return []
        from rdkit import Chem

        mol = Chem.MolFromSmiles(gas_smiles)
        if mol is None:
            return []
        mol = Chem.AddHs(mol)

        seen_smi: Set[str] = set()
        results: List[Tuple[str, Dict[str, int], str]] = []
        queue: List[Any] = [mol]

        while queue:
            m = queue.pop(0)
            if m.GetNumAtoms() < 1:
                continue
            # Register this fragment
            ec: Dict[str, int] = {}
            for atom in m.GetAtoms():
                sym = atom.GetSymbol()
                ec[sym] = ec.get(sym, 0) + 1
            if not ec:
                continue
            label = _build_label_from_elements(ec)
            try:
                frag_smi = Chem.MolToSmiles(Chem.RemoveHs(m))
            except Exception:
                frag_smi = ""
            key = (label, frag_smi)
            key_str = f"{label}|{frag_smi}"
            if key_str in seen_smi:
                continue
            seen_smi.add(key_str)
            results.append((label, ec, frag_smi))

            # Expand further by breaking every bond in this fragment
            for bond in m.GetBonds():
                a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                emol = Chem.RWMol(m)
                try:
                    emol.RemoveBond(a1, a2)
                    frags = Chem.GetMolFrags(emol, asMols=True, sanitizeFrags=False)
                except Exception:
                    continue
                for f in frags:
                    if f.GetNumAtoms() > 0:
                        queue.append(f)

        return results

    def _seed_gas_partners(self, relevant_elements: Optional[Set[str]] = None) -> None:
        """Dissociate every default gas partner recursively and register
        every fragment (*H, *O, *OH, *CO, *N, ...) into the species
        registry. After seeding, all "+H / +O / +OH / +CO / Eley-Rideal"
        operations emerge from ``_generate_association_candidates``.

        When ``relevant_elements`` is given (typically the union of the
        reactant's and target's element sets), any gas partner whose
        atoms fall entirely outside that set is skipped. This keeps the
        registry focused on the elements actually present in the task
        while remaining fully generic (no hardcoded chemistry — the
        filter is just set intersection on element symbols).
        """
        if self._gas_partners_seeded:
            return
        self._gas_partners_seeded = True

        def _gas_elements(smi: str) -> Set[str]:
            return set(parse_species_elements(_smiles_to_ads_label(smi)).keys())

        if relevant_elements is not None:
            # Always keep H2 (needed for any hydrogenation/dehydrogenation)
            relevant_elements = set(relevant_elements) | {"H"}
            filtered = [
                smi for smi in self._gas_partner_smiles
                if _gas_elements(smi).issubset(relevant_elements)
            ]
        else:
            filtered = list(self._gas_partner_smiles)

        for gas_smi in filtered:
            gas_label = _smiles_to_ads_label(gas_smi)
            gas_elements = parse_species_elements(gas_label)
            if gas_label not in self._species_registry:
                self._species_registry[gas_label] = {
                    "elements": gas_elements, "smiles": gas_smi,
                }
                formula_key = _build_label_from_elements(gas_elements, adsorbed=False)
                self._formula_registry.setdefault(formula_key, [])
                if gas_label not in self._formula_registry[formula_key]:
                    self._formula_registry[formula_key].append(gas_label)

            for frag_label, frag_elements, frag_smi in self._fragment_gas_molecule(gas_smi):
                if frag_label in self._species_registry:
                    continue
                self._species_registry[frag_label] = {
                    "elements": frag_elements, "smiles": frag_smi or None,
                }
                formula_key = _build_label_from_elements(frag_elements, adsorbed=False)
                self._formula_registry.setdefault(formula_key, [])
                if frag_label not in self._formula_registry[formula_key]:
                    self._formula_registry[formula_key].append(frag_label)

    def expand_from_state(
        self,
        species_label: str,
        elements: Optional[Dict[str, int]] = None,
    ) -> List[Dict[str, Any]]:
        """Generate all candidate next-states from a given adsorbate.

        All operations are expressed as one of two generic primitives:
        ``_enumerate_bond_breaking`` (dissociation) and
        ``_enumerate_bond_formation`` + registry partners (association,
        hydrogenation, oxidation, hydroxylation, Eley-Rideal, PCET,
        coupling — all element-agnostic). Gas-phase partners that
        typical catalysis implicitly assumes (H2, O2, H2O, CO, CO2, N2)
        are seeded into the registry once at first call, so users
        provide only the initial and target states.
        """
        self._seed_gas_partners(self._relevant_elements)

        if elements is None:
            elements = parse_species_elements(species_label)
        candidates: List[Dict[str, Any]] = []
        is_adsorbed = species_label.startswith("*") and "(g)" not in species_label

        if not is_adsorbed:
            # Gas-phase: molecular adsorption + dissociative adsorption
            candidates.extend(self._gas_phase_reactions(species_label, elements))
            for c in candidates:
                c.setdefault("source_backend", "systematic_crn")
            return candidates

        smiles = self._label_to_smiles(species_label)

        # 1–2. Bond dissociation + Hydrogenation
        if smiles and _check_rdkit():
            candidates.extend(self._rdkit_expand(species_label, smiles, elements))
        else:
            candidates.extend(self._formula_expand(species_label, elements))

        # Universal bond-formation: one primitive handles hydrogenation,
        # oxidation, hydroxylation, coupling, Eley-Rideal (gas partners
        # are seeded into the registry at init) and PCET (same physics
        # as +H; the electrochem label is applied downstream).
        candidates.extend(self._generate_association_candidates(species_label, elements))

        # Formula-level dissociation (inverse of association) — required so
        # atoms can also be *removed*, e.g. *CHO → *CO + *H, *C6H14O →
        # *C6H14 + *O. Without this, the search can only monotonically
        # add atoms and cannot reach targets that have fewer of some
        # element than the reactant (e.g. *CO to *C6H14 loses one O).
        candidates.extend(self._generate_dissociation_candidates(species_label, elements))

        # Associative desorption (2*A → A₂(g)) is a composition change at
        # the gas interface; not expressible as a single bond-break on a
        # single core, so we keep a small generic helper for it.
        candidates.extend(self._associative_desorption_candidates(species_label, elements))

        # Isomerization (same formula, different connectivity) — uses
        # RDKit H-migration, element-agnostic.
        candidates.extend(self._isomerization_candidates(species_label, elements, smiles))

        # Register this species for future association reactions
        if is_adsorbed:
            self._species_registry[species_label] = {"elements": elements, "smiles": smiles}
            formula_key = _build_label_from_elements(elements, adsorbed=False)
            self._formula_registry.setdefault(formula_key, [])
            if species_label not in self._formula_registry[formula_key]:
                self._formula_registry[formula_key].append(species_label)

        for c in candidates:
            c.setdefault("source_backend", "systematic_crn")
        return candidates

    def _rdkit_expand(
        self, label: str, smiles: str, elements: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """Generate candidates via the two generic primitives.

        Only emits operations that do NOT require a partner species:
        bond dissociation (scission within the core) and desorption
        (surface → gas). All bond-formation operations (hydrogenation,
        oxidation, hydroxylation, coupling, Eley-Rideal, PCET) live in
        ``_generate_association_candidates`` and share one primitive
        with partners drawn from the species registry.
        """
        candidates: List[Dict[str, Any]] = []

        # --- Bond dissociation ---
        for _parent, frag_labels, frag_elem_list, bond_desc, frag_smi_list in _enumerate_bond_breaks(smiles):
            for i, (fl, fe) in enumerate(zip(frag_labels, frag_elem_list)):
                other = frag_labels[1 - i] if len(frag_labels) == 2 else ""
                # Cache each fragment's SMILES so that when it later
                # becomes an expansion core, it can participate in
                # association reactions (without this, any fragment
                # label not in the hardcoded _KNOWN dict would have no
                # SMILES on its next visit).
                if i < len(frag_smi_list) and frag_smi_list[i]:
                    self._smiles_cache[fl] = frag_smi_list[i]
                candidates.append({
                    "product_label": fl,
                    "product_elements": fe,
                    "operation_type": "dissociation",
                    "bond_changes": f"{bond_desc} → {fl}" + (f" + {other}" if other else ""),
                })

        # --- Desorption ---
        if sum(elements.values()) <= 6:
            gas_label = label.replace("*", "") + "(g)"
            candidates.append({
                "product_label": gas_label,
                "product_elements": dict(elements),
                "operation_type": "desorption",
                "bond_changes": "desorb from surface",
            })

        return candidates

    def _formula_expand(self, label: str, elements: Dict[str, int]) -> List[Dict[str, Any]]:
        """Minimal fallback when RDKit is unavailable.

        Without RDKit we cannot validate product valence, so we only emit
        the two operations that do not require creating a new chemical
        structure: dehydrogenation (drop one H) and desorption (phase
        change, same composition). Bond-formation candidates are
        deliberately NOT emitted, since every production-quality fix
        for valence-invalid formulas lives in RDKit sanitize.
        """
        candidates: List[Dict[str, Any]] = []

        # Dehydrogenation — purely subtractive, cannot over-valence
        if elements.get("H", 0) > 0:
            dh = dict(elements)
            dh["H"] -= 1
            if dh["H"] == 0:
                del dh["H"]
            candidates.append({
                "product_label": _build_label_from_elements(dh),
                "product_elements": dh,
                "operation_type": "dehydrogenation",
                "bond_changes": "break X-H",
            })

        # Desorption — pure phase change, same composition
        if sum(elements.values()) <= 6:
            candidates.append({
                "product_label": label.replace("*", "") + "(g)",
                "product_elements": dict(elements),
                "operation_type": "desorption",
                "bond_changes": "desorb from surface",
            })

        return candidates

    # --- 3. Universal surface association (*A + *B → *AB) ---

    def _generate_association_candidates(
        self, species_label: str, elements: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """Generate all "core + partner → product" candidates generically.

        Pure formula-level: sum element counts of core and partner,
        validate via the generic ``_is_plausible_formula`` valence check,
        emit a Hill-format label. No SMILES, no ``_KNOWN`` dict, no
        element-specific code — the primitive works identically for any
        combination of elements RDKit knows valences for.

        Every partner in the species registry (gas partners seeded at
        init + species discovered during search) participates, so
        +H / +O / +OH / +CO / Eley-Rideal / PCET / C-C / O-O / N-N
        coupling / self-dimerization all emerge from this single loop —
        no named operations.
        """
        candidates: List[Dict[str, Any]] = []
        seen_products: Set[str] = set()
        my_heavy = sum(v for k, v in elements.items() if k != "H")

        def _emit_for_partner(partner_label: str, partner_elem: Dict[str, int],
                              is_self: bool):
            other_heavy = sum(v for k, v in partner_elem.items() if k != "H")
            if other_heavy > self.max_partner_heavy_atoms:
                return
            total_heavy = my_heavy + other_heavy
            if total_heavy > self.max_heavy_atoms:
                return
            total_c = elements.get("C", 0) + partner_elem.get("C", 0)
            if total_c > self.max_carbon:
                return

            product_elem: Dict[str, int] = {}
            for e in set(list(elements.keys()) + list(partner_elem.keys())):
                product_elem[e] = elements.get(e, 0) + partner_elem.get(e, 0)
            if not _is_plausible_formula(product_elem):
                return

            product_label = _build_label_from_elements(product_elem)
            if product_label in seen_products or product_label == species_label:
                return
            seen_products.add(product_label)
            op_type = _infer_operation_type(
                elements, product_elem, species_label, product_label,
            )
            rxn_desc = (
                f"2×{species_label} → {product_label}" if is_self
                else f"{species_label} + {partner_label} → {product_label}"
            )
            candidates.append({
                "product_label": product_label,
                "product_elements": product_elem,
                "operation_type": op_type,
                "bond_changes": rxn_desc,
            })

        for partner_label, partner_info in self._species_registry.items():
            if partner_label == species_label:
                continue
            _emit_for_partner(
                partner_label,
                partner_info.get("elements", {}),
                is_self=False,
            )

        _emit_for_partner(species_label, elements, is_self=True)

        return candidates

    def _generate_dissociation_candidates(
        self, species_label: str, elements: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """Generic formula-level dissociation: the inverse of association.

        For each registered species R, check whether ``elements - R`` is a
        plausible formula too. If yes, emit ``species → R + complement``
        as a dissociation candidate with the complement as the primary
        product. This lets the search remove atoms (needed when the
        target formula has fewer of some element than the reactant),
        without requiring a SMILES for the core.

        Element-agnostic; the only chemistry in use is the valence-sum
        plausibility check on both fragments.
        """
        candidates: List[Dict[str, Any]] = []
        seen_products: Set[str] = set()

        for partner_label, partner_info in self._species_registry.items():
            partner_elem = partner_info.get("elements", {})
            # Compute complement = elements - partner
            complement: Dict[str, int] = {}
            valid = True
            for e, c in elements.items():
                remaining = c - partner_elem.get(e, 0)
                if remaining < 0:
                    valid = False
                    break
                if remaining > 0:
                    complement[e] = remaining
            if not valid:
                continue
            # Need partner to exactly overlap elements; any extra partner
            # element means it's not a real fragment of this core
            if any(partner_elem.get(e, 0) > elements.get(e, 0) for e in partner_elem):
                continue
            if not complement:
                continue  # partner == core, not dissociation
            if sum(partner_elem.values()) == 0:
                continue  # empty partner
            if not _is_plausible_formula(partner_elem):
                continue
            if not _is_plausible_formula(complement):
                continue

            # Emit both fragments as primary-product candidates
            for prod_elem, other_elem, other_label in (
                (complement, partner_elem, partner_label),
                (partner_elem, complement, _build_label_from_elements(complement)),
            ):
                prod_label = _build_label_from_elements(prod_elem)
                if prod_label in seen_products or prod_label == species_label:
                    continue
                seen_products.add(prod_label)
                candidates.append({
                    "product_label": prod_label,
                    "product_elements": prod_elem,
                    "operation_type": "dissociation",
                    "bond_changes": f"{species_label} → {prod_label} + {other_label}",
                })

        return candidates

    # --- Dissociative adsorption + Associative desorption ---
    # (Eley-Rideal and PCET are expressed as association with gas
    # partners in the registry and are generated by
    # _generate_association_candidates — no dedicated helpers needed.)

    def _gas_phase_reactions(
        self, species_label: str, elements: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """For gas-phase input: molecular adsorption + dissociative adsorption."""
        candidates: List[Dict[str, Any]] = []
        ads_label = "*" + species_label.replace("(g)", "").strip().lstrip("*")

        # Molecular adsorption: A(g) → *A
        candidates.append({
            "product_label": ads_label,
            "product_elements": dict(elements),
            "operation_type": "adsorption",
            "bond_changes": f"{species_label} → {ads_label}",
        })

        # Dissociative adsorption: A₂(g) → 2*A (for H2, O2, N2)
        _DISSOCIATIVE = {
            "H2": ("*H", {"H": 1}),
            "O2": ("*O", {"O": 1}),
            "N2": ("*N", {"N": 1}),
        }
        cleaned = species_label.replace("(g)", "").replace("*", "").strip()
        if cleaned in _DISSOCIATIVE:
            frag_label, frag_elem = _DISSOCIATIVE[cleaned]
            candidates.append({
                "product_label": frag_label,
                "product_elements": frag_elem,
                "operation_type": "dissociative_adsorption",
                "bond_changes": f"{species_label} → 2{frag_label}",
            })

        return candidates

    def _associative_desorption_candidates(
        self, species_label: str, elements: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """Associative desorption: 2*A → A₂(g) (for *H, *O, *N)."""
        candidates: List[Dict[str, Any]] = []
        _ASSOC = {
            "*H": "H2(g)", "*O": "O2(g)", "*N": "N2(g)",
        }
        if species_label in _ASSOC:
            candidates.append({
                "product_label": _ASSOC[species_label],
                "product_elements": {k: v * 2 for k, v in elements.items()},
                "operation_type": "associative_desorption",
                "bond_changes": f"2{species_label} → {_ASSOC[species_label]}",
            })
        return candidates

    # --- Isomerization (1,2-H migration) ---

    def _isomerization_candidates(
        self, species_label: str, elements: Dict[str, int],
        smiles: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Generate isomerization candidates: same formula, different structure.

        Uses RDKit to enumerate H-migration products (move H from one heavy
        atom to an adjacent one).
        """
        candidates: List[Dict[str, Any]] = []
        if not smiles or not _check_rdkit():
            return candidates

        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return candidates

        original_formula = rdMolDescriptors.CalcMolFormula(mol)
        original_canonical = Chem.MolToSmiles(mol)
        seen: Set[str] = {original_canonical}

        mol_h = Chem.AddHs(mol)
        # For each H, try moving it to each adjacent heavy atom's neighbor
        for atom in mol_h.GetAtoms():
            if atom.GetSymbol() != "H":
                continue
            h_idx = atom.GetIdx()
            # Find what the H is bonded to
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            if not neighbors:
                continue
            source_idx = neighbors[0]
            source_atom = mol_h.GetAtomWithIdx(source_idx)
            # Find heavy-atom neighbors of source (potential H destinations)
            for dest in source_atom.GetNeighbors():
                dest_idx = dest.GetIdx()
                if dest_idx == h_idx or dest.GetSymbol() == "H":
                    continue
                # Try moving H: remove H-source bond, add H-dest bond
                emol = Chem.RWMol(Chem.RWMol(mol_h))
                try:
                    emol.RemoveBond(h_idx, source_idx)
                    emol.AddBond(h_idx, dest_idx, Chem.BondType.SINGLE)
                    Chem.SanitizeMol(emol)
                    product = Chem.RemoveHs(emol)
                    psmi = Chem.MolToSmiles(product)
                except Exception:
                    continue
                if not psmi or psmi in seen:
                    continue
                # Verify same formula
                pm = Chem.MolFromSmiles(psmi)
                if pm and rdMolDescriptors.CalcMolFormula(pm) == original_formula:
                    seen.add(psmi)
                    product_label = _smiles_to_ads_label(psmi)
                    candidates.append({
                        "product_label": product_label,
                        "product_elements": dict(elements),
                        "operation_type": "isomerization",
                        "bond_changes": f"H-migration: {species_label} → {product_label}",
                    })

        # Also check formula registry for known isomers
        formula_key = _build_label_from_elements(elements, adsorbed=False)
        for known in self._formula_registry.get(formula_key, []):
            if known != species_label:
                candidates.append({
                    "product_label": known,
                    "product_elements": dict(elements),
                    "operation_type": "isomerization",
                    "bond_changes": f"isomerize {species_label} ↔ {known}",
                })

        return candidates

    @staticmethod
    def _try_couple_smiles(smi_a: Optional[str], smi_b: Optional[str]) -> Optional[str]:
        """Try to create a coupled product SMILES by joining two fragments with a C-C bond."""
        if not smi_a or not smi_b or not _check_rdkit():
            return None
        from rdkit import Chem
        mol_a = Chem.MolFromSmiles(smi_a)
        mol_b = Chem.MolFromSmiles(smi_b)
        if mol_a is None or mol_b is None:
            return None
        try:
            # Combine and add C-C bond between first C atoms
            combo = Chem.CombineMols(mol_a, mol_b)
            emol = Chem.RWMol(combo)
            # Find first C in each fragment
            c_a = next((a.GetIdx() for a in mol_a.GetAtoms() if a.GetSymbol() == "C"), None)
            c_b_offset = mol_a.GetNumAtoms()
            c_b = next((a.GetIdx() + c_b_offset for a in mol_b.GetAtoms() if a.GetSymbol() == "C"), None)
            if c_a is not None and c_b is not None:
                emol.AddBond(c_a, c_b, Chem.BondType.SINGLE)
                try:
                    Chem.SanitizeMol(emol)
                    return Chem.MolToSmiles(emol)
                except Exception:
                    pass
        except Exception:
            pass
        return None

    def _label_to_smiles(self, label: str) -> Optional[str]:
        """Convert adsorbate label to SMILES."""
        if label in self._smiles_cache:
            return self._smiles_cache[label]

        cleaned = label.replace("*", "").strip()
        result = _try_parse_smiles(cleaned)
        if result:
            self._smiles_cache[label] = result
        return result


# ═══════════════════════════════════════════════════════════════════════════
# Router: dispatches to Module 1 and/or Module 2
# ═══════════════════════════════════════════════════════════════════════════

class CandidateGeneratorRouter:
    """Routes to agent-guided and/or systematic exploration, merges results."""

    def __init__(
        self,
        enable_care: bool = True,
        enable_generic: bool = True,
        care_bridge: Any = None,
    ):
        self.enable_care = enable_care
        self.enable_generic = enable_generic
        self.care_bridge = care_bridge
        self._systematic = SystematicCRNExplorer(care_bridge=care_bridge)

    def configure_task(
        self,
        relevant_elements: Optional[Set[str]] = None,
        reactant_elements: Optional[Dict[str, int]] = None,
        target_elements: Optional[Dict[str, int]] = None,
        margin: int = 2,
    ) -> None:
        """Forward task scope to the underlying generator.

        When ``reactant_elements`` and ``target_elements`` are given, the
        generator also auto-sizes ``max_carbon`` / ``max_heavy_atoms`` /
        ``max_partner_heavy_atoms`` from the formulas so each new
        reaction gets an appropriate search window without manual tuning.
        """
        self._systematic.configure_task(
            relevant_elements=relevant_elements,
            reactant_elements=reactant_elements,
            target_elements=target_elements,
            margin=margin,
        )

    def generate_candidates(
        self,
        species_label: str,
        elements: Dict[str, int],
        *,
        care_network: Any = None,
    ) -> List[Dict[str, Any]]:
        """Generate candidate next-steps (used by systematic search engine)."""
        all_candidates: List[Dict[str, Any]] = []

        if self.enable_generic:
            all_candidates.extend(self._systematic.expand_from_state(species_label, elements))

        if self.enable_care and self.care_bridge is not None and care_network is not None:
            try:
                _, care_steps = self.care_bridge.reaction_network_to_candidates(care_network)
                for step in care_steps:
                    step["source_backend"] = "care"
                all_candidates.extend(care_steps)
            except Exception as exc:
                logger.warning("CARE candidate generation failed: %s", exc)

        return self.merge_and_deduplicate(all_candidates)

    @staticmethod
    def merge_and_deduplicate(candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        seen: Set[str] = set()
        deduped: List[Dict[str, Any]] = []
        for c in candidates:
            key = c.get("product_label", "")
            if key and key not in seen:
                seen.add(key)
                deduped.append(c)
        return deduped


# ═══════════════════════════════════════════════════════════════════════════
# RDKit helpers (shared by SystematicCRNExplorer)
# ═══════════════════════════════════════════════════════════════════════════

def _try_parse_smiles(cleaned: str) -> Optional[str]:
    """Try to convert a species name to SMILES."""
    if not _check_rdkit():
        return None
    from rdkit import Chem

    # Check catalysis-label → SMILES mapping FIRST. Labels like "CO",
    # "NO", "CH3" are not valid SMILES representations of the intended
    # catalytic species (e.g. SMILES "CO" is methanol, not carbon
    # monoxide). The KNOWN table holds the chemistry-correct mapping.
    _KNOWN: Dict[str, str] = {
        "CO2": "O=C=O", "CO": "[C-]#[O+]", "COOH": "OC=O",
        "CHO": "[CH]=O", "COH": "[C]=O",
        "CHOH": "O[CH]", "CH2O": "C=O",
        "CH3OH": "CO", "CH2OH": "[CH2]O", "CH3O": "CO", "OCH3": "CO",
        "CH4": "C", "CH3": "[CH3]", "CH2": "[CH2]", "CH": "[CH]",
        "OH": "[OH]", "H2O": "O", "H2": "[H][H]",
        "C2H4": "C=C", "C2H5OH": "CCO", "HCOOH": "OC=O",
        "NH3": "N", "NO": "[N]=O", "N2O": "N=NO", "NO2": "O=[N]=O",
        "SO2": "O=S=O", "H2S": "S",
    }
    known = _KNOWN.get(cleaned)
    if known is not None:
        return known

    # Fall back to treating the label as SMILES directly (works for
    # labels that happen to also be valid SMILES, e.g. formatted by
    # RDKit itself).
    mol = Chem.MolFromSmiles(cleaned)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    return None


# SMILES → readable adsorbate name mapping
# Maps canonical SMILES to commonly used catalysis names
_SMILES_TO_NAME: Dict[str, str] = {
    "O=CO": "*HCOOH",     # formic acid (CH2O2)
    "[O]C=O": "*COOH",    # carboxyl radical (CHO2)
    "O=C=O": "*CO2",
    "[C-]#[O+]": "*CO",
    "C=O": "*CH2O",       # formaldehyde
    "[CH]=O": "*CHO",     # formyl
    "CO": "*CH3OH",       # methanol
    "[CH2]O": "*CH2OH",
    "C": "*CH4",          # methane
    "[CH3]": "*CH3",
    "[CH2]": "*CH2",
    "[CH]": "*CH",
    "O": "*H2O",
    "[OH]": "*OH",
    "N": "*NH3",
    "[NH2]": "*NH2",
    "[NH]": "*NH",
}


def _smiles_to_ads_label(smiles: str) -> str:
    """Convert SMILES to readable adsorbate label like ``*HCOOH``, ``*CHO``."""
    if not _check_rdkit():
        return f"*{smiles}"

    # Check known name first
    if smiles in _SMILES_TO_NAME:
        return _SMILES_TO_NAME[smiles]

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is not None:
        canonical = Chem.MolToSmiles(mol)
        if canonical in _SMILES_TO_NAME:
            return _SMILES_TO_NAME[canonical]
        return f"*{rdMolDescriptors.CalcMolFormula(mol)}"

    # Unsanitized fallback (radical fragments)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is not None:
        ec: Dict[str, int] = {}
        for atom in mol.GetAtoms():
            s = atom.GetSymbol()
            if s != "H":
                ec[s] = ec.get(s, 0) + 1
        return f"*{_build_label_from_elements(ec, adsorbed=False)}"
    return f"*{smiles}"


def _enumerate_bond_breaks(smiles: str) -> List[Tuple[str, List[str], List[Dict[str, int]], str]]:
    """Break every non-equivalent bond; return fragment labels + element dicts."""
    if not _check_rdkit():
        return []
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    mol = Chem.AddHs(mol)

    try:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        ranks = list(range(mol.GetNumAtoms()))

    seen: Set[Tuple[int, int]] = set()
    results = []

    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        rp = (min(ranks[a1], ranks[a2]), max(ranks[a1], ranks[a2]))
        if rp in seen:
            continue
        seen.add(rp)

        sym1 = mol.GetAtomWithIdx(a1).GetSymbol()
        sym2 = mol.GetAtomWithIdx(a2).GetSymbol()
        desc = f"break {sym1}-{sym2}"

        emol = Chem.RWMol(mol)
        emol.RemoveBond(a1, a2)
        try:
            frags = Chem.GetMolFrags(emol, asMols=True, sanitizeFrags=False)
        except Exception:
            continue

        frag_labels = []
        frag_elem_list = []
        frag_smi_list = []
        for frag in frags:
            ec: Dict[str, int] = {}
            for atom in frag.GetAtoms():
                s = atom.GetSymbol()
                ec[s] = ec.get(s, 0) + 1
            if ec:
                frag_labels.append(_build_label_from_elements(ec))
                frag_elem_list.append(ec)
                try:
                    frag_smi = Chem.MolToSmiles(Chem.RemoveHs(frag))
                except Exception:
                    frag_smi = ""
                frag_smi_list.append(frag_smi)

        if frag_labels:
            results.append((smiles, frag_labels, frag_elem_list, desc, frag_smi_list))

    return results


_BOND_FORMATION_CACHE: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}


def _element_default_valence(symbol: str) -> Optional[int]:
    """Default valence for any element via RDKit's periodic table.

    Element-agnostic: covers every element RDKit knows. Atoms without a
    defined default valence (most transition metals) return None so the
    caller can treat them permissively.
    """
    if not _check_rdkit():
        return None
    from rdkit import Chem
    try:
        pt = Chem.GetPeriodicTable()
        atomic_num = pt.GetAtomicNumber(symbol)
        if atomic_num <= 0:
            return None
        v = pt.GetDefaultValence(atomic_num)
        return v if v > 0 else None
    except Exception:
        return None


def _is_plausible_formula(elements: Dict[str, int]) -> bool:
    """Generic valence-sum connectivity check.

    Rejects formulas that cannot form a connected molecular graph given
    standard valences (e.g. *CH5, *C2H7, *NH4). Unknown-valence atoms
    (transition metals, exotics) bypass the check and are accepted.
    No reaction- or surface-specific chemistry hardcoded.
    """
    elements = {e: c for e, c in elements.items() if c and c > 0}
    if not elements:
        return False
    valences: Dict[str, int] = {}
    for e in elements:
        v = _element_default_valence(e)
        if v is None:
            return True
        valences[e] = v
    n_atoms = sum(elements.values())
    if n_atoms == 1:
        return True
    total_valence = sum(valences[e] * c for e, c in elements.items())
    return total_valence >= 2 * (n_atoms - 1)


def _enumerate_bond_formation(
    core_smiles: str,
    partner_smiles: str,
) -> List[Tuple[str, str]]:
    """Generic "make one bond between core and partner" primitive.

    Element-agnostic. The partner can be a single atom (``[H]``, ``[O]``),
    a small radical (``[OH]``, ``[CH3]``), a full adsorbate, or a gas
    molecule — the primitive treats them identically. Used to express
    hydrogenation, oxidation, hydroxylation, C-C coupling, Eley-Rideal,
    and (in electrochem context) PCET as one operation: every flavour is
    just "form a bond with whichever partner the registry happens to hold".

    Two strategies are tried:
        (1) Direct single bond between each core atom and each partner
            atom. RDKit sanitize rejects valence-invalid products, so e.g.
            ``*CH4 + [H]`` produces nothing (no *CH5 can escape).
        (2) If strategy 1 yields nothing (fully saturated core), reduce
            one adjacent π-bond in the core by one order and then form
            the new bond. This generalises the "reduce triple bond + add
            H to *CO" pattern to any partner without H-specific code.

    Returns a list of ``(product_smiles, bond_description)`` pairs. The
    caller is responsible for translating SMILES → adsorbate label.
    """
    if not _check_rdkit():
        return []
    cache_key = (core_smiles, partner_smiles)
    cached = _BOND_FORMATION_CACHE.get(cache_key)
    if cached is not None:
        return cached

    from rdkit import Chem

    core_raw = Chem.MolFromSmiles(core_smiles)
    partner_raw = Chem.MolFromSmiles(partner_smiles)
    if core_raw is None or partner_raw is None:
        _BOND_FORMATION_CACHE[cache_key] = []
        return []

    # Make every H explicit before manipulating bonds; otherwise RDKit
    # silently rewrites implicit H counts and a valence-overflowing
    # "addition" round-trips to the core molecule unchanged.
    core = Chem.AddHs(core_raw)
    partner = Chem.AddHs(partner_raw)

    # Expected product element counts = core + partner. Double-check
    # formula after sanitize so no-op absorptions are caught.
    def _atom_counts(mol) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for atom in Chem.AddHs(mol).GetAtoms():
            sym = atom.GetSymbol()
            counts[sym] = counts.get(sym, 0) + 1
        return counts

    core_counts = _atom_counts(core_raw)
    partner_counts = _atom_counts(partner_raw)
    expected: Dict[str, int] = {}
    for sym in set(core_counts) | set(partner_counts):
        expected[sym] = core_counts.get(sym, 0) + partner_counts.get(sym, 0)

    results: List[Tuple[str, str]] = []
    seen: Set[str] = set()

    def _try_combine(core_mol, partner_mol, core_idx, partner_idx, desc):
        combined = Chem.RWMol(Chem.CombineMols(core_mol, partner_mol))
        p_offset = core_mol.GetNumAtoms()
        try:
            combined.AddBond(
                core_idx, partner_idx + p_offset, Chem.BondType.SINGLE,
            )
            Chem.SanitizeMol(combined)
            prod = Chem.RemoveHs(combined)
            psmi = Chem.MolToSmiles(prod)
        except Exception:
            return
        if not psmi or psmi in seen:
            return
        prod_mol = Chem.MolFromSmiles(psmi)
        if prod_mol is None:
            return
        if _atom_counts(prod_mol) != expected:
            return
        seen.add(psmi)
        results.append((psmi, desc))

    # Only try combining via non-hydrogen atoms of the core (H atoms on
    # core are already explicit; bonding a new partner TO a hydrogen is
    # not a meaningful single-step elementary operation). Partner side
    # can attach via any atom.
    def _non_h_indices(mol) -> List[int]:
        return [
            a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1
        ]

    core_sites = _non_h_indices(core)
    partner_sites = list(range(partner.GetNumAtoms()))

    # --- Strategy 1: direct bond between each core site and each partner site ---
    for ci in core_sites:
        c_sym = core.GetAtomWithIdx(ci).GetSymbol()
        for pi in partner_sites:
            p_sym = partner.GetAtomWithIdx(pi).GetSymbol()
            _try_combine(core, partner, ci, pi, f"form {c_sym}-{p_sym} bond")

    if results:
        _BOND_FORMATION_CACHE[cache_key] = results
        return results

    # --- Strategy 2: reduce a multi-bond in core, then form the new bond ---
    # Only runs when strategy 1 produced nothing (core is fully saturated).
    for bond in core.GetBonds():
        bt = bond.GetBondType()
        if bt == Chem.BondType.DOUBLE:
            new_bt = Chem.BondType.SINGLE
        elif bt == Chem.BondType.TRIPLE:
            new_bt = Chem.BondType.DOUBLE
        else:
            continue
        a1_idx = bond.GetBeginAtomIdx()
        a2_idx = bond.GetEndAtomIdx()
        if core.GetAtomWithIdx(a1_idx).GetAtomicNum() == 1:
            continue
        if core.GetAtomWithIdx(a2_idx).GetAtomicNum() == 1:
            continue
        a1_sym = core.GetAtomWithIdx(a1_idx).GetSymbol()
        a2_sym = core.GetAtomWithIdx(a2_idx).GetSymbol()
        for attach_idx, attach_sym in ((a1_idx, a1_sym), (a2_idx, a2_sym)):
            core_copy = Chem.RWMol(core)
            core_copy.GetBondBetweenAtoms(a1_idx, a2_idx).SetBondType(new_bt)
            for pi in partner_sites:
                p_sym = partner.GetAtomWithIdx(pi).GetSymbol()
                _try_combine(
                    core_copy, partner, attach_idx, pi,
                    f"reduce {a1_sym}={a2_sym}, form {attach_sym}-{p_sym}",
                )

    _BOND_FORMATION_CACHE[cache_key] = results
    return results


def _enumerate_h_additions(smiles: str) -> List[Tuple[str, str]]:
    """Add one H to each unique site; return (product_smiles, description).

    Two strategies:
    1. Direct H addition (works for unsaturated atoms like CH2, NH2).
    2. Bond-order reduction + H addition: for fully saturated molecules like
       O=C=O, reduce one double/triple bond to single then add H on both ends.
       This produces e.g. O=C=O + H → OC(=O)[H] (COOH) and O=C(-H)=O (CHO2).
    """
    if not _check_rdkit():
        return []
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    results: List[Tuple[str, str]] = []
    seen: Set[str] = set()

    def _try_add(product_mol, desc_prefix: str):
        try:
            Chem.SanitizeMol(product_mol)
            product = Chem.RemoveHs(product_mol)
            psmi = Chem.MolToSmiles(product)
        except Exception:
            return
        if psmi and psmi not in seen:
            seen.add(psmi)
            pm = Chem.MolFromSmiles(psmi)
            if pm:
                pf = rdMolDescriptors.CalcMolFormula(pm)
                results.append((psmi, f"{desc_prefix} → *{pf}"))

    # Strategy 1: direct H addition
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        sym = atom.GetSymbol()
        emol = Chem.RWMol(Chem.AddHs(mol))
        h_idx = emol.AddAtom(Chem.Atom(1))
        emol.AddBond(idx, h_idx, Chem.BondType.SINGLE)
        _try_add(emol, f"add H to {sym}")

    # Strategy 2: reduce bond order + add H (for saturated molecules like CO2)
    # For each double/triple bond, reduce it by 1 order and add H to one of the atoms
    if not results:
        for bond in mol.GetBonds():
            bt = bond.GetBondType()
            a1_idx = bond.GetBeginAtomIdx()
            a2_idx = bond.GetEndAtomIdx()
            a1_sym = mol.GetAtomWithIdx(a1_idx).GetSymbol()
            a2_sym = mol.GetAtomWithIdx(a2_idx).GetSymbol()

            if bt == Chem.BondType.DOUBLE:
                new_bt = Chem.BondType.SINGLE
            elif bt == Chem.BondType.TRIPLE:
                new_bt = Chem.BondType.DOUBLE
            else:
                continue

            # Add H to atom 1 (reduce bond, H goes to a1)
            emol = Chem.RWMol(Chem.AddHs(mol))
            emol.GetBondBetweenAtoms(a1_idx, a2_idx).SetBondType(new_bt)
            h_idx = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(a1_idx, h_idx, Chem.BondType.SINGLE)
            _try_add(emol, f"add H to {a1_sym} (reduce {a1_sym}={a2_sym})")

            # Add H to atom 2 (reduce bond, H goes to a2)
            emol2 = Chem.RWMol(Chem.AddHs(mol))
            emol2.GetBondBetweenAtoms(a1_idx, a2_idx).SetBondType(new_bt)
            h_idx2 = emol2.AddAtom(Chem.Atom(1))
            emol2.AddBond(a2_idx, h_idx2, Chem.BondType.SINGLE)
            _try_add(emol2, f"add H to {a2_sym} (reduce {a1_sym}={a2_sym})")

    return results


def _lookup_compositions_for_elements(elements: Dict[str, int]) -> List[str]:
    """Fallback for RDKit H-addition failures (e.g. CO triple bond)."""
    if not _check_rdkit():
        return []
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    key = frozenset((k, v) for k, v in elements.items() if v > 0)
    _DB: Dict[FrozenSet[Tuple[str, int]], List[str]] = {
        frozenset({("C", 1), ("H", 1), ("O", 1)}): ["[CH]=O", "[C]=O"],
        frozenset({("C", 1), ("H", 1), ("O", 2)}): ["OC=O"],
        frozenset({("C", 1), ("H", 2), ("O", 1)}): ["C=O", "[CH2]O"],
    }
    smiles_list = _DB.get(key, [])
    results = []
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m:
            f = rdMolDescriptors.CalcMolFormula(m)
            label = f"*{f}"
            if label not in results:
                results.append(label)
    return results


# ═══════════════════════════════════════════════════════════════════════════
# Common helpers
# ═══════════════════════════════════════════════════════════════════════════

def _build_label_from_elements(elements: Dict[str, int], adsorbed: bool = True) -> str:
    """Build ``*CHO`` style label from element counts (Hill order)."""
    parts = []
    for elem in ["C", "H"]:
        count = elements.get(elem, 0)
        if count > 0:
            parts.append(elem if count == 1 else f"{elem}{count}")
    for elem in sorted(elements.keys()):
        if elem in ("C", "H"):
            continue
        count = elements[elem]
        if count > 0:
            parts.append(elem if count == 1 else f"{elem}{count}")
    formula = "".join(parts)
    return f"*{formula}" if adsorbed else formula


def parse_species_elements(label: str) -> Dict[str, int]:
    """Parse element counts from ``*CHO``, ``CH3OH(g)``, etc."""
    cleaned = label.replace("*", "").replace("(g)", "").replace("(s)", "").replace("(l)", "")
    if "+" in cleaned:
        cleaned = cleaned.split("+")[0].strip()
    elements: Dict[str, int] = {}
    for match in re.finditer(r"([A-Z][a-z]?)(\d*)", cleaned):
        elem = match.group(1)
        count = int(match.group(2) or 1)
        if elem:
            elements[elem] = elements.get(elem, 0) + count
    return elements


__all__ = [
    "AgentGuidedPathwayGenerator",
    "SystematicCRNExplorer",
    "CandidateGeneratorRouter",
    "parse_species_elements",
]
