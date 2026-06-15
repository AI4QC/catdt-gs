"""CatDT universal mechanism search schema definitions.

These are CatDT's own representations for multi-pathway mechanism search,
independent of any particular backend (CARE, etc.).  They live between
Agent3 output and Agent4/5 input.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------

class ReactionOperationType(str, Enum):
    """Types of elementary reaction operations."""
    ADSORPTION = "adsorption"
    DESORPTION = "desorption"
    DISSOCIATION = "dissociation"
    ASSOCIATION = "association"
    HYDROGENATION = "hydrogenation"
    DEHYDROGENATION = "dehydrogenation"
    PCET = "pcet"                      # proton-coupled electron transfer
    COUPLING = "coupling"              # C-C, C-N, etc.
    REARRANGEMENT = "rearrangement"
    ELEY_RIDEAL = "eley_rideal"
    OTHER = "other"


class StatePhase(str, Enum):
    """Phase of a catalytic state species."""
    GAS = "gas"
    ADSORBED = "ads"
    SURFACE = "surf"
    SOLVATED = "solv"


class EnergyBackendType(str, Enum):
    """Backend used for quick free-energy estimation."""
    THERMAL_UMA = "thermal_uma"        # UMA adsorption energy
    ELECTRO_CHE = "electro_che"        # computational hydrogen electrode
    HEURISTIC = "heuristic"            # rule-based estimate
    EXTERNAL = "external"              # user-supplied value


# ---------------------------------------------------------------------------
# State-level representations
# ---------------------------------------------------------------------------

class CatalyticStateRecord(BaseModel):
    """Lightweight, hashable record for one catalytic micro-state.

    This is intentionally richer than CARE's ``Intermediate``: it can encode
    surface identity, adsorption site type, coverage, and co-adsorbate
    configuration.
    """
    state_id: str = Field(
        ..., description="Unique canonical identifier, e.g. 'Cu111_hollow_CO*'",
    )
    species_label: str = Field(
        ..., description="Adsorbate / gas label, e.g. '*CO', 'H2(g)'",
    )
    phase: StatePhase = StatePhase.ADSORBED
    elements: Dict[str, int] = Field(
        default_factory=dict,
        description="Element counts of the species, e.g. {'C': 1, 'O': 1}",
    )
    surface_id: str = Field(
        default="",
        description="Surface identity tag, e.g. 'Cu111' or 'SrTiO3_001'",
    )
    site_type: str = Field(
        default="",
        description="Adsorption site type, e.g. 'hollow', 'bridge', 'atop'",
    )
    co_adsorbates: List[str] = Field(
        default_factory=list,
        description="Other species sharing the surface, e.g. ['*OH']",
    )
    coverage: float = Field(
        default=0.0, ge=0.0,
        description="Fractional coverage (0-1) for this species",
    )
    structure_path: Optional[str] = Field(
        default=None,
        description="Path to the materialized ASE structure file (.vasp/.traj)",
    )
    extra: Dict[str, Any] = Field(
        default_factory=dict,
        description="Arbitrary metadata (CARE intermediate ref, etc.)",
    )


class StateEnergyEstimate(BaseModel):
    """Quick free-energy estimate for a single state."""
    state_id: str
    free_energy_eV: Optional[float] = None
    adsorption_energy_eV: Optional[float] = None
    backend: EnergyBackendType = EnergyBackendType.HEURISTIC
    confidence: float = Field(default=0.0, ge=0.0, le=1.0)
    details: Dict[str, Any] = Field(default_factory=dict)


# ---------------------------------------------------------------------------
# Step-level representations
# ---------------------------------------------------------------------------

class ElementaryStepCandidate(BaseModel):
    """One candidate elementary reaction step between two catalytic states."""
    step_id: str = Field(
        ..., description="Unique step identifier",
    )
    reactant_state_id: str
    product_state_id: str
    operation_type: ReactionOperationType = ReactionOperationType.OTHER
    bond_changes: str = Field(
        default="",
        description="Concise bond change description, e.g. 'break C-O', 'form C-H'",
    )
    estimated_barrier_eV: Optional[float] = None
    estimated_reaction_energy_eV: Optional[float] = None
    source_backend: str = Field(
        default="",
        description="Which backend generated this candidate, e.g. 'care', 'generic_ops', 'llm'",
    )
    confidence: float = Field(default=0.0, ge=0.0, le=1.0)
    extra: Dict[str, Any] = Field(default_factory=dict)


# ---------------------------------------------------------------------------
# Pathway-level representations
# ---------------------------------------------------------------------------

class PathwayStateEntry(BaseModel):
    """One state node in a candidate pathway."""
    index: int
    state: CatalyticStateRecord
    energy_estimate: Optional[StateEnergyEstimate] = None


class PathwayStepEntry(BaseModel):
    """One step edge in a candidate pathway."""
    index: int
    step: ElementaryStepCandidate
    reactant_structure_path: Optional[str] = None
    product_structure_path: Optional[str] = None


class CandidatePathway(BaseModel):
    """A single candidate reaction pathway (sequence of states + steps)."""
    pathway_id: str
    description: str = ""
    states: List[PathwayStateEntry] = Field(default_factory=list)
    steps: List[PathwayStepEntry] = Field(default_factory=list)
    total_free_energy_change_eV: Optional[float] = None
    max_step_energy_eV: Optional[float] = None
    is_retained: bool = True
    prune_reason: str = ""
    export_dir: Optional[str] = None


class MultiPathwaySearchResult(BaseModel):
    """Top-level output of the mechanism search stage."""
    context_summary: str = ""
    all_pathways: List[CandidatePathway] = Field(default_factory=list)
    retained_pathways: List[str] = Field(
        default_factory=list,
        description="pathway_id list of pathways that survived pruning",
    )
    pruned_pathways: List[str] = Field(
        default_factory=list,
        description="pathway_id list of pathways removed by pruning",
    )
    search_stats: Dict[str, Any] = Field(default_factory=dict)
    export_manifest: Dict[str, Any] = Field(default_factory=dict)


# ---------------------------------------------------------------------------
# Context / Config
# ---------------------------------------------------------------------------

class MechanismContext(BaseModel):
    """Structured context for the mechanism search stage.

    Created from Agent3 output + user-supplied reaction specification.
    """
    bulk_formula: str = ""
    surface_id: str = ""
    surface_facet: str = ""
    initial_state: str = Field(
        ..., description="Starting adsorbate/species label",
    )
    target_state: str = Field(
        ..., description="Desired final state label",
    )
    environment: Dict[str, Any] = Field(
        default_factory=dict,
        description="T, P, U, pH, solvent, etc.",
    )
    known_intermediates: List[str] = Field(
        default_factory=list,
        description="User-supplied hints about intermediate species",
    )
    constraints: str = ""
    surface_structure_path: Optional[str] = None
    clean_slab_path: Optional[str] = None
    surface_indices: List[int] = Field(default_factory=list)
    adsorbate_indices: List[int] = Field(default_factory=list)
    fixed_adsorption_site: Optional[List[float]] = Field(
        default=None,
        description=(
            "If set, all intermediates are evaluated at this fixed [x,y,z] "
            "adsorption site (from Agent2/3 best site), skipping site search. "
            "Dramatically speeds up UMA evaluation."
        ),
    )


class ExplorationMode(str, Enum):
    """Which exploration module(s) to run."""
    AGENT_GUIDED = "agent_guided"      # LLM recommends pathway → fast evaluate
    SYSTEMATIC = "systematic"          # CARE-style CRN build from scratch
    BOTH_SEQUENTIAL = "both_sequential"  # agent-guided first, then systematic
    BOTH_PARALLEL = "both_parallel"    # run both, merge results


class MechanismSearchConfig(BaseModel):
    """Tunable parameters for the search engine."""
    exploration_mode: ExplorationMode = Field(
        default=ExplorationMode.AGENT_GUIDED,
        description=(
            "agent_guided: LLM recommends pathway, fast evaluation. "
            "systematic: CRN build from scratch (CARE-style). "
            "both_sequential: agent-guided first, then systematic for uncovered branches. "
            "both_parallel: run both, merge and deduplicate."
        ),
    )
    max_depth: int = Field(default=8, ge=1)
    beam_width: int = Field(default=8, ge=1)
    delta_keep_eV: float = Field(
        default=0.10, ge=0.0,
        description="Sibling energy gap <= this → keep both branches",
    )
    delta_prune_eV: float = Field(
        default=0.30, ge=0.0,
        description="Sibling energy gap >= this → prune the higher branch",
    )
    max_state_evaluations: int = Field(default=40, ge=1)
    energy_backend: EnergyBackendType = EnergyBackendType.THERMAL_UMA
    enable_care_backend: bool = True
    enable_generic_ops_backend: bool = True
    enable_llm_expansion: bool = False


# ---------------------------------------------------------------------------
# Module exports
# ---------------------------------------------------------------------------

__all__ = [
    "ReactionOperationType",
    "StatePhase",
    "EnergyBackendType",
    "CatalyticStateRecord",
    "StateEnergyEstimate",
    "ElementaryStepCandidate",
    "PathwayStateEntry",
    "PathwayStepEntry",
    "CandidatePathway",
    "MultiPathwaySearchResult",
    "MechanismContext",
    "MechanismSearchConfig",
    "ExplorationMode",
]
