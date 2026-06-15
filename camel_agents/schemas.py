"""Shared schema definitions for CAMEL CatDT workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


class ReactionContext(BaseModel):
    reaction_type: str
    initial_adsorbate: str
    intermediates: List[str]
    constraints: str = ""


class AtomAddSpec(BaseModel):
    species: str = ""
    element: str = ""
    symbol: str = ""
    position: Optional[List[float]] = None
    reactant_position: Optional[List[float]] = None
    product_position: Optional[List[float]] = None
    reason: str = ""


class AtomModifySpec(BaseModel):
    index: int
    reactant_position: Optional[List[float]] = None
    product_position: Optional[List[float]] = None


class PathwayStepSpec(BaseModel):
    step_name: str
    reactant_formula: str
    product_formula: str
    reaction_type: str
    atoms_to_add: List[AtomAddSpec] = Field(default_factory=list)
    atoms_to_modify: List[AtomModifySpec] = Field(default_factory=list)
    atoms_to_remove: List[Any] = Field(default_factory=list)
    confidence: float = 0.0


class PathwayDesign(BaseModel):
    pathway_name: str
    description: str
    overall_reaction: str
    steps: List[PathwayStepSpec]
    confidence: float = 0.0


class ValidationReport(BaseModel):
    status: str
    issues: List[str] = Field(default_factory=list)
    feedback: str = ""


@dataclass
class WorkflowState:
    """LLM iteration / conversation state for a pipeline run.

    After the Phase-1 refactor (2026-04), physical surface state and file paths
    live on CatalystSurface — not here. WorkflowState only holds LLM iteration
    state, retrieval context, and orchestration bookkeeping. Access the surface
    via ``workflow._get_surface()`` (or ``surface.state`` in reverse).
    """

    run_id: str = ""
    output_base_dir: str = ""
    reaction_description: str = ""
    reaction_context: Optional[ReactionContext] = None
    intermediates: List[str] = field(default_factory=list)
    pathway_design: Optional[PathwayDesign] = None
    validation_report: Optional[ValidationReport] = None
    step_structures: List[Dict[str, Any]] = field(default_factory=list)
    tool_baseline_steps: List[Dict[str, Any]] = field(default_factory=list)
    memento_retrieval: Dict[str, Any] = field(default_factory=dict)
    knowledge_retrieval: Dict[str, Any] = field(default_factory=dict)
    skill_retrieval: Dict[str, Any] = field(default_factory=dict)
    memento_query_text: str = ""
    energy_gate_report: Dict[str, Any] = field(default_factory=dict)
    adsorption_energies: Dict[str, float] = field(default_factory=dict)
    adsorbate_energies: Dict[str, float] = field(default_factory=dict)
    # Vacuum-phase vibrational frequencies (cm^-1) per adsorbate — required
    # for CatMAP ideal_gas thermo mode to compute TΔS(T) correctly.
    gas_frequencies: Dict[str, List[float]] = field(default_factory=dict)
    neb_results: Dict[str, Any] = field(default_factory=dict)
    last_feedback: str = ""
    iteration_history: List[Dict[str, Any]] = field(default_factory=list)

    # Multi-pathway mechanism search fields (populated when enable_mechanism_search=True)
    mechanism_context: Optional[Dict[str, Any]] = None
    candidate_pathways: List[Dict[str, Any]] = field(default_factory=list)
    retained_pathways: List[Dict[str, Any]] = field(default_factory=list)
    pruned_pathways: List[Dict[str, Any]] = field(default_factory=list)
    pathway_manifests: Dict[str, Any] = field(default_factory=dict)
    mechanism_search_result: Optional[Dict[str, Any]] = None

    # Multi-facet fields (populated when multi_facet=True)
    all_surface_paths: Dict[str, str] = field(default_factory=dict)
    all_area_fractions: Dict[str, float] = field(default_factory=dict)
    all_miller_indices: Dict[str, tuple] = field(default_factory=dict)
    facet_results: List[Any] = field(default_factory=list)  # List[FacetResult]
    aggregated_tof: Optional[float] = None
    aggregated_production_rates: Dict[str, float] = field(default_factory=dict)
    current_facet_id: Optional[str] = None


@dataclass
class FacetResult:
    """Per-facet pipeline output for multi-facet Wulff-weighted aggregation."""
    facet_id: str = ""
    miller_index: Optional[tuple] = None
    area_fraction: float = 0.0
    surface_path: Optional[str] = None
    pathway_result_pkl: Optional[str] = None
    kmc_result_pkl: Optional[str] = None
    tof: float = 0.0
    production_rates: Dict[str, float] = field(default_factory=dict)
    neb_summary: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None


class CatDTConfig(BaseModel):
    run_id: str = Field(
        default_factory=lambda: datetime.now().strftime("%Y%m%d_%H%M%S"),
        description="Unique identifier for this workflow run.",
    )
    output_base_dir: str = Field(default="output/catdt_workflow")

    bulk_structure_path: Optional[str] = None
    initial_surface_path: Optional[str] = None
    reaction_description: Optional[str] = None

    top_n_surfaces: int = Field(3, ge=1)
    num_adsorption_sites: int = Field(5, ge=1)
    enable_llm_adsorption_review: bool = False
    llm_adsorption_review_context: Optional[str] = None

    adsorbate_elements_for_mc: Optional[List[str]] = None
    mc_temperature_k: Optional[float] = Field(None, gt=0)
    mc_total_sweeps: int = Field(100, ge=1)
    clean_slab_path_for_mc: Optional[str] = None

    # Two-step reconstruction: Agent3a(clean) → Agent2 → Agent3b(with adsorbate)
    mc_two_step_reconstruction: bool = Field(
        True,
        description="Enable two-step reconstruction: first reconstruct clean slab, "
                    "then adsorb, then reconstruct again with adsorbate present.",
    )
    mc_clean_sweeps: Optional[int] = Field(
        None, ge=1,
        description="MC sweeps for clean-slab reconstruction (step 3a). "
                    "Defaults to mc_total_sweeps // 2 if not set.",
    )

    # Energy model for VSSR-MC surface reconstruction
    mc_energy_model: str = Field(
        "CHGNet",
        description="Energy model for MC surface energy evaluation. "
                    "'CHGNet' (NFF framework, default) or 'UMA' (FairChem framework).",
    )
    mc_chem_pots: Optional[Dict[str, float]] = Field(
        None,
        description="Optional per-element chemical potentials for VSSR-MC.",
    )
    mc_num_adsorbates_override: Optional[int] = Field(
        None,
        ge=0,
        description="Optional canonical adsorbate count override for VSSR-MC.",
    )
    mc_use_seed_for_virtual_sites: bool = Field(
        False,
        description="Generate VSSR virtual sites on the seed structure instead of "
                    "an adsorbate-free slab.",
    )
    mc_adsorbate_exclusion_radius_A: Optional[float] = Field(
        None,
        ge=0,
        description="Optional radius in A for excluding VSSR virtual Ti/O sites "
                    "around molecular Agent2 adsorbates before Agent3b.",
    )
    mc_existing_atom_exclusion_radius_A: Optional[float] = Field(
        None,
        ge=0,
        description="Optional radius in A for excluding VSSR virtual sites "
                    "around pre-existing non-molecular atoms before Agent3b.",
    )
    stop_after_agent3b: bool = Field(
        False,
        description="Run only Agent1/reaction parsing/Agent2/Agent3b and return "
                    "before Agent4/5, NEB, adsorption-energy, and KMC stages.",
    )

    # Agent2 adsorption backend — AdsorbML (default) or AdsorbDiff (legacy)
    adsorption_backend: str = Field(
        "adsorbml",
        description="Agent2 backend for initial adsorbate placement. "
                    "'adsorbml' (heuristic enumeration + UMA relax + chemisorption "
                    "filter — default, produces real chemisorbed minima) or "
                    "'adsorbdiff' (legacy diffusion model using proxy energy, "
                    "often returns physisorbed / unbound configs).",
    )
    adsorbml_num_sites: int = Field(
        20, ge=1,
        description="AdsorbML: number of candidate sites to enumerate + relax.",
    )
    adsorbml_placement_mode: str = Field(
        "random_site_heuristic_placement",
        description="AdsorbML placement mode: 'heuristic' / 'random' / "
                    "'random_site_heuristic_placement'.",
    )
    adsorbml_interstitial_gap: float = Field(
        0.1, ge=0,
        description="AdsorbML: minimum initial adsorbate-slab gap in A. "
                    "For saturated PDH molecules on oxide overlayers, larger "
                    "values such as 2.0 can avoid artificial C-H dissociation "
                    "during Agent2 relaxation.",
    )
    adsorbml_relax_steps: int = Field(
        200, ge=0,
        description="AdsorbML: ASE optimizer max steps per enumerated candidate.",
    )
    adsorbml_relax_fmax: float = Field(
        0.05, gt=0,
        description="AdsorbML: force convergence target in eV/A.",
    )

    # Electrochemical parameters (liquid-solid interface)
    potential_she: Optional[float] = Field(
        None,
        description="Electrode potential vs. SHE (V). If set, enables electrochemical "
                    "Pourbaix mode for surface reconstruction.",
    )
    ph: Optional[float] = Field(
        None, ge=0, le=14,
        description="Solution pH for electrochemical reconstruction.",
    )

    # Materials Project API key (for Pourbaix diagram data)
    mp_api_key: Optional[str] = Field(
        None,
        description="Materials Project API key for Pourbaix data. "
                    "Falls back to MP_API_KEY env var if not set.",
    )

    # Multi-facet Wulff-weighted microkinetics
    multi_facet: bool = Field(
        False,
        description="Run the full Agent2→NEB→CatMAP pipeline independently on each "
                    "of the top-N Wulff-exposed facets, then aggregate TOF by area fraction: "
                    "TOF_total = Σ (area_fraction_i × TOF_i). "
                    "When False (default), only the most exposed facet is processed.",
    )
    multi_facet_top_n: int = Field(
        3, ge=1, le=10,
        description="Number of top Wulff facets to process when multi_facet=True.",
    )

    calculate_barriers: bool = True
    num_pathway_sites: int = Field(3, ge=1)
    neb_n_frames: int = Field(5, ge=3)
    neb_max_steps: int = Field(80, ge=1)
    adsorption_relax_max_steps: int = Field(120, ge=1)
    stagea_relax_steps: int = Field(
        200, ge=0,
        description="Stage A endpoint relaxation max steps after atom add/remove.",
    )
    stagea_relax_fmax: float = Field(
        0.05, gt=0,
        description="Stage A endpoint relaxation force convergence target in eV/A.",
    )

    kmc_temperature_k: Optional[float] = Field(None, gt=0)
    gas_pressures: Optional[Dict[str, float]] = None

    # Mechanism search stage (Agent3 → Agent4/5 upstream)
    enable_mechanism_search: bool = Field(
        False,
        description="Enable multi-pathway mechanism search between Agent3 and Agent4/5.",
    )
    mechanism_exploration_mode: str = Field(
        "agent_guided",
        description=(
            "Exploration module: "
            "'agent_guided' = LLM recommends pathway (fast), "
            "'systematic' = CRN from scratch (thorough), "
            "'both_sequential' = agent-guided first then systematic, "
            "'both_parallel' = run both and merge."
        ),
    )
    mechanism_max_depth: int = Field(8, ge=1)
    mechanism_beam_width: int = Field(8, ge=1)
    mechanism_delta_keep: float = Field(
        0.10, ge=0.0,
        description="Sibling energy gap <= this → keep both branches (eV).",
    )
    mechanism_delta_prune: float = Field(
        0.30, ge=0.0,
        description="Sibling energy gap >= this → prune the higher branch (eV).",
    )
    mechanism_max_state_evaluations: int = Field(40, ge=1)

    # Systematic search strategy: 'quick' (per-path eval), 'pruning'
    # (layer-by-layer beam BFS on CRN), or 'mcts' (UCB1 tree search).
    mechanism_systematic_strategy: str = Field(
        "pruning",
        description=(
            "Phase-3 algorithm inside 'systematic' exploration mode: "
            "'quick' = score Stage-1 paths by max-intermediate UMA energy; "
            "'pruning' = layer-by-layer beam BFS on the CRN graph with "
            "per-parent sibling pruning + distance-bucketed round-robin "
            "beam cap; 'mcts' = Monte Carlo Tree Search with UCB1."
        ),
    )
    mechanism_mcts_iterations: int = Field(
        100, ge=1,
        description="Number of MCTS iterations (selection+expansion+rollout+backprop).",
    )
    mechanism_mcts_c: float = Field(
        1.4, ge=0.0,
        description="UCB1 exploration constant. Higher = more exploration.",
    )
    mechanism_mcts_rollout_depth: int = Field(
        4, ge=0,
        description="Max steps taken during random rollout from an expanded leaf.",
    )
    mechanism_mcts_rollout_bias_beta: float = Field(
        1.0, ge=0.0,
        description=(
            "Softmax inverse-temperature for target-aware rollout sampling. "
            "Candidates are picked ∝ exp(-beta * element_count_distance_to_target). "
            "0 = pure random (ignore target), higher = greedier toward target formula. "
            "Generic: only depends on element-count distance, not on reaction chemistry."
        ),
    )
    mechanism_mcts_expansion_bias_beta: float = Field(
        0.5, ge=0.0,
        description=(
            "Same softmax weight, applied at tree-expansion time. Keeps some "
            "exploration (beta smaller than rollout_bias) so UCB1 still gets "
            "diverse branches to evaluate."
        ),
    )
    mechanism_mcts_distance_penalty_eV: float = Field(
        0.3, ge=0.0,
        description=(
            "eV-equivalent cost per unit of element-count distance between a "
            "rollout endpoint and the target composition. Small positive "
            "value provides a progress signal for UCB1 when state energies "
            "are flat or uncomputed; keep much smaller than typical "
            "adsorption energies so real energetics still dominate."
        ),
    )

    class Config:
        arbitrary_types_allowed = True

    @classmethod
    def create_from_args(cls, **kwargs):
        return cls(**kwargs)

    @classmethod
    def from_kwargs(cls, **kwargs):
        return cls(**kwargs)
