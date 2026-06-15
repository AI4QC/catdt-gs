"""Typed data model for CatDT catalyst surface pipeline.

Central classes:
- StructureEndpoint: ASE Atoms + index metadata, always in sync
- ReactionStep: typed replacement for step_structures Dict[str, Any]
- CatalystSurface: one surface through the entire pipeline — SINGLE SOURCE OF TRUTH
- IntermediateAdsorbate: one adsorbate intermediate (formula + indices + positions)
- Agent*Result: per-agent typed result containers

Post-Phase-1 (2026-04) refactor: all physical state (structures, indices, file paths)
lives on CatalystSurface. WorkflowState (schemas.py) only holds LLM iteration state
and is reached via ``surface.state``. After ``surface.freeze()`` (called at the end
of agent3), ``surface_indices``, ``first_adsorbate_indices``, and ``top_layer_indices``
are immutable — writes raise RuntimeError.
"""

from __future__ import annotations

import json
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from ase import Atoms
from ase.io import write as ase_write

from camel_agents.schemas import PathwayDesign, ReactionContext, ValidationReport


# Per-agent subdirectory layout under surface.output_dir/
AGENT_DIR_LAYOUT: Dict[str, str] = {
    "agent1": "01_surface",
    "agent2": "02_adsorption",
    "agent3": "03_reconstruct",
    "agent3a": "03_reconstruct",
    "agent3b": "03_reconstruct",
    "agent45": "04_pathway",
    "agent4": "04_pathway",
    "agent5": "04_pathway",
    "agent6": "05_neb",
    "agent7": "06_kmc",
    "checkpoints": "checkpoints",
    "visualization": "visualization",
}


# Fields that become immutable after CatalystSurface.freeze()
_FROZEN_FIELDS = frozenset({
    "surface_indices",
    "first_adsorbate_indices",
    "top_layer_indices",
    "n_surface_atoms",
})


# ── Structure Endpoint ───────────────────────────────────────────────────


@dataclass
class StructureEndpoint:
    """ASE Atoms with index metadata that survives copy/relax/sort."""

    atoms: Atoms
    surface_indices: List[int]
    adsorbate_indices: List[int]
    staged_indices: List[int] = field(default_factory=list)

    @property
    def n_surface(self) -> int:
        return len(self.surface_indices)

    @property
    def n_adsorbate(self) -> int:
        return len(self.adsorbate_indices)

    @property
    def n_staged(self) -> int:
        return len(self.staged_indices)

    def copy(self) -> StructureEndpoint:
        ep = StructureEndpoint(
            atoms=self.atoms.copy(),
            surface_indices=list(self.surface_indices),
            adsorbate_indices=list(self.adsorbate_indices),
            staged_indices=list(self.staged_indices),
        )
        ep.sync_to_atoms_info()
        return ep

    def sync_to_atoms_info(self) -> None:
        self.atoms.info["surface_indices"] = list(self.surface_indices)
        self.atoms.info["adsorbate_indices"] = list(self.adsorbate_indices)
        if self.staged_indices:
            self.atoms.info["staged_indices"] = list(self.staged_indices)

    @classmethod
    def from_atoms(cls, atoms: Atoms) -> StructureEndpoint:
        return cls(
            atoms=atoms.copy(),
            surface_indices=list(atoms.info.get("surface_indices", [])),
            adsorbate_indices=list(atoms.info.get("adsorbate_indices", [])),
            staged_indices=list(atoms.info.get("staged_indices", [])),
        )

    @classmethod
    def from_parts(
        cls,
        atoms: Atoms,
        surface_indices: List[int],
        adsorbate_indices: List[int],
        staged_indices: Optional[List[int]] = None,
    ) -> StructureEndpoint:
        ep = cls(
            atoms=atoms.copy(),
            surface_indices=list(surface_indices),
            adsorbate_indices=list(adsorbate_indices),
            staged_indices=list(staged_indices or []),
        )
        ep.sync_to_atoms_info()
        return ep

    def formula(self) -> str:
        return self.atoms.get_chemical_formula("hill")

    def __len__(self) -> int:
        return len(self.atoms)

    @property
    def pymatgen_structure(self):
        from pymatgen.io.ase import AseAtomsAdaptor
        return AseAtomsAdaptor.get_structure(self.atoms)


# ── Reaction Step ────────────────────────────────────────────────────────


@dataclass
class ReactionStep:
    """Typed replacement for step_structures Dict[str, Any]."""

    name: str
    reactant: StructureEndpoint
    product: StructureEndpoint
    reactant_formula: str = ""
    product_formula: str = ""
    reaction_type: str = ""
    staging_hints: List[Dict[str, Any]] = field(default_factory=list)

    # Filled after NEB
    neb_result: Any = None
    ea_forward: Optional[float] = None
    ea_reverse: Optional[float] = None
    reaction_energy: Optional[float] = None
    neb_converged: Optional[bool] = None
    neb_trajectory_path: Optional[str] = None
    neb_interpolation_gif: Optional[str] = None
    neb_energy_profile: Optional[List[float]] = None

    # Filled after pre-NEB relax
    pre_relax_reactant_path: Optional[str] = None
    pre_relax_product_path: Optional[str] = None
    pre_relax_reactant_energy: Optional[float] = None
    pre_relax_product_energy: Optional[float] = None

    # Filled after energy gate
    energy_gate_status: Optional[str] = None
    energy_gate_report: Optional[Dict[str, Any]] = None

    # Visualization paths
    reactant_vasp_path: Optional[str] = None
    product_vasp_path: Optional[str] = None
    reactant_image_path: Optional[str] = None
    product_image_path: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        self.reactant.sync_to_atoms_info()
        self.product.sync_to_atoms_info()
        return {
            "name": self.name,
            "reactant": self.reactant.atoms,
            "product": self.product.atoms,
            "reactant_formula": self.reactant_formula,
            "product_formula": self.product_formula,
            "reaction_type": self.reaction_type,
            "reactant_adsorbate_indices": list(self.reactant.adsorbate_indices),
            "product_adsorbate_indices": list(self.product.adsorbate_indices),
            "reactant_fixed_atom_count": self.reactant.n_surface,
            "product_fixed_atom_count": self.product.n_surface,
            "reactant_staged_indices": list(self.reactant.staged_indices),
            "product_staged_indices": list(self.product.staged_indices),
            "reactant_fixed_reference_signature": "",
            "product_fixed_reference_signature": "",
            "staging_hints": list(self.staging_hints),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> ReactionStep:
        r_atoms = d["reactant"]
        p_atoms = d["product"]
        if not isinstance(r_atoms, Atoms) or not isinstance(p_atoms, Atoms):
            raise TypeError("reactant and product must be ASE Atoms")

        r_ads = list(d.get("reactant_adsorbate_indices", []))
        p_ads = list(d.get("product_adsorbate_indices", []))
        r_staged = list(d.get("reactant_staged_indices", []))
        p_staged = list(d.get("product_staged_indices", []))

        r_surf = [i for i in range(len(r_atoms)) if i not in set(r_ads) and i not in set(r_staged)]
        p_surf = [i for i in range(len(p_atoms)) if i not in set(p_ads) and i not in set(p_staged)]

        reactant = StructureEndpoint(
            atoms=r_atoms.copy() if hasattr(r_atoms, "copy") else r_atoms,
            surface_indices=r_surf,
            adsorbate_indices=r_ads,
            staged_indices=r_staged,
        )
        product = StructureEndpoint(
            atoms=p_atoms.copy() if hasattr(p_atoms, "copy") else p_atoms,
            surface_indices=p_surf,
            adsorbate_indices=p_ads,
            staged_indices=p_staged,
        )

        return cls(
            name=str(d.get("name", "")),
            reactant=reactant,
            product=product,
            reactant_formula=str(d.get("reactant_formula", "")),
            product_formula=str(d.get("product_formula", "")),
            reaction_type=str(d.get("reaction_type", "")),
            staging_hints=list(d.get("staging_hints", [])),
        )


# ── Agent Result Containers ──────────────────────────────────────────────


@dataclass
class Agent1Result:
    """Surface generation (SurFF) result."""
    clean_slab: Atoms
    surface_indices: List[int] = field(default_factory=list)
    bulk_path: Optional[str] = None
    miller_index: Optional[Tuple[int, ...]] = None
    surface_energy: Optional[float] = None
    slab_path: Optional[str] = None
    result_pkl_path: Optional[str] = None
    wulff_surfaces: List[Dict[str, Any]] = field(default_factory=list)
    visualization_paths: Dict[str, str] = field(default_factory=dict)


@dataclass
class Agent2Result:
    """Adsorption site prediction (AdsorbDiff) result."""
    best_structure: Atoms
    adsorption_energy: float
    surface_indices: List[int] = field(default_factory=list)
    adsorbate_indices: List[int] = field(default_factory=list)
    num_sites_tried: int = 0
    best_site_index: int = 0
    all_energies: List[float] = field(default_factory=list)
    all_configurations: List[Atoms] = field(default_factory=list)
    result_pkl_path: Optional[str] = None
    best_config_path: Optional[str] = None
    llm_review: Optional[Dict[str, Any]] = None
    visualization_paths: Dict[str, str] = field(default_factory=dict)


@dataclass
class Agent3aResult:
    """Clean slab UMA relaxation result."""
    relaxed_slab: Atoms
    energy: float = 0.0
    converged: bool = False
    optimization_steps: int = 0
    surface_indices: List[int] = field(default_factory=list)
    relaxed_slab_path: Optional[str] = None


@dataclass
class Agent3bResult:
    """MC surface reconstruction result."""
    reconstructed_surface: Atoms
    mc_energy: Optional[float] = None
    mc_sweeps: int = 0
    mc_converged: bool = False

    # Index tracking through MC
    pre_mc_surface_indices: List[int] = field(default_factory=list)
    pre_mc_adsorbate_indices: List[int] = field(default_factory=list)
    post_mc_surface_indices: List[int] = field(default_factory=list)
    post_mc_adsorbate_indices: List[int] = field(default_factory=list)
    mc_added_adatom_indices: List[int] = field(default_factory=list)

    # File paths
    result_pkl_path: Optional[str] = None
    reconstructed_path: Optional[str] = None
    mc_log_path: Optional[str] = None
    trajectory_path: Optional[str] = None

    # Visualization
    mc_trajectory_gif: Optional[str] = None
    visualization_paths: Dict[str, str] = field(default_factory=dict)


@dataclass
class Agent3Result:
    """Combined Agent3a + Agent3b results."""
    agent3a: Optional[Agent3aResult] = None
    agent3b: Optional[Agent3bResult] = None

    @property
    def reconstructed_surface(self) -> Optional[Atoms]:
        if self.agent3b is not None:
            return self.agent3b.reconstructed_surface
        if self.agent3a is not None:
            return self.agent3a.relaxed_slab
        return None

    @property
    def post_mc_surface_indices(self) -> List[int]:
        if self.agent3b is not None:
            return self.agent3b.post_mc_surface_indices
        if self.agent3a is not None:
            return self.agent3a.surface_indices
        return []

    @property
    def post_mc_adsorbate_indices(self) -> List[int]:
        if self.agent3b is not None:
            return self.agent3b.post_mc_adsorbate_indices
        return []


@dataclass
class PreplacedIntermediate:
    """A single preplaced intermediate structure from tool-based adsorption."""
    label: str
    structure: Atoms
    energy: Optional[float] = None
    surface_indices: List[int] = field(default_factory=list)
    adsorbate_indices: List[int] = field(default_factory=list)
    source: str = ""
    vasp_path: Optional[str] = None


@dataclass
class IntermediateAdsorbate:
    """One adsorbate intermediate on a frozen catalyst surface.

    The catalyst surface itself (its atom count, order, and canonical indices)
    is frozen after agent3. Each reaction-pathway intermediate is stored as a
    snapshot of the full structure (frozen surface atoms + this intermediate's
    adsorbate atoms), with ``adsorbate_indices`` pointing into that snapshot.

    ``parent_formula`` establishes the serial chain: each subsequent intermediate
    is derived from the previous step's clean product.
    """

    formula: str                              # e.g. "*CO", "*CHO", "CH4(g)"
    atoms: Atoms                              # full snapshot (surface + this adsorbate)
    adsorbate_indices: List[int]              # indices of adsorbate atoms in ``atoms``
    positions: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    energy: Optional[float] = None
    source: str = ""                          # "agent2" | "agent45_stageA" | ...
    parent_formula: Optional[str] = None      # previous intermediate in the chain
    file_path: Optional[Path] = None          # on-disk reference (VASP / xyz)

    def __post_init__(self) -> None:
        if self.positions.size == 0 and self.adsorbate_indices:
            pos = self.atoms.get_positions()
            self.positions = np.array([pos[i] for i in self.adsorbate_indices], dtype=float)

    @property
    def n_adsorbate(self) -> int:
        return len(self.adsorbate_indices)

    def symbols(self) -> List[str]:
        syms = self.atoms.get_chemical_symbols()
        return [syms[i] for i in self.adsorbate_indices if 0 <= i < len(syms)]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "formula": self.formula,
            "adsorbate_indices": list(self.adsorbate_indices),
            "symbols": self.symbols(),
            "positions": self.positions.tolist(),
            "energy": self.energy,
            "source": self.source,
            "parent_formula": self.parent_formula,
            "file_path": str(self.file_path) if self.file_path else None,
            "n_atoms_total": len(self.atoms),
        }


@dataclass
class Agent45Result:
    """Pathway design (Agent4) + validation (Agent5) result."""
    pathway_design: Optional[PathwayDesign] = None
    validation_report: Optional[ValidationReport] = None
    steps: List[ReactionStep] = field(default_factory=list)
    iteration_count: int = 0

    # Preplaced intermediates
    preplaced_intermediates: Dict[str, PreplacedIntermediate] = field(default_factory=dict)

    # Per-step energy gate reports
    energy_gate_reports: List[Dict[str, Any]] = field(default_factory=list)

    # Prompt snapshots (agent_name -> prompt text)
    prompt_snapshots: Dict[str, str] = field(default_factory=dict)

    # Iteration history
    iteration_history: List[Dict[str, Any]] = field(default_factory=list)

    # Visualization
    baseline_export_dir: Optional[str] = None
    step_visualization_dir: Optional[str] = None
    visualization_paths: Dict[str, str] = field(default_factory=dict)


@dataclass
class Agent6Result:
    """NEB barrier calculation results."""
    neb_results: Dict[str, Any] = field(default_factory=dict)
    barriers_forward: Dict[str, float] = field(default_factory=dict)
    barriers_reverse: Dict[str, float] = field(default_factory=dict)
    reaction_energies: Dict[str, float] = field(default_factory=dict)
    converged: Dict[str, bool] = field(default_factory=dict)
    rate_determining_step: Optional[str] = None
    max_barrier: Optional[float] = None

    # File paths
    neb_summary_json_path: Optional[str] = None
    neb_output_dir: Optional[str] = None
    trajectory_paths: Dict[str, str] = field(default_factory=dict)
    interpolation_gif_paths: Dict[str, str] = field(default_factory=dict)
    energy_profile_paths: Dict[str, str] = field(default_factory=dict)
    energy_profiles: Dict[str, List[float]] = field(default_factory=dict)

    # Pre-NEB relaxation
    pre_relax_output_dir: Optional[str] = None
    pre_relax_paths: Dict[str, Dict[str, str]] = field(default_factory=dict)


@dataclass
class AdsorptionEnergiesResult:
    """Adsorption energy computation for all intermediates."""
    adsorbate_energies: Dict[str, float] = field(default_factory=dict)
    adsorption_energies: Dict[str, float] = field(default_factory=dict)
    best_configurations: Dict[str, Atoms] = field(default_factory=dict)
    surface_energy: Optional[float] = None
    output_dir: Optional[str] = None
    energy_diagram_path: Optional[str] = None
    visualization_paths: Dict[str, str] = field(default_factory=dict)


@dataclass
class CompletePathwayInfo:
    """Summary of the full reaction pathway energetics."""
    surface_formula: str = ""
    intermediates: List[str] = field(default_factory=list)
    adsorbate_energies: Dict[str, float] = field(default_factory=dict)
    adsorption_energies: Dict[str, float] = field(default_factory=dict)
    rate_determining_step: Optional[str] = None
    max_barrier: Optional[float] = None
    overall_reaction_energy: Optional[float] = None
    result_pkl_path: Optional[str] = None


@dataclass
class Agent7Result:
    """Microkinetic modeling (CatMAP/KMC) result."""
    tof: Optional[float] = None
    production_rates: Dict[str, float] = field(default_factory=dict)
    selectivity: Dict[str, float] = field(default_factory=dict)
    coverages: Dict[str, float] = field(default_factory=dict)

    # CatMAP specific
    catmap_result_path: Optional[str] = None
    mkm_path: Optional[str] = None
    energies_txt_path: Optional[str] = None

    # KMC trajectory
    state_history: List[str] = field(default_factory=list)
    time_history: List[float] = field(default_factory=list)

    # Visualization
    kmc_dynamics_gif: Optional[str] = None
    energy_diagram_path: Optional[str] = None
    visualization_paths: Dict[str, str] = field(default_factory=dict)


# ── Central Surface Object ───────────────────────────────────────────────


@dataclass
class CatalystSurface:
    """One catalyst surface — THE single source of truth through the entire pipeline.

    Lifecycle:
      1. Created at pipeline start from a clean slab.
      2. Agent1 records facet/miller/surface-energy.
      3. Agent2 attaches the first adsorbate → stored as an IntermediateAdsorbate.
      4. Agent3 (VSSR-MC reconstruction + UMA relax) produces the canonical surface
         geometry. At the end of agent3 we MUST call ``freeze(top_layer_indices=...)``:
         - ``surface_indices``: immutable catalyst atoms (indices 0..N_surface-1)
         - ``first_adsorbate_indices``: immutable initial adsorbate indices
         - ``top_layer_indices``: subset of surface_indices that stay relaxable in
           downstream optimizations (the two topmost physical layers). Computed
           from the relaxed slab by z-layer clustering.
      5. Agent4/5 extend ``intermediate_adsorbates`` sequentially — each new
         intermediate references its parent.
      6. Agent6 runs NEB using pairs of intermediates + staging.
      7. Agent7 runs kinetics.

    Any write to a frozen field after ``freeze()`` raises RuntimeError.
    """

    surface_id: str
    endpoint: StructureEndpoint
    clean_slab: Optional[Atoms] = None
    output_dir: Optional[str] = None

    # ── File path pointers (mutable — track canonical structure files on disk) ──
    # These were previously in WorkflowState but are physical artefacts, so they
    # belong to the surface. Set by agents as they run.
    surface_path: Optional[str] = None                     # Agent1 slab input
    clean_slab_path: Optional[str] = None                  # Agent3a relaxed clean slab
    surface_with_adsorbate_path: Optional[str] = None      # Agent2 best adsorbed config
    reconstructed_surface_path: Optional[str] = None       # Agent3b reconstructed surface

    # ── Frozen after agent3 ──
    first_adsorbate_indices: List[int] = field(default_factory=list)
    top_layer_indices: List[int] = field(default_factory=list)
    n_surface_atoms: int = 0
    _frozen: bool = field(default=False, repr=False)

    # ── Intermediate adsorbates (appended sequentially by agent2, agent4/5) ──
    # key = intermediate formula (e.g. "*CO"), value = IntermediateAdsorbate
    intermediate_adsorbates: Dict[str, IntermediateAdsorbate] = field(default_factory=dict)
    intermediate_order: List[str] = field(default_factory=list)

    # Reaction context
    reaction_context: Optional[ReactionContext] = None
    reaction_description: str = ""
    intermediates: List[str] = field(default_factory=list)

    # Reaction steps (built sequentially by Agent4/5)
    reaction_steps: List[ReactionStep] = field(default_factory=list)

    # Pathway energetics summary
    pathway_info: Optional[CompletePathwayInfo] = None
    adsorption_energies_result: Optional[AdsorptionEnergiesResult] = None

    # Agent results — populated as pipeline progresses
    agent1: Optional[Agent1Result] = None
    agent2: Optional[Agent2Result] = None
    agent3: Optional[Agent3Result] = None
    agent45: Optional[Agent45Result] = None
    agent6: Optional[Agent6Result] = None
    agent7: Optional[Agent7Result] = None

    # Workflow state — LLM iteration state only (WorkflowState from schemas.py)
    state: Any = None

    # ── Frozen-field enforcement ──

    def __setattr__(self, name: str, value: Any) -> None:
        if name in _FROZEN_FIELDS and getattr(self, "_frozen", False):
            current = getattr(self, name, None)
            if current != value:
                raise RuntimeError(
                    f"CatalystSurface.{name} is frozen after agent3; "
                    f"attempted to change {current!r} → {value!r}"
                )
        object.__setattr__(self, name, value)

    # ── Convenience accessors ──

    @property
    def atoms(self) -> Atoms:
        return self.endpoint.atoms

    @property
    def surface_indices(self) -> List[int]:
        return self.endpoint.surface_indices

    @property
    def adsorbate_indices(self) -> List[int]:
        """Indices of the CURRENT endpoint's adsorbate.

        Note: after freeze(), ``first_adsorbate_indices`` holds the canonical
        initial adsorbate. Use that for "what was agent2's adsorbate".
        """
        return self.endpoint.adsorbate_indices

    @property
    def atoms_frozen_reference(self) -> Atoms:
        """The canonical post-agent3 atoms snapshot (surface + first adsorbate)."""
        return self.endpoint.atoms

    @property
    def is_frozen(self) -> bool:
        return self._frozen

    @property
    def current_intermediate(self) -> Optional[IntermediateAdsorbate]:
        """Most recently appended intermediate (the rolling "current" state)."""
        if not self.intermediate_order:
            return None
        return self.intermediate_adsorbates.get(self.intermediate_order[-1])

    @property
    def pymatgen_structure(self):
        return self.endpoint.pymatgen_structure

    @property
    def barriers(self) -> Dict[str, float]:
        if self.agent6 is not None:
            return dict(self.agent6.barriers_forward)
        return {}

    @property
    def tof(self) -> Optional[float]:
        if self.agent7 is not None:
            return self.agent7.tof
        return None

    @property
    def production_rates(self) -> Dict[str, float]:
        if self.agent7 is not None:
            return dict(self.agent7.production_rates)
        return {}

    # ── Freeze (called at end of agent3) ──

    def freeze(
        self,
        top_layer_indices: List[int],
        first_adsorbate_indices: Optional[List[int]] = None,
    ) -> None:
        """Lock surface_indices, first_adsorbate_indices, and top_layer_indices.

        After this call any attempt to reassign these fields raises RuntimeError.
        The frozen endpoint atoms are also copied into ``clean_slab`` if absent.

        Args:
            top_layer_indices: indices of atoms to keep relaxable during
                downstream optimizations (typically the two topmost physical
                layers, computed via z-layer clustering on the post-MC slab).
            first_adsorbate_indices: initial adsorbate indices. If None, uses
                the current endpoint.adsorbate_indices.
        """
        if self._frozen:
            return  # idempotent

        ads = list(first_adsorbate_indices) if first_adsorbate_indices is not None \
            else list(self.endpoint.adsorbate_indices)

        object.__setattr__(self, "first_adsorbate_indices", list(ads))
        object.__setattr__(self, "top_layer_indices", list(top_layer_indices))
        object.__setattr__(self, "n_surface_atoms", len(self.endpoint.surface_indices))
        object.__setattr__(self, "_frozen", True)

        if self.clean_slab is None:
            object.__setattr__(self, "clean_slab", self.endpoint.atoms.copy())

    def assert_frozen(self) -> None:
        if not self._frozen:
            raise RuntimeError(
                "CatalystSurface has not been frozen yet. Call surface.freeze(...) "
                "at the end of agent3 before any pathway / NEB operation."
            )

    # ── Intermediate management ──

    def add_intermediate(
        self,
        formula: str,
        atoms: Atoms,
        adsorbate_indices: List[int],
        source: str,
        parent_formula: Optional[str] = None,
        energy: Optional[float] = None,
        file_path: Optional[Union[str, Path]] = None,
    ) -> IntermediateAdsorbate:
        """Register a new intermediate adsorbate in the serial chain.

        Any subsequent re-add of the same formula overwrites the entry (but
        preserves its position in ``intermediate_order``).
        """
        if self._frozen and len(atoms) < self.n_surface_atoms:
            raise RuntimeError(
                f"Intermediate '{formula}': structure has {len(atoms)} atoms, "
                f"fewer than the frozen catalyst atom count ({self.n_surface_atoms})."
            )

        fp: Optional[Path] = Path(file_path) if file_path else None
        entry = IntermediateAdsorbate(
            formula=formula,
            atoms=atoms.copy() if hasattr(atoms, "copy") else atoms,
            adsorbate_indices=list(adsorbate_indices),
            energy=energy,
            source=source,
            parent_formula=parent_formula,
            file_path=fp,
        )

        if formula not in self.intermediate_adsorbates:
            self.intermediate_order.append(formula)
        self.intermediate_adsorbates[formula] = entry
        return entry

    def get_intermediate(self, formula: str) -> IntermediateAdsorbate:
        if formula not in self.intermediate_adsorbates:
            raise KeyError(
                f"Intermediate '{formula}' not registered. "
                f"Available: {list(self.intermediate_adsorbates)}"
            )
        return self.intermediate_adsorbates[formula]

    def has_intermediate(self, formula: str) -> bool:
        return formula in self.intermediate_adsorbates

    def list_intermediates(self) -> List[str]:
        return list(self.intermediate_order)

    def iter_intermediate_chain(self) -> List[IntermediateAdsorbate]:
        return [self.intermediate_adsorbates[f] for f in self.intermediate_order]

    # ── File I/O methods (agent-indexed directory layout) ──

    def agent_dir(self, agent_key: str, subdir: Optional[str] = None) -> Path:
        """Return the output directory for a given agent, creating it if needed.

        Layout under ``self.output_dir``:
            01_surface/     agent1
            02_adsorption/  agent2
            03_reconstruct/ agent3 (agent3a + agent3b)
            04_pathway/     agent4 / agent5 / agent45
            05_neb/         agent6
            06_kmc/         agent7
            checkpoints/    pipeline checkpoints
            visualization/  cross-agent visualizations
        """
        if self.output_dir is None:
            raise RuntimeError(
                "CatalystSurface.output_dir is not set; cannot resolve agent_dir."
            )
        name = AGENT_DIR_LAYOUT.get(agent_key, agent_key)
        path = Path(self.output_dir) / name
        if subdir:
            path = path / subdir
        path.mkdir(parents=True, exist_ok=True)
        return path

    def write_structure(
        self,
        agent_key: str,
        filename: str,
        atoms: Optional[Atoms] = None,
        subdir: Optional[str] = None,
        fmt: Optional[str] = None,
    ) -> Path:
        """Write an ASE Atoms object under the agent's output directory."""
        target_dir = self.agent_dir(agent_key, subdir=subdir)
        path = target_dir / filename
        atoms_to_write = atoms if atoms is not None else self.endpoint.atoms
        ase_write(str(path), atoms_to_write, format=fmt) if fmt else ase_write(str(path), atoms_to_write)
        return path

    def write_json(
        self,
        agent_key: str,
        filename: str,
        payload: Dict[str, Any],
        subdir: Optional[str] = None,
    ) -> Path:
        target_dir = self.agent_dir(agent_key, subdir=subdir)
        path = target_dir / filename
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=_json_default)
        return path

    def register_intermediate_file(
        self,
        formula: str,
        atoms: Atoms,
        adsorbate_indices: List[int],
        source: str,
        parent_formula: Optional[str] = None,
        agent_key: str = "agent45",
        filename: Optional[str] = None,
    ) -> IntermediateAdsorbate:
        """Convenience: add intermediate + write its VASP file under agent_key dir."""
        safe = _sanitize_filename(formula)
        fname = filename or f"intermediate_{len(self.intermediate_order):02d}_{safe}.vasp"
        path = self.write_structure(agent_key, fname, atoms=atoms, subdir="intermediates")
        return self.add_intermediate(
            formula=formula,
            atoms=atoms,
            adsorbate_indices=adsorbate_indices,
            source=source,
            parent_formula=parent_formula,
            file_path=path,
        )

    def intermediates_summary(self) -> Dict[str, Any]:
        return {
            "order": list(self.intermediate_order),
            "entries": {
                f: self.intermediate_adsorbates[f].to_dict()
                for f in self.intermediate_order
            },
        }

    # ── Lifecycle methods ──

    @classmethod
    def from_slab(
        cls,
        slab: Atoms,
        surface_id: str,
        output_dir: Optional[str] = None,
    ) -> CatalystSurface:
        surface_indices = list(range(len(slab)))
        ep = StructureEndpoint.from_parts(slab, surface_indices, [])
        return cls(
            surface_id=surface_id,
            endpoint=ep,
            clean_slab=slab.copy(),
            output_dir=output_dir,
        )

    def copy(self) -> CatalystSurface:
        new = CatalystSurface(
            surface_id=self.surface_id,
            endpoint=self.endpoint.copy(),
            clean_slab=self.clean_slab.copy() if self.clean_slab is not None else None,
            output_dir=self.output_dir,
            surface_path=self.surface_path,
            clean_slab_path=self.clean_slab_path,
            surface_with_adsorbate_path=self.surface_with_adsorbate_path,
            reconstructed_surface_path=self.reconstructed_surface_path,
            first_adsorbate_indices=list(self.first_adsorbate_indices),
            top_layer_indices=list(self.top_layer_indices),
            n_surface_atoms=self.n_surface_atoms,
            intermediate_adsorbates={
                k: IntermediateAdsorbate(
                    formula=v.formula,
                    atoms=v.atoms.copy(),
                    adsorbate_indices=list(v.adsorbate_indices),
                    positions=np.array(v.positions, dtype=float),
                    energy=v.energy,
                    source=v.source,
                    parent_formula=v.parent_formula,
                    file_path=v.file_path,
                )
                for k, v in self.intermediate_adsorbates.items()
            },
            intermediate_order=list(self.intermediate_order),
            reaction_context=self.reaction_context,
            reaction_description=self.reaction_description,
            intermediates=list(self.intermediates),
            reaction_steps=list(self.reaction_steps),
            pathway_info=self.pathway_info,
            adsorption_energies_result=self.adsorption_energies_result,
            agent1=self.agent1,
            agent2=self.agent2,
            agent3=self.agent3,
            agent45=self.agent45,
            agent6=self.agent6,
            agent7=self.agent7,
            state=self.state,
        )
        object.__setattr__(new, "_frozen", self._frozen)
        return new

    def update_endpoint(
        self,
        atoms: Atoms,
        surface_indices: List[int],
        adsorbate_indices: List[int],
        staged_indices: Optional[List[int]] = None,
    ) -> None:
        """Replace the current endpoint atoms + index partition.

        After ``freeze()`` the **composition** of the surface (multiset of
        atomic symbols marked as surface) must stay identical across calls.
        Individual index values may shift when adsorbate atoms are inserted
        or deleted from the middle of the structure — that's a legitimate
        consequence of Stage A add/remove operations and does not violate the
        Phase-1 invariant that the physical catalyst atoms are fixed.

        Raises RuntimeError if the surface composition changes after freeze.
        """
        if self._frozen:
            from collections import Counter
            old_surf_syms = [
                self.endpoint.atoms.get_chemical_symbols()[i]
                for i in self.endpoint.surface_indices
                if 0 <= i < len(self.endpoint.atoms)
            ]
            new_surf_syms = [
                atoms.get_chemical_symbols()[i]
                for i in surface_indices
                if 0 <= i < len(atoms)
            ]
            old_cnt = Counter(old_surf_syms)
            new_cnt = Counter(new_surf_syms)
            if old_cnt != new_cnt:
                raise RuntimeError(
                    "update_endpoint: surface composition changed after freeze() "
                    f"(was {dict(old_cnt)}, got {dict(new_cnt)}). "
                    "Phase-1 invariant: catalyst atoms are fixed after agent3."
                )
        self.endpoint = StructureEndpoint.from_parts(
            atoms, surface_indices, adsorbate_indices, staged_indices,
        )

    # ── Checkpoint ──

    def save_checkpoint(self, path: Optional[Union[str, Path]] = None) -> Path:
        if path is None:
            path = self.agent_dir("checkpoints") / f"{self.surface_id}.pkl"
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)
        return path

    @classmethod
    def load_checkpoint(cls, path: Union[str, Path]) -> CatalystSurface:
        with open(path, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, cls):
            raise TypeError(f"Expected CatalystSurface, got {type(obj).__name__}")
        return obj

    # ── Summary / export ──

    def summary(self) -> Dict[str, Any]:
        s: Dict[str, Any] = {
            "surface_id": self.surface_id,
            "formula": self.endpoint.formula(),
            "n_atoms": len(self.endpoint),
            "n_surface": self.endpoint.n_surface,
            "n_adsorbate": self.endpoint.n_adsorbate,
            "n_first_adsorbate": len(self.first_adsorbate_indices),
            "n_top_layer": len(self.top_layer_indices),
            "is_frozen": self._frozen,
            "intermediates": list(self.intermediates),
            "intermediate_chain": list(self.intermediate_order),
            "n_reaction_steps": len(self.reaction_steps),
        }
        if self.agent2 is not None:
            s["adsorption_energy"] = self.agent2.adsorption_energy
        if self.agent6 is not None:
            s["barriers"] = dict(self.agent6.barriers_forward)
            s["rate_determining_step"] = self.agent6.rate_determining_step
            s["max_barrier"] = self.agent6.max_barrier
        if self.agent7 is not None:
            s["tof"] = self.agent7.tof
            s["production_rates"] = dict(self.agent7.production_rates)
        for i, step in enumerate(self.reaction_steps):
            key = f"step_{i+1}"
            s[key] = {
                "name": step.name,
                "ea_forward": step.ea_forward,
                "ea_reverse": step.ea_reverse,
                "reaction_energy": step.reaction_energy,
            }
        return s


# ── Module-level helpers ─────────────────────────────────────────────────


def _sanitize_filename(name: str) -> str:
    """Convert an adsorbate formula (e.g. '*CHO', 'CH4(g)') into a filename-safe token."""
    bad = {"*": "s", "(": "_", ")": "", " ": "_", "/": "_", "\\": "_", ":": "_"}
    out = name
    for k, v in bad.items():
        out = out.replace(k, v)
    return out or "intermediate"


def _json_default(obj: Any) -> Any:
    """Default serializer for write_json — handles numpy + Path + Atoms."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, Atoms):
        return {
            "symbols": obj.get_chemical_symbols(),
            "positions": obj.get_positions().tolist(),
            "cell": obj.get_cell().array.tolist(),
            "pbc": list(obj.get_pbc()),
        }
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
