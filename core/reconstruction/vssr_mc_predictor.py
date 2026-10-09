"""
VSSR-MC (Virtual Surface Site Relaxation - Monte Carlo) Predictor

This module provides a wrapper class for the surface-sampling (VSSR-MC) algorithm
to sample surface reconstructions across compositional and configurational spaces.

The VSSR-MC algorithm generates virtual adsorption sites above a pristine surface
and uses Monte Carlo sampling to explore different surface reconstructions by
adding/removing atoms at these sites.

Usage:
    from vssr_mc_predictor import VSSRMCPredictor

    # Initialize predictor
    predictor = VSSRMCPredictor(
        surface_sampling_root="/path/to/surface-sampling",
        model_type="CHGNetNFF",
    )

    # Run sampling - offset_data is automatically generated
    result = predictor.sample(
        surface=slab,  # ASE Atoms or path to structure file
        adsorbates=["O", "Sr", "Ti"],  # elements for surface reconstruction
        total_sweeps=100,
        temperature=1.0,
    )

    print(result.summary())
"""

import os
import sys
import pickle
import tempfile
import shutil
import logging
import json
from typing import Union, Optional, List, Dict, Any, Literal
from dataclasses import dataclass, field
from pathlib import Path
from collections import Counter

import numpy as np


def _patch_potcar_spec_accessor() -> bool:
    """Let pymatgen's PotcarCorrection read emmet's PotcarSpec objects.

    ``PotcarCorrection.get_correction`` does ``{dct.get("titel") ... for dct in potcar_spec}``,
    which assumes the Materials Project returns ``potcar_spec`` as a list of dicts. Current emmet
    returns pydantic ``PotcarSpec`` models instead, so the call raises
    ``AttributeError: 'PotcarSpec' object has no attribute 'get'`` and every Pourbaix reference
    lookup fails. This adds a dict-style accessor to that model; no value is changed, and the
    function is a no-op once the upstream versions agree again.

    Returns True when the accessor is present afterwards.
    """
    try:
        from emmet.core.vasp.calculation import PotcarSpec
    except Exception:
        return False
    if not hasattr(PotcarSpec, "get"):
        PotcarSpec.get = lambda self, key, default=None: getattr(self, key, default)
    return True

# 势函数的唯一配置入口（CATDT_FAIRCHEM_MODEL / _TASK / _MODEL_PATH）。
# 该模块也会被 run_vssr_mc_subprocess.py 以脚本方式导入，故留一条路径兜底。
try:
    from core.fairchem_config import (
        local_checkpoint_candidates,
        resolve_fairchem_model,
    )
except ImportError:  # pragma: no cover - fallback for direct-script execution
    _REPO_ROOT = str(Path(__file__).resolve().parents[2])
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)
    from core.fairchem_config import (
        local_checkpoint_candidates,
        resolve_fairchem_model,
    )


@dataclass
class SampledStructure:
    """Information about a sampled structure"""
    sweep_number: int
    atoms: "ase.Atoms"
    energy: float  # Surface excess energy (eV) - relative to bulk reference
    num_adsorbates: int
    acceptance_rate: float
    total_energy: Optional[float] = None  # CHGNet absolute total energy (eV)

    def __repr__(self):
        formula = self.atoms.get_chemical_formula() if self.atoms else "N/A"
        return (f"SampledStructure(sweep={self.sweep_number}, "
                f"formula={formula}, E={self.energy:.4f} eV, "
                f"n_ads={self.num_adsorbates})")


@dataclass
class VSSRMCResult:
    """Complete VSSR-MC sampling result"""
    surface_name: str
    model_type: str
    total_sweeps: int
    sweep_size: int
    temperature: float
    canonical: bool
    structures: List[SampledStructure]
    energy_history: List[float]
    acceptance_history: List[float]
    adsorption_count_history: List[int]
    best_structure: SampledStructure
    lowest_energy_structure: SampledStructure
    offset_data: Optional[Dict] = None
    output_dir: Optional[str] = None

    def get_structures_by_energy(self, top_n: int = 10) -> List[SampledStructure]:
        """Get structures sorted by energy (lowest first)"""
        sorted_structures = sorted(self.structures, key=lambda x: x.energy)
        return sorted_structures[:top_n]

    def get_unique_structures(self, energy_threshold: float = 0.01) -> List[SampledStructure]:
        """Get unique structures based on energy threshold"""
        if not self.structures:
            return []
        sorted_structures = sorted(self.structures, key=lambda x: x.energy)
        unique = [sorted_structures[0]]
        for s in sorted_structures[1:]:
            if abs(s.energy - unique[-1].energy) > energy_threshold:
                unique.append(s)
        return unique

    def summary(self, top_n: int = 5) -> str:
        """Generate summary of sampling results"""
        lines = [
            "VSSR-MC Surface Reconstruction Sampling Results",
            "=" * 60,
            f"Surface: {self.surface_name}",
            f"Model: {self.model_type}",
            f"Sampling mode: {'Canonical' if self.canonical else 'Semi-grand canonical'}",
            f"Total sweeps: {self.total_sweeps}",
            f"Sweep size: {self.sweep_size}",
            f"Temperature: {self.temperature} kT",
            "",
            f"Lowest energy structure:",
            f"  Surface energy: {self.lowest_energy_structure.energy:.4f} eV",
            f"  Total energy: {self.lowest_energy_structure.total_energy:.4f} eV" if self.lowest_energy_structure.total_energy is not None else "  Total energy: N/A",
            f"  Formula: {self.lowest_energy_structure.atoms.get_chemical_formula()}",
            f"  Adsorbates: {self.lowest_energy_structure.num_adsorbates}",
            "",
            f"Energy statistics:",
            f"  Min: {float(np.min(self.energy_history)):.4f} eV",
            f"  Max: {float(np.max(self.energy_history)):.4f} eV",
            f"  Mean: {float(np.mean(self.energy_history)):.4f} eV",
            f"  Std: {float(np.std(self.energy_history)):.4f} eV",
            "",
            f"Average acceptance rate: {float(np.mean(self.acceptance_history)):.2%}",
            "",
            f"Top {min(top_n, len(self.structures))} structures by energy:",
            "-" * 60,
        ]
        for i, s in enumerate(self.get_structures_by_energy(top_n), 1):
            e_total_str = f", E_total={s.total_energy:.4f} eV" if s.total_energy is not None else ""
            lines.append(f"  {i}. E_surf={s.energy:.4f} eV{e_total_str}, "
                        f"{s.atoms.get_chemical_formula()}, "
                        f"n_ads={s.num_adsorbates}")
        return "\n".join(lines)


class VSSRMCPredictor:
    """
    VSSR-MC Surface Reconstruction Sampler

    This class wraps the surface-sampling (VSSR-MC) algorithm to provide
    a simple interface for sampling surface reconstructions across
    compositional and configurational spaces.

    The predictor automatically generates offset_data from the input structure
    using bulk energy calculations from the NFF model.

    Parameters
    ----------
    surface_sampling_root : str
        Path to the surface-sampling repository root
    model_type : str, default="CHGNetNFF"
        Type of NFF model: "CHGNetNFF", "PaiNN", or "NffScaleMACE"
    model_paths : list of str, optional
        Paths to trained model checkpoints. If not provided and model_type
        is "CHGNetNFF", uses pre-trained CHGNet model.
    device : str, default="cuda"
        Device to use: "cuda" or "cpu"
    cutoff : float, optional
        Neighbor list cutoff. Uses model-specific default if not provided.
    work_dir : str, optional
        Working directory. Creates temp dir if not specified.
    keep_files : bool, default=True
        Whether to keep output files after sampling
    verbose : bool, default=True
        Whether to print progress information
    bulk_energy_overrides : dict, optional
        Explicit bulk energies used instead of MP/model lookups during
        offset_data auto-generation. Keys are formulas/elements; values are
        eV/atom for elements and eV per formula unit for the reference
        compound. (Alternatively, pass a complete ``offset_data`` to
        ``sample()`` to skip auto-generation entirely.)
    """

    DEFAULT_CUTOFFS = {
        "CHGNetNFF": 6.0,
        "NffScaleMACE": 5.0,
        "PaiNN": 5.0,
        "UMA": 12.0,
    }

    # Max distance (Å) from existing surface atoms for locally expanded
    # virtual sites in _expand_virtual_sites_locally.
    _VIRTUAL_SITE_NEIGHBOR_CUTOFF = 4.5

    def __init__(
        self,
        surface_sampling_root: str,
        model_type: Literal["CHGNetNFF", "PaiNN", "NffScaleMACE", "UMA"] = "CHGNetNFF",
        model_paths: Optional[List[str]] = None,
        device: Literal["cuda", "cpu"] = "cuda",
        cutoff: Optional[float] = None,
        work_dir: Optional[str] = None,
        keep_files: bool = True,
        verbose: bool = True,
        # Electrochemical (Pourbaix) parameters
        potential_she: Optional[float] = None,
        ph: Optional[float] = None,
        bulk_energy_overrides: Optional[Dict[str, float]] = None,
    ):
        self.surface_sampling_root = os.path.abspath(surface_sampling_root)
        self.model_type = model_type
        self.model_paths = model_paths or []
        self.device = device
        self.cutoff = cutoff or self.DEFAULT_CUTOFFS.get(model_type, 6.0)
        self.keep_files = keep_files
        self.verbose = verbose
        self.potential_she = potential_she
        self.ph = ph
        self.bulk_energy_overrides = bulk_energy_overrides or {}

        # Set up work directory
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="vssr_mc_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # Add surface-sampling to path
        if self.surface_sampling_root not in sys.path:
            sys.path.insert(0, self.surface_sampling_root)

        # Validate environment
        self._validate_environment()

        # Lazy loading of models and calculators
        self._models = None
        self._calculator = None
        self._logger = None
        self._uma_predictor = None

        if self.verbose:
            print(f"VSSRMCPredictor initialized:")
            print(f"  Surface-sampling root: {self.surface_sampling_root}")
            print(f"  Model type: {self.model_type}")
            print(f"  Device: {self.device}")
            print(f"  Cutoff: {self.cutoff}")
            print(f"  Work dir: {self.work_dir}")

    def _validate_environment(self):
        """Validate the environment and dependencies"""
        if not os.path.exists(self.surface_sampling_root):
            raise FileNotFoundError(
                f"Surface-sampling root not found: {self.surface_sampling_root}"
            )

        mcmc_init = os.path.join(self.surface_sampling_root, "mcmc", "__init__.py")
        if not os.path.exists(mcmc_init):
            raise FileNotFoundError(
                f"mcmc package not found in {self.surface_sampling_root}"
            )

    def _log(self, message: str):
        """Log a message"""
        if self.verbose:
            print(message)

    def _setup_logger(self, log_file: Path) -> logging.Logger:
        """Set up logger for MCMC"""
        from mcmc.utils import setup_logger
        return setup_logger(
            "mcmc",
            log_file,
            level=logging.INFO if self.verbose else logging.WARNING
        )

    def _load_models(self):
        """Load NFF models"""
        if self._models is None:
            self._log("Loading NFF models...")

            from nff.train.builders.model import load_model
            from nff.utils.cuda import get_final_device

            device = get_final_device(self.device)

            # Both CHGNetNFF and NffScaleMACE support loading pre-trained models without paths
            if self.model_type in ["CHGNetNFF", "NffScaleMACE"] and not self.model_paths:
                self._log(f"  Using pre-trained {self.model_type} model")
                # Use higher cutoff (8.0) to avoid isolated atom issues when adsorbates
                # are slightly far from the surface (e.g., after diffusion)
                model = load_model(
                    "",
                    model_type=self.model_type,
                    map_location=device,
                    cutoff=8.0,  # Increased from default 5.0 to handle edge cases
                )
                self._models = [model]
            elif self.model_paths:
                self._models = []
                for path in self.model_paths:
                    model = load_model(
                        path,
                        model_type=self.model_type,
                        map_location=device
                    )
                    self._models.append(model)
            else:
                raise ValueError(
                    f"Model type {self.model_type} requires model_paths to be specified"
                )
            self._log(f"  Loaded {len(self._models)} model(s)")

        return self._models

    def _load_calculator(
        self,
        calc_settings: Optional[Dict] = None,
        potential_she: Optional[float] = None,
        ph: Optional[float] = None,
    ) -> "Calculator":
        """Load surface energy calculator.

        Selects from 4 calculator types based on model_type and electrochemical params:
          CHGNet + thermal  → EnsembleNFFSurface
          CHGNet + Pourbaix → NFFPourbaix
          UMA    + thermal  → UMASurfaceCalculator
          UMA    + Pourbaix → UMAPourbaixCalculator
        """
        is_pourbaix = (potential_she is not None and ph is not None)

        if self.model_type == "UMA":
            calc = self._load_uma_calculator(calc_settings, is_pourbaix, potential_she, ph)
        else:
            calc = self._load_nff_calculator(calc_settings, is_pourbaix, potential_she, ph)

        return calc

    def _load_nff_calculator(self, calc_settings, is_pourbaix, potential_she, ph):
        """Load NFF-based calculator (CHGNetNFF, NffScaleMACE, PaiNN)."""
        from nff.utils.cuda import get_final_device

        device = get_final_device(self.device)
        models = self._load_models()
        if not models:
            raise ValueError("No models loaded.")

        model_units = models[0].units if self.model_type != "PaiNN" else "kcal/mol"

        if is_pourbaix:
            from mcmc.calculators import NFFPourbaix
            from mcmc.pourbaix.atoms import PourbaixAtom
            calc = NFFPourbaix(
                models[0],
                device=device,
                model_units=model_units,
                prediction_units="eV",
            )
            pourbaix_settings = dict(calc_settings or {})
            pourbaix_settings["phi"] = potential_she
            pourbaix_settings["pH"] = ph
            # Auto-generate pourbaix_atoms if not provided:
            # each element gets a default PourbaixAtom with bulk energy from offset_data
            if "pourbaix_atoms" not in pourbaix_settings:
                # Use VSSR-MC's generate_pourbaix_atoms() to derive thermodynamically
                # consistent Pourbaix parameters from pymatgen phase/Pourbaix diagrams.
                pourbaix_settings["pourbaix_atoms"] = self._generate_pourbaix_atoms(
                    elements=list(pourbaix_settings.get("chem_pots", {}).keys()),
                    phi=potential_she,
                    pH=ph,
                )
            calc.set(**pourbaix_settings)
            self._log(f"  Loaded NFFPourbaix (U={potential_she} V, pH={ph})")
        else:
            from mcmc.calculators import EnsembleNFFSurface
            calc = EnsembleNFFSurface(
                models,
                device=device,
                model_units=model_units,
                prediction_units="eV",
                offset_units="eV" if self.model_type != "PaiNN" else "atomic",
            )
            if calc_settings:
                calc.set(**calc_settings)
            self._log(f"  Loaded EnsembleNFFSurface (thermal)")

        return calc

    def _load_uma_calculator(self, calc_settings, is_pourbaix, potential_she, ph):
        """Load UMA (FairChem) based calculator.

        Uses the same import pattern as core/pathway/fairchem_predictor.py:
        local checkpoint → HuggingFace fallback.
        """
        try:
            from core.reconstruction.uma_surface_calculator import (
                UMAPourbaixCalculator,
                UMASurfaceCalculator,
            )
        except ModuleNotFoundError:
            from uma_surface_calculator import (
                UMAPourbaixCalculator,
                UMASurfaceCalculator,
            )

        self._log("Loading UMA (FairChem) model...")

        try:
            from fairchem.core.calculate import FAIRChemCalculator
        except ImportError:
            from fairchem.core import FAIRChemCalculator

        predictor = self._get_uma_predictor()

        base_calc = FAIRChemCalculator(
            predictor, task_name=self._get_fairchem_config().task_name
        )

        offset_data = (calc_settings or {}).get("offset_data")
        chem_pots = (calc_settings or {}).get("chem_pots")

        if is_pourbaix:
            pourbaix_atoms = (calc_settings or {}).get("pourbaix_atoms", {})
            if not pourbaix_atoms:
                # Use VSSR-MC's generate_pourbaix_atoms() for thermodynamically
                # consistent parameters, then convert to dict format for UMA adapter
                all_elems = set()
                if offset_data:
                    all_elems.update(k for k in offset_data.get("bulk_energies", {}) if len(k) <= 2)
                if chem_pots:
                    all_elems.update(k for k in chem_pots if len(k) <= 2)
                # UMAPourbaixCalculator accepts both PourbaixAtom objects and dicts
                pourbaix_atoms = self._generate_pourbaix_atoms(
                    elements=list(all_elems), phi=potential_she, pH=ph,
                )
            temperature_k = 300.0
            if calc_settings and "temperature" in calc_settings:
                temperature_k = calc_settings["temperature"] / 8.617333262e-5
            calc = UMAPourbaixCalculator(
                base_calculator=base_calc,
                pourbaix_atoms=pourbaix_atoms,
                potential_she=potential_she,
                ph=ph,
                temperature_k=temperature_k,
            )
            self._log(f"  UMA Pourbaix (U={potential_she} V, pH={ph})")
        else:
            calc = UMASurfaceCalculator(
                base_calculator=base_calc,
                offset_data=offset_data,
                chem_pots=chem_pots,
            )
            if calc_settings:
                calc.set(**{k: v for k, v in calc_settings.items()
                           if k not in ("offset_data", "chem_pots")})
            self._log(f"  UMA Surface (thermal)")

        return calc

    def _get_uma_predictor(self):
        """Load the UMA (FairChem) predict unit once and cache it on self.

        Used by both _load_uma_calculator() and _compute_energy_uma() so the
        model weights are only loaded a single time per predictor instance.
        """
        if self._uma_predictor is not None:
            return self._uma_predictor

        import torch
        if hasattr(torch.serialization, "add_safe_globals"):
            torch.serialization.add_safe_globals([slice])

        # Import fairchem components (same pattern as fairchem_predictor.py)
        try:
            from fairchem.core.calculate import pretrained_mlip
        except ImportError:
            from fairchem.core import pretrained_mlip
        from fairchem.core.units.mlip_unit import load_predict_unit

        # 势函数选择统一由 core.fairchem_config 解析（单一配置入口）。
        # 没有 CATDT_FAIRCHEM_* 时结果仍是 uma-s-1p1 / oc20 与本地 checkpoint。
        cfg = self._get_fairchem_config()
        uma_model_name = cfg.model_name

        # Try local checkpoints first. Only checkpoints that really belong to
        # `uma_model_name` are considered: otherwise asking for another
        # potential would silently load uma-s-1p1.pt.
        local_candidates = [cfg.model_path] if cfg.model_path else []
        local_candidates.extend(
            local_checkpoint_candidates(uma_model_name, include_legacy=False)
        )
        for path in local_candidates:
            if path and os.path.exists(path):
                try:
                    self._uma_predictor = load_predict_unit(
                        path=path, device=self.device,
                        overrides={"backbone": {"always_use_pbc": False}},
                    )
                    self._log(f"  MLIP from local: {path}")
                    break
                except Exception as e:
                    self._log(f"  Local load failed ({path}): {e}")

        # Fallback to the fairchem pretrained registry / HuggingFace cache.
        # If that fails it raises - it never substitutes another potential.
        if self._uma_predictor is None:
            self._uma_predictor = pretrained_mlip.get_predict_unit(
                uma_model_name, device=self.device,
            )
            self._log(f"  MLIP from fairchem registry: {uma_model_name}")

        return self._uma_predictor

    def _get_fairchem_config(self):
        """Resolved potential selection for this predictor (cached)."""
        cfg = getattr(self, "_fairchem_config", None)
        if cfg is None:
            cfg = resolve_fairchem_model(
                model_path=self.model_paths[0] if self.model_paths else None,
                context="VSSRMCPredictor",
            )
            self._fairchem_config = cfg
        return cfg

    def _generate_pourbaix_atoms(self, elements, phi, pH):
        """Generate PourbaixAtom objects using VSSR-MC's native Pourbaix module.

        Uses pymatgen PhaseDiagram + PourbaixDiagram to derive thermodynamically
        consistent parameters (num_e, num_H, delta_G2_std) for each element.
        """
        from mcmc.pourbaix.atoms import generate_pourbaix_atoms

        # Build phase diagram and Pourbaix diagram from pymatgen/MP data
        from pymatgen.analysis.phase_diagram import PhaseDiagram
        from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixEntry
        from pymatgen.entries.computed_entries import ComputedEntry
        from pymatgen.core import Composition, Element

        clean_elements = sorted(set(e for e in elements if len(e) <= 2 and e not in ("X",)))
        self._log(f"  Generating Pourbaix atoms for elements: {clean_elements}")

        # emmet/pymatgen disagree on the shape of potcar_spec; without this the lookup below
        # raises AttributeError and no Pourbaix reference can be built (see the function's docstring)
        _patch_potcar_spec_accessor()

        try:
            # Try Materials Project API
            from mp_api.client import MPRester
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mpr = MPRester()

                # Get phase diagram entries
                entries = mpr.get_entries_in_chemsys(
                    clean_elements + ["O", "H"],
                )
                pd = PhaseDiagram(entries)

                # Get Pourbaix entries
                pbx_entries = mpr.get_pourbaix_entries(clean_elements)
                pbx = PourbaixDiagram(pbx_entries)

            pourbaix_atoms = generate_pourbaix_atoms(pd, pbx, phi, pH, clean_elements)
            self._log(f"  Generated from MP: { {k: f'n_e={v.num_e}' for k,v in pourbaix_atoms.items()} }")
            return pourbaix_atoms

        except Exception as e:
            # No silent fallback: PourbaixAtom with delta_G2_std=0.0 and
            # atom_std_state_energy=0.0 would degrade the electrochemistry
            # to "-n_e*U only" without any warning to the user.
            raise RuntimeError(
                f"Could not generate Pourbaix atoms from Materials Project "
                f"for elements {clean_elements}: {e}. Provide a Materials "
                f"Project API key (MP_API_KEY) or pass explicit "
                f"'pourbaix_atoms' via calc_settings."
            ) from e

    def _convert_surface(self, surface) -> "ase.Atoms":
        """Convert input surface to ASE Atoms"""
        from ase import Atoms
        from ase.io import read

        if isinstance(surface, str):
            if not os.path.exists(surface):
                raise FileNotFoundError(f"Surface file not found: {surface}")
            return read(surface)
        elif isinstance(surface, Atoms):
            return surface.copy()
        else:
            # Try pymatgen Structure
            try:
                from pymatgen.io.ase import AseAtomsAdaptor
                from pymatgen.core import Structure
                if isinstance(surface, Structure):
                    return AseAtomsAdaptor.get_atoms(surface)
            except ImportError:
                pass
            raise TypeError(
                f"Unsupported surface type: {type(surface)}. "
                "Expected str (file path), ase.Atoms, or pymatgen.Structure"
            )

    def _prepare_slab_batch(
        self,
        atoms: "ase.Atoms",
        device: str,
    ) -> "AtomsBatch":
        """Prepare AtomsBatch from ASE Atoms"""
        from mcmc.utils.misc import get_atoms_batch

        return get_atoms_batch(
            atoms,
            nff_cutoff=self.cutoff,
            device=device,
            props={"energy": 0, "energy_grad": []},
        )

    def _compute_bulk_energy(self, formula: str, per_atom: bool = True) -> float:
        """
        Compute bulk energy for a given formula using the current model
        (NFF for CHGNet/PaiNN/MACE, FairChem for UMA).

        Parameters
        ----------
        formula : str
            Chemical formula (e.g., "Cu", "O2", "SrTiO3")
        per_atom : bool
            If True, return energy per atom

        Returns
        -------
        float
            Bulk energy in eV

        Raises
        ------
        RuntimeError
            If no bulk structure can be obtained or the energy calculation
            fails. A silent 0.0 fallback would corrupt offset_data (every
            surface excess energy would be shifted by N x E_bulk), so the
            caller must either fix the environment (e.g. MP API key) or
            supply explicit values via ``bulk_energy_overrides``.
        """
        from ase.build import bulk as ase_bulk
        import warnings

        atoms = None

        # Try to get bulk structure from Materials Project
        try:
            from mp_api.client import MPRester
            from pymatgen.io.ase import AseAtomsAdaptor
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mpr = MPRester()
                docs = mpr.materials.summary.search(
                    formula=formula,
                    fields=["material_id", "structure", "energy_per_atom"],
                )
                if docs:
                    docs_sorted = sorted(docs, key=lambda x: x.energy_per_atom if x.energy_per_atom else 0)
                    structure = docs_sorted[0].structure
                    atoms = AseAtomsAdaptor.get_atoms(structure)
                    self._log(f"  Got {formula} structure from Materials Project")
        except Exception as e:
            self._log(f"  Materials Project lookup failed for {formula}: {e}")

        # Fall back to simple bulk structures for elements
        if atoms is None:
            try:
                from pymatgen.core import Composition
                comp = Composition(formula)
                if len(comp.elements) == 1:
                    atoms = ase_bulk(formula)
                    self._log(f"  Using ASE bulk structure for {formula}")
                else:
                    raise RuntimeError(
                        f"Could not obtain a bulk structure for compound "
                        f"{formula}: Materials Project lookup failed and no "
                        f"elemental ASE bulk fallback exists for compounds. "
                        f"Provide an MP API key (MP_API_KEY) or pass explicit "
                        f"bulk energies via bulk_energy_overrides / offset_data."
                    )
            except RuntimeError:
                raise
            except Exception as e2:
                raise RuntimeError(
                    f"Could not obtain a bulk structure for {formula}: {e2}. "
                    f"Provide an MP API key (MP_API_KEY) or pass explicit "
                    f"bulk energies via bulk_energy_overrides / offset_data."
                ) from e2

        atoms.pbc = [True, True, True]

        try:
            if self.model_type == "UMA":
                energy = self._compute_energy_uma(atoms)
            else:
                energy = self._compute_energy_nff(atoms)
        except Exception as e:
            raise RuntimeError(
                f"Bulk energy calculation failed for {formula} with model "
                f"{self.model_type}: {e}. Pass explicit bulk energies via "
                f"bulk_energy_overrides / offset_data to bypass the "
                f"calculation."
            ) from e

        if per_atom:
            return energy / len(atoms)
        return energy

    def _compute_energy_nff(self, atoms) -> float:
        """Compute energy using NFF model (CHGNet/PaiNN/MACE)."""
        from nff.utils.cuda import get_final_device
        from mcmc.calculators import EnsembleNFFSurface

        device = get_final_device(self.device)
        models = self._load_models()
        model_units = models[0].units if self.model_type != "PaiNN" else "kcal/mol"

        calc = EnsembleNFFSurface(
            models, device=device, model_units=model_units,
            prediction_units="eV",
            offset_units="eV" if self.model_type != "PaiNN" else "atomic",
        )

        slab_batch = self._prepare_slab_batch(atoms, device)
        energy = calc.get_potential_energy(atoms=slab_batch)
        if hasattr(energy, '__iter__') and not isinstance(energy, str):
            energy = float(energy.item() if hasattr(energy, 'item') else energy[0])
        return float(energy)

    def _compute_energy_uma(self, atoms) -> float:
        """Compute energy using UMA (FairChem) model."""
        try:
            from fairchem.core.calculate import FAIRChemCalculator
        except ImportError:
            from fairchem.core import FAIRChemCalculator

        calc = FAIRChemCalculator(
            self._get_uma_predictor(),
            task_name=self._get_fairchem_config().task_name,
        )
        atoms_copy = atoms.copy()
        atoms_copy.calc = calc
        return float(atoms_copy.get_potential_energy())

    def _auto_generate_offset_data(
        self,
        slab: "ase.Atoms",
        adsorbates: Optional[List[str]] = None,
    ) -> Dict:
        """
        Automatically generate offset_data from the slab composition.

        The offset_data contains:
        - bulk_energies: bulk energy per formula unit for each element and reference compound
        - stoics: stoichiometry of each element in the reference formula
        - ref_formula: reference compound formula (derived from slab)
        - ref_element: reference element for normalization (usually the most common cation)

        Parameters
        ----------
        slab : ase.Atoms
            The pristine slab structure
        adsorbates : list of str, optional
            Additional adsorbate elements

        Returns
        -------
        dict
            offset_data dictionary
        """
        from collections import Counter
        from pymatgen.core import Composition
        from ase.data import chemical_symbols

        self._log("Auto-generating offset_data from slab composition...")

        # Get slab composition
        slab_symbols = slab.get_chemical_symbols()
        composition = Counter(slab_symbols)
        elements = list(composition.keys())

        # Include adsorbate elements if provided
        if adsorbates:
            for ads in adsorbates:
                if ads not in elements:
                    elements.append(ads)

        self._log(f"  Elements: {elements}")

        # Determine reference formula from slab composition
        # Use pymatgen to get reduced formula
        comp = Composition(composition)
        ref_formula = comp.reduced_formula
        # Get stoichiometry from reduced formula, not original composition
        reduced_comp = Composition(ref_formula)
        el_amt = reduced_comp.get_el_amt_dict()
        stoics = {str(el): int(amt) for el, amt in el_amt.items()}
        for el in elements:
            stoics.setdefault(str(el), 0)

        self._log(f"  Reference formula: {ref_formula}")
        self._log(f"  Stoichiometry: {stoics}")

        # Determine reference element (usually a cation that's not O/N/S/F/Cl)
        # Priority: transition metals > main group metals > non-metals
        anions = {"O", "N", "S", "F", "Cl", "Br", "I"}
        cations = [el for el in stoics.keys() if el not in anions]
        if cations:
            ref_element = cations[0]  # First non-anion element
        else:
            ref_element = list(stoics.keys())[0]

        self._log(f"  Reference element: {ref_element}")

        # Calculate bulk energies
        bulk_energies = {}

        # First, compute bulk energy for the reference compound.
        # The consumer (UMASurfaceCalculator._compute_surface_energy /
        # EnsembleNFFSurface) multiplies bulk_energies[ref_formula] by the
        # number of reference-element atoms, so it must be stored per
        # formula unit — not as the total energy of whatever MP cell matched.
        if ref_formula in self.bulk_energy_overrides:
            bulk_energies[ref_formula] = float(self.bulk_energy_overrides[ref_formula])
            self._log(f"    {ref_formula}: {bulk_energies[ref_formula]:.4f} eV/f.u. (override)")
        else:
            self._log(f"  Computing bulk energy for {ref_formula}...")
            ref_energy_per_atom = self._compute_bulk_energy(ref_formula, per_atom=True)
            n_atoms_per_fu = sum(stoics.values())
            bulk_energies[ref_formula] = ref_energy_per_atom * n_atoms_per_fu
            self._log(f"    {ref_formula}: {bulk_energies[ref_formula]:.4f} eV/f.u.")

        # Compute bulk energies for individual elements
        for el in elements:
            if el not in bulk_energies:
                if el in self.bulk_energy_overrides:
                    bulk_energies[el] = float(self.bulk_energy_overrides[el])
                    self._log(f"    {el}: {bulk_energies[el]:.4f} eV/atom (override)")
                    continue
                self._log(f"  Computing bulk energy for {el}...")
                el_energy = self._compute_bulk_energy(el, per_atom=True)
                bulk_energies[el] = el_energy
                self._log(f"    {el}: {el_energy:.4f} eV/atom")

        offset_data = {
            "bulk_energies": bulk_energies,
            "stoics": stoics,
            "ref_formula": ref_formula,
            "ref_element": ref_element,
        }

        self._log(f"  Generated offset_data: {json.dumps(offset_data, indent=2)}")

        return offset_data

    def _auto_generate_chem_pots(
        self,
        slab: "ase.Atoms",
        adsorbates: Optional[List[str]] = None,
        ref_element: Optional[str] = None,
    ) -> Dict[str, float]:
        """
        Generate default chemical potentials.

        By default, sets all chemical potentials to 0, which corresponds to
        elemental reference states.

        Parameters
        ----------
        slab : ase.Atoms
            The slab structure
        adsorbates : list of str, optional
            Adsorbate elements
        ref_element : str, optional
            Reference element (will have chem_pot = 0)

        Returns
        -------
        dict
            Chemical potentials for each element
        """
        composition = Counter(slab.get_chemical_symbols())
        elements = list(composition.keys())

        if adsorbates:
            for ads in adsorbates:
                if ads not in elements:
                    elements.append(ads)

        # Default: all chemical potentials = 0
        chem_pots = {el: 0.0 for el in elements}

        return chem_pots

    def sample(
        self,
        surface: Union[str, "ase.Atoms", "pymatgen.core.Structure"],
        adsorbates: Optional[List[str]] = None,
        canonical: bool = False,
        num_adsorbates: int = 0,
        adsorbate_counts: Optional[Dict[str, int]] = None,
        total_sweeps: int = 100,
        sweep_size: int = 20,
        temperature: float = 1.0,
        perform_annealing: bool = False,
        annealing_alpha: float = 0.99,
        chem_pots: Optional[Dict[str, float]] = None,
        offset_data: Optional[Dict] = None,
        auto_offset: bool = True,
        calc_settings: Optional[Dict] = None,
        system_settings: Optional[Dict] = None,
        output_dir: Optional[str] = None,
        run_name: Optional[str] = None,
        clean_slab: Optional["ase.Atoms"] = None,
        surface_indices: Optional[List[int]] = None,
        adsorbate_indices: Optional[List[int]] = None,
        potential_she: Optional[float] = None,
        ph: Optional[float] = None,
    ) -> VSSRMCResult:
        """
        Perform VSSR-MC sampling on a surface.

        This method samples surface reconstructions across compositional and
        configurational spaces using Monte Carlo with virtual adsorption sites.

        Parameters
        ----------
        surface : str or ase.Atoms or pymatgen.Structure
            Input surface structure (with adsorbate if doing canonical sampling)
        adsorbates : list of str, optional
            List of elements that can be added/removed at virtual sites
        canonical : bool, default=False
            If True, perform canonical (NVT) sampling with fixed composition.
            If False, perform semi-grand canonical sampling with variable composition.
        num_adsorbates : int, default=0
            Number of adsorbate atoms (required if canonical=True)
        adsorbate_counts : dict, optional
            Fixed canonical composition for virtual-site adsorbates. When
            provided, the total count is inferred from this map.
        total_sweeps : int, default=100
            Number of MC sweeps to perform
        sweep_size : int, default=20
            Number of MC steps per sweep
        temperature : float, default=1.0
            Temperature in kT units
        perform_annealing : bool, default=False
            Whether to perform simulated annealing
        annealing_alpha : float, default=0.99
            Annealing schedule parameter (T_new = alpha * T_old)
        chem_pots : dict, optional
            Chemical potentials for each element. If None, auto-generated with all 0.
        offset_data : dict, optional
            Offset data for surface energy calculation. If None and auto_offset=True,
            automatically generated from slab composition.
        auto_offset : bool, default=True
            If True and offset_data is None, automatically generate offset_data.
        calc_settings : dict, optional
            Additional calculator settings
        system_settings : dict, optional
            Additional system settings (surface_depth, near_reduce, etc.)
        output_dir : str, optional
            Directory to save results
        run_name : str, optional
            Name for this sampling run
        clean_slab : ase.Atoms, optional
            Clean slab without adsorbates, used for virtual site generation.
            If provided, virtual sites are generated on this clean surface
            (preventing sites from being placed above adsorbates).

        Returns
        -------
        VSSRMCResult
            Sampling results including all structures and statistics
        """
        from nff.utils.cuda import get_final_device
        from mcmc import MCMC
        from mcmc.system import SurfaceSystem

        # Convert and prepare surface
        self._log("\n" + "=" * 60)
        self._log("VSSR-MC Surface Reconstruction Sampling")
        self._log("=" * 60)

        surface_atoms = self._convert_surface(surface)
        surface_formula = surface_atoms.get_chemical_formula()
        self._log(f"Surface: {surface_formula}")
        self._log(f"Number of atoms: {len(surface_atoms)}")

        # Set up output directory
        if output_dir is None:
            output_dir = os.path.join(
                self.work_dir,
                run_name or f"vssr_mc_{surface_formula}"
            )
        os.makedirs(output_dir, exist_ok=True)

        # Set up logger
        logger = self._setup_logger(Path(output_dir) / "mc.log")

        # Auto-generate offset_data if needed
        if offset_data is None and auto_offset:
            offset_data = self._auto_generate_offset_data(surface_atoms, adsorbates)

        # Auto-generate chemical potentials if needed
        if chem_pots is None:
            ref_element = offset_data.get("ref_element") if offset_data else None
            chem_pots = self._auto_generate_chem_pots(surface_atoms, adsorbates, ref_element)
            self._log(f"Chemical potentials: {chem_pots}")

        # Prepare calculator settings
        _calc_settings = calc_settings.copy() if calc_settings else {}
        _calc_settings["chem_pots"] = chem_pots
        if offset_data:
            _calc_settings["offset_data"] = offset_data

        # Prepare system settings
        _system_settings = {
            "surface_name": run_name or surface_formula,
            "cutoff": self.cutoff,
            "surface_depth": 2,  # Allow top 2 layers to move during reconstruction
            "near_reduce": 0.01,
            "planar_distance": 1.5,  # 虚拟位点在表面上方的高度
            "no_obtuse_hollow": True,
            "ads_site_type": "all",
        }
        if system_settings:
            _system_settings.update(system_settings)

        # Inject external surface/adsorbate indices (overrides z-based detection)
        if surface_indices is not None and adsorbate_indices is not None:
            _system_settings["external_surface_idx"] = list(surface_indices)
            _system_settings["external_adsorbate_idx"] = list(adsorbate_indices)
            self._log(
                f"Using external indices: {len(surface_indices)} surface, "
                f"{len(adsorbate_indices)} adsorbate atoms"
            )

        # 如果指定了max_site_height，使用它来覆盖planar_distance
        # 这可以防止在吸附物上方生成虚拟位点
        if "max_site_height" in _system_settings:
            max_height = _system_settings.pop("max_site_height")
            _system_settings["planar_distance"] = max_height
            self._log(f"Limiting virtual site height to {max_height} Å above surface")

        # Get device
        device = get_final_device(self.device)

        # Load calculator — select mode based on electrochemical params
        eff_potential = potential_she if potential_she is not None else self.potential_she
        eff_ph = ph if ph is not None else self.ph
        self._log(f"Loading {self.model_type} calculator...")
        calculator = self._load_calculator(
            _calc_settings,
            potential_she=eff_potential,
            ph=eff_ph,
        )

        # 预先生成虚拟位点坐标（如果提供了clean_slab）
        # 使用clean_slab（不含吸附物）生成虚拟位点，避免在吸附物上方生成
        ads_coords = None
        if clean_slab is not None:
            self._log(f"Generating virtual sites on clean slab: {clean_slab.get_chemical_formula()}")

            from pymatgen.io.ase import AseAtomsAdaptor
            from pymatgen.analysis.adsorption import AdsorbateSiteFinder

            pmg_clean = AseAtomsAdaptor.get_structure(clean_slab)
            site_finder = AdsorbateSiteFinder(pmg_clean)
            self._log(f"  planar_distance: {_system_settings.get('planar_distance', 1.5)} Å")

            ads_site_type = _system_settings.get("ads_site_type", "all")
            ads_positions = site_finder.find_adsorption_sites(
                put_inside=True,
                symm_reduce=_system_settings.get("symm_reduce", False),
                near_reduce=_system_settings.get("near_reduce", 0.01),
                distance=_system_settings.get("planar_distance", 1.5),
                no_obtuse_hollow=_system_settings.get("no_obtuse_hollow", True),
            )[ads_site_type]

            ads_coords = ads_positions
            overlayer_radius = float(
                _system_settings.get("existing_atom_exclusion_radius_A", 0.0) or 0.0
            )
            if overlayer_radius > 0:
                excluded = set(int(i) for i in adsorbate_indices or [])
                reference_indices = [
                    idx for idx in range(len(surface_atoms))
                    if idx not in excluded
                ]
                n_before = len(ads_coords)
                ads_coords = self._filter_virtual_sites_near_atoms(
                    np.asarray(ads_coords, dtype=float),
                    reference_atoms=surface_atoms,
                    reference_indices=reference_indices,
                    exclusion_radius=overlayer_radius,
                )
                self._log(
                    f"  Filtered virtual sites near existing surface atoms: "
                    f"{n_before} -> {len(ads_coords)} "
                    f"(radius={overlayer_radius:.2f} Å)"
                )
            min_site_distance = float(
                _system_settings.get("min_virtual_site_distance_A", 0.0) or 0.0
            )
            if min_site_distance > 0:
                n_before = len(ads_coords)
                ads_coords = self._filter_virtual_sites_by_site_spacing(
                    np.asarray(ads_coords, dtype=float),
                    cell=surface_atoms.cell,
                    pbc=surface_atoms.pbc,
                    min_site_distance=min_site_distance,
                )
                self._log(
                    f"  Filtered virtual sites by site-site spacing: "
                    f"{n_before} -> {len(ads_coords)} "
                    f"(min={min_site_distance:.2f} Å)"
                )
            max_surface_distance = float(
                _system_settings.get("max_virtual_site_distance_to_surface_A", 0.0) or 0.0
            )
            if max_surface_distance > 0:
                excluded = set(int(i) for i in adsorbate_indices or [])
                reference_indices = [
                    idx for idx in range(len(surface_atoms))
                    if idx not in excluded
                ]
                n_before = len(ads_coords)
                ads_coords = self._filter_virtual_sites_by_max_distance_to_atoms(
                    np.asarray(ads_coords, dtype=float),
                    reference_atoms=surface_atoms,
                    reference_indices=reference_indices,
                    max_distance=max_surface_distance,
                )
                self._log(
                    f"  Filtered virtual sites far from existing surface: "
                    f"{n_before} -> {len(ads_coords)} "
                    f"(max={max_surface_distance:.2f} Å)"
                )
            min_virtual_site_count = int(
                _system_settings.get("virtual_site_min_count", 0) or 0
            )
            local_expansion_radius = float(
                _system_settings.get("virtual_site_local_expansion_radius_A", 0.0)
                or 0.0
            )
            if (
                min_virtual_site_count > 0
                and local_expansion_radius > 0
                and len(ads_coords) < min_virtual_site_count
            ):
                excluded = set(int(i) for i in adsorbate_indices or [])
                reference_indices = [
                    idx for idx in range(len(surface_atoms))
                    if idx not in excluded
                ]
                n_before = len(ads_coords)
                ads_coords = self._expand_virtual_sites_locally(
                    np.asarray(ads_coords, dtype=float),
                    reference_atoms=surface_atoms,
                    reference_indices=reference_indices,
                    radius=local_expansion_radius,
                    min_count=min_virtual_site_count,
                )
                self._log(
                    f"  Expanded local virtual sites: {n_before} -> {len(ads_coords)} "
                    f"(radius={local_expansion_radius:.2f} Å, min={min_virtual_site_count})"
                )
            exclusion_radius = float(
                _system_settings.get("adsorbate_exclusion_radius_A", 0.0) or 0.0
            )
            if exclusion_radius > 0 and adsorbate_indices:
                n_before = len(ads_coords)
                ads_coords = self._filter_virtual_sites_near_adsorbate(
                    np.asarray(ads_coords, dtype=float),
                    surface_atoms=surface_atoms,
                    adsorbate_indices=adsorbate_indices,
                    exclusion_radius=exclusion_radius,
                )
                self._log(
                    f"  Filtered virtual sites near adsorbate: "
                    f"{n_before} -> {len(ads_coords)} "
                    f"(radius={exclusion_radius:.2f} Å)"
                )
            self._log(f"  Generated {len(ads_coords)} virtual sites on clean surface")
            if len(ads_coords) > 0:
                self._log(f"  First virtual site at: {ads_coords[0]}")

        # Prepare AtomsBatch (with adsorbate)
        self._log("Preparing surface system...")
        slab_batch = self._prepare_slab_batch(surface_atoms, device)

        # Create SurfaceSystem
        # 如果提供了ads_coords，SurfaceSystem将使用预先生成的虚拟位点
        surface_system = SurfaceSystem(
            slab_batch,
            calc=calculator,
            ads_coords=ads_coords,  # 使用在clean_slab上生成的虚拟位点
            system_settings=_system_settings,
            save_folder=output_dir,
        )

        # Save initial structure with virtual adsorption sites
        surface_system.all_atoms.write(Path(output_dir) / "all_virtual_ads.cif")

        initial_energy = surface_system.get_surface_energy()
        self._log(f"Initial energy: {float(initial_energy):.4f} eV")
        self._log(f"Number of adsorption sites: {len(surface_system.ads_coords)}")

        # Set up MCMC
        sampling_settings = {
            "canonical": canonical,
            "total_sweeps": total_sweeps,
            "sweep_size": sweep_size,
            "start_temp": temperature,
            "perform_annealing": perform_annealing,
            "alpha": annealing_alpha,
            "run_folder": output_dir,
        }

        mcmc = MCMC(
            adsorbates=adsorbates,
            canonical=canonical,
            num_ads_atoms=num_adsorbates if canonical else 0,
            adsorbate_counts=adsorbate_counts if canonical else None,
        )

        # Run sampling
        self._log(f"\nRunning VSSR-MC sampling...")
        self._log(f"  Mode: {'Canonical' if canonical else 'Semi-grand canonical'}")
        self._log(f"  Total sweeps: {total_sweeps}")
        self._log(f"  Sweep size: {sweep_size}")
        self._log(f"  Temperature: {temperature} kT")
        if adsorbates:
            self._log(f"  Adsorbates: {adsorbates}")

        results = mcmc.run(
            surface=surface_system,
            logger=logger,
            **sampling_settings,
        )

        self._log("Sampling complete!")

        # Process results
        structures = []
        for i, (surf_sys, energy, accept_rate, n_ads) in enumerate(zip(
            results["history"],
            results["energy_hist"],
            results["frac_accept_hist"],
            results["adsorption_count_hist"],
        )):
            if hasattr(surf_sys, 'real_atoms'):
                atoms = surf_sys.real_atoms.copy()
            elif hasattr(surf_sys, 'atoms'):
                atoms = surf_sys.atoms.copy()
            else:
                atoms = surf_sys

            # Compute absolute total energy from surface energy + bulk reference
            total_e = None
            if offset_data:
                bulk_energies = offset_data.get("bulk_energies", {})
                stoics = offset_data.get("stoics", {})
                ref_element = offset_data.get("ref_element")
                if ref_element and ref_element in bulk_energies:
                    atom_counts = Counter(atoms.get_chemical_symbols())
                    n_ref = atom_counts.get(ref_element, 0)
                    e_bulk_ref = bulk_energies[ref_element]
                    total_e = float(energy) + n_ref * e_bulk_ref

            sampled = SampledStructure(
                sweep_number=i + 1,
                atoms=atoms,
                energy=float(energy),
                num_adsorbates=n_ads,
                acceptance_rate=accept_rate,
                total_energy=total_e,
            )
            structures.append(sampled)

        # Find best structures
        lowest_energy_idx = np.argmin(results["energy_hist"])
        lowest_energy_structure = structures[lowest_energy_idx]

        # Final structure
        best_structure = structures[-1]

        # Save structures
        with open(Path(output_dir) / f"{len(structures)}_mcmc_structures.pkl", "wb") as f:
            pickle.dump([s.atoms for s in structures], f)

        # Save offset_data
        if offset_data:
            with open(Path(output_dir) / "offset_data.json", "w") as f:
                json.dump(offset_data, f, indent=2)

        # Save settings
        all_settings = {
            "system_settings": _system_settings,
            "sampling_settings": sampling_settings,
            "calc_settings": {k: v for k, v in _calc_settings.items() if k != "offset_data"},
            "offset_data": offset_data,
        }
        with open(Path(output_dir) / "settings.json", "w") as f:
            json.dump(all_settings, f, indent=2, default=str)

        result = VSSRMCResult(
            surface_name=_system_settings["surface_name"],
            model_type=self.model_type,
            total_sweeps=total_sweeps,
            sweep_size=sweep_size,
            temperature=temperature,
            canonical=canonical,
            structures=structures,
            energy_history=results["energy_hist"],
            acceptance_history=results["frac_accept_hist"],
            adsorption_count_history=results["adsorption_count_hist"],
            best_structure=best_structure,
            lowest_energy_structure=lowest_energy_structure,
            offset_data=offset_data,
            output_dir=output_dir,
        )

        self._log(f"\n{result.summary()}")
        self._log(f"\nResults saved to: {output_dir}")

        return result

    @staticmethod
    def _filter_virtual_sites_near_adsorbate(
        ads_coords: np.ndarray,
        surface_atoms: "ase.Atoms",
        adsorbate_indices: List[int],
        exclusion_radius: float,
    ) -> np.ndarray:
        """Remove virtual sites closer than exclusion_radius to molecular adsorbate atoms."""
        return VSSRMCPredictor._filter_virtual_sites_near_atoms(
            ads_coords,
            reference_atoms=surface_atoms,
            reference_indices=adsorbate_indices,
            exclusion_radius=exclusion_radius,
        )

    @staticmethod
    def _filter_virtual_sites_near_atoms(
        ads_coords: np.ndarray,
        reference_atoms: "ase.Atoms",
        reference_indices: List[int],
        exclusion_radius: float,
    ) -> np.ndarray:
        """Remove virtual sites closer than exclusion_radius to selected existing atoms."""
        coords = np.asarray(ads_coords, dtype=float)
        if len(coords) == 0 or exclusion_radius <= 0 or not reference_indices:
            return coords

        valid_ref = [
            int(i) for i in reference_indices
            if 0 <= int(i) < len(reference_atoms)
        ]
        if not valid_ref:
            return coords

        from ase.geometry import get_distances

        ref_positions = reference_atoms.get_positions()[valid_ref]
        _, dists = get_distances(
            coords,
            ref_positions,
            cell=reference_atoms.cell,
            pbc=reference_atoms.pbc,
        )
        min_dists = dists.min(axis=1)
        return coords[min_dists >= exclusion_radius]

    @staticmethod
    def _filter_virtual_sites_by_site_spacing(
        ads_coords: np.ndarray,
        cell: Any,
        pbc: Any,
        min_site_distance: float,
    ) -> np.ndarray:
        """Greedily keep virtual sites separated by at least min_site_distance."""
        coords = np.asarray(ads_coords, dtype=float)
        if len(coords) == 0 or min_site_distance <= 0:
            return coords

        from ase import Atoms

        kept: list[np.ndarray] = []
        for coord in coords:
            if not kept:
                kept.append(coord)
                continue
            trial = Atoms(
                "H" * (len(kept) + 1),
                positions=np.vstack([kept, coord]),
                cell=cell,
                pbc=pbc,
            )
            site_idx = len(trial) - 1
            min_dist = min(
                float(trial.get_distance(site_idx, idx, mic=True))
                for idx in range(len(kept))
            )
            if min_dist >= min_site_distance:
                kept.append(coord)

        return np.asarray(kept, dtype=float)

    @staticmethod
    def _filter_virtual_sites_by_max_distance_to_atoms(
        ads_coords: np.ndarray,
        reference_atoms: "ase.Atoms",
        reference_indices: List[int],
        max_distance: float,
    ) -> np.ndarray:
        """Remove virtual sites farther than max_distance from selected atoms."""
        coords = np.asarray(ads_coords, dtype=float)
        if len(coords) == 0 or max_distance <= 0 or not reference_indices:
            return coords

        valid_ref = [
            int(i) for i in reference_indices
            if 0 <= int(i) < len(reference_atoms)
        ]
        if not valid_ref:
            return coords

        from ase.geometry import get_distances

        ref_positions = reference_atoms.get_positions()[valid_ref]
        _, dists = get_distances(
            coords,
            ref_positions,
            cell=reference_atoms.cell,
            pbc=reference_atoms.pbc,
        )
        min_dists = dists.min(axis=1)
        return coords[min_dists <= max_distance]

    @staticmethod
    def _expand_virtual_sites_locally(
        ads_coords: np.ndarray,
        reference_atoms: "ase.Atoms",
        reference_indices: List[int],
        radius: float,
        min_count: int,
    ) -> np.ndarray:
        """Add local tangential neighbors around real near-surface virtual sites."""
        coords = np.asarray(ads_coords, dtype=float)
        if len(coords) == 0 or radius <= 0 or min_count <= 0:
            return coords
        if len(coords) >= min_count:
            return coords

        valid_ref = [
            int(i) for i in reference_indices
            if 0 <= int(i) < len(reference_atoms)
        ]
        if not valid_ref:
            return coords

        ref_positions = reference_atoms.get_positions()[valid_ref]
        center = ref_positions.mean(axis=0)
        expanded: list[np.ndarray] = []
        for coord in coords:
            radial = np.asarray(coord, dtype=float) - center
            norm = float(np.linalg.norm(radial))
            normal = radial / norm if norm > 1e-8 else np.array([0.0, 0.0, 1.0])
            ref_axis = np.array([0.0, 0.0, 1.0])
            if abs(float(np.dot(normal, ref_axis))) > 0.9:
                ref_axis = np.array([1.0, 0.0, 0.0])
            tangent_u = np.cross(normal, ref_axis)
            tangent_u /= max(float(np.linalg.norm(tangent_u)), 1e-12)
            tangent_v = np.cross(normal, tangent_u)
            tangent_v /= max(float(np.linalg.norm(tangent_v)), 1e-12)
            # Build offsets fresh for each coordinate (a shared list would
            # leak tangents from the previous iteration).
            offsets = [
                np.array([0.0, 0.0, 0.0]),
                radius * tangent_u,
                -radius * tangent_u,
                radius * tangent_v,
                -radius * tangent_v,
                radius * (tangent_u + tangent_v) / np.sqrt(2.0),
                radius * (tangent_u - tangent_v) / np.sqrt(2.0),
            ]
            for offset in offsets:
                candidate = np.asarray(coord, dtype=float) + np.asarray(offset, dtype=float)
                if any(
                    float(np.linalg.norm(candidate - existing)) < 1e-8
                    for existing in expanded
                ):
                    continue
                trial = reference_atoms.copy()
                trial.append("H")
                site_idx = len(trial) - 1
                positions = trial.get_positions()
                positions[site_idx] = candidate
                trial.set_positions(positions)
                min_dist = min(
                    float(trial.get_distance(site_idx, ref_idx, mic=True))
                    for ref_idx in valid_ref
                )
                if min_dist <= VSSRMCPredictor._VIRTUAL_SITE_NEIGHBOR_CUTOFF:
                    expanded.append(candidate)

        return np.asarray(expanded, dtype=float) if expanded else coords

    def __del__(self):
        """Clean up temporary directory"""
        if hasattr(self, '_temp_dir') and self._temp_dir and not self.keep_files:
            try:
                shutil.rmtree(self._temp_dir)
            except:
                pass


def sample_surface_reconstruction(
    surface: Union[str, "ase.Atoms"],
    surface_sampling_root: str,
    model_type: str = "CHGNetNFF",
    adsorbates: Optional[List[str]] = None,
    total_sweeps: int = 100,
    temperature: float = 1.0,
    output_dir: Optional[str] = None,
    device: str = "cuda",
    **kwargs,
) -> VSSRMCResult:
    """
    Convenience function for VSSR-MC surface reconstruction sampling.

    Parameters
    ----------
    surface : str or ase.Atoms
        Input surface structure
    surface_sampling_root : str
        Path to surface-sampling repository
    model_type : str, default="CHGNetNFF"
        Model type to use
    adsorbates : list of str, optional
        List of adsorbate elements
    total_sweeps : int, default=100
        Number of MC sweeps
    temperature : float, default=1.0
        Temperature in kT
    output_dir : str, optional
        Output directory
    device : str, default="cuda"
        Device to use
    **kwargs
        Additional arguments

    Returns
    -------
    VSSRMCResult
        Sampling results
    """
    predictor = VSSRMCPredictor(
        surface_sampling_root=surface_sampling_root,
        model_type=model_type,
        device=device,
        keep_files=True if output_dir else False,
    )
    return predictor.sample(
        surface=surface,
        adsorbates=adsorbates,
        total_sweeps=total_sweeps,
        temperature=temperature,
        output_dir=output_dir,
        **kwargs,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="VSSR-MC Surface Reconstruction Sampler")
    parser.add_argument("surface", help="Path to surface structure file")
    parser.add_argument("--surface-sampling-root",
                       default="deps/surface-sampling",
                       help="Path to surface-sampling repository")
    parser.add_argument("--model-type", default="CHGNetNFF",
                       choices=["CHGNetNFF", "PaiNN", "NffScaleMACE", "UMA"],
                       help="NFF model type")
    parser.add_argument("--adsorbates", nargs="*", default=None,
                       help="Adsorbate elements (e.g., O H)")
    parser.add_argument("--total-sweeps", type=int, default=50,
                       help="Number of MC sweeps")
    parser.add_argument("--temperature", type=float, default=1.0,
                       help="Temperature in kT")
    parser.add_argument("--canonical", action="store_true",
                       help="Use canonical (NVT) sampling")
    parser.add_argument("--num-adsorbates", type=int, default=0,
                       help="Number of adsorbates (for canonical)")
    parser.add_argument("--cpu", action="store_true",
                       help="Use CPU instead of GPU")
    parser.add_argument("--output-dir", default=None,
                       help="Output directory")
    parser.add_argument("--no-auto-offset", action="store_true",
                       help="Disable automatic offset_data generation")

    args = parser.parse_args()

    predictor = VSSRMCPredictor(
        surface_sampling_root=args.surface_sampling_root,
        model_type=args.model_type,
        device="cpu" if args.cpu else "cuda",
        keep_files=True,
    )

    result = predictor.sample(
        surface=args.surface,
        adsorbates=args.adsorbates,
        canonical=args.canonical,
        num_adsorbates=args.num_adsorbates,
        total_sweeps=args.total_sweeps,
        temperature=args.temperature,
        output_dir=args.output_dir,
        auto_offset=not args.no_auto_offset,
    )

    print("\n" + "=" * 60)
    print("SAMPLING COMPLETE")
    print("=" * 60)
    print(result.summary())
