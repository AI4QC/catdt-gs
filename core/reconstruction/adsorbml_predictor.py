"""
AdsorbML Predictor — 基于 fairchem AdsorbML 算法的吸附位点预测封装

与 ``adsorbdiff_predictor.AdsorbDiffPredictor`` API 兼容，可作为 agent2
(``predict_adsorption_sites``) 的可选后端。

算法流程 (Lan et al., npj Comp Mater 2023)：
    1. 用 fairchem ``AdsorbateSlabConfig`` 按 heuristic/random 模式枚举候选位点
    2. 用 fairchem ``FAIRChemCalculator`` (UMA 等预训练 MLFF) 做全结构弛豫
    3. 取能量最低的配置作为全局最优吸附 (GMAE)

与 AdsorbDiff 的区别：
    * 不依赖 diffusion 模型，直接用 heuristic (ontop/bridge/hollow) + random 采样
    * 弛豫用 UMA (catdt 已在 NEB 中使用)，不需要额外 checkpoint
    * 结果 dataclass 与 AdsorbDiffPredictor 完全相同，可直接替换

使用方式:
    from core.reconstruction.adsorbml_predictor import AdsorbMLPredictor

    predictor = AdsorbMLPredictor(
        fairchem_model="uma-s-1p1",
        use_gpu=True,
        num_sites=20,
        placement_mode="random_site_heuristic_placement",
    )
    result = predictor.predict(
        surface="slab.vasp",
        adsorbate="*CO",
    )

注：dataclass 复用 adsorbdiff_predictor 中的 AdsorptionSite / AdsorptionResult /
PredictionOutput，以便 agent2 两个 backend 接口完全一致。
"""

import os
import shutil
import tempfile
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from core.reconstruction.adsorbdiff_predictor import (
    AdsorptionResult,
    AdsorptionSite,
    PredictionOutput,
)


def _cuda_available() -> bool:
    """Return True only when torch can actually use CUDA."""
    try:
        import torch

        return bool(torch.cuda.is_available())
    except Exception:
        return False


def _resolve_device(use_gpu: bool, requested_device: Optional[str] = None) -> str:
    """Resolve AdsorbML/Fairchem device with env override and CPU fallback."""
    env_device = (
        os.getenv("CATDT_ADSORBML_DEVICE", "") or os.getenv("CATDT_FAIRCHEM_DEVICE", "")
    ).strip().lower()
    device = env_device or (requested_device or ("cuda" if use_gpu else "cpu")).strip().lower()
    if device not in {"cuda", "cpu"}:
        warnings.warn(f"Unsupported AdsorbML device '{device}', falling back to cpu.")
        device = "cpu"
    if device == "cuda" and not _cuda_available():
        warnings.warn("CUDA is unavailable for AdsorbML; falling back to CPU.")
        device = "cpu"
    return device


class AdsorbMLPredictor:
    """
    AdsorbML 吸附位点预测器

    Parameters
    ----------
    fairchem_model : str, default="uma-s-1p1"
        fairchem 预训练 MLFF 模型名 (``uma-s-1p1``, ``uma-m-1p1`` 等)。
    fairchem_model_path : str, optional
        本地 checkpoint 路径；若提供则优先于 ``fairchem_model``。
    task_name : str, default="oc20"
        UMA 多任务头选择。气固界面用 ``oc20``；其它任务参考 fairchem 文档。
    use_gpu : bool, default=True
    num_sites : int, default=20
        采样位点数量。heuristic 模式下实际数量由 slab 对称性决定。
    num_augmentations_per_site : int, default=1
        每个位点随机旋转增强次数。
    placement_mode : str, default="random_site_heuristic_placement"
        fairchem ``AdsorbateSlabConfig.mode``：
        ``random`` / ``heuristic`` / ``random_site_heuristic_placement``
    interstitial_gap : float, default=0.1
        吸附物与 slab 的最小距离 (Å)。
    relax_steps : int, default=200
    relax_fmax : float, default=0.05
    seed : int, optional
    work_dir : str, optional
        工作目录；不指定则创建临时目录。
    keep_files : bool, default=False
    verbose : bool, default=True
    """

    def __init__(
        self,
        fairchem_model: str = "uma-s-1p1",
        fairchem_model_path: Optional[str] = None,
        task_name: str = "oc20",
        use_gpu: bool = True,
        device: Optional[str] = None,
        num_sites: int = 20,
        num_augmentations_per_site: int = 1,
        placement_mode: str = "random_site_heuristic_placement",
        interstitial_gap: float = 0.1,
        relax_steps: int = 200,
        relax_fmax: float = 0.05,
        seed: Optional[int] = 0,
        work_dir: Optional[str] = None,
        keep_files: bool = False,
        verbose: bool = True,
    ):
        self.fairchem_model = fairchem_model
        self.fairchem_model_path = (
            os.path.abspath(fairchem_model_path) if fairchem_model_path else None
        )
        self.task_name = task_name
        self.device = _resolve_device(use_gpu=use_gpu, requested_device=device)
        self.use_gpu = self.device == "cuda"
        self.num_sites = num_sites
        self.num_augmentations_per_site = num_augmentations_per_site
        self.placement_mode = placement_mode
        env_gap = os.getenv("CATDT_ADSORBML_INTERSTITIAL_GAP", "").strip()
        if env_gap:
            try:
                interstitial_gap = float(env_gap)
            except ValueError:
                warnings.warn(
                    f"Invalid CATDT_ADSORBML_INTERSTITIAL_GAP='{env_gap}', "
                    f"using constructor value {interstitial_gap}."
                )
        self.interstitial_gap = interstitial_gap
        self.relax_steps = relax_steps
        self.relax_fmax = relax_fmax
        self.seed = seed
        self.keep_files = keep_files
        self.verbose = verbose

        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="adsorbml_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        if seed is not None:
            np.random.seed(seed)

        self._calculator = None

        if self.verbose:
            print("AdsorbMLPredictor initialized:")
            print(f"  Model: {self.fairchem_model_path or self.fairchem_model}")
            print(f"  Task: {self.task_name}")
            print(f"  Device: {self.device}")
            print(f"  Placement mode: {self.placement_mode}")
            print(f"  interstitial_gap: {self.interstitial_gap}")
            print(f"  num_sites: {self.num_sites}")
            print(f"  Work dir: {self.work_dir}")

    # ------------------------------------------------------------------
    # Model loading
    # ------------------------------------------------------------------
    def _load_calculator(self):
        if self._calculator is not None:
            return self._calculator

        self._log("Loading fairchem MLFF calculator...")
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        if self.fairchem_model_path:
            predict_unit = pretrained_mlip.load_predict_unit(
                path=self.fairchem_model_path, device=self.device
            )
        else:
            predict_unit = pretrained_mlip.get_predict_unit(
                self.fairchem_model, device=self.device
            )

        self._calculator = FAIRChemCalculator(
            predict_unit=predict_unit,
            task_name=self.task_name,
            seed=self.seed,
        )
        self._log(f"  Loaded {self.fairchem_model} on {self.device}")
        return self._calculator

    def _log(self, message: str):
        if self.verbose:
            print(message)

    # ------------------------------------------------------------------
    # Input conversion
    # ------------------------------------------------------------------
    def _convert_surface(self, surface) -> "ase.Atoms":
        from ase import Atoms
        from ase.io import read

        if isinstance(surface, str):
            if not os.path.exists(surface):
                raise FileNotFoundError(f"Surface file not found: {surface}")
            return read(surface)
        if isinstance(surface, Atoms):
            return surface.copy()
        try:
            from pymatgen.core import Structure
            from pymatgen.io.ase import AseAtomsAdaptor

            if isinstance(surface, Structure):
                return AseAtomsAdaptor.get_atoms(surface)
        except ImportError:
            pass
        raise TypeError(
            f"Unsupported surface type: {type(surface)}. "
            "Expected str / ase.Atoms / pymatgen.Structure"
        )

    def _convert_adsorbate(self, adsorbate) -> Any:
        """
        Return a fairchem Adsorbate instance. Accepts:
          * fairchem Adsorbate (returned as-is — the typed path callers
            should prefer, so binding_indices are carried by the type)
          * SMILES string (via fairchem's built-in db — fuzzy match,
            unreliable for non-catalog labels)
          * file path to a structure
          * integer db id
          * ase.Atoms (binding defaults to atom 0)
        """
        from ase import Atoms
        from fairchem.data.oc.core import Adsorbate

        if isinstance(adsorbate, Adsorbate):
            return adsorbate
        if isinstance(adsorbate, str):
            if os.path.exists(adsorbate):
                from ase.io import read

                atoms = read(adsorbate)
                return Adsorbate(
                    adsorbate_atoms=atoms, adsorbate_binding_indices=[0]
                )
            return Adsorbate(adsorbate_smiles_from_db=adsorbate)
        if isinstance(adsorbate, int):
            return Adsorbate(adsorbate_id_from_db=adsorbate)
        if isinstance(adsorbate, Atoms):
            return Adsorbate(
                adsorbate_atoms=adsorbate.copy(), adsorbate_binding_indices=[0]
            )
        raise TypeError(
            f"Unsupported adsorbate type: {type(adsorbate)}. "
            "Expected fairchem Adsorbate, str (SMILES or path), int (db id), "
            "or ase.Atoms."
        )

    def _tag_surface_atoms(self, atoms: "ase.Atoms") -> "ase.Atoms":
        """Ensure slab has tags: 0=bulk, 1=surface."""
        tags = np.asarray(atoms.get_tags())
        if np.any(tags == 1):
            return atoms

        z = atoms.get_positions()[:, 2]
        z_max = z.max()
        new_tags = np.where(z > z_max - 2.5, 1, 0)
        atoms.set_tags(new_tags.tolist())
        return atoms

    def _build_slab(self, slab_atoms: "ase.Atoms"):
        """Bypass fairchem Slab validation to accept arbitrary ASE slabs."""
        from fairchem.data.oc.core import Slab

        slab = object.__new__(Slab)
        slab.atoms = slab_atoms
        slab.bulk = None
        slab.millers = (1, 1, 1)
        slab.shift = 0.0
        slab.top = True
        return slab

    def _extract_adsorbate_from_surface(
        self, atoms: "ase.Atoms", adsorbate_tag: int = 2
    ) -> Tuple["ase.Atoms", "ase.Atoms"]:
        tags = np.asarray(atoms.get_tags())
        ads_mask = tags == adsorbate_tag
        slab = atoms[~ads_mask]
        adsorbate = atoms[ads_mask]
        slab.set_tags([1 if t == 1 else 0 for t in slab.get_tags()])
        return slab, adsorbate

    # ------------------------------------------------------------------
    # Placement + relaxation
    # ------------------------------------------------------------------
    def _create_adslab_config(
        self,
        slab_atoms: "ase.Atoms",
        adsorbate: Any,
    ) -> Tuple[List["ase.Atoms"], List[Dict]]:
        from fairchem.data.oc.core import AdsorbateSlabConfig

        slab = self._build_slab(slab_atoms)
        cfg = AdsorbateSlabConfig(
            slab=slab,
            adsorbate=adsorbate,
            num_sites=self.num_sites,
            num_augmentations_per_site=self.num_augmentations_per_site,
            interstitial_gap=self.interstitial_gap,
            mode=self.placement_mode,
        )
        return cfg.atoms_list, cfg.metadata_list

    def _run_relaxation(
        self, adslab: "ase.Atoms", traj_path: Optional[str] = None
    ) -> Tuple["ase.Atoms", float, float]:
        from ase.constraints import FixAtoms
        from ase.optimize import BFGS

        calc = self._load_calculator()
        adslab = adslab.copy()

        # fairchem MLFFs require full 3D PBC; ASE slab builders leave z=False.
        adslab.set_pbc([True, True, True])

        # Freeze bulk atoms (tag==0) during relaxation, per OC20 convention.
        tags = np.asarray(adslab.get_tags())
        bulk_idx = np.where(tags == 0)[0].tolist()
        if bulk_idx and not any(
            isinstance(c, FixAtoms) for c in adslab.constraints
        ):
            adslab.set_constraint(FixAtoms(indices=bulk_idx))

        adslab.calc = calc
        opt = BFGS(adslab, trajectory=traj_path, logfile=None)
        try:
            opt.run(fmax=self.relax_fmax, steps=self.relax_steps)
        except Exception as e:  # noqa: BLE001
            self._log(f"    Warning: relaxation error: {e}")

        energy = adslab.get_potential_energy()
        forces = adslab.get_forces()
        fmax = float(np.max(np.linalg.norm(forces, axis=1)))
        return adslab, float(energy), fmax

    def _get_adsorbate_info(
        self, atoms: "ase.Atoms"
    ) -> Tuple[np.ndarray, np.ndarray]:
        tags = np.asarray(atoms.get_tags())
        ads_mask = tags == 2
        ads_positions = atoms.get_positions()[ads_mask]
        if len(ads_positions) == 0:
            return np.zeros(3), np.array([0.0, 0.0, 1.0])

        com = ads_positions.mean(axis=0)
        if len(ads_positions) > 1:
            direction = ads_positions[-1] - ads_positions[0]
            norm = np.linalg.norm(direction)
            orient = direction / norm if norm > 1e-6 else np.array([0.0, 0.0, 1.0])
        else:
            orient = np.array([0.0, 0.0, 1.0])
        return com, orient

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def predict(
        self,
        surface: Union[str, "ase.Atoms"],
        adsorbate: Union[str, int, "ase.Atoms", None] = None,
        has_adsorbate: bool = False,
        adsorbate_tag: int = 2,
        binding_indices: Optional[List[int]] = None,
        num_samples: Optional[int] = None,
        output_dir: Optional[str] = None,
        save_trajectory: bool = True,
    ) -> PredictionOutput:
        """
        预测最优吸附位点 (API 与 ``AdsorbDiffPredictor.predict`` 一致)。
        """
        if num_samples is not None:
            self.num_sites = num_samples

        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "prediction")
        os.makedirs(output_dir, exist_ok=True)

        self._log(f"\n{'=' * 60}")
        self._log("AdsorbML Prediction")
        self._log(f"{'=' * 60}")

        surface_atoms = self._convert_surface(surface)
        self._log(f"Surface: {surface_atoms.get_chemical_formula()}")

        if has_adsorbate:
            slab_atoms, adsorbate_atoms = self._extract_adsorbate_from_surface(
                surface_atoms, adsorbate_tag
            )
            self._log(
                f"  Extracted adsorbate: {adsorbate_atoms.get_chemical_formula()}"
            )
            adslab_list = [surface_atoms.copy()]
            metadata_list = [
                {
                    "site": self._get_adsorbate_info(surface_atoms)[0],
                    "xyz_angles": (0.0, 0.0, 0.0),
                }
            ]
            ads_formula = adsorbate_atoms.get_chemical_formula()
        else:
            if adsorbate is None:
                raise ValueError(
                    "`adsorbate` is required when has_adsorbate=False"
                )
            slab_atoms = self._tag_surface_atoms(surface_atoms)
            adsorbate_obj = self._convert_adsorbate(adsorbate)
            if binding_indices is not None:
                adsorbate_obj.binding_indices = binding_indices
            ads_formula = adsorbate_obj.atoms.get_chemical_formula()
            self._log(f"Adsorbate: {ads_formula}")

            self._log(
                f"\nGenerating up to {self.num_sites} placements "
                f"(mode={self.placement_mode})..."
            )
            adslab_list, metadata_list = self._create_adslab_config(
                slab_atoms, adsorbate_obj
            )
            self._log(f"  Created {len(adslab_list)} configurations")

        if not adslab_list:
            raise RuntimeError(
                "AdsorbML could not enumerate any adsorbate placements. "
                "Try increasing num_sites or switching placement_mode."
            )

        traj_dir = os.path.join(output_dir, "trajectories")
        os.makedirs(traj_dir, exist_ok=True)

        self._log("\nRelaxing candidates with fairchem MLFF...")
        results: List[AdsorptionResult] = []
        for i, (adslab, metadata) in enumerate(zip(adslab_list, metadata_list)):
            self._log(f"  Configuration {i + 1}/{len(adslab_list)}...")
            initial_atoms = adslab.copy()

            relax_traj = (
                os.path.join(traj_dir, f"relax_{i}.traj") if save_trajectory else None
            )
            try:
                relaxed, energy, fmax = self._run_relaxation(adslab, relax_traj)
            except Exception as e:  # noqa: BLE001
                self._log(f"    ✗ Relaxation failed: {e}")
                continue

            final_com, final_orient = self._get_adsorbate_info(relaxed)

            site_pos = np.asarray(metadata.get("site", np.zeros(3)))
            results.append(
                AdsorptionResult(
                    initial_structure=initial_atoms,
                    final_structure=relaxed,
                    adsorbate_position=final_com,
                    adsorbate_orientation=final_orient,
                    energy=energy,
                    forces_max=fmax,
                    trajectory_path=relax_traj,
                    site_info=AdsorptionSite(
                        site_id=i,
                        position=site_pos,
                        binding_type=self.placement_mode,
                    ),
                )
            )

        if not results:
            raise RuntimeError("AdsorbML: all candidate relaxations failed.")

        # ── Chemisorption- and integrity-aware selection ──────────────
        # UMA's PES can have "physisorbed / vacuum" basins slightly lower
        # than the "chemisorbed" basin for weakly-binding adsorbates (e.g.
        # CO2 on Ga-rich Ni5Ga3(211) where bent CO2^δ- is ~0.1 eV ABOVE
        # linear floating CO2). For building a reaction pathway we need
        # the chemically-bound starting state, not the thermodynamically
        # lowest one. For saturated molecules such as C3H8 on TiOx we must
        # also reject candidates where adsorption relaxation already broke
        # internal adsorbate bonds, otherwise Agent2 silently starts from a
        # hydroxylated/dehydrogenated state instead of molecular propane.
        bound_results = [
            r for r in results
            if self._is_chemisorbed(r.final_structure)
        ]
        intact_by_id = {
            id(r): self._adsorbate_molecule_intact(
                r.initial_structure, r.final_structure
            )
            for r in results
        }
        intact_results = [r for r in results if intact_by_id[id(r)]]
        bound_intact_results = [
            r for r in bound_results
            if intact_by_id[id(r)]
        ]
        require_chemisorbed = os.getenv(
            "CATDT_REQUIRE_CHEMISORBED_INITIAL", "0"
        ).strip().lower() in {"1", "true", "yes", "on"}
        require_intact = os.getenv(
            "CATDT_REQUIRE_INTACT_ADSORBATE_INITIAL", "0"
        ).strip().lower() in {"1", "true", "yes", "on"}

        if bound_intact_results:
            best_result = min(bound_intact_results, key=lambda r: r.energy)
            selection_note = (
                f"Selected CHEMISORBED + INTACT candidate "
                f"(best of {len(bound_intact_results)} bound/intact, "
                f"{len(bound_results)} bound, {len(intact_results)} intact "
                f"out of {len(results)} total): E={best_result.energy:.3f} eV"
            )
            chemisorbed_selected = True
            intact_selected = True
        elif intact_results:
            if require_chemisorbed:
                raise RuntimeError(
                    f"No chemisorbed candidate found among molecularly intact "
                    f"candidates in "
                    f"{len(results)} relaxed configurations for adsorbate "
                    f"{ads_formula}. Refusing to continue from an intact but "
                    "physisorbed fallback because "
                    "CATDT_REQUIRE_CHEMISORBED_INITIAL=1. Inspect AdsorbML "
                    f"outputs in {output_dir} or adjust placement controls "
                    "(for PDH/C3H8 on TiOx, try "
                    "CATDT_ADSORBML_INTERSTITIAL_GAP=2.0)."
                )
            best_result = min(intact_results, key=lambda r: r.energy)
            chemisorbed_selected = self._is_chemisorbed(best_result.final_structure)
            selection_note = (
                f"⚠ No chemisorbed molecularly intact candidate found "
                f"({len(bound_results)} bound, {len(intact_results)} intact "
                f"out of {len(results)} total). Selected lowest-energy INTACT "
                f"candidate instead: E={best_result.energy:.3f} eV"
            )
            intact_selected = True
        else:
            if require_intact:
                raise RuntimeError(
                    f"No molecularly intact candidate found in {len(results)} "
                    f"relaxed configurations for adsorbate {ads_formula}. "
                    "Refusing to continue because "
                    "CATDT_REQUIRE_INTACT_ADSORBATE_INITIAL=1. Inspect "
                    f"AdsorbML outputs in {output_dir} and adjust placement "
                    "controls."
                )
            fallback_pool = bound_results or results
            if require_chemisorbed:
                raise RuntimeError(
                    f"No chemisorbed candidate found in {len(results)} relaxed "
                    f"configurations for adsorbate {ads_formula}. Refusing to "
                    "continue from a physisorbed fallback because "
                    "CATDT_REQUIRE_CHEMISORBED_INITIAL=1. Inspect AdsorbML "
                    f"outputs in {output_dir} and provide a chemically bound "
                    "initial structure or adjusted placement."
                )
            best_result = min(fallback_pool, key=lambda r: r.energy)
            selection_note = (
                f"⚠ No molecularly intact candidate found in {len(results)} "
                f"relaxed configurations. Falling back to lowest-energy "
                f"{'chemisorbed' if bound_results else 'overall'} candidate, "
                f"which may be dissociated: E={best_result.energy:.3f} eV"
            )
            chemisorbed_selected = bool(bound_results)
            intact_selected = False

        self._log(f"\n{selection_note}")

        output = PredictionOutput(
            surface_formula=slab_atoms.get_chemical_formula(),
            adsorbate_formula=ads_formula,
            num_sites_sampled=len(results),
            results=results,
            best_result=best_result,
            output_dir=output_dir,
            use_proxy_energy=False,
        )
        # Stash selection metadata on the output for downstream logging
        output.chemisorbed_selected = chemisorbed_selected
        output.n_bound_candidates = len(bound_results)
        output.intact_adsorbate_selected = intact_selected
        output.n_intact_candidates = len(intact_results)
        output.n_bound_intact_candidates = len(bound_intact_results)
        output.selection_note = selection_note

        self._log(f"\n{output.summary()}")

        from ase.io import write

        best_path = os.path.join(output_dir, "best_config.vasp")
        write(best_path, best_result.final_structure, format="vasp")
        self._log(f"\nBest configuration saved to: {best_path}")
        return output

    @staticmethod
    def _adsorbate_internal_bonds(
        atoms: "ase.Atoms",
        cov_scale: float = 1.25,
    ) -> set[tuple[int, int]]:
        """Return covalent bonds within tag=2 adsorbate atoms."""
        from ase.data import atomic_numbers, covalent_radii

        tags = np.asarray(atoms.get_tags())
        ads_idx = np.where(tags == 2)[0].tolist()
        if len(ads_idx) < 2:
            return set()

        symbols = atoms.get_chemical_symbols()
        bonds: set[tuple[int, int]] = set()
        for pos_i, i in enumerate(ads_idx):
            ri = covalent_radii[atomic_numbers[symbols[i]]]
            for j in ads_idx[pos_i + 1:]:
                rj = covalent_radii[atomic_numbers[symbols[j]]]
                cutoff = cov_scale * (ri + rj)
                d = float(atoms.get_distance(i, j, mic=True))
                if d < cutoff:
                    bonds.add((i, j))
        return bonds

    @classmethod
    def _adsorbate_molecule_intact(
        cls,
        initial_atoms: "ase.Atoms",
        final_atoms: "ase.Atoms",
    ) -> bool:
        """True when all initial internal adsorbate bonds survive relaxation."""
        initial_tags = np.asarray(initial_atoms.get_tags())
        final_tags = np.asarray(final_atoms.get_tags())
        if len(initial_atoms) != len(final_atoms):
            return False
        if np.count_nonzero(initial_tags == 2) != np.count_nonzero(final_tags == 2):
            return False

        initial_bonds = cls._adsorbate_internal_bonds(initial_atoms)
        if not initial_bonds:
            return True

        final_bonds = cls._adsorbate_internal_bonds(final_atoms)
        return initial_bonds.issubset(final_bonds)

    @staticmethod
    def _is_chemisorbed(atoms: "ase.Atoms", cov_scale: float = 1.3) -> bool:
        """True iff at least one adsorbate atom (tag=2) is within covalent
        bond range of any non-adsorbate atom (surface)."""
        import numpy as np
        from ase.data import covalent_radii, atomic_numbers

        tags = atoms.get_tags()
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()

        ads_mask = tags == 2
        surf_mask = ~ads_mask
        ads_idx = np.where(ads_mask)[0]
        surf_idx = np.where(surf_mask)[0]
        if len(ads_idx) == 0 or len(surf_idx) == 0:
            return False

        for a in ads_idx:
            ra = covalent_radii[atomic_numbers[symbols[a]]]
            for s in surf_idx:
                rs = covalent_radii[atomic_numbers[symbols[s]]]
                cutoff = cov_scale * (ra + rs)
                d = float(np.linalg.norm(positions[a] - positions[s]))
                if d < cutoff:
                    return True
        return False

    def predict_single(
        self,
        adslab: "ase.Atoms",
        output_dir: Optional[str] = None,
    ) -> AdsorptionResult:
        """对单个已放置好的 adslab 进行 MLFF 弛豫。"""
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "single_prediction")
        os.makedirs(output_dir, exist_ok=True)

        initial_atoms = adslab.copy()
        relax_traj = os.path.join(output_dir, "relax.traj")
        relaxed, energy, fmax = self._run_relaxation(adslab, relax_traj)
        final_com, final_orient = self._get_adsorbate_info(relaxed)

        return AdsorptionResult(
            initial_structure=initial_atoms,
            final_structure=relaxed,
            adsorbate_position=final_com,
            adsorbate_orientation=final_orient,
            energy=energy,
            forces_max=fmax,
            trajectory_path=relax_traj,
        )

    def __del__(self):
        if getattr(self, "_temp_dir", None) and not self.keep_files:
            try:
                shutil.rmtree(self._temp_dir)
            except Exception:
                pass


# --------------------------------------------------------------------------
# Convenience wrapper
# --------------------------------------------------------------------------
def predict_adsorption_site_adsorbml(
    surface: Union[str, "ase.Atoms"],
    adsorbate: Union[str, int, "ase.Atoms"],
    fairchem_model: str = "uma-s-1p1",
    use_gpu: bool = True,
    num_samples: int = 20,
    output_dir: Optional[str] = None,
    **kwargs,
) -> PredictionOutput:
    predictor = AdsorbMLPredictor(
        fairchem_model=fairchem_model,
        use_gpu=use_gpu,
        num_sites=num_samples,
        keep_files=bool(output_dir),
        **kwargs,
    )
    return predictor.predict(
        surface=surface,
        adsorbate=adsorbate,
        has_adsorbate=False,
        output_dir=output_dir,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="AdsorbML Adsorption Site Predictor")
    parser.add_argument("surface", help="Path to surface structure file")
    parser.add_argument("--adsorbate", default="*CO", help="Adsorbate SMILES")
    parser.add_argument("--model", default="uma-s-1p1", help="fairchem model name")
    parser.add_argument("--num-samples", type=int, default=20)
    parser.add_argument(
        "--mode",
        default="random_site_heuristic_placement",
        choices=[
            "random",
            "heuristic",
            "random_site_heuristic_placement",
        ],
    )
    parser.add_argument("--cpu", action="store_true")
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    predictor = AdsorbMLPredictor(
        fairchem_model=args.model,
        use_gpu=not args.cpu,
        num_sites=args.num_samples,
        placement_mode=args.mode,
        keep_files=True,
    )
    result = predictor.predict(
        surface=args.surface,
        adsorbate=args.adsorbate,
        output_dir=args.output_dir,
    )
    print("\n" + "=" * 60)
    print("PREDICTION COMPLETE")
    print("=" * 60)
    print(result.summary())
