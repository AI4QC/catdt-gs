"""
Fairchem Predictor - 催化剂自由能预测封装类

这个模块提供了一个完整的封装类，用于使用 Fairchem 模型预测催化剂表面的吸附能和自由能。

Fairchem 提供了多个预训练模型（如 UMA、eSCN、GemNet 等）用于催化剂系统的能量预测。

使用方式:
    from fairchem_predictor import FairchemPredictor

    # 初始化预测器
    predictor = FairchemPredictor(
        fairchem_root="/path/to/fairchem",
        model_name="uma-s-1p1",  # 或 "uma-m-1p1"
        use_gpu=True,
    )

    # 方式1: 预测单个结构的能量
    result = predictor.predict_energy(
        structure="surface_with_adsorbate.vasp",
        relax=True,  # 是否进行结构优化
    )
    print(f"Energy: {result.energy:.3f} eV")

    # 方式2: 预测吸附能
    ads_result = predictor.predict_adsorption_energy(
        surface="clean_surface.vasp",
        adsorbate="*CO",  # 或提供分子结构
        num_sites=10,  # 尝试多少个吸附位点
    )
    print(f"Adsorption energy: {ads_result.adsorption_energy:.3f} eV")

    # 方式3: 批量预测多个吸附质的能量
    pathway_result = predictor.predict_pathway_energies(
        surface="surface.vasp",
        adsorbates=["*O", "*OH", "*OOH"],
        num_sites=5,
    )
    for ads, energy in pathway_result.adsorbate_energies.items():
        print(f"{ads}: {energy:.3f} eV")

作者: Claude
"""

import os
import sys
import tempfile
import logging
import importlib
import traceback
from typing import Union, Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
import warnings

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import LBFGS, BFGS


@dataclass
class EnergyResult:
    """单个结构的能量预测结果"""
    structure_name: str
    initial_structure: Atoms
    final_structure: Atoms
    energy: float  # eV
    forces: np.ndarray  # eV/Å
    relaxed: bool
    converged: bool
    optimization_steps: int
    output_file: Optional[str] = None

    def __repr__(self):
        return (f"EnergyResult(structure={self.structure_name}, "
                f"energy={self.energy:.4f} eV, relaxed={self.relaxed})")


@dataclass
class AdsorptionResult:
    """吸附能预测结果"""
    adsorbate: str
    surface_energy: float  # eV
    adsorbate_surface_energy: float  # eV
    adsorbate_reference_energy: float  # eV — atomic-reference based formation energy
    adsorption_energy: float  # eV — E_ads in atomic-reference frame
    num_sites_tried: int
    best_site_index: int
    best_configuration: Atoms
    all_configurations: List[Atoms] = field(default_factory=list)
    all_energies: List[float] = field(default_factory=list)
    is_anomalous: bool = False

    # Physical adsorption energy referenced to the adsorbate relaxed in a
    # vacuum box (same level of theory). This is what microkinetic models
    # expect when gas species are pinned at formation_energy = 0. If vacuum
    # relaxation was not run, `gas_phase_energy` is None and physical_E_ads
    # falls back to the atomic-reference adsorption_energy.
    gas_phase_energy: Optional[float] = None  # eV — E(adsorbate) relaxed in vacuum
    physical_adsorption_energy: Optional[float] = None  # eV — E(A@S) - E(S) - E(A_gas_DFT)
    gas_phase_frequencies_ev: Optional[List[float]] = None  # positive vib energies, eV
    gas_phase_frequencies_cm: Optional[List[float]] = None  # positive vib frequencies, cm^-1

    def __repr__(self):
        eads_str = (f"E_ads_phys={self.physical_adsorption_energy:.4f}"
                    if self.physical_adsorption_energy is not None
                    else f"E_ads={self.adsorption_energy:.4f}")
        return (f"AdsorptionResult(adsorbate={self.adsorbate}, "
                f"{eads_str} eV, sites={self.num_sites_tried})")


@dataclass
class PathwayResult:
    """反应路径能量预测结果"""
    surface_formula: str
    surface_energy: float  # eV
    adsorbates: List[str]
    adsorbate_energies: Dict[str, float]  # adsorbate -> energy (eV)
    adsorption_energies: Dict[str, float]  # adsorbate -> E_ads (eV)
    best_configurations: Dict[str, Atoms]
    relative_energies: Dict[str, float]  # relative to most stable
    output_dir: str
    # Vacuum-phase vibrational frequencies in cm^-1, keyed by the same labels
    # as adsorbate_energies. Populated when predict_pathway_energies completes
    # a gas-phase relaxation + vibration calc for each adsorbate. Required for
    # CatMAP's ideal_gas thermo mode.
    gas_frequencies: Dict[str, List[float]] = field(default_factory=dict)

    def summary(self) -> str:
        """生成结果摘要"""
        lines = [
            "Fairchem Pathway Energy Prediction Results",
            "=" * 60,
            f"Surface: {self.surface_formula}",
            f"Surface energy: {self.surface_energy:.4f} eV",
            f"Number of adsorbates: {len(self.adsorbates)}",
            "",
            "Adsorption energies:",
        ]

        # 按吸附能排序
        sorted_ads = sorted(self.adsorption_energies.items(),
                          key=lambda x: x[1])

        for ads, e_ads in sorted_ads:
            total_e = self.adsorbate_energies[ads]
            rel_e = self.relative_energies[ads]
            lines.append(f"  {ads:10s}: E_ads = {e_ads:7.3f} eV, "
                        f"E_total = {total_e:10.3f} eV, "
                        f"E_rel = {rel_e:6.3f} eV")

        lines.append("")
        lines.append(f"Output directory: {self.output_dir}")

        return "\n".join(lines)


class FairchemPredictor:
    """
    Fairchem 自由能预测器

    这个类封装了 Fairchem 项目的功能，提供简洁的 API 用于预测催化剂表面的
    吸附能和自由能。

    Parameters
    ----------
    fairchem_root : str
        Fairchem 项目的根目录路径
    model_name : str, default="uma-s-1p1"
        使用的模型名称。可选:
        - "uma-s-1p1": UMA small model (快速)
        - "uma-m-1p1": UMA medium model (更准确)
        - 或提供本地模型文件的完整路径
    model_path : str, optional
        本地模型文件路径。如果提供，将使用本地模型而不是从 Hugging Face 下载
    task_name : str, default="oc20"
        任务类型。催化应用使用 "oc20"
    use_gpu : bool, default=True
        是否使用 GPU
    device : str, optional
        指定设备 ("cuda" 或 "cpu")
    inference_mode : str, default="default"
        推理模式: "default" 或 "turbo" (更快)
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """

    # 原子参考能量 (eV) - 用于计算吸附能
    ATOMIC_REFERENCE_ENERGIES = {
        "H": -3.477,
        "C": -7.282,
        "N": -8.083,
        "O": -7.204,
        "F": -4.891,
        "S": -4.659,
    }

    def __init__(
        self,
        fairchem_root: str,
        model_name: str = "uma-s-1p1",
        model_path: Optional[str] = None,
        task_name: str = "oc20",
        use_gpu: bool = True,
        device: Optional[str] = None,
        inference_mode: str = "default",
        work_dir: Optional[str] = None,
        keep_files: bool = True,
        verbose: bool = True,
    ):
        self.fairchem_root = os.path.abspath(fairchem_root)
        self.model_name = model_name
        self.model_path = os.path.abspath(model_path) if model_path else None
        self.task_name = task_name
        env_device = str(os.getenv("CATDT_FAIRCHEM_DEVICE", "")).strip().lower()
        requested_device = str(device or ("cuda" if use_gpu else "cpu")).strip().lower()
        resolved_device = env_device or requested_device or "cuda"
        if resolved_device not in {"cuda", "cpu"}:
            warnings.warn(
                f"Unsupported Fairchem device '{resolved_device}', falling back to "
                f"{requested_device or 'cuda'}'."
            )
            resolved_device = requested_device or "cuda"
        self.device = resolved_device
        self.use_gpu = self.device == "cuda"
        self.inference_mode = inference_mode
        self.keep_files = keep_files
        self.verbose = verbose

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="fairchem_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 添加 Fairchem 路径到系统路径
        fairchem_src = os.path.join(self.fairchem_root, "src")
        if fairchem_src not in sys.path:
            sys.path.insert(0, fairchem_src)

        # 初始化日志
        self.logger = logging.getLogger(__name__)
        if self.verbose:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.WARNING)

        # 延迟加载模型（在第一次使用时）
        self._predictor = None
        self._calculator = None

    def _load_model(self):
        """延迟加载模型"""
        if self._predictor is not None:
            return

        # PyTorch 2.8+ defaults weights_only=True; UMA checkpoints need `slice` whitelisted
        import torch
        if hasattr(torch.serialization, "add_safe_globals"):
            torch.serialization.add_safe_globals([slice])

        # 🔧 确保使用 deps/fairchem 的版本，而不是环境中已安装的旧版本
        fairchem_src_path = os.path.join(self.fairchem_root, "src")
        fairchem_pkg_path = os.path.join(fairchem_src_path, "fairchem")
        if not os.path.isdir(fairchem_pkg_path):
            raise FileNotFoundError(
                f"fairchem source directory not found: {fairchem_pkg_path}"
            )

        os.environ.setdefault("CATDT_DISABLE_WANDB", "1")

        if fairchem_src_path in sys.path:
            sys.path.remove(fairchem_src_path)
        sys.path.insert(0, fairchem_src_path)

        # 避免在长时间多 episode 运行中反复清空 fairchem 相关模块。
        # 这种“硬重置”会触发不可预测的导入态损坏（出现随机 sys/sys、tuple callable 等异常）。
        # 仅在显式开启时执行。
        hard_reset = str(os.getenv("CATDT_FAIRCHEM_HARD_RESET", "0")).strip().lower() in {
            "1",
            "true",
            "yes",
        }
        def _reset_fairchem_modules(full: bool = False) -> None:
            local_src_abs = os.path.abspath(fairchem_src_path)
            for mod_name, mod in list(sys.modules.items()):
                if not mod_name.startswith("fairchem"):
                    continue
                if full:
                    sys.modules.pop(mod_name, None)
                    continue
                mod_file = getattr(mod, "__file__", "") or ""
                mod_file_abs = os.path.abspath(mod_file) if mod_file else ""
                if mod_file_abs.startswith(local_src_abs):
                    continue
                sys.modules.pop(mod_name, None)
            importlib.invalidate_caches()

        def _import_fairchem_components() -> Tuple[Any, Any, Any]:
            try:
                try:
                    from fairchem.core.calculate import pretrained_mlip, FAIRChemCalculator
                except ImportError:
                    from fairchem.core import pretrained_mlip, FAIRChemCalculator
                from fairchem.core.units.mlip_unit import load_predict_unit
                return pretrained_mlip, FAIRChemCalculator, load_predict_unit
            except ImportError as e:
                raise ImportError(
                    f"Failed to import fairchem modules from {fairchem_src_path}. "
                    f"Make sure fairchem is properly installed. Error: {e}"
                ) from e

        if hard_reset:
            _reset_fairchem_modules(full=False)

        pretrained_mlip, FAIRChemCalculator, load_predict_unit = _import_fairchem_components()

        # 检查是否使用本地模型
        if self.model_path:
            if not os.path.exists(self.model_path):
                raise FileNotFoundError(f"Model file not found: {self.model_path}")

            self.logger.info(f"Loading model from local path: {self.model_path}")

            # 使用本地模型文件
            if self.inference_mode == "turbo":
                self._predictor = load_predict_unit(
                    path=self.model_path,
                    device=self.device,
                    inference_settings="turbo"
                )
            else:
                self._predictor = load_predict_unit(
                    path=self.model_path,
                    device=self.device
                )
        else:
            # 优先使用本地模型，避免网络下载阻塞；失败后再尝试 HuggingFace。
            local_checkpoints = [
                f"deps/fairchem_models/{self.model_name}.pt",  # 本地UMA模型
                "deps/fairchem_models/uma-s-1p1.pt",
                "deps/fairchem_models/uma-m-1p1.pt",
                "data/ocp/checkpoints/2024-04-12-11-03-28/best_checkpoint.pt",
                "data/ocp/checkpoints/2024-03-03-23-57-52/best_checkpoint.pt",
            ]
            existing_checkpoints = []
            for checkpoint in local_checkpoints:
                if os.path.exists(checkpoint):
                    existing_checkpoints.append(checkpoint)

            fallback_model = existing_checkpoints[0] if existing_checkpoints else None
            local_errors: List[str] = []
            if existing_checkpoints:
                max_local_attempts = max(
                    1,
                    int(os.getenv("CATDT_LOCAL_MODEL_LOAD_ATTEMPTS", "3")),
                )
                for checkpoint in existing_checkpoints:
                    self.logger.info(f"Attempting to use local checkpoint first: {checkpoint}")
                    loaded = False
                    last_err = None
                    for attempt in range(1, max_local_attempts + 1):
                        try:
                            self._predictor = load_predict_unit(
                                path=checkpoint,
                                device=self.device,
                                overrides={"backbone": {"always_use_pbc": False}},
                            )
                            self.logger.info(
                                "Successfully loaded local checkpoint %s (attempt %d/%d)",
                                checkpoint,
                                attempt,
                                max_local_attempts,
                            )
                            loaded = True
                            break
                        except Exception as err:
                            last_err = err
                            self.logger.warning(
                                "Failed to load local checkpoint %s (attempt %d/%d): %s",
                                checkpoint,
                                attempt,
                                max_local_attempts,
                                err,
                            )
                            self.logger.debug(
                                "Local checkpoint traceback:\n%s",
                                traceback.format_exc(),
                            )
                            if attempt < max_local_attempts:
                                self.logger.warning(
                                    "Retrying local checkpoint load after full fairchem module reset..."
                                )
                                try:
                                    _reset_fairchem_modules(full=True)
                                    pretrained_mlip, FAIRChemCalculator, load_predict_unit = (
                                        _import_fairchem_components()
                                    )
                                except Exception as reset_err:
                                    self.logger.warning(
                                        "Fairchem reset/import failed before retry: %s",
                                        reset_err,
                                    )
                                    break
                    if loaded:
                        break
                    local_errors.append(f"{checkpoint}: {last_err}")
                    self._predictor = None

            local_error = "; ".join(local_errors) if local_errors else None

            if self._predictor is None:
                self.logger.info(f"Attempting to load model from HuggingFace: {self.model_name}")
                try:
                    if self.inference_mode == "turbo":
                        self._predictor = pretrained_mlip.get_predict_unit(
                            self.model_name,
                            device=self.device,
                            inference_settings="turbo"
                        )
                    else:
                        self._predictor = pretrained_mlip.get_predict_unit(
                            self.model_name,
                            device=self.device
                        )
                    self.logger.info(f"Successfully loaded {self.model_name} from HuggingFace")
                except Exception as hf_error:
                    if fallback_model is None:
                        raise RuntimeError(
                            f"Failed to download {self.model_name} from HuggingFace and no local fallback model found.\n"
                            f"HuggingFace error: {hf_error}\n"
                            f"Please either:\n"
                            f"  (1) Log in to HuggingFace: huggingface-cli login\n"
                            f"  (2) Provide a local compatible model via model_path parameter\n"
                            f"  (3) Set calculate_barriers=False to skip NEB calculations"
                        )

                    raise RuntimeError(
                        f"Failed to load both local checkpoint and HuggingFace model.\n"
                        f"Local checkpoint: {fallback_model}\n"
                        f"Local error: {local_error}\n"
                        f"HuggingFace error: {hf_error}\n\n"
                        f"Please either:\n"
                        f"  (1) Provide a compatible local model via model_path parameter\n"
                        f"  (2) Log in to HuggingFace: huggingface-cli login\n"
                        f"  (3) Set calculate_barriers=False to skip barrier calculations"
                    )

        # 创建计算器
        self._calculator = FAIRChemCalculator(
            self._predictor,
            task_name=self.task_name
        )

        self.logger.info("Model loaded successfully")

    def _relax_adsorbate_in_vacuum(
        self,
        adsorbate_atoms: Atoms,
        fmax: float = 0.05,
        max_steps: int = 200,
        vacuum_box: float = 20.0,
        compute_vibrations: bool = False,
    ) -> Union[float, Tuple[float, Optional[List[float]], Optional[Atoms]]]:
        """Relax `adsorbate_atoms` in a cubic vacuum box.

        Returns the relaxed energy (eV). If ``compute_vibrations`` is True,
        returns ``(energy, freq_cm1_list, relaxed_atoms)`` where ``freq_cm1_list``
        is the list of real vibrational frequencies in cm^-1 (imaginary modes
        are filtered out, translational/rotational zero modes removed).
        """
        self._load_model()
        if not hasattr(self, "_vacuum_energy_cache"):
            self._vacuum_energy_cache: Dict[str, Tuple[float, Optional[List[float]], Optional[Atoms]]] = {}

        formula_key = adsorbate_atoms.get_chemical_formula() or "_empty"
        cached = self._vacuum_energy_cache.get(formula_key)
        if cached is not None:
            gas_energy, freq_list, relaxed_atoms = cached
            # If caller asked for vibrations and we didn't cache them, fall
            # through to recompute; otherwise reuse.
            if (not compute_vibrations) or (freq_list is not None):
                self.logger.info(
                    "  Gas-phase energy reused from cache: %s → %.4f eV%s",
                    formula_key, gas_energy,
                    f" ({len(freq_list)} freqs)" if freq_list else "",
                )
                if compute_vibrations:
                    return gas_energy, freq_list, relaxed_atoms
                return gas_energy

        # Build an isolated gas-phase structure: center in vacuum cubic cell.
        gas_atoms = adsorbate_atoms.copy()
        gas_atoms.set_pbc(False)
        gas_atoms.set_cell([vacuum_box, vacuum_box, vacuum_box])
        gas_atoms.center()
        try:
            gas_atoms.set_constraint()
        except Exception:
            pass

        try:
            result = self.predict_energy(
                gas_atoms, relax=True, fmax=fmax, max_steps=max_steps,
            )
            gas_energy = float(result.energy)
            relaxed_atoms = result.final_structure.copy()
            self.logger.info(
                "  Relaxed %s in vacuum (%d atoms): E = %.4f eV",
                formula_key, len(gas_atoms), gas_energy,
            )
        except Exception as exc:
            self.logger.warning(
                "Vacuum relaxation failed for %s (%s); using single-point fallback",
                formula_key, exc,
            )
            try:
                sp_atoms = gas_atoms.copy()
                sp_atoms.calc = self._calculator
                gas_energy = float(sp_atoms.get_potential_energy())
                relaxed_atoms = sp_atoms.copy()
            except Exception as exc2:
                self.logger.error(
                    "Single-point vacuum calc also failed for %s: %s",
                    formula_key, exc2,
                )
                self._vacuum_energy_cache[formula_key] = (float("nan"), None, None)
                if compute_vibrations:
                    return float("nan"), None, None
                return float("nan")

        freq_list_cm: Optional[List[float]] = None
        if compute_vibrations and len(relaxed_atoms) >= 2:
            try:
                freq_list_cm = self._vibrations_in_vacuum(relaxed_atoms)
                self.logger.info(
                    "  Gas-phase vibrations %s: %d real modes (%.0f..%.0f cm^-1)",
                    formula_key, len(freq_list_cm),
                    min(freq_list_cm) if freq_list_cm else 0.0,
                    max(freq_list_cm) if freq_list_cm else 0.0,
                )
            except Exception as exc:
                self.logger.warning(
                    "Vibrational frequency calc failed for %s: %s — "
                    "falling back to no-frequency (ideal_gas will fail).",
                    formula_key, exc,
                )
                freq_list_cm = None

        self._vacuum_energy_cache[formula_key] = (gas_energy, freq_list_cm, relaxed_atoms)
        if compute_vibrations:
            return gas_energy, freq_list_cm, relaxed_atoms
        return gas_energy

    def _vibrations_in_vacuum(
        self,
        atoms: Atoms,
        delta: float = 0.02,
    ) -> List[float]:
        """Run ASE Vibrations on an isolated gas-phase molecule.

        Returns a sorted list of real (non-imaginary) vibrational frequencies
        in cm^-1. Rigid-body translational/rotational modes are already near
        zero and are filtered by the ``freq > 20 cm^-1`` cutoff.
        """
        import tempfile
        from ase.vibrations import Vibrations
        self._load_model()
        work = atoms.copy()
        work.set_pbc(False)
        work.set_constraint()
        work.calc = self._calculator
        with tempfile.TemporaryDirectory() as tmpdir:
            name = os.path.join(tmpdir, "vib")
            vib = Vibrations(work, name=name, delta=delta)
            vib.run()
            energies = vib.get_energies()  # eV, complex
            # Convert to cm^-1, keep only real positive above threshold.
            out: List[float] = []
            for e in energies:
                val = float(e.real) if hasattr(e, "real") else float(e)
                if abs(val) < 1e-6:
                    continue
                # Negative or imaginary values indicate saddle modes; drop.
                if val <= 0:
                    continue
                # Convert eV → cm^-1 (1 eV = 8065.54 cm^-1)
                freq_cm = val * 8065.54429
                if freq_cm < 20.0:  # drop rigid-body near-zero modes
                    continue
                out.append(freq_cm)
        return sorted(out)

    def _prepare_slab_for_adsorbate_placement(self, atoms: Atoms) -> Atoms:
        """
        准备表面用于吸附质放置

        确保表面满足 Fairchem Slab 的要求：
        - Tagged（表面原子 tag=1，体相原子 tag=0）
        - Tiled（xy 方向 >= 8 埃）
        - 添加了 FixAtoms constraints

        Parameters
        ----------
        atoms : Atoms
            原始表面结构

        Returns
        -------
        Atoms
            准备好的表面结构
        """
        from ase.constraints import FixAtoms

        prepared_atoms = atoms.copy()

        # 1. Tag 原子
        # 简单策略：基于 z 坐标区分表面和体相
        # 最上面 2 层为表面原子（tag=1），其余为体相原子（tag=0）
        z_coords = prepared_atoms.positions[:, 2]
        z_max = z_coords.max()

        # 识别层（通过 z 坐标分组）
        z_unique = np.unique(np.round(z_coords, decimals=2))
        z_sorted = np.sort(z_unique)[::-1]  # 从上到下排序

        # 标记前 2-3 层为表面原子
        n_surface_layers = min(3, len(z_sorted))
        surface_z_threshold = z_sorted[n_surface_layers-1] if len(z_sorted) > 0 else z_max - 2.0

        tags = np.zeros(len(prepared_atoms), dtype=int)
        for i, z in enumerate(z_coords):
            if z >= surface_z_threshold:
                tags[i] = 1  # 表面原子
            else:
                tags[i] = 0  # 体相原子

        prepared_atoms.set_tags(tags)

        # 确保 PBC 正确设置（表面应该在 xy 方向周期）
        prepared_atoms.pbc = True

        # 2. 检查是否需要 tile
        cell = prepared_atoms.cell
        min_ab = 8.0

        if cell[0, 0] < min_ab or cell[1, 1] < min_ab:
            # 需要 tile
            nx = int(np.ceil(min_ab / cell[0, 0]))
            ny = int(np.ceil(min_ab / cell[1, 1]))
            if nx > 1 or ny > 1:
                prepared_atoms = prepared_atoms * (nx, ny, 1)
                # 确保 PBC 仍然正确
                prepared_atoms.pbc = True
                self.logger.info(f"Tiled surface to {nx}x{ny}x1")

        # 3. 添加 constraints（固定体相原子）
        indices_to_fix = [i for i, tag in enumerate(prepared_atoms.get_tags()) if tag == 0]
        if indices_to_fix:
            constraint = FixAtoms(indices=indices_to_fix)
            prepared_atoms.set_constraint(constraint)
            self.logger.info(f"Fixed {len(indices_to_fix)} subsurface atoms")

        return prepared_atoms

    def _get_binding_atom_info(
        self,
        config: Atoms,
        n_surface: int
    ) -> Tuple[int, np.ndarray]:
        """
        获取吸附物的结合原子信息（z坐标最低的原子）

        Parameters
        ----------
        config : Atoms
            表面+吸附物的结构
        n_surface : int
            表面原子数量

        Returns
        -------
        Tuple[int, np.ndarray]
            (结合原子在吸附物中的索引, 结合原子的位置)
        """
        adsorbate_indices = list(range(n_surface, len(config)))
        if not adsorbate_indices:
            return -1, np.zeros(3)

        adsorbate_positions = config.positions[adsorbate_indices]
        binding_atom_local_idx = np.argmin(adsorbate_positions[:, 2])
        binding_atom_global_idx = adsorbate_indices[binding_atom_local_idx]
        binding_atom_pos = config.positions[binding_atom_global_idx].copy()

        return binding_atom_local_idx, binding_atom_pos

    def _check_binding_atom_drift(
        self,
        actual_pos: np.ndarray,
        target_pos: np.ndarray,
        xy_threshold: float = 1.0,
    ) -> Tuple[bool, float]:
        """
        检查结合原子是否从目标位点漂移

        Parameters
        ----------
        actual_pos : np.ndarray
            实际位置
        target_pos : np.ndarray
            目标位置
        xy_threshold : float
            xy平面内的漂移阈值（埃）

        Returns
        -------
        Tuple[bool, float]
            (是否漂移超出阈值, xy平面内的漂移距离)
        """
        xy_drift = np.sqrt((actual_pos[0] - target_pos[0])**2 +
                          (actual_pos[1] - target_pos[1])**2)
        is_drifted = xy_drift > xy_threshold
        return is_drifted, xy_drift

    def _optimize_with_xy_constraint(
        self,
        config: Atoms,
        n_surface: int,
        target_xy: np.ndarray,
        fmax: float,
        max_steps: int,
    ) -> Atoms:
        """
        使用xy方向软约束优化结构

        通过在每几步后重置结合原子的xy位置来实现软约束

        Parameters
        ----------
        config : Atoms
            初始结构
        n_surface : int
            表面原子数量
        target_xy : np.ndarray
            目标xy位置 [x, y]
        fmax : float
            力收敛标准
        max_steps : int
            最大优化步数

        Returns
        -------
        Atoms
            优化后的结构
        """
        from ase.constraints import FixCartesian

        atoms = config.copy()
        atoms.calc = self._calculator

        # 找到结合原子索引
        adsorbate_indices = list(range(n_surface, len(atoms)))
        if not adsorbate_indices:
            return atoms

        adsorbate_positions = atoms.positions[adsorbate_indices]
        binding_atom_local_idx = np.argmin(adsorbate_positions[:, 2])
        binding_atom_global_idx = adsorbate_indices[binding_atom_local_idx]

        # 保存原始约束
        original_constraints = atoms.constraints.copy() if atoms.constraints else []

        # 添加xy方向约束（只允许z方向移动）
        xy_constraint = FixCartesian(binding_atom_global_idx, mask=[True, True, False])
        atoms.set_constraint(original_constraints + [xy_constraint])

        # 先将结合原子移动到目标xy位置
        current_pos = atoms.positions[binding_atom_global_idx].copy()
        atoms.positions[binding_atom_global_idx, 0] = target_xy[0]
        atoms.positions[binding_atom_global_idx, 1] = target_xy[1]

        self.logger.info(f"  Applying xy-constraint at ({target_xy[0]:.3f}, {target_xy[1]:.3f})")
        self.logger.info(f"  Binding atom moved from ({current_pos[0]:.3f}, {current_pos[1]:.3f}) to target")

        # 优化
        opt = LBFGS(atoms)
        try:
            opt.run(fmax=fmax, steps=max_steps)
            self.logger.info(f"  Constrained optimization converged in {opt.get_number_of_steps()} steps")
        except Exception as e:
            self.logger.warning(f"  Constrained optimization warning: {e}")

        # 恢复原始约束
        atoms.set_constraint(original_constraints)

        return atoms

    def _predict_adsorption_energy_manual(
        self,
        surface_atoms: Atoms,
        surface_energy: float,
        adsorbate: Union[str, Atoms],
        adsorbate_atoms: Atoms,
        adsorbate_ref_energy: float,
        num_sites: int,
        fmax: float,
        max_steps: int,
        fixed_site_position: Optional[np.ndarray] = None,
    ) -> AdsorptionResult:
        """
        手动放置吸附质并预测吸附能（后备方法）

        当 Fairchem 的 AdsorbateSlabConfig 失败时使用此方法
        或当提供了固定吸附位点时直接使用
        """
        self.logger.info("Using manual adsorbate placement")

        # 如果提供了固定吸附位点，只使用这个位点
        if fixed_site_position is not None:
            self.logger.info(f"Using fixed adsorption site at {fixed_site_position}")
            sites = [fixed_site_position]
        else:
            # 在表面最高的几个原子周围生成位点
            z_coords = surface_atoms.positions[:, 2]
            z_max = z_coords.max()

            # 找到表面原子（z > z_max - 1.5 埃）
            surface_atom_indices = np.where(z_coords > z_max - 1.5)[0]

            if len(surface_atom_indices) == 0:
                raise RuntimeError("No surface atoms found")

            # 生成吸附位点：在表面原子上方
            sites = []
            for idx in surface_atom_indices[:num_sites]:
                site = surface_atoms.positions[idx].copy()
                site[2] = z_max + 2.0  # 在表面上方 2 埃
                sites.append(site)

            # 如果表面原子数量不足，在它们周围生成更多位点
            if len(sites) < num_sites:
                # 在现有位点周围添加偏移位点
                for i in range(len(sites), num_sites):
                    base_site = sites[i % len(sites)].copy()
                    # 添加随机偏移
                    base_site[0] += np.random.uniform(-1.5, 1.5)
                    base_site[1] += np.random.uniform(-1.5, 1.5)
                    sites.append(base_site)

        # 优化每个配置
        all_energies = []
        all_configs = []

        for i, site in enumerate(sites):
            self.logger.info(f"Relaxing configuration {i+1}/{len(sites)}...")

            try:
                # 创建配置
                config = surface_atoms.copy()
                ads_copy = adsorbate_atoms.copy()

                # 将吸附质的结合原子（z坐标最低的原子）移动到位点
                # 这确保了不同吸附质使用相同的结合位点
                ads_positions = ads_copy.get_positions()
                binding_atom_idx = np.argmin(ads_positions[:, 2])  # z坐标最低的原子
                binding_atom_pos = ads_positions[binding_atom_idx]

                # 计算需要的平移量，使结合原子位于site位置
                translation = site - binding_atom_pos
                ads_copy.translate(translation)

                self.logger.info(f"  Placed binding atom ({ads_copy.get_chemical_symbols()[binding_atom_idx]}) at site")

                # 合并（确保 PBC 一致）
                n_surface = len(config)
                config.extend(ads_copy)
                # 确保整个系统的 PBC 正确设置
                config.pbc = surface_atoms.pbc

                # 正确设置tags：保留表面原有的tags，只为吸附质设置tag=2
                if surface_atoms.has('tags'):
                    # 使用表面原有的tags
                    surface_tags = list(surface_atoms.get_tags())
                    self.logger.debug(f"Using original surface tags: {np.unique(surface_tags)}")
                else:
                    # 如果没有tags，从constraints推断
                    if config.constraints:
                        fixed_indices = []
                        for constraint in config.constraints:
                            if hasattr(constraint, 'index'):
                                fixed_indices.extend(constraint.index)
                        surface_tags = [0 if i in fixed_indices else 1 for i in range(n_surface)]
                    else:
                        # 没有constraints和tags，假设底部原子是bulk
                        surface_tags = [1] * n_surface

                # 吸附质原子tag=2
                tags = surface_tags + [2] * len(ads_copy)
                config.set_tags(tags)

                # 验证tags设置
                set_tags = config.get_tags()
                self.logger.info(f"Configuration {i+1} tags before optimization: unique={np.unique(set_tags)}, counts={np.bincount(set_tags)}")

                # 优化
                result = self.predict_energy(
                    config,
                    relax=True,
                    fmax=fmax,
                    max_steps=max_steps,
                )

                # 验证优化后tags
                final_tags = result.final_structure.get_tags()
                self.logger.info(f"Configuration {i+1} tags after optimization: unique={np.unique(final_tags)}, counts={np.bincount(final_tags)}")

                final_config = result.final_structure
                final_energy = result.energy

                # 检查结合原子是否从目标位点漂移（仅在使用固定位点时检查）
                if fixed_site_position is not None:
                    _, actual_binding_pos = self._get_binding_atom_info(final_config, n_surface)
                    is_drifted, xy_drift = self._check_binding_atom_drift(
                        actual_binding_pos, site, xy_threshold=1.0
                    )

                    self.logger.info(f"  Binding atom position after relaxation: ({actual_binding_pos[0]:.3f}, {actual_binding_pos[1]:.3f}, {actual_binding_pos[2]:.3f})")
                    self.logger.info(f"  Target site: ({site[0]:.3f}, {site[1]:.3f}, {site[2]:.3f})")
                    self.logger.info(f"  XY drift: {xy_drift:.3f} Å")

                    if is_drifted:
                        self.logger.warning(f"  ⚠️ Binding atom drifted {xy_drift:.2f} Å from target site! Re-optimizing with constraint...")

                        # 重新优化：使用xy约束
                        constrained_config = self._optimize_with_xy_constraint(
                            config=config,  # 使用初始配置（吸附物在目标位点）
                            n_surface=n_surface,
                            target_xy=site[:2],  # 只用xy坐标
                            fmax=fmax,
                            max_steps=max_steps,
                        )

                        # 计算约束优化后的能量
                        constrained_config.calc = self._calculator
                        constrained_energy = constrained_config.get_potential_energy()

                        # 检查约束后的位置
                        _, constrained_binding_pos = self._get_binding_atom_info(constrained_config, n_surface)
                        _, constrained_drift = self._check_binding_atom_drift(
                            constrained_binding_pos, site, xy_threshold=1.0
                        )

                        self.logger.info(f"  After constrained optimization:")
                        self.logger.info(f"    Position: ({constrained_binding_pos[0]:.3f}, {constrained_binding_pos[1]:.3f}, {constrained_binding_pos[2]:.3f})")
                        self.logger.info(f"    XY drift: {constrained_drift:.3f} Å")
                        self.logger.info(f"    Energy: {constrained_energy:.3f} eV (unconstrained: {final_energy:.3f} eV)")

                        # 使用约束优化的结果（保持位点一致性比能量更重要）
                        final_config = constrained_config
                        final_energy = constrained_energy

                        # 恢复tags
                        if final_config.has('tags'):
                            pass  # tags应该已经在约束优化中保留
                        else:
                            final_config.set_tags(tags)

                all_energies.append(final_energy)
                all_configs.append(final_config)

            except Exception as e:
                self.logger.warning(f"Failed to relax configuration {i+1}: {e}")
                continue

        if not all_energies:
            raise RuntimeError("All configurations failed to relax")

        # 找到最低能量配置
        best_idx = np.argmin(all_energies)
        best_energy = all_energies[best_idx]
        best_config = all_configs[best_idx]

        # 计算吸附能
        adsorption_energy = best_energy - surface_energy - adsorbate_ref_energy

        # Physical adsorption energy via vacuum-phase reference (see
        # predict_adsorption_energy for rationale). Keeps the automatic
        # (AdsorbateSlabConfig) and manual (fixed-site) paths consistent.
        gas_phase_energy = None
        physical_adsorption_energy = None
        gas_phase_freq_cm: Optional[List[float]] = None
        gas_phase_freq_ev: Optional[List[float]] = None
        try:
            gas_out = self._relax_adsorbate_in_vacuum(
                adsorbate_atoms, fmax=fmax, max_steps=max_steps,
                compute_vibrations=True,
            )
            if isinstance(gas_out, tuple):
                gas_phase_energy, gas_phase_freq_cm, _ = gas_out
            else:
                gas_phase_energy = gas_out
            if gas_phase_energy is not None and not (
                isinstance(gas_phase_energy, float) and (gas_phase_energy != gas_phase_energy)
            ):
                physical_adsorption_energy = float(
                    best_energy - surface_energy - gas_phase_energy
                )
                self.logger.info(
                    "  Physical E_ads (manual path) = %.3f eV  "
                    "(atomic-ref E_ads = %.3f eV)",
                    physical_adsorption_energy, adsorption_energy,
                )
            if gas_phase_freq_cm:
                gas_phase_freq_ev = [f / 8065.54429 for f in gas_phase_freq_cm]
        except Exception as exc:
            self.logger.warning(
                "Gas-phase reference unavailable for %s (manual): %s",
                adsorbate_atoms.get_chemical_formula(), exc,
            )

        return AdsorptionResult(
            adsorbate=str(adsorbate),
            surface_energy=surface_energy,
            adsorbate_surface_energy=best_energy,
            adsorbate_reference_energy=adsorbate_ref_energy,
            adsorption_energy=adsorption_energy,
            num_sites_tried=len(sites),
            best_site_index=best_idx,
            best_configuration=best_config,
            all_configurations=all_configs,
            all_energies=all_energies,
            is_anomalous=False,
            gas_phase_energy=gas_phase_energy,
            physical_adsorption_energy=physical_adsorption_energy,
            gas_phase_frequencies_ev=gas_phase_freq_ev,
            gas_phase_frequencies_cm=gas_phase_freq_cm,
        )

    def predict_energy(
        self,
        structure: Union[str, Atoms],
        relax: bool = True,
        fmax: float = 0.05,
        max_steps: int = 300,
        optimizer: str = "LBFGS",
        output_file: Optional[str] = None,
    ) -> EnergyResult:
        """
        预测单个结构的能量

        Parameters
        ----------
        structure : str or Atoms
            输入结构（文件路径或 ASE Atoms 对象）
        relax : bool, default=True
            是否进行结构优化
        fmax : float, default=0.05
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        optimizer : str, default="LBFGS"
            优化器类型 ("LBFGS" 或 "BFGS")
        output_file : str, optional
            输出文件路径

        Returns
        -------
        EnergyResult
            能量预测结果
        """
        self._load_model()

        # 读取结构
        if isinstance(structure, str):
            atoms = read(structure)
            structure_name = os.path.basename(structure)
        else:
            atoms = structure.copy()
            structure_name = atoms.get_chemical_formula()

        # 确保 PBC 设置一致（FAIRChemCalculator 要求）
        # PBC 必须是全 True 或全 False
        pbc = atoms.pbc
        if any(pbc) and not all(pbc):
            # 如果有部分为 True，则全部设为 True
            self.logger.warning(f"Inconsistent PBC {pbc}, setting all to True")
            atoms.pbc = True

        initial_structure = atoms.copy()
        atoms.calc = self._calculator

        # 保存原始tags（fairchem计算器可能会修改tags）
        original_tags = atoms.get_tags().copy() if atoms.has('tags') else None
        self.logger.debug(f"Before optimization - tags: {original_tags if original_tags is not None else 'NO TAGS'}")

        # 结构优化
        converged = False
        optimization_steps = 0

        if relax:
            self.logger.info(f"Relaxing structure: {structure_name}")

            if optimizer == "LBFGS":
                opt = LBFGS(atoms)
            else:
                opt = BFGS(atoms)

            try:
                opt.run(fmax=fmax, steps=max_steps)
                converged = True
                optimization_steps = opt.get_number_of_steps()
                self.logger.info(f"Converged in {optimization_steps} steps")
            except Exception as e:
                self.logger.warning(f"Optimization failed: {e}")
                optimization_steps = max_steps

        # 检查优化后的tags
        self.logger.debug(f"After optimization - tags: {atoms.get_tags() if atoms.has('tags') else 'NO TAGS'}")

        # 恢复原始tags（必须在获取能量前恢复，因为fairchem可能依赖tags）
        if original_tags is not None and len(original_tags) == len(atoms):
            atoms.set_tags(original_tags)
            self.logger.debug(f"After restoring tags - tags: {atoms.get_tags()}")
        else:
            if original_tags is not None:
                self.logger.warning(f"Cannot restore tags: length mismatch ({len(original_tags)} vs {len(atoms)})")

        # 获取最终能量和力
        try:
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
        except Exception as e:
            self.logger.error(f"Failed to get energy/forces: {e}")
            # 如果失败，尝试重新设置calculator
            atoms.calc = self._calculator
            if original_tags is not None:
                atoms.set_tags(original_tags)
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()

        # 保存结果
        if output_file:
            write(output_file, atoms)
            self.logger.info(f"Saved relaxed structure to {output_file}")

        return EnergyResult(
            structure_name=structure_name,
            initial_structure=initial_structure,
            final_structure=atoms,
            energy=energy,
            forces=forces,
            relaxed=relax,
            converged=converged,
            optimization_steps=optimization_steps,
            output_file=output_file,
        )

    def predict_adsorption_energy(
        self,
        surface: Union[str, Atoms],
        adsorbate: Union[str, Atoms],
        num_sites: int = 10,
        mode: str = "random_site_heuristic_placement",
        fmax: float = 0.05,
        max_steps: int = 300,
        detect_anomalies: bool = True,
        fixed_site_position: Optional[np.ndarray] = None,
        skip_surface_relaxation: bool = False,
    ) -> AdsorptionResult:
        """
        预测吸附能

        Parameters
        ----------
        surface : str or Atoms
            清洁表面结构
        adsorbate : str or Atoms
            吸附质（SMILES 字符串或 Atoms 对象）
        num_sites : int, default=10
            尝试的吸附位点数量
        mode : str, default="random_site_heuristic_placement"
            位点生成模式
        fmax : float, default=0.05
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        detect_anomalies : bool, default=True
            是否检测异常（解吸、解离等）
        fixed_site_position : np.ndarray, optional
            固定的吸附位点坐标 (x, y, z)，如果提供则只在此位点添加吸附物
        skip_surface_relaxation : bool, default=False
            是否跳过表面优化（当表面已经优化过时使用）

        Returns
        -------
        AdsorptionResult
            吸附能预测结果
        """
        self._load_model()

        # 读取表面结构
        if isinstance(surface, str):
            surface_atoms = read(surface)
        else:
            surface_atoms = surface.copy()

        # 1. 优化清洁表面（如果需要）
        if skip_surface_relaxation:
            self.logger.info("Skipping surface relaxation (using pre-optimized surface)")
            # 保留原始surface，包括tags
            surface_atoms_with_tags = surface_atoms.copy()

            # 计算表面能量
            surface_atoms_copy = surface_atoms.copy()
            surface_atoms_copy.calc = self._calculator
            surface_energy = surface_atoms_copy.get_potential_energy()

            # 创建result，确保final_structure保留tags
            surface_result = EnergyResult(
                structure_name=surface_atoms.get_chemical_formula(),
                initial_structure=surface_atoms,
                final_structure=surface_atoms_with_tags,  # 保留tags的surface
                energy=surface_energy,
                forces=np.zeros((len(surface_atoms), 3)),
                relaxed=False,
                converged=True,
                optimization_steps=0,
            )

            self.logger.info(f"  Surface energy: {surface_energy:.3f} eV")
            if surface_atoms_with_tags.has('tags'):
                tags = surface_atoms_with_tags.get_tags()
                self.logger.info(f"  Surface tags: {np.unique(tags)} (counts: {np.bincount(tags)})")
        else:
            self.logger.info("Relaxing clean surface...")
            surface_result = self.predict_energy(
                surface_atoms,
                relax=True,
                fmax=fmax,
                max_steps=max_steps,
            )
            surface_energy = surface_result.energy

        # 2. 创建吸附质对象
        adsorbate_obj = None
        if isinstance(adsorbate, str):
            from fairchem.data.oc.core.adsorbate import Adsorbate
            from core.pathway.adsorbate_resolver import resolve_adsorbate

            resolved_atoms, resolved_binding = resolve_adsorbate(adsorbate)
            adsorbate_obj = Adsorbate(
                adsorbate_atoms=resolved_atoms,
                adsorbate_binding_indices=resolved_binding,
            )
            adsorbate_atoms = adsorbate_obj.atoms
        else:
            adsorbate_atoms = adsorbate.copy()

        # 3. 计算吸附质参考能量
        adsorbate_symbols = adsorbate_atoms.get_chemical_symbols()
        adsorbate_ref_energy = sum([
            self.ATOMIC_REFERENCE_ENERGIES.get(sym, 0.0)
            for sym in adsorbate_symbols
        ])

        # 如果提供了固定吸附位点，直接使用手动放置方法
        if fixed_site_position is not None:
            self.logger.info(f"Using fixed adsorption site at {fixed_site_position}")
            return self._predict_adsorption_energy_manual(
                surface_atoms=surface_result.final_structure,
                surface_energy=surface_energy,
                adsorbate=adsorbate,
                adsorbate_atoms=adsorbate_atoms,
                adsorbate_ref_energy=adsorbate_ref_energy,
                num_sites=1,  # 只用一个位点
                fmax=fmax,
                max_steps=max_steps,
                fixed_site_position=fixed_site_position,
            )

        try:
            from fairchem.data.oc.core import (
                Adsorbate, Slab, AdsorbateSlabConfig,
            )
        except ImportError as e:
            raise ImportError(
                f"Failed to import fairchem.data.oc modules. Error: {e}"
            )

        # 为了使用 AdsorbateSlabConfig，我们需要确保表面满足要求：
        # - 被 tagged（表面原子 tag=1，体相原子 tag=0）
        # - 被 tiled（xy 方向 >= 8 埃）
        # - 添加了 constraints
        slab_atoms = self._prepare_slab_for_adsorbate_placement(
            surface_result.final_structure
        )

        # 创建 Slab 对象
        try:
            slab = Slab(slab_atoms=slab_atoms)
        except Exception as e:
            self.logger.warning(f"Failed to create Slab object: {e}")
            self.logger.info("Falling back to manual adsorbate placement")
            # 如果 Slab 创建失败，使用手动放置方法
            return self._predict_adsorption_energy_manual(
                surface_atoms=surface_result.final_structure,
                surface_energy=surface_energy,
                adsorbate=adsorbate,
                adsorbate_atoms=adsorbate_atoms,
                adsorbate_ref_energy=adsorbate_ref_energy,
                num_sites=num_sites,
                fmax=fmax,
                max_steps=max_steps,
            )

        self.logger.info(f"Generating {num_sites} adsorbate configurations...")

        # 创建吸附质对象（如果还没有）
        if isinstance(adsorbate, str):
            ads_for_config = adsorbate_obj
        else:
            # 需要创建 Adsorbate 对象
            try:
                ads_for_config = Adsorbate(
                    adsorbate_atoms=adsorbate_atoms,
                    adsorbate_binding_indices=[0]  # 假设第一个原子是绑定原子
                )
            except Exception as e:
                self.logger.warning(f"Failed to create Adsorbate object: {e}")
                # 使用手动方法
                return self._predict_adsorption_energy_manual(
                    surface_atoms=surface_result.final_structure,
                    surface_energy=surface_energy,
                    adsorbate=adsorbate,
                    adsorbate_atoms=adsorbate_atoms,
                    adsorbate_ref_energy=adsorbate_ref_energy,
                    num_sites=num_sites,
                    fmax=fmax,
                    max_steps=max_steps,
                )

        try:
            configs = AdsorbateSlabConfig(
                slab=slab,
                adsorbate=ads_for_config,
                mode=mode,
                num_sites=num_sites,
                num_augmentations_per_site=1,
                interstitial_gap=0.1,
            ).atoms_list
        except Exception as e:
            self.logger.warning(f"AdsorbateSlabConfig failed: {e}")
            self.logger.info("Falling back to manual adsorbate placement")
            # 使用手动放置方法
            return self._predict_adsorption_energy_manual(
                surface_atoms=surface_result.final_structure,
                surface_energy=surface_energy,
                adsorbate=adsorbate,
                adsorbate_atoms=adsorbate_atoms,
                adsorbate_ref_energy=adsorbate_ref_energy,
                num_sites=num_sites,
                fmax=fmax,
                max_steps=max_steps,
            )

        # 5. 优化所有配置
        all_energies = []
        all_configs = []

        for i, config in enumerate(configs):
            self.logger.info(f"Relaxing configuration {i+1}/{len(configs)}...")

            try:
                result = self.predict_energy(
                    config,
                    relax=True,
                    fmax=fmax,
                    max_steps=max_steps,
                )
                all_energies.append(result.energy)
                all_configs.append(result.final_structure)
            except Exception as e:
                self.logger.warning(f"Failed to relax configuration {i+1}: {e}")
                continue

        if not all_energies:
            raise RuntimeError("All configurations failed to relax")

        # 6. 找到最低能量配置
        best_idx = np.argmin(all_energies)
        best_energy = all_energies[best_idx]
        best_config = all_configs[best_idx]

        # 7. 计算吸附能: E_ads = E(ads+surf) - E(surf) - E(ads_ref)
        adsorption_energy = best_energy - surface_energy - adsorbate_ref_energy

        # 7b. Physical adsorption energy: relax the adsorbate in a vacuum box
        # and reference against that. Yields E_ads consistent with the "gas
        # species at formation_energy = 0" convention CatMAP expects.
        gas_phase_energy = None
        physical_adsorption_energy = None
        gas_phase_freq_cm: Optional[List[float]] = None
        gas_phase_freq_ev: Optional[List[float]] = None
        try:
            gas_out = self._relax_adsorbate_in_vacuum(
                adsorbate_atoms, fmax=fmax, max_steps=max_steps,
                compute_vibrations=True,
            )
            if isinstance(gas_out, tuple):
                gas_phase_energy, gas_phase_freq_cm, _ = gas_out
            else:
                gas_phase_energy = gas_out
            if gas_phase_energy is not None and not (
                isinstance(gas_phase_energy, float) and (gas_phase_energy != gas_phase_energy)
            ):
                physical_adsorption_energy = float(
                    best_energy - surface_energy - gas_phase_energy
                )
                self.logger.info(
                    "  Physical E_ads = %.3f eV  (atomic-ref E_ads = %.3f eV)",
                    physical_adsorption_energy, adsorption_energy,
                )
            if gas_phase_freq_cm:
                gas_phase_freq_ev = [f / 8065.54429 for f in gas_phase_freq_cm]
        except Exception as exc:
            self.logger.warning(
                "Gas-phase reference unavailable for %s: %s — microkinetic "
                "model will fall back to atomic-reference formation energy.",
                adsorbate_atoms.get_chemical_formula(), exc,
            )

        # 8. 异常检测
        is_anomalous = False
        if detect_anomalies:
            try:
                from fairchem.data.oc.utils import DetectTrajAnomaly

                detector = DetectTrajAnomaly(
                    configs[best_idx],
                    best_config,
                    best_config.get_tags()
                )

                is_anomalous = (
                    detector.is_adsorbate_dissociated() or
                    detector.is_adsorbate_desorbed() or
                    detector.has_surface_changed() or
                    detector.is_adsorbate_intercalated()
                )

                if is_anomalous:
                    self.logger.warning("Anomaly detected in relaxation!")
            except Exception as e:
                self.logger.warning(f"Anomaly detection failed: {e}")

        return AdsorptionResult(
            adsorbate=str(adsorbate),
            surface_energy=surface_energy,
            adsorbate_surface_energy=best_energy,
            adsorbate_reference_energy=adsorbate_ref_energy,
            adsorption_energy=adsorption_energy,
            num_sites_tried=len(configs),
            best_site_index=best_idx,
            best_configuration=best_config,
            all_configurations=all_configs,
            all_energies=all_energies,
            is_anomalous=is_anomalous,
            gas_phase_energy=gas_phase_energy,
            physical_adsorption_energy=physical_adsorption_energy,
            gas_phase_frequencies_ev=gas_phase_freq_ev,
            gas_phase_frequencies_cm=gas_phase_freq_cm,
        )

    def predict_pathway_energies(
        self,
        surface: Union[str, Atoms],
        adsorbates: List[Union[str, Atoms]],
        num_sites: int = 5,
        fmax: float = 0.05,
        max_steps: int = 300,
        output_dir: Optional[str] = None,
        fixed_adsorption_site: Optional[np.ndarray] = None,
        skip_surface_relaxation: bool = False,
    ) -> PathwayResult:
        """
        预测反应路径上多个吸附质的能量

        Parameters
        ----------
        surface : str or Atoms
            催化剂表面结构（clean表面，不含吸附物）
        adsorbates : list of str or Atoms
            吸附质列表
        num_sites : int, default=5
            每个吸附质尝试的位点数（如果提供fixed_adsorption_site则忽略）
        fmax : float, default=0.05
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录
        fixed_adsorption_site : np.ndarray, optional
            固定的吸附位点坐标 (x, y, z)，如果提供则所有吸附物都使用此位点
        skip_surface_relaxation : bool, default=False
            是否跳过该函数内的 clean-surface relax（已优化表面可启用以提速/稳健）

        Returns
        -------
        PathwayResult
            反应路径能量预测结果
        """
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "pathway_results")
        os.makedirs(output_dir, exist_ok=True)
        self._load_model()

        # 读取表面
        if isinstance(surface, str):
            surface_atoms = read(surface)
        else:
            surface_atoms = surface.copy()

        surface_formula = surface_atoms.get_chemical_formula()

        # 1. 处理清洁表面
        self.logger.info("=" * 60)
        if skip_surface_relaxation:
            self.logger.info("Skipping clean surface relaxation in pathway prediction.")
            self._load_model()
            surface_relaxed = surface_atoms.copy()
            surface_calc = surface_atoms.copy()
            surface_calc.calc = self._calculator
            surface_energy = float(surface_calc.get_potential_energy())
        else:
            self.logger.info("Relaxing clean surface...")
            surface_result = self.predict_energy(
                surface_atoms,
                relax=True,
                fmax=fmax,
                max_steps=max_steps,
            )
            surface_relaxed = surface_result.final_structure
            surface_energy = surface_result.energy

        # 保存清洁表面
        surface_file = os.path.join(output_dir, "surface_relaxed.vasp")
        write(surface_file, surface_relaxed, format="vasp")

        # 2. 预测每个吸附质
        adsorbate_energies = {}
        adsorption_energies = {}
        best_configurations = {}
        # gas-phase (vacuum) frequencies in cm^-1, keyed by adsorbate label
        gas_frequencies: Dict[str, List[float]] = {}
        failed_adsorbates: List[Tuple[str, str]] = []

        # 跟踪实际的参考吸附位点（从第一个非"*"吸附物的结果中获取）
        # 这确保了所有后续吸附物使用相同的位点，即使优化过程中发生了微小漂移
        actual_reference_site = None
        n_surface_atoms = len(surface_relaxed)

        for i, ads in enumerate(adsorbates):
            self.logger.info("=" * 60)
            self.logger.info(f"Processing adsorbate {i+1}/{len(adsorbates)}: {ads}")

            try:
                # 特殊处理：如果吸附物是 "*"
                ads_str = str(ads).strip()
                if ads_str == "*" or ads_str == "":
                    ads_name = "*"

                    # 判断是吸附的起始态还是脱附的产物态
                    if i == 0:
                        # 第一个中间体是"*"：吸附反应的起始态（clean surface）
                        self.logger.info("  Initial state (*) - clean surface")
                        adsorbate_energies[ads_name] = surface_energy
                        adsorption_energies[ads_name] = 0.0
                        best_configurations[ads_name] = surface_relaxed.copy()

                        ads_file = os.path.join(output_dir, f"{ads_name}_relaxed.vasp")
                        write(ads_file, surface_relaxed, format="vasp")
                        self.logger.info(f"  E_surface = {surface_energy:.3f} eV")
                        continue
                    else:
                        # 中间或末尾的"*"：脱附反应的产物态
                        # 需要将前一个中间体的分子放在表面上方3-4埃（代表气相/弱吸附态）
                        self.logger.info("  Desorption product (*) - molecule above surface")

                        # 获取前一个中间体
                        prev_ads_name = str(adsorbates[i-1]).strip()
                        if prev_ads_name in best_configurations:
                            prev_config = best_configurations[prev_ads_name]

                            # 识别前一个配置中的吸附分子原子
                            n_surface = len(surface_relaxed)
                            adsorbate_indices = list(range(n_surface, len(prev_config)))

                            if adsorbate_indices:
                                # 创建脱附态：分子在表面上方3.5埃
                                desorbed_config = surface_relaxed.copy()

                                # 提取吸附分子（需要先复制再修改位置）
                                ads_atoms_list = []
                                for idx in adsorbate_indices:
                                    atom = prev_config[idx]
                                    ads_atoms_list.append((atom.symbol, atom.position.copy()))

                                # 计算需要的z方向移动
                                surface_top_z = surface_relaxed.positions[:, 2].max()
                                ads_positions = np.array([pos for _, pos in ads_atoms_list])
                                ads_center_z = ads_positions[:, 2].mean()
                                z_shift = (surface_top_z + 3.5) - ads_center_z

                                # 合并表面和移位后的分子
                                from ase import Atom
                                for symbol, pos in ads_atoms_list:
                                    new_pos = pos.copy()
                                    new_pos[2] += z_shift
                                    desorbed_config.append(Atom(symbol, position=new_pos))

                                # 计算脱附态能量（不优化，因为这是gas phase）
                                self.logger.info("  Computing energy of desorbed state (molecule at 3.5 Å above surface)...")
                                desorbed_result = self.predict_energy(
                                    desorbed_config,
                                    relax=False,  # 不优化，保持gas phase几何
                                )

                                adsorbate_energies[ads_name] = desorbed_result.energy
                                # 脱附能 = E_desorbed - E_surface (应该是正值)
                                adsorption_energies[ads_name] = desorbed_result.energy - surface_energy
                                best_configurations[ads_name] = desorbed_result.final_structure

                                # 保存脱附态结构
                                ads_file = os.path.join(output_dir, f"{ads_name}_relaxed.vasp")
                                write(ads_file, desorbed_result.final_structure, format="vasp")

                                self.logger.info(f"  E_desorbed = {desorbed_result.energy:.3f} eV")
                                self.logger.info(f"  E_des = {adsorption_energies[ads_name]:.3f} eV (desorption energy)")
                                continue
                            else:
                                self.logger.warning(f"  Previous configuration {prev_ads_name} has no adsorbate atoms")
                        else:
                                self.logger.warning(f"  Previous adsorbate {prev_ads_name} not found in configurations")

                        # Fallback: use clean surface
                        self.logger.warning("  Falling back to clean surface for '*'")
                        adsorbate_energies[ads_name] = surface_energy
                        adsorption_energies[ads_name] = 0.0
                        best_configurations[ads_name] = surface_relaxed.copy()

                        ads_file = os.path.join(output_dir, f"{ads_name}_relaxed.vasp")
                        write(ads_file, surface_relaxed, format="vasp")
                        continue

                if fixed_adsorption_site is not None:
                    # 使用固定吸附位点
                    # 优先使用actual_reference_site（如果已从之前的吸附物获取）
                    # 这确保了所有吸附物使用完全相同的位点
                    site_to_use = actual_reference_site if actual_reference_site is not None else fixed_adsorption_site
                    self.logger.info(f"  Using fixed adsorption site at: {site_to_use}")
                    if actual_reference_site is not None:
                        self.logger.info(f"  (Updated from first intermediate's actual position)")

                    ads_result = self.predict_adsorption_energy(
                        surface=surface_relaxed,
                        adsorbate=ads,
                        num_sites=1,  # 只试一个位点
                        fmax=fmax,
                        max_steps=max_steps,
                        fixed_site_position=site_to_use,  # 传递固定位点
                        skip_surface_relaxation=True,  # 跳过表面优化（已经优化过了）
                    )
                else:
                    # 正常模式 - 尝试多个吸附位点
                    ads_result = self.predict_adsorption_energy(
                        surface=surface_relaxed,
                        adsorbate=ads,
                        num_sites=num_sites,
                        fmax=fmax,
                        max_steps=max_steps,
                        skip_surface_relaxation=True,  # 跳过表面优化（已经优化过了）
                    )

                ads_name = str(ads)
                adsorbate_energies[ads_name] = ads_result.adsorbate_surface_energy
                # Prefer the physical E_ads (gas-phase-referenced) when the
                # vacuum relaxation succeeded. The atomic-reference fallback
                # is kept if vacuum relaxation failed.
                if ads_result.physical_adsorption_energy is not None:
                    adsorption_energies[ads_name] = float(
                        ads_result.physical_adsorption_energy
                    )
                    self.logger.info(
                        "  E_ads (gas-phase ref) = %.3f eV  "
                        "(atomic-ref: %.3f eV)",
                        ads_result.physical_adsorption_energy,
                        ads_result.adsorption_energy,
                    )
                else:
                    adsorption_energies[ads_name] = ads_result.adsorption_energy
                    self.logger.info(
                        "  E_ads (atomic ref, no gas-phase calc) = %.3f eV",
                        ads_result.adsorption_energy,
                    )
                best_configurations[ads_name] = ads_result.best_configuration
                if ads_result.gas_phase_frequencies_cm:
                    gas_frequencies[ads_name] = list(ads_result.gas_phase_frequencies_cm)

                # 保存最佳配置
                ads_file = os.path.join(output_dir, f"{ads_name}_relaxed.vasp")
                write(ads_file, ads_result.best_configuration, format="vasp")

                # 更新参考吸附位点（从第一个成功处理的非"*"吸附物获取）
                # 这个位点将用于所有后续吸附物，确保位点一致性
                if actual_reference_site is None and fixed_adsorption_site is not None:
                    _, binding_pos = self._get_binding_atom_info(
                        ads_result.best_configuration, n_surface_atoms
                    )
                    actual_reference_site = binding_pos.copy()
                    self.logger.info(f"  ✓ Reference site updated from {ads_name}: ({actual_reference_site[0]:.3f}, {actual_reference_site[1]:.3f}, {actual_reference_site[2]:.3f})")
                    self.logger.info(f"    All subsequent intermediates will use this site for consistency")

            except Exception as e:
                self.logger.exception(f"Failed to process {ads}: {e}")
                failed_adsorbates.append((str(ads), str(e)))
                continue

        # 3. 计算相对能量（相对于最稳定的吸附质）
        if not adsorbate_energies:
            detail = "; ".join([f"{name}: {err}" for name, err in failed_adsorbates]) or "unknown"
            raise RuntimeError(
                "No adsorbate energies were produced in predict_pathway_energies. "
                f"All adsorbates failed. details={detail}"
            )
        min_energy = min(adsorbate_energies.values())
        relative_energies = {
            ads: e - min_energy
            for ads, e in adsorbate_energies.items()
        }

        return PathwayResult(
            surface_formula=surface_formula,
            surface_energy=surface_energy,
            adsorbates=list(adsorbate_energies.keys()),
            adsorbate_energies=adsorbate_energies,
            adsorption_energies=adsorption_energies,
            best_configurations=best_configurations,
            relative_energies=relative_energies,
            output_dir=output_dir,
            gas_frequencies=gas_frequencies,
        )

    def __del__(self):
        """清理临时文件"""
        if not self.keep_files and self._temp_dir and os.path.exists(self._temp_dir):
            import shutil
            shutil.rmtree(self._temp_dir)


if __name__ == "__main__":
    # 简单测试
    print("Fairchem Predictor module loaded successfully")
    print("Available models: uma-s-1p1, uma-m-1p1")
    print("See documentation for usage examples")
