"""
AdsorbDiff Predictor - 吸附位点预测封装类

这个模块提供了一个完整的封装类，用于预测最优吸附位点和取向。

使用方式:
    from adsorbdiff_predictor import AdsorbDiffPredictor

    # 初始化预测器
    predictor = AdsorbDiffPredictor(
        adsorbdiff_root="/path/to/AdsorbDiff",
        checkpoint_path="/path/to/checkpoint.pt",  # 可选，默认使用PaiNN模型
        use_gpu=True,
    )

    # 方式1: 表面+吸附物分别提供
    result = predictor.predict(
        surface="slab.vasp",
        adsorbate="*CO",  # 可以是SMILES字符串或ASE Atoms
        has_adsorbate=False,
    )

    # 方式2: 已有吸附物的表面结构
    result = predictor.predict(
        surface="ads_slab.vasp",
        has_adsorbate=True,  # 吸附物已在结构中
        adsorbate_tag=2,  # 吸附物原子的tag
    )

作者: Claude
"""

import os
import sys
import shutil
import tempfile
import pickle
import warnings
from typing import Union, Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field

import numpy as np


@dataclass
class AdsorptionSite:
    """单个吸附位点信息"""
    site_id: int
    position: np.ndarray  # 吸附位点坐标
    binding_type: str  # "ontop", "bridge", "hollow", "random"
    energy: Optional[float] = None  # 预测能量 (eV)


@dataclass
class AdsorptionResult:
    """吸附预测结果"""
    initial_structure: "ase.Atoms"  # 初始结构
    final_structure: "ase.Atoms"  # 最终弛豫结构
    adsorbate_position: np.ndarray  # 最终吸附位置 (COM)
    adsorbate_orientation: np.ndarray  # 最终取向
    energy: float  # 预测能量 (eV)
    forces_max: float  # 最大残余力 (eV/Å)
    trajectory_path: Optional[str] = None  # 轨迹文件路径
    site_info: Optional[AdsorptionSite] = None  # 初始位点信息

    def __repr__(self):
        return (f"AdsorptionResult(E={self.energy:.4f} eV, "
                f"Fmax={self.forces_max:.4f} eV/Å, "
                f"ads_pos={self.adsorbate_position})")


@dataclass
class PredictionOutput:
    """完整的预测输出"""
    surface_formula: str
    adsorbate_formula: str
    num_sites_sampled: int
    results: List[AdsorptionResult]
    best_result: AdsorptionResult
    output_dir: Optional[str] = None
    use_proxy_energy: bool = False  # 是否使用了代理能量

    def get_top_n(self, n: int) -> List[AdsorptionResult]:
        """获取能量最低的前N个结果"""
        sorted_results = sorted(self.results, key=lambda x: x.energy)
        return sorted_results[:n]

    def summary(self, top_n: int = 5) -> str:
        """生成预测结果摘要"""
        energy_label = "Proxy Energy" if self.use_proxy_energy else "Energy"
        lines = [
            f"Adsorption Prediction Results",
            f"="*60,
            f"Surface: {self.surface_formula}",
            f"Adsorbate: {self.adsorbate_formula}",
            f"Sites sampled: {self.num_sites_sampled}",
        ]
        if self.use_proxy_energy:
            lines.append(f"Note: Using proxy energy (no relaxation model provided)")
        lines.extend([
            f"",
            f"Best configuration:",
            f"  {energy_label}: {self.best_result.energy:.4f}" + (" (proxy)" if self.use_proxy_energy else " eV"),
            f"  Position: {self.best_result.adsorbate_position}",
            f"  Max force: {self.best_result.forces_max:.4f} eV/Å",
            f"",
            f"Top {min(top_n, len(self.results))} configurations:",
            f"-"*60,
        ])
        for i, r in enumerate(self.get_top_n(top_n), 1):
            lines.append(f"  {i}. {energy_label}={r.energy:.4f}, Fmax={r.forces_max:.4f} eV/Å")
        return "\n".join(lines)


class AdsorbDiffPredictor:
    """
    AdsorbDiff吸附位点预测器

    这个类封装了AdsorbDiff项目的所有功能，提供简洁的API用于预测最优吸附位点和取向。

    Parameters
    ----------
    adsorbdiff_root : str
        AdsorbDiff项目的根目录路径
    checkpoint_path : str, optional
        扩散模型checkpoint文件路径。如果不指定，将使用默认的PaiNN模型
    relax_checkpoint_path : str, optional
        弛豫模型checkpoint文件路径（用于能量和力的计算）。
        如果不指定，将使用代理能量进行配置排序（基于吸附位置而非实际能量计算）
    use_gpu : bool, default=True
        是否使用GPU进行推理
    cutoff : float, default=12.0
        截断半径 (Å)
    max_neighbors : int, default=50
        最大邻居数
    num_sites : int, default=10
        采样的吸附位点数量
    num_augmentations_per_site : int, default=1
        每个位点的增强数量
    placement_mode : str, default="heuristic"
        吸附物放置模式: "random", "heuristic", "random_site_heuristic_placement"
    interstitial_gap : float, default=2.0
        吸附物与表面的最小距离 (Å)
    diffusion_steps : int, default=100
        扩散步数
    ads_std_low : float, default=0.1
        吸附物位置噪声的最小标准差
    ads_std_high : float, default=10.0
        吸附物位置噪声的最大标准差
    rot_std_low : float, default=0.01
        旋转噪声的最小标准差
    rot_std_high : float, default=1.55
        旋转噪声的最大标准差
    relax_steps : int, default=100
        弛豫优化最大步数
    relax_fmax : float, default=0.05
        弛豫收敛力阈值 (eV/Å)
    seed : int, optional
        随机种子
    work_dir : str, optional
        工作目录。如果不指定，将创建临时目录
    keep_files : bool, default=False
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """

    def __init__(
        self,
        adsorbdiff_root: str,
        checkpoint_path: Optional[str] = None,
        relax_checkpoint_path: Optional[str] = None,
        use_gpu: bool = True,
        cutoff: float = 12.0,
        max_neighbors: int = 50,
        num_sites: int = 10,
        num_augmentations_per_site: int = 1,
        placement_mode: str = "heuristic",
        interstitial_gap: float = 0.5,  # Gap between adsorbate and surface (was 2.0, too high)
        diffusion_steps: int = 100,
        ads_std_low: float = 0.1,
        ads_std_high: float = 10.0,
        rot_std_low: float = 0.01,
        rot_std_high: float = 1.55,
        relax_steps: int = 100,
        relax_fmax: float = 0.05,
        seed: Optional[int] = 0,
        work_dir: Optional[str] = None,
        keep_files: bool = False,
        verbose: bool = True,
    ):
        self.adsorbdiff_root = os.path.abspath(adsorbdiff_root)
        self.use_gpu = use_gpu
        self.cutoff = cutoff
        self.max_neighbors = max_neighbors
        self.num_sites = num_sites
        self.num_augmentations_per_site = num_augmentations_per_site
        self.placement_mode = placement_mode
        self.interstitial_gap = interstitial_gap
        self.diffusion_steps = diffusion_steps
        self.ads_std_low = ads_std_low
        self.ads_std_high = ads_std_high
        self.rot_std_low = rot_std_low
        self.rot_std_high = rot_std_high
        self.relax_steps = relax_steps
        self.relax_fmax = relax_fmax
        self.seed = seed
        self.keep_files = keep_files
        self.verbose = verbose

        # 设置checkpoint路径
        if checkpoint_path is None:
            default_painn = os.path.join(self.adsorbdiff_root, "PT_zeroshot_painn.pt")
            default_eqv2 = os.path.join(self.adsorbdiff_root, "PT_fewshot_eqv2_cond.pt")
            if os.path.exists(default_painn):
                self.checkpoint_path = default_painn
                self.model_name = "PaiNN"
            elif os.path.exists(default_eqv2):
                self.checkpoint_path = default_eqv2
                self.model_name = "EquiformerV2"
            else:
                raise FileNotFoundError(
                    "未找到默认checkpoint文件。请指定checkpoint_path参数。"
                )
        else:
            self.checkpoint_path = os.path.abspath(checkpoint_path)
            self.model_name = os.path.basename(checkpoint_path)

        self.relax_checkpoint_path = relax_checkpoint_path

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="adsorbdiff_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 添加AdsorbDiff路径到系统路径
        if self.adsorbdiff_root not in sys.path:
            sys.path.insert(0, self.adsorbdiff_root)

        # 设置随机种子
        if seed is not None:
            np.random.seed(seed)

        # 验证环境
        self._validate_environment()

        # 延迟加载模型（第一次预测时加载）
        self._diff_calculator = None
        self._relax_calculator = None

        if self.verbose:
            print(f"AdsorbDiffPredictor initialized:")
            print(f"  AdsorbDiff root: {self.adsorbdiff_root}")
            print(f"  Model: {self.model_name}")
            print(f"  Checkpoint: {self.checkpoint_path}")
            print(f"  Use GPU: {self.use_gpu}")
            print(f"  Work dir: {self.work_dir}")

    def _validate_environment(self):
        """验证运行环境"""
        if not os.path.exists(self.adsorbdiff_root):
            raise FileNotFoundError(f"AdsorbDiff root not found: {self.adsorbdiff_root}")
        if not os.path.exists(self.checkpoint_path):
            raise FileNotFoundError(f"Checkpoint not found: {self.checkpoint_path}")

        # 检查必要的数据库文件
        ads_db = os.path.join(self.adsorbdiff_root, "adsorbdiff/placement/pkls/adsorbates.pkl")
        if not os.path.exists(ads_db):
            raise FileNotFoundError(f"Adsorbate database not found: {ads_db}")

    def _log(self, message: str):
        """输出日志"""
        if self.verbose:
            print(message)

    def _load_diffusion_calculator(self):
        """加载扩散模型"""
        if self._diff_calculator is None:
            self._log("Loading diffusion model...")
            from adsorbdiff import AdsorbDiffCalculator
            self._diff_calculator = AdsorbDiffCalculator(
                checkpoint_path=self.checkpoint_path,
                cutoff=self.cutoff,
                max_neighbors=self.max_neighbors,
                cpu=not self.use_gpu,
                seed=self.seed,
            )
            self._log(f"  Loaded {self.model_name} diffusion model")
        return self._diff_calculator

    def _load_relax_calculator(self):
        """加载弛豫模型"""
        if self._relax_calculator is None:
            if self.relax_checkpoint_path:
                self._log("Loading relaxation model...")
                # 尝试使用fairchem的OCPCalculator
                try:
                    from fairchem.core.common.relaxation.ase_utils import OCPCalculator
                    self._relax_calculator = OCPCalculator(
                        checkpoint_path=self.relax_checkpoint_path,
                        cpu=not self.use_gpu,
                    )
                    self._log("  Loaded relaxation model (fairchem OCPCalculator)")
                except ImportError:
                    # 回退到AdsorbDiff的calculator
                    from adsorbdiff import AdsorbDiffCalculator
                    self._relax_calculator = AdsorbDiffCalculator(
                        checkpoint_path=self.relax_checkpoint_path,
                        cutoff=self.cutoff,
                        max_neighbors=self.max_neighbors,
                        cpu=not self.use_gpu,
                        seed=self.seed,
                    )
                    self._log("  Loaded relaxation model (AdsorbDiff)")
            else:
                # 没有弛豫模型时返回None，使用代理能量
                self._log("  No relaxation model provided, using proxy energy")
                return None
        return self._relax_calculator

    def _convert_surface(self, surface) -> "ase.Atoms":
        """将输入表面转换为ASE Atoms"""
        from ase import Atoms
        from ase.io import read

        if isinstance(surface, str):
            if not os.path.exists(surface):
                raise FileNotFoundError(f"Surface file not found: {surface}")
            return read(surface)
        elif isinstance(surface, Atoms):
            return surface.copy()
        else:
            # 尝试从pymatgen Structure转换
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

    def _convert_adsorbate(self, adsorbate) -> "ase.Atoms":
        """将输入吸附物转换为ASE Atoms"""
        from ase import Atoms

        if isinstance(adsorbate, str):
            # 尝试作为SMILES解析
            if adsorbate.startswith("*") or len(adsorbate) <= 10:
                # 从数据库查找
                return self._get_adsorbate_from_smiles(adsorbate)
            elif os.path.exists(adsorbate):
                from ase.io import read
                return read(adsorbate)
            else:
                return self._get_adsorbate_from_smiles(adsorbate)
        elif isinstance(adsorbate, Atoms):
            return adsorbate.copy()
        elif isinstance(adsorbate, int):
            # 从数据库按ID查找
            return self._get_adsorbate_from_db(adsorbate)
        else:
            raise TypeError(
                f"Unsupported adsorbate type: {type(adsorbate)}. "
                "Expected str (SMILES or file path), int (database ID), or ase.Atoms"
            )

    def _get_adsorbate_from_smiles(self, smiles: str) -> "ase.Atoms":
        """从数据库获取吸附物（按SMILES）"""
        ads_db_path = os.path.join(
            self.adsorbdiff_root, "adsorbdiff/placement/pkls/adsorbates.pkl"
        )
        with open(ads_db_path, 'rb') as f:
            ads_db = pickle.load(f)

        # 搜索匹配的SMILES
        for idx, ads_info in ads_db.items():
            if ads_info[1] == smiles:
                self._log(f"  Found adsorbate in database: {smiles} (id={idx})")
                return ads_info[0].copy()

        # 如果没找到，列出可用的吸附物
        available = [ads_info[1] for ads_info in ads_db.values()][:20]
        raise ValueError(
            f"Adsorbate '{smiles}' not found in database. "
            f"Available adsorbates (first 20): {available}"
        )

    def _get_adsorbate_from_db(self, idx: int) -> "ase.Atoms":
        """从数据库获取吸附物（按ID）"""
        ads_db_path = os.path.join(
            self.adsorbdiff_root, "adsorbdiff/placement/pkls/adsorbates.pkl"
        )
        with open(ads_db_path, 'rb') as f:
            ads_db = pickle.load(f)

        if idx not in ads_db:
            raise ValueError(f"Adsorbate ID {idx} not found in database")

        return ads_db[idx][0].copy()

    def _extract_adsorbate_from_surface(
        self,
        atoms: "ase.Atoms",
        adsorbate_tag: int = 2
    ) -> Tuple["ase.Atoms", "ase.Atoms"]:
        """从表面结构中提取吸附物"""
        tags = atoms.get_tags()
        ads_mask = tags == adsorbate_tag
        slab_mask = ~ads_mask

        slab = atoms[slab_mask]
        adsorbate = atoms[ads_mask]

        # 设置正确的tags
        slab.set_tags([1 if t == 1 else 0 for t in slab.get_tags()])

        return slab, adsorbate

    def _tag_surface_atoms(self, atoms: "ase.Atoms") -> "ase.Atoms":
        """为表面原子添加tags (1=表面, 0=subsurface)"""
        tags = atoms.get_tags()
        if np.any(tags == 1):
            # 已经有tags
            return atoms

        # 使用简单的高度判断
        positions = atoms.get_positions()
        z_coords = positions[:, 2]
        z_max = z_coords.max()
        z_threshold = z_max - 2.5  # 表面2.5 Å以内的原子

        new_tags = np.where(z_coords > z_threshold, 1, 0)
        atoms.set_tags(new_tags)
        return atoms

    def _create_adslab_config(
        self,
        slab_atoms: "ase.Atoms",
        adsorbate_atoms: "ase.Atoms",
        binding_indices: List[int] = None,
    ) -> List["ase.Atoms"]:
        """创建吸附物-表面配置"""
        from adsorbdiff.placement import Adsorbate, AdsorbateSlabConfig, Slab, Bulk

        # 创建临时Bulk对象
        class DummyBulk:
            def __init__(self, atoms):
                self.atoms = atoms
                self.src_id = "custom"

        bulk = DummyBulk(slab_atoms)

        # 创建Slab对象（绕过验证）
        slab = object.__new__(Slab)
        slab.bulk = bulk
        slab.atoms = slab_atoms
        slab.millers = (1, 1, 1)
        slab.shift = 0.0
        slab.top = True

        # 确保slab有surface tags
        if not np.any(slab_atoms.get_tags() == 1):
            slab.atoms = self._tag_surface_atoms(slab_atoms)

        # 创建Adsorbate对象
        adsorbate = Adsorbate(
            adsorbate_atoms=adsorbate_atoms,
            adsorbate_binding_indices=binding_indices or [0],
        )

        # 创建吸附配置
        config = AdsorbateSlabConfig(
            slab=slab,
            adsorbate=adsorbate,
            num_sites=self.num_sites,
            num_augmentations_per_site=self.num_augmentations_per_site,
            interstitial_gap=self.interstitial_gap,
            mode=self.placement_mode,
        )

        return config.atoms_list, config.metadata_list

    def _run_diffusion(
        self,
        adslab: "ase.Atoms",
        traj_dir: str,
    ) -> "ase.Atoms":
        """运行扩散采样"""
        calc = self._load_diffusion_calculator()
        return calc.run_diffusion(adslab, trajectory=traj_dir)

    def _run_relaxation(
        self,
        adslab: "ase.Atoms",
        traj_path: str = None,
    ) -> Tuple["ase.Atoms", float, float]:
        """运行结构弛豫"""
        # 如果没有弛豫模型，使用代理能量
        if self.relax_checkpoint_path is None:
            return self._compute_proxy_energy(adslab)

        from ase.optimize import BFGS

        calc = self._load_relax_calculator()
        adslab.calc = calc

        if traj_path:
            opt = BFGS(adslab, trajectory=traj_path, logfile=None)
        else:
            opt = BFGS(adslab, logfile=None)

        try:
            opt.run(fmax=self.relax_fmax, steps=self.relax_steps)
        except Exception as e:
            self._log(f"  Warning: Relaxation error: {e}")

        energy = adslab.get_potential_energy()
        forces = adslab.get_forces()
        fmax = np.max(np.linalg.norm(forces, axis=1))

        return adslab, energy, fmax

    def _compute_proxy_energy(
        self,
        adslab: "ase.Atoms",
    ) -> Tuple["ase.Atoms", float, float]:
        """
        计算代理能量（当没有弛豫模型时使用）

        使用吸附物与表面的距离和位置作为代理能量指标。
        较低的能量表示更稳定的配置。
        """
        tags = adslab.get_tags()
        ads_mask = tags == 2
        surf_mask = tags == 1

        ads_positions = adslab.positions[ads_mask]
        surf_positions = adslab.positions[surf_mask]

        if len(ads_positions) == 0 or len(surf_positions) == 0:
            return adslab, 0.0, 0.0

        # 计算吸附物质心
        ads_com = ads_positions.mean(axis=0)

        # 计算到表面原子的最小距离
        distances = []
        for surf_pos in surf_positions:
            dist = np.linalg.norm(ads_com[:2] - surf_pos[:2])  # xy平面距离
            distances.append(dist)
        min_xy_dist = min(distances)

        # 计算吸附物到表面的高度
        surf_z_max = surf_positions[:, 2].max()
        ads_z_min = ads_positions[:, 2].min()
        height = ads_z_min - surf_z_max

        # 代理能量：距离表面原子越近、高度在合理范围内越好
        # 理想高度约1.5-2.5 Å
        ideal_height = 2.0
        height_penalty = (height - ideal_height) ** 2

        # 代理能量 (越低越好)
        proxy_energy = min_xy_dist * 0.5 + height_penalty * 0.3

        # 代理力 (使用一个固定的小值)
        proxy_fmax = 0.1

        return adslab, proxy_energy, proxy_fmax

    def _get_adsorbate_info(self, atoms: "ase.Atoms") -> Tuple[np.ndarray, np.ndarray]:
        """获取吸附物的位置和取向信息"""
        tags = atoms.get_tags()
        ads_mask = tags == 2
        ads_positions = atoms.positions[ads_mask]

        if len(ads_positions) == 0:
            return np.array([0, 0, 0]), np.array([0, 0, 1])

        # 计算质心
        com = ads_positions.mean(axis=0)

        # 计算主轴方向（简化）
        if len(ads_positions) > 1:
            # 使用第一个到最后一个原子的方向
            orientation = ads_positions[-1] - ads_positions[0]
            norm = np.linalg.norm(orientation)
            if norm > 1e-6:
                orientation = orientation / norm
            else:
                orientation = np.array([0, 0, 1])
        else:
            orientation = np.array([0, 0, 1])

        return com, orientation

    def predict(
        self,
        surface: Union[str, "ase.Atoms", "pymatgen.core.Structure"],
        adsorbate: Union[str, int, "ase.Atoms"] = None,
        has_adsorbate: bool = False,
        adsorbate_tag: int = 2,
        binding_indices: List[int] = None,
        num_samples: int = None,
        output_dir: Optional[str] = None,
        save_trajectory: bool = True,
    ) -> PredictionOutput:
        """
        预测最优吸附位点和取向

        Parameters
        ----------
        surface : str or ase.Atoms or pymatgen.Structure
            表面结构。可以是文件路径、ASE Atoms对象或pymatgen Structure
        adsorbate : str or int or ase.Atoms, optional
            吸附物。可以是:
            - SMILES字符串 (如 "*CO", "*O", "*OH")
            - 数据库ID (整数)
            - ASE Atoms对象
            - 文件路径
            如果has_adsorbate=True，则忽略此参数
        has_adsorbate : bool, default=False
            表面结构中是否已包含吸附物
        adsorbate_tag : int, default=2
            当has_adsorbate=True时，用于识别吸附物原子的tag值
        binding_indices : list of int, optional
            吸附物的结合原子索引列表
        num_samples : int, optional
            采样的配置数量。如果不指定，使用初始化时的num_sites
        output_dir : str, optional
            输出目录
        save_trajectory : bool, default=True
            是否保存扩散轨迹

        Returns
        -------
        PredictionOutput
            预测结果，包含最优配置和所有采样配置
        """
        if num_samples is not None:
            self.num_sites = num_samples

        # 设置输出目录
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "prediction")
        os.makedirs(output_dir, exist_ok=True)

        self._log(f"\n{'='*60}")
        self._log(f"AdsorbDiff Prediction")
        self._log(f"{'='*60}")

        # 转换表面结构
        surface_atoms = self._convert_surface(surface)
        self._log(f"Surface: {surface_atoms.get_chemical_formula()}")

        # 处理吸附物
        if has_adsorbate:
            # 从表面提取吸附物
            slab_atoms, adsorbate_atoms = self._extract_adsorbate_from_surface(
                surface_atoms, adsorbate_tag
            )
            self._log(f"  Extracted adsorbate from surface: {adsorbate_atoms.get_chemical_formula()}")
        else:
            if adsorbate is None:
                raise ValueError("adsorbate must be provided when has_adsorbate=False")
            slab_atoms = surface_atoms
            adsorbate_atoms = self._convert_adsorbate(adsorbate)
            self._log(f"Adsorbate: {adsorbate_atoms.get_chemical_formula()}")

        # 确保表面有正确的tags
        slab_atoms = self._tag_surface_atoms(slab_atoms)

        # 创建吸附配置
        self._log(f"\nGenerating {self.num_sites} adsorption configurations...")
        if has_adsorbate:
            # 使用已有的吸附物位置作为初始配置
            adslab_list = [surface_atoms.copy()]
            metadata_list = [{"site": self._get_adsorbate_info(surface_atoms)[0]}]
        else:
            adslab_list, metadata_list = self._create_adslab_config(
                slab_atoms, adsorbate_atoms, binding_indices
            )

        self._log(f"  Created {len(adslab_list)} configurations")

        # 运行扩散和弛豫
        results = []
        traj_dir = os.path.join(output_dir, "trajectories")
        os.makedirs(traj_dir, exist_ok=True)

        self._log(f"\nRunning diffusion and relaxation...")
        for i, (adslab, metadata) in enumerate(zip(adslab_list, metadata_list)):
            self._log(f"  Processing configuration {i+1}/{len(adslab_list)}...")

            initial_atoms = adslab.copy()

            # 运行扩散
            diff_traj_path = os.path.join(traj_dir, f"diff_{i}")
            os.makedirs(diff_traj_path, exist_ok=True)

            try:
                diffused_adslab = self._run_diffusion(adslab, diff_traj_path)
            except Exception as e:
                self._log(f"    Warning: Diffusion failed: {e}")
                diffused_adslab = adslab.copy()

            # 获取扩散后的吸附位置
            ads_com, ads_orient = self._get_adsorbate_info(diffused_adslab)

            # 运行弛豫
            relax_traj_path = os.path.join(traj_dir, f"relax_{i}.traj") if save_trajectory else None
            relaxed_adslab, energy, fmax = self._run_relaxation(
                diffused_adslab, relax_traj_path
            )

            # 获取最终位置
            final_com, final_orient = self._get_adsorbate_info(relaxed_adslab)

            result = AdsorptionResult(
                initial_structure=initial_atoms,
                final_structure=relaxed_adslab,
                adsorbate_position=final_com,
                adsorbate_orientation=final_orient,
                energy=energy,
                forces_max=fmax,
                trajectory_path=relax_traj_path,
                site_info=AdsorptionSite(
                    site_id=i,
                    position=metadata.get("site", np.array([0, 0, 0])),
                    binding_type=self.placement_mode,
                ),
            )
            results.append(result)

        # 找到最优配置
        best_result = min(results, key=lambda x: x.energy)

        output = PredictionOutput(
            surface_formula=slab_atoms.get_chemical_formula(),
            adsorbate_formula=adsorbate_atoms.get_chemical_formula(),
            num_sites_sampled=len(results),
            results=results,
            best_result=best_result,
            output_dir=output_dir,
            use_proxy_energy=(self.relax_checkpoint_path is None),
        )

        self._log(f"\n{output.summary()}")

        # 保存最优结构
        from ase.io import write
        best_path = os.path.join(output_dir, "best_config.vasp")
        write(best_path, best_result.final_structure, format="vasp")
        self._log(f"\nBest configuration saved to: {best_path}")

        return output

    def predict_single(
        self,
        adslab: "ase.Atoms",
        output_dir: Optional[str] = None,
    ) -> AdsorptionResult:
        """
        对单个吸附物-表面配置进行预测

        Parameters
        ----------
        adslab : ase.Atoms
            已经放置好吸附物的表面结构（吸附物tag=2）
        output_dir : str, optional
            输出目录

        Returns
        -------
        AdsorptionResult
            预测结果
        """
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "single_prediction")
        os.makedirs(output_dir, exist_ok=True)

        initial_atoms = adslab.copy()

        # 运行扩散
        traj_dir = os.path.join(output_dir, "trajectory")
        os.makedirs(traj_dir, exist_ok=True)

        diffused_adslab = self._run_diffusion(adslab, traj_dir)

        # 运行弛豫
        relax_traj_path = os.path.join(output_dir, "relax.traj")
        relaxed_adslab, energy, fmax = self._run_relaxation(
            diffused_adslab, relax_traj_path
        )

        # 获取最终位置
        final_com, final_orient = self._get_adsorbate_info(relaxed_adslab)

        return AdsorptionResult(
            initial_structure=initial_atoms,
            final_structure=relaxed_adslab,
            adsorbate_position=final_com,
            adsorbate_orientation=final_orient,
            energy=energy,
            forces_max=fmax,
            trajectory_path=relax_traj_path,
        )

    def list_available_adsorbates(self) -> List[Dict]:
        """列出数据库中所有可用的吸附物"""
        ads_db_path = os.path.join(
            self.adsorbdiff_root, "adsorbdiff/placement/pkls/adsorbates.pkl"
        )
        with open(ads_db_path, 'rb') as f:
            ads_db = pickle.load(f)

        adsorbates = []
        for idx, ads_info in ads_db.items():
            adsorbates.append({
                "id": idx,
                "formula": ads_info[0].get_chemical_formula(),
                "smiles": ads_info[1],
                "num_atoms": len(ads_info[0]),
            })
        return adsorbates

    def __del__(self):
        """清理临时目录"""
        if hasattr(self, '_temp_dir') and self._temp_dir and not self.keep_files:
            try:
                shutil.rmtree(self._temp_dir)
            except:
                pass


# 便捷函数
def predict_adsorption_site(
    surface: Union[str, "ase.Atoms"],
    adsorbate: Union[str, int, "ase.Atoms"],
    adsorbdiff_root: str,
    checkpoint_path: Optional[str] = None,
    use_gpu: bool = True,
    num_samples: int = 10,
    output_dir: Optional[str] = None,
    **kwargs,
) -> PredictionOutput:
    """
    便捷函数：预测最优吸附位点

    Parameters
    ----------
    surface : str or ase.Atoms
        表面结构
    adsorbate : str or int or ase.Atoms
        吸附物 (SMILES, 数据库ID, 或ASE Atoms)
    adsorbdiff_root : str
        AdsorbDiff项目根目录
    checkpoint_path : str, optional
        模型checkpoint路径
    use_gpu : bool, default=True
        是否使用GPU
    num_samples : int, default=10
        采样配置数量
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数传递给AdsorbDiffPredictor

    Returns
    -------
    PredictionOutput
        预测结果
    """
    predictor = AdsorbDiffPredictor(
        adsorbdiff_root=adsorbdiff_root,
        checkpoint_path=checkpoint_path,
        use_gpu=use_gpu,
        num_sites=num_samples,
        keep_files=True if output_dir else False,
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

    parser = argparse.ArgumentParser(description="AdsorbDiff Adsorption Site Predictor")
    parser.add_argument("surface", help="Path to surface structure file")
    parser.add_argument("--adsorbate", default="*CO", help="Adsorbate (SMILES or file path)")
    parser.add_argument("--adsorbdiff-root", default="deps/AdsorbDiff",
                       help="Path to AdsorbDiff root directory")
    parser.add_argument("--checkpoint", default=None, help="Path to model checkpoint")
    parser.add_argument("--num-samples", type=int, default=5, help="Number of sites to sample")
    parser.add_argument("--cpu", action="store_true", help="Use CPU instead of GPU")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--has-adsorbate", action="store_true",
                       help="Surface already contains adsorbate")

    args = parser.parse_args()

    predictor = AdsorbDiffPredictor(
        adsorbdiff_root=args.adsorbdiff_root,
        checkpoint_path=args.checkpoint,
        use_gpu=not args.cpu,
        num_sites=args.num_samples,
        keep_files=True,
    )

    result = predictor.predict(
        surface=args.surface,
        adsorbate=args.adsorbate if not args.has_adsorbate else None,
        has_adsorbate=args.has_adsorbate,
        output_dir=args.output_dir,
    )

    print("\n" + "="*60)
    print("PREDICTION COMPLETE")
    print("="*60)
    print(result.summary())
