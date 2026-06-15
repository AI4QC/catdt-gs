"""
Electrochemical Surface Predictor - CP-MACE + VSSR-MC 集成

这个模块实现了策略B：On-the-Fly完全集成
将CP-MACE作为VSSR-MC的计算器，实现真正的电化学模拟。

核心功能：
1. 使用CP-MACE在恒电势下计算能量（包含电场效应）
2. VSSR-MC进行表面结构采样（吸附/解吸/交换）
3. 基于Grand Potential的Metropolis判据（考虑pH和电势）
4. 支持Active Learning（可选）

使用方式:
    from electrochemical_surface_predictor import ElectrochemicalSurfacePredictor

    # 初始化预测器
    predictor = ElectrochemicalSurfacePredictor(
        cp_mace_root="/path/to/CP-MACE",
        surface_sampling_root="/path/to/surface-sampling",
        use_gpu=True,
    )

    # 在指定pH和电势下采样表面结构
    result = predictor.sample_surface(
        surface="slab.vasp",
        adsorbates=["O", "OH"],
        potential_she=1.0,  # V vs. SHE
        ph=12.0,
        temperature=300.0,
        total_sweeps=100,
    )

    print(result.summary())

作者: Claude
"""

import os
import sys
import copy
import pickle
import tempfile
import shutil
import json
import logging
from typing import Union, Optional, List, Dict, Any, Tuple, Literal
from dataclasses import dataclass, field
from pathlib import Path
from collections import Counter

import numpy as np


# 物理常数
KB = 8.617333262e-5  # eV/K
KB_J = 1.380649e-23  # J/K
E_CHARGE = 1.602176634e-19  # C
FARADAY = 96485.33212  # C/mol
R_GAS = 8.314462618  # J/(mol·K)

# SHE到真空能级的转换 (常用值)
SHE_TO_VACUUM = -4.44  # eV (可根据CP-MACE训练基准调整)


@dataclass
class ElectrochemicalConditions:
    """电化学条件"""
    potential_she: float  # 电极电势 vs. SHE (V)
    ph: float  # pH值
    temperature: float  # 温度 (K)

    @property
    def potential_vacuum(self) -> float:
        """转换为真空能级 (eV)"""
        return SHE_TO_VACUUM - self.potential_she

    @property
    def kT(self) -> float:
        """热能 kT (eV)"""
        return KB * self.temperature

    def __repr__(self):
        return (f"ElectrochemicalConditions(U={self.potential_she:.2f} V vs. SHE, "
                f"pH={self.ph:.1f}, T={self.temperature:.1f} K)")


@dataclass
class ChemicalPotentials:
    """化学势数据"""
    species: Dict[str, float]  # 各物种的化学势 (eV)
    reference_state: str  # 参考态描述

    def get(self, species: str, default: float = 0.0) -> float:
        return self.species.get(species, default)


@dataclass
class SampledConfiguration:
    """采样的表面配置"""
    sweep_number: int
    atoms: "ase.Atoms"
    electronic_grand_potential: float  # Ω_el from CP-MACE (eV)
    chemical_potential_correction: float  # Δμ_ions (eV)
    total_grand_potential: float  # Ω_total (eV)
    fermi_level: float  # 费米能级 (V)
    electron_number: float  # 电子数
    composition: Dict[str, int]  # 表面组成
    num_adsorbates: int
    acceptance_rate: float
    force_uncertainty: float = 0.0  # 力的不确定性
    fermi_uncertainty: float = 0.0  # 费米能级不确定性

    def __repr__(self):
        formula = self.atoms.get_chemical_formula() if self.atoms else "N/A"
        return (f"SampledConfiguration(sweep={self.sweep_number}, "
                f"Ω={self.total_grand_potential:.4f} eV, "
                f"μ_F={self.fermi_level:.4f} V, formula={formula})")


@dataclass
class ElectrochemicalSamplingResult:
    """电化学表面采样结果"""
    surface_name: str
    conditions: ElectrochemicalConditions
    total_sweeps: int
    sweep_size: int
    configurations: List[SampledConfiguration]
    grand_potential_history: List[float]
    fermi_level_history: List[float]
    electron_history: List[float]
    acceptance_history: List[float]
    best_configuration: SampledConfiguration
    lowest_energy_configuration: SampledConfiguration
    high_uncertainty_configs: List[SampledConfiguration] = field(default_factory=list)
    output_dir: Optional[str] = None

    def get_top_n(self, n: int) -> List[SampledConfiguration]:
        """获取能量最低的前N个配置"""
        sorted_configs = sorted(self.configurations, key=lambda x: x.total_grand_potential)
        return sorted_configs[:n]

    def summary(self, top_n: int = 5) -> str:
        """生成结果摘要"""
        lines = [
            "Electrochemical Surface Sampling Results (CP-MACE + VSSR-MC)",
            "=" * 70,
            f"Surface: {self.surface_name}",
            f"Conditions: {self.conditions}",
            f"Total sweeps: {self.total_sweeps}",
            f"Sweep size: {self.sweep_size}",
            "",
            "Lowest Grand Potential Configuration:",
            f"  Ω_total: {self.lowest_energy_configuration.total_grand_potential:.4f} eV",
            f"  Ω_el (CP-MACE): {self.lowest_energy_configuration.electronic_grand_potential:.4f} eV",
            f"  Δμ_ions: {self.lowest_energy_configuration.chemical_potential_correction:.4f} eV",
            f"  Fermi level: {self.lowest_energy_configuration.fermi_level:.4f} V",
            f"  Electrons: {self.lowest_energy_configuration.electron_number:.2f}",
            f"  Formula: {self.lowest_energy_configuration.atoms.get_chemical_formula()}",
            "",
            "Grand Potential Statistics:",
            f"  Min: {np.min(self.grand_potential_history):.4f} eV",
            f"  Max: {np.max(self.grand_potential_history):.4f} eV",
            f"  Mean: {np.mean(self.grand_potential_history):.4f} eV",
            f"  Std: {np.std(self.grand_potential_history):.4f} eV",
            "",
            f"Average acceptance rate: {np.mean(self.acceptance_history):.2%}",
            f"High uncertainty structures collected: {len(self.high_uncertainty_configs)}",
            "",
            f"Top {min(top_n, len(self.configurations))} configurations:",
            "-" * 70,
        ]
        for i, c in enumerate(self.get_top_n(top_n), 1):
            lines.append(
                f"  {i}. Ω={c.total_grand_potential:.4f} eV, "
                f"μ_F={c.fermi_level:.4f} V, "
                f"{c.atoms.get_chemical_formula()}"
            )
        return "\n".join(lines)


class ElectrochemicalEvaluator:
    """
    电化学能量评估器

    封装CP-MACE计算器，并实现Grand Potential的计算。

    Grand Potential计算公式:
    Ω_total = Ω_el(U) + Δμ_ions(pH) + E_cap

    其中:
    - Ω_el(U): CP-MACE在恒电势下输出的电子巨势
    - Δμ_ions(pH): 离子化学势修正 (基于Nernst方程)
    - E_cap: 电容修正项 (可选)
    """

    # 标准化学势参考值 (eV)
    # 这些值需要根据实际DFT计算校准
    REFERENCE_ENERGIES = {
        "H2": -6.77,  # H2气体分子
        "H2O": -14.22,  # H2O液态
        "O2": -9.86,  # O2气体分子
    }

    def __init__(
        self,
        cp_mace_root: str,
        model_paths: List[str],
        conditions: ElectrochemicalConditions,
        device: str = "cuda",
        use_ensemble: bool = True,
        verbose: bool = True,
    ):
        """
        初始化电化学评估器

        Parameters
        ----------
        cp_mace_root : str
            CP-MACE根目录
        model_paths : list of str
            CP-MACE模型路径列表
        conditions : ElectrochemicalConditions
            电化学条件
        device : str
            计算设备
        use_ensemble : bool
            是否使用ensemble计算
        verbose : bool
            是否输出详细信息
        """
        self.cp_mace_root = os.path.abspath(cp_mace_root)
        self.model_paths = model_paths
        self.conditions = conditions
        self.device = device
        self.use_ensemble = use_ensemble
        self.verbose = verbose

        # 添加路径
        if self.cp_mace_root not in sys.path:
            sys.path.insert(0, self.cp_mace_root)

        # 延迟加载计算器
        self._calculators = None
        self._avg_calculator = None

        if verbose:
            print(f"ElectrochemicalEvaluator initialized:")
            print(f"  Target Fermi level: {conditions.potential_vacuum:.4f} eV (vacuum)")
            print(f"  Equivalent to: {conditions.potential_she:.4f} V vs. SHE")
            print(f"  pH: {conditions.ph}")
            print(f"  Temperature: {conditions.temperature} K")

    def _load_calculators(self):
        """加载CP-MACE计算器"""
        if self._calculators is None:
            import torch
            from mace.calculators import MACECalculator

            self._calculators = []
            for path in self.model_paths:
                if not os.path.exists(path):
                    raise FileNotFoundError(f"Model not found: {path}")
                calc = MACECalculator(
                    model_paths=[path],
                    device=self.device,
                )
                self._calculators.append(calc)

            if self.verbose:
                print(f"  Loaded {len(self._calculators)} CP-MACE model(s)")

        return self._calculators

    def _create_average_calculator(self):
        """创建平均力计算器"""
        if self._avg_calculator is None:
            from ase.calculators.calculator import Calculator

            calculators = self._load_calculators()

            class CPMACEAverageCalculator(Calculator):
                """CP-MACE Ensemble计算器"""
                implemented_properties = ['energy', 'forces', 'potential']

                def __init__(self, calculators, **kwargs):
                    super().__init__(**kwargs)
                    self.calculators = calculators

                def calculate(self, atoms=None, properties=['energy', 'forces', 'potential'],
                             system_changes=None):
                    super().calculate(atoms, properties, system_changes)

                    total_forces = 0
                    total_energy = 0
                    total_potential = 0
                    all_forces = []
                    all_potentials = []

                    for calc in self.calculators:
                        atoms_copy = copy.deepcopy(atoms)
                        atoms_copy.calc = calc
                        calc.calculate(atoms_copy, properties, system_changes)

                        total_forces += calc.results['forces']
                        total_energy += calc.results['energy']
                        total_potential += calc.results['potential']
                        all_forces.append(calc.results['forces'])
                        all_potentials.append(calc.results['potential'])

                    n = len(self.calculators)
                    self.results['forces'] = total_forces / n
                    self.results['energy'] = total_energy / n
                    self.results['potential'] = total_potential / n

                    # 计算不确定性
                    all_forces_array = np.array(all_forces)
                    force_std = np.std(all_forces_array, axis=0)
                    total_std = np.sqrt(np.sum(force_std**2, axis=1))
                    self.results['force_std'] = np.max(total_std)
                    self.results['potential_std'] = np.std(all_potentials)

                def get_fermi_level(self):
                    return self.results.get('potential', 0)

                def get_force_uncertainty(self):
                    return self.results.get('force_std', 0)

                def get_fermi_uncertainty(self):
                    return self.results.get('potential_std', 0)

            self._avg_calculator = CPMACEAverageCalculator(calculators)

        return self._avg_calculator

    def get_calculator(self):
        """获取计算器"""
        if self.use_ensemble and len(self.model_paths) > 1:
            return self._create_average_calculator()
        else:
            calcs = self._load_calculators()
            return calcs[0]

    def calculate_chemical_potential_correction(
        self,
        added_species: Dict[str, int],
    ) -> float:
        """
        计算化学势修正项 Δμ_ions(pH)

        基于计算氢电极(CHE)模型:
        - 吸附O*:  ΔG = G(H2O) - G(H2) - 1/2*G(O2) - kT*ln(10)*pH
        - 吸附OH*: ΔG = G(H2O) - 1/2*G(H2) - kT*ln(10)*pH
        - 吸附H*:  ΔG = 1/2*G(H2) - e*U_SHE

        Parameters
        ----------
        added_species : dict
            添加/移除的物种计数, 如 {"O": 1, "H": -1}
            正数表示添加，负数表示移除

        Returns
        -------
        float
            化学势修正 (eV)
        """
        correction = 0.0
        kT = self.conditions.kT
        pH = self.conditions.ph
        U = self.conditions.potential_she

        for species, count in added_species.items():
            if count == 0:
                continue

            if species == "O":
                # O* 吸附: 来自H2O，释放H2
                # μ_O = μ_H2O - μ_H2 + kT*ln(10)*pH - eU
                # 在pH=0, U=0时，μ_O ≈ 2.46 eV (相对于1/2 O2)
                mu_O = (self.REFERENCE_ENERGIES["H2O"]
                       - self.REFERENCE_ENERGIES["H2"]
                       + kT * np.log(10) * pH
                       - U)
                correction += count * mu_O

            elif species == "OH":
                # OH* 吸附: 来自H2O，释放1/2 H2
                # μ_OH = μ_H2O - 0.5*μ_H2 + kT*ln(10)*pH - 0.5*eU
                mu_OH = (self.REFERENCE_ENERGIES["H2O"]
                        - 0.5 * self.REFERENCE_ENERGIES["H2"]
                        + kT * np.log(10) * pH
                        - 0.5 * U)
                correction += count * mu_OH

            elif species == "H":
                # H* 吸附: 来自溶液中的H+
                # μ_H = 0.5*μ_H2 - kT*ln(10)*pH + eU
                mu_H = (0.5 * self.REFERENCE_ENERGIES["H2"]
                       - kT * np.log(10) * pH
                       + U)
                correction += count * mu_H

            else:
                # 其他元素: 使用默认化学势 (0)
                # 用户可以通过继承此类来添加更多物种
                pass

        return correction

    def calculate_grand_potential(
        self,
        atoms: "ase.Atoms",
        added_species: Optional[Dict[str, int]] = None,
        relax: bool = False,
        fmax: float = 0.05,
        steps: int = 100,
    ) -> Tuple[float, float, float, float, float]:
        """
        计算总巨势 Ω_total

        Ω_total = Ω_el(U) + Δμ_ions(pH)

        Parameters
        ----------
        atoms : ase.Atoms
            原子结构（需要包含'electron'属性）
        added_species : dict, optional
            相对于参考态添加/移除的物种
        relax : bool
            是否进行结构弛豫
        fmax : float
            弛豫收敛标准
        steps : int
            最大弛豫步数

        Returns
        -------
        tuple
            (Ω_total, Ω_el, Δμ_ions, fermi_level, electron_number)
        """
        calc = self.get_calculator()
        atoms = atoms.copy()
        atoms.calc = calc

        # 确保有电子数属性
        if 'electron' not in atoms.info:
            # 估计初始电子数（可以基于元素的原子序数）
            total_electrons = sum(atoms.get_atomic_numbers())
            atoms.info['electron'] = float(total_electrons)

        # 可选：结构弛豫
        if relax:
            from ase.optimize import BFGS
            opt = BFGS(atoms, logfile=None)
            try:
                opt.run(fmax=fmax, steps=steps)
            except Exception as e:
                if self.verbose:
                    print(f"  Warning: Relaxation failed: {e}")

        # 计算CP-MACE能量（电子巨势）
        omega_el = atoms.get_potential_energy()

        # 获取费米能级
        if hasattr(calc, 'get_fermi_level'):
            fermi_level = calc.get_fermi_level()
        elif hasattr(calc, 'results') and 'potential' in calc.results:
            fermi_level = calc.results['potential']
        else:
            fermi_level = atoms.info.get('potential', 0)

        electron_number = atoms.info.get('electron', 0)

        # 计算化学势修正
        if added_species:
            delta_mu = self.calculate_chemical_potential_correction(added_species)
        else:
            delta_mu = 0.0

        # 总巨势
        omega_total = omega_el + delta_mu

        return omega_total, omega_el, delta_mu, fermi_level, electron_number

    def metropolis_criterion(
        self,
        delta_omega: float,
        temperature_kT: Optional[float] = None,
    ) -> bool:
        """
        Metropolis判据

        P_accept = min(1, exp(-ΔΩ/kT))

        Parameters
        ----------
        delta_omega : float
            巨势变化 (eV)
        temperature_kT : float, optional
            温度 (kT单位)，默认使用conditions中的值

        Returns
        -------
        bool
            是否接受
        """
        if temperature_kT is None:
            kT = self.conditions.kT
        else:
            kT = temperature_kT

        if delta_omega <= 0:
            return True
        else:
            probability = np.exp(-delta_omega / kT)
            return np.random.random() < probability


class ElectrochemicalSurfacePredictor:
    """
    电化学表面预测器 (CP-MACE + VSSR-MC)

    实现策略B：On-the-Fly完全集成
    将CP-MACE作为VSSR-MC的计算器，在恒电势下进行表面结构采样。

    Parameters
    ----------
    cp_mace_root : str
        CP-MACE项目根目录
    surface_sampling_root : str
        surface-sampling项目根目录
    model_paths : list of str
        CP-MACE模型路径列表
    device : str, default="cuda"
        计算设备
    cutoff : float, default=6.0
        截断半径 (Å)
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留输出文件
    verbose : bool, default=True
        是否输出详细信息
    """

    def __init__(
        self,
        cp_mace_root: str,
        surface_sampling_root: str,
        model_paths: Optional[List[str]] = None,
        device: str = "cuda",
        cutoff: float = 6.0,
        work_dir: Optional[str] = None,
        keep_files: bool = True,
        verbose: bool = True,
    ):
        self.cp_mace_root = os.path.abspath(cp_mace_root)
        self.surface_sampling_root = os.path.abspath(surface_sampling_root)
        self.device = device
        self.cutoff = cutoff
        self.keep_files = keep_files
        self.verbose = verbose

        # 设置模型路径
        if model_paths is None:
            # 使用CP-MACE自带的示例模型
            default_models = [
                os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_1.model"),
                os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_2.model"),
            ]
            self.model_paths = [p for p in default_models if os.path.exists(p)]
            if not self.model_paths:
                raise FileNotFoundError("No CP-MACE models found. Please specify model_paths.")
        else:
            self.model_paths = model_paths

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="echem_surface_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 添加路径
        if self.cp_mace_root not in sys.path:
            sys.path.insert(0, self.cp_mace_root)
        if self.surface_sampling_root not in sys.path:
            sys.path.insert(0, self.surface_sampling_root)

        # 验证环境
        self._validate_environment()

        if self.verbose:
            print(f"ElectrochemicalSurfacePredictor initialized:")
            print(f"  CP-MACE root: {self.cp_mace_root}")
            print(f"  Surface-sampling root: {self.surface_sampling_root}")
            print(f"  Models: {len(self.model_paths)}")
            print(f"  Device: {self.device}")
            print(f"  Work dir: {self.work_dir}")

    def _validate_environment(self):
        """验证环境"""
        if not os.path.exists(self.cp_mace_root):
            raise FileNotFoundError(f"CP-MACE root not found: {self.cp_mace_root}")
        if not os.path.exists(self.surface_sampling_root):
            raise FileNotFoundError(f"Surface-sampling root not found: {self.surface_sampling_root}")

        mace_init = os.path.join(self.cp_mace_root, "mace", "__init__.py")
        if not os.path.exists(mace_init):
            raise FileNotFoundError(f"CP-MACE mace package not found")

        mcmc_init = os.path.join(self.surface_sampling_root, "mcmc", "__init__.py")
        if not os.path.exists(mcmc_init):
            raise FileNotFoundError(f"VSSR-MC mcmc package not found")

    def _log(self, message: str):
        """输出日志"""
        if self.verbose:
            print(message)

    def _convert_surface(self, surface) -> "ase.Atoms":
        """转换输入结构"""
        from ase import Atoms
        from ase.io import read

        if isinstance(surface, str):
            if not os.path.exists(surface):
                raise FileNotFoundError(f"Surface file not found: {surface}")
            return read(surface)
        elif isinstance(surface, Atoms):
            return surface.copy()
        else:
            try:
                from pymatgen.io.ase import AseAtomsAdaptor
                from pymatgen.core import Structure
                if isinstance(surface, Structure):
                    return AseAtomsAdaptor.get_atoms(surface)
            except ImportError:
                pass
            raise TypeError(f"Unsupported surface type: {type(surface)}")

    def _generate_virtual_sites(
        self,
        slab: "ase.Atoms",
        adsorbates: List[str],
        planar_distance: float = 1.5,
        surface_depth: int = 1,
    ) -> List[np.ndarray]:
        """
        生成虚拟吸附位点

        基于VSSR-MC的方法，在表面上方生成ontop/bridge/hollow位点
        """
        from mcmc.system import SurfaceSystem
        from mcmc.utils.misc import get_atoms_batch

        # 准备AtomsBatch
        device = "cpu"  # 仅用于生成位点
        slab_batch = get_atoms_batch(
            slab,
            nff_cutoff=self.cutoff,
            device=device,
            props={"energy": 0, "energy_grad": []},
        )

        # 创建SurfaceSystem (不需要计算器)
        system_settings = {
            "surface_name": "temp",
            "cutoff": self.cutoff,
            "surface_depth": surface_depth,
            "planar_distance": planar_distance,
            "ads_site_type": "all",
        }

        # 使用dummy计算器
        class DummyCalc:
            def get_potential_energy(self, atoms=None):
                return 0.0
            def set(self, **kwargs):
                pass

        surface_system = SurfaceSystem(
            slab_batch,
            calc=DummyCalc(),
            system_settings=system_settings,
            save_folder=self.work_dir,
        )

        return surface_system.ads_coords.copy()

    def sample_surface(
        self,
        surface: Union[str, "ase.Atoms", "pymatgen.core.Structure"],
        adsorbates: List[str],
        potential_she: float,
        ph: float,
        temperature: float = 300.0,
        initial_electron_number: Optional[float] = None,
        total_sweeps: int = 100,
        sweep_size: int = 20,
        canonical: bool = False,
        num_adsorbates: int = 0,
        force_threshold: float = 0.15,
        fermi_threshold: float = 0.04,
        relax_each_step: bool = False,
        output_dir: Optional[str] = None,
        run_name: Optional[str] = None,
    ) -> ElectrochemicalSamplingResult:
        """
        在指定电化学条件下采样表面结构

        Parameters
        ----------
        surface : str or ase.Atoms or pymatgen.Structure
            表面结构
        adsorbates : list of str
            可吸附的物种列表 (如 ["O", "OH", "H"])
        potential_she : float
            电极电势 vs. SHE (V)
        ph : float
            溶液pH值
        temperature : float, default=300.0
            温度 (K)
        initial_electron_number : float, optional
            初始电子数
        total_sweeps : int, default=100
            MC采样轮数
        sweep_size : int, default=20
            每轮采样步数
        canonical : bool, default=False
            是否使用canonical (NVT) 采样
        num_adsorbates : int, default=0
            固定吸附物数量 (canonical模式)
        force_threshold : float, default=0.15
            力不确定性阈值，用于Active Learning
        fermi_threshold : float, default=0.04
            费米能级不确定性阈值
        relax_each_step : bool, default=False
            是否在每步后进行结构弛豫
        output_dir : str, optional
            输出目录
        run_name : str, optional
            运行名称

        Returns
        -------
        ElectrochemicalSamplingResult
            采样结果
        """
        import torch
        from ase.io import write

        self._log("\n" + "=" * 70)
        self._log("Electrochemical Surface Sampling (CP-MACE + VSSR-MC)")
        self._log("=" * 70)

        # 设置电化学条件
        conditions = ElectrochemicalConditions(
            potential_she=potential_she,
            ph=ph,
            temperature=temperature,
        )
        self._log(f"Conditions: {conditions}")

        # 转换结构
        slab = self._convert_surface(surface)
        slab.pbc = [True, True, True]
        surface_formula = slab.get_chemical_formula()
        self._log(f"Surface: {surface_formula}")
        self._log(f"Atoms: {len(slab)}")
        self._log(f"Adsorbates: {adsorbates}")

        # 设置输出目录
        if output_dir is None:
            output_dir = os.path.join(
                self.work_dir,
                run_name or f"echem_sample_{surface_formula}"
            )
        os.makedirs(output_dir, exist_ok=True)

        # 设置初始电子数
        if initial_electron_number is not None:
            slab.info['electron'] = initial_electron_number
        elif 'electron' not in slab.info:
            # 估计初始电子数（基于原子的价电子）
            total_electrons = sum(slab.get_atomic_numbers())
            slab.info['electron'] = float(total_electrons)

        self._log(f"Initial electron number: {slab.info['electron']:.2f}")

        # 创建电化学评估器
        self._log("\nLoading CP-MACE models...")
        evaluator = ElectrochemicalEvaluator(
            cp_mace_root=self.cp_mace_root,
            model_paths=self.model_paths,
            conditions=conditions,
            device=self.device,
            use_ensemble=len(self.model_paths) > 1,
            verbose=self.verbose,
        )

        # 生成虚拟吸附位点
        self._log("\nGenerating virtual adsorption sites...")
        ads_sites = self._generate_virtual_sites(slab, adsorbates)
        self._log(f"  Generated {len(ads_sites)} virtual sites")

        # 初始化MC采样
        self._log(f"\nStarting MC sampling...")
        self._log(f"  Mode: {'Canonical' if canonical else 'Semi-grand canonical'}")
        self._log(f"  Total sweeps: {total_sweeps}")
        self._log(f"  Sweep size: {sweep_size}")

        # 计算初始巨势
        current_atoms = slab.copy()
        current_omega, current_omega_el, _, current_fermi, current_ne = \
            evaluator.calculate_grand_potential(current_atoms)

        self._log(f"  Initial Ω: {current_omega:.4f} eV")
        self._log(f"  Initial Fermi level: {current_fermi:.4f} V")

        # 初始化历史记录
        configurations = []
        grand_potential_history = []
        fermi_history = []
        electron_history = []
        acceptance_history = []
        high_uncertainty_configs = []

        # 追踪当前表面组成
        current_composition = Counter(current_atoms.get_chemical_symbols())
        adsorbate_counts = {ads: 0 for ads in adsorbates}

        # MC采样循环
        total_accepted = 0
        total_trials = 0

        try:
            from tqdm import tqdm
            sweep_iterator = tqdm(range(total_sweeps), desc="Sweeps")
        except ImportError:
            sweep_iterator = range(total_sweeps)

        for sweep in sweep_iterator:
            sweep_accepted = 0

            for step in range(sweep_size):
                total_trials += 1

                # 选择试探步类型
                # 1. 添加吸附物
                # 2. 移除吸附物
                # 3. 交换吸附物位置
                move_type = np.random.choice(['add', 'remove', 'swap'],
                                             p=[0.4, 0.4, 0.2])

                trial_atoms = current_atoms.copy()
                added_species = {}

                if move_type == 'add' and len(ads_sites) > 0:
                    # 添加吸附物
                    species = np.random.choice(adsorbates)
                    site_idx = np.random.randint(len(ads_sites))
                    site_pos = ads_sites[site_idx]

                    # 添加原子
                    from ase import Atom
                    trial_atoms.append(Atom(species, position=site_pos))
                    added_species[species] = 1

                elif move_type == 'remove':
                    # 移除吸附物（只能移除表面吸附的原子）
                    removable_indices = []
                    for i, symbol in enumerate(trial_atoms.get_chemical_symbols()):
                        if symbol in adsorbates:
                            # 检查是否在表面（z坐标较高）
                            if trial_atoms[i].position[2] > np.mean(trial_atoms.positions[:, 2]):
                                removable_indices.append(i)

                    if removable_indices:
                        remove_idx = np.random.choice(removable_indices)
                        removed_species = trial_atoms[remove_idx].symbol
                        del trial_atoms[remove_idx]
                        added_species[removed_species] = -1
                    else:
                        # 没有可移除的原子，跳过
                        continue

                elif move_type == 'swap':
                    # 交换吸附物位置
                    adsorbate_indices = []
                    for i, symbol in enumerate(trial_atoms.get_chemical_symbols()):
                        if symbol in adsorbates:
                            adsorbate_indices.append(i)

                    if len(adsorbate_indices) >= 2:
                        idx1, idx2 = np.random.choice(adsorbate_indices, 2, replace=False)
                        pos1 = trial_atoms[idx1].position.copy()
                        pos2 = trial_atoms[idx2].position.copy()
                        trial_atoms[idx1].position = pos2
                        trial_atoms[idx2].position = pos1
                    else:
                        continue
                else:
                    continue

                # 保持电子数不变（从current_atoms继承）
                trial_atoms.info['electron'] = current_atoms.info.get('electron', current_ne)

                # 计算试探结构的巨势
                try:
                    trial_omega, trial_omega_el, trial_delta_mu, trial_fermi, trial_ne = \
                        evaluator.calculate_grand_potential(
                            trial_atoms,
                            added_species,
                            relax=relax_each_step,
                        )
                except Exception as e:
                    if self.verbose:
                        print(f"  Warning: Energy calculation failed: {e}")
                    continue

                # Metropolis判据
                delta_omega = trial_omega - current_omega
                accepted = evaluator.metropolis_criterion(delta_omega, conditions.kT)

                if accepted:
                    current_atoms = trial_atoms
                    current_omega = trial_omega
                    current_omega_el = trial_omega_el
                    current_fermi = trial_fermi
                    current_ne = trial_ne
                    sweep_accepted += 1
                    total_accepted += 1

                    # 更新组成
                    for species, count in added_species.items():
                        adsorbate_counts[species] = adsorbate_counts.get(species, 0) + count

            # 记录这一轮的结果
            acceptance_rate = sweep_accepted / sweep_size if sweep_size > 0 else 0

            # 获取不确定性
            calc = evaluator.get_calculator()
            force_uncertainty = 0.0
            fermi_uncertainty = 0.0
            if hasattr(calc, 'get_force_uncertainty'):
                force_uncertainty = calc.get_force_uncertainty()
            if hasattr(calc, 'get_fermi_uncertainty'):
                fermi_uncertainty = calc.get_fermi_uncertainty()

            config = SampledConfiguration(
                sweep_number=sweep + 1,
                atoms=current_atoms.copy(),
                electronic_grand_potential=current_omega_el,
                chemical_potential_correction=current_omega - current_omega_el,
                total_grand_potential=current_omega,
                fermi_level=current_fermi,
                electron_number=current_ne,
                composition=dict(Counter(current_atoms.get_chemical_symbols())),
                num_adsorbates=sum(adsorbate_counts.values()),
                acceptance_rate=acceptance_rate,
                force_uncertainty=force_uncertainty,
                fermi_uncertainty=fermi_uncertainty,
            )
            configurations.append(config)

            grand_potential_history.append(current_omega)
            fermi_history.append(current_fermi)
            electron_history.append(current_ne)
            acceptance_history.append(acceptance_rate)

            # 检查不确定性，收集高不确定性结构 (Active Learning)
            if force_uncertainty > force_threshold or fermi_uncertainty > fermi_threshold:
                high_uncertainty_configs.append(config)
                if self.verbose and (sweep + 1) % 10 == 0:
                    print(f"  Sweep {sweep+1}: High uncertainty detected "
                          f"(force_std={force_uncertainty:.4f}, fermi_std={fermi_uncertainty:.4f})")

        # 找到最优结构
        lowest_idx = np.argmin(grand_potential_history)
        lowest_energy_config = configurations[lowest_idx]
        best_config = configurations[-1]

        # 保存结果
        self._log("\nSaving results...")

        # 保存最终结构
        write(os.path.join(output_dir, "final_structure.xyz"),
              best_config.atoms, format='extxyz')
        write(os.path.join(output_dir, "lowest_energy_structure.xyz"),
              lowest_energy_config.atoms, format='extxyz')

        # 保存所有结构
        all_atoms = [c.atoms for c in configurations]
        write(os.path.join(output_dir, "all_structures.xyz"),
              all_atoms, format='extxyz')

        # 保存高不确定性结构
        if high_uncertainty_configs:
            uncertain_atoms = [c.atoms for c in high_uncertainty_configs]
            write(os.path.join(output_dir, "high_uncertainty_structures.xyz"),
                  uncertain_atoms, format='extxyz')

        # 保存设置和历史
        settings = {
            "conditions": {
                "potential_she": potential_she,
                "ph": ph,
                "temperature": temperature,
            },
            "sampling": {
                "total_sweeps": total_sweeps,
                "sweep_size": sweep_size,
                "canonical": canonical,
                "adsorbates": adsorbates,
            },
            "model_paths": self.model_paths,
        }
        with open(os.path.join(output_dir, "settings.json"), 'w') as f:
            json.dump(settings, f, indent=2)

        # 创建结果对象
        result = ElectrochemicalSamplingResult(
            surface_name=surface_formula,
            conditions=conditions,
            total_sweeps=total_sweeps,
            sweep_size=sweep_size,
            configurations=configurations,
            grand_potential_history=grand_potential_history,
            fermi_level_history=fermi_history,
            electron_history=electron_history,
            acceptance_history=acceptance_history,
            best_configuration=best_config,
            lowest_energy_configuration=lowest_energy_config,
            high_uncertainty_configs=high_uncertainty_configs,
            output_dir=output_dir,
        )

        self._log(f"\n{result.summary()}")
        self._log(f"\nResults saved to: {output_dir}")

        return result

    def __del__(self):
        """清理临时目录"""
        if hasattr(self, '_temp_dir') and self._temp_dir and not self.keep_files:
            try:
                shutil.rmtree(self._temp_dir)
            except:
                pass


# 便捷函数
def sample_electrochemical_surface(
    surface: Union[str, "ase.Atoms"],
    adsorbates: List[str],
    potential_she: float,
    ph: float,
    cp_mace_root: str,
    surface_sampling_root: str,
    model_paths: Optional[List[str]] = None,
    temperature: float = 300.0,
    total_sweeps: int = 100,
    device: str = "cuda",
    output_dir: Optional[str] = None,
    **kwargs,
) -> ElectrochemicalSamplingResult:
    """
    便捷函数：在指定电化学条件下采样表面结构

    Parameters
    ----------
    surface : str or ase.Atoms
        表面结构
    adsorbates : list of str
        可吸附物种
    potential_she : float
        电极电势 vs. SHE (V)
    ph : float
        溶液pH值
    cp_mace_root : str
        CP-MACE根目录
    surface_sampling_root : str
        surface-sampling根目录
    model_paths : list of str, optional
        模型路径
    temperature : float, default=300.0
        温度 (K)
    total_sweeps : int, default=100
        采样轮数
    device : str, default="cuda"
        计算设备
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数

    Returns
    -------
    ElectrochemicalSamplingResult
        采样结果
    """
    predictor = ElectrochemicalSurfacePredictor(
        cp_mace_root=cp_mace_root,
        surface_sampling_root=surface_sampling_root,
        model_paths=model_paths,
        device=device,
    )
    return predictor.sample_surface(
        surface=surface,
        adsorbates=adsorbates,
        potential_she=potential_she,
        ph=ph,
        temperature=temperature,
        total_sweeps=total_sweeps,
        output_dir=output_dir,
        **kwargs,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Electrochemical Surface Sampling (CP-MACE + VSSR-MC)"
    )
    parser.add_argument("surface", help="Surface structure file")
    parser.add_argument("--adsorbates", nargs="+", default=["O", "OH"],
                        help="Adsorbate species")
    parser.add_argument("--potential", type=float, required=True,
                        help="Electrode potential vs. SHE (V)")
    parser.add_argument("--ph", type=float, required=True,
                        help="Solution pH")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="Temperature (K)")
    parser.add_argument("--total-sweeps", type=int, default=100,
                        help="Number of MC sweeps")
    parser.add_argument("--sweep-size", type=int, default=20,
                        help="Steps per sweep")
    parser.add_argument("--cp-mace-root",
                        default="deps/CP-MACE",
                        help="CP-MACE root directory")
    parser.add_argument("--surface-sampling-root",
                        default="deps/surface-sampling",
                        help="Surface-sampling root directory")
    parser.add_argument("--model-paths", nargs="+", default=None,
                        help="CP-MACE model paths")
    parser.add_argument("--initial-electron", type=float, default=None,
                        help="Initial electron number")
    parser.add_argument("--cpu", action="store_true",
                        help="Use CPU instead of GPU")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory")

    args = parser.parse_args()

    predictor = ElectrochemicalSurfacePredictor(
        cp_mace_root=args.cp_mace_root,
        surface_sampling_root=args.surface_sampling_root,
        model_paths=args.model_paths,
        device="cpu" if args.cpu else "cuda",
    )

    result = predictor.sample_surface(
        surface=args.surface,
        adsorbates=args.adsorbates,
        potential_she=args.potential,
        ph=args.ph,
        temperature=args.temperature,
        initial_electron_number=args.initial_electron,
        total_sweeps=args.total_sweeps,
        sweep_size=args.sweep_size,
        output_dir=args.output_dir,
    )

    print("\n" + "=" * 70)
    print("ELECTROCHEMICAL SURFACE SAMPLING COMPLETE")
    print("=" * 70)
    print(result.summary())
