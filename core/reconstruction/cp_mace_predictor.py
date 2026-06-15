"""
CP-MACE Predictor - 常电位分子动力学模拟封装类

这个模块提供了一个完整的封装类，用于训练CP-MACE模型并进行常电位分子动力学模拟。

CP-MACE (Constant-Potential MACE) 是MACE框架的扩展，支持在巨正则系综条件下
对电化学界面进行常电位分子模拟。

使用方式:
    from cp_mace_predictor import CPMACEPredictor

    # 初始化预测器
    predictor = CPMACEPredictor(
        cp_mace_root="/path/to/CP-MACE",
        use_gpu=True,
    )

    # 方式1: 使用预训练模型直接模拟
    result = predictor.simulate(
        structure="init.xyz",  # 带有吸附物的表面结构
        model_paths=["model1.model", "model2.model"],
        target_potential=-3.36,  # 目标电极电位 (V)
        temperature=300.0,
        steps=1000,
    )

    # 方式2: 先训练模型，再模拟
    train_result = predictor.train(
        train_file="train.xyz",  # 训练数据文件
        model_name="my_model",
        max_num_epochs=300,
    )

    sim_result = predictor.simulate(
        structure="init.xyz",
        model_paths=[train_result.model_path],
        target_potential=-3.36,
    )

作者: Claude
"""

import os
import sys
import shutil
import tempfile
import subprocess
import copy
import yaml
import gc
import logging
from typing import Union, Optional, List, Dict, Any, Tuple, Literal
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class TrainingResult:
    """训练结果"""
    model_name: str
    model_path: str
    checkpoint_dir: str
    train_file: str
    valid_file: Optional[str]
    num_epochs: int
    final_loss: Optional[float]
    log_file: str

    def __repr__(self):
        return (f"TrainingResult(model={self.model_name}, "
                f"path={self.model_path}, epochs={self.num_epochs})")


@dataclass
class SimulationStep:
    """单步模拟信息"""
    step: int
    time_fs: float
    energy: float  # eV
    kinetic_energy: float  # eV
    potential_energy: float  # eV
    temperature: float  # K
    fermi_level: float  # V
    force_std: float  # eV/Å
    fermi_std: float  # V
    electron_number: float


@dataclass
class SimulationResult:
    """模拟结果"""
    surface_formula: str
    model_paths: List[str]
    target_potential: float
    temperature: float
    total_steps: int
    completed_steps: int
    timestep_fs: float
    trajectory: List["ase.Atoms"]
    energy_history: List[float]
    fermi_history: List[float]
    electron_history: List[float]
    force_std_history: List[float]
    fermi_std_history: List[float]
    final_structure: "ase.Atoms"
    output_dir: str
    high_uncertainty_structures: List["ase.Atoms"] = field(default_factory=list)

    def summary(self, show_last_n: int = 10) -> str:
        """生成模拟结果摘要"""
        lines = [
            "CP-MACE Constant-Potential MD Simulation Results",
            "=" * 60,
            f"Surface: {self.surface_formula}",
            f"Models: {len(self.model_paths)}",
            f"Target potential: {self.target_potential:.4f} V",
            f"Temperature: {self.temperature:.1f} K",
            f"Timestep: {self.timestep_fs:.1f} fs",
            f"Steps: {self.completed_steps}/{self.total_steps}",
            "",
            "Energy statistics (eV):",
            f"  Final: {self.energy_history[-1]:.4f}" if self.energy_history else "  N/A",
            f"  Mean: {np.mean(self.energy_history):.4f}" if self.energy_history else "  N/A",
            f"  Std: {np.std(self.energy_history):.4f}" if self.energy_history else "  N/A",
            "",
            "Fermi level statistics (V):",
            f"  Final: {self.fermi_history[-1]:.4f}" if self.fermi_history else "  N/A",
            f"  Mean: {np.mean(self.fermi_history):.4f}" if self.fermi_history else "  N/A",
            f"  Std: {np.std(self.fermi_history):.4f}" if self.fermi_history else "  N/A",
            "",
            f"Electron number (final): {self.electron_history[-1]:.2f}" if self.electron_history else "",
            f"High uncertainty structures collected: {len(self.high_uncertainty_structures)}",
            "",
            f"Output directory: {self.output_dir}",
        ]
        return "\n".join(lines)


class CPMACEPredictor:
    """
    CP-MACE 常电位分子动力学预测器

    这个类封装了CP-MACE项目的所有功能，提供简洁的API用于训练CP-MACE模型
    并进行常电位分子动力学模拟。

    Parameters
    ----------
    cp_mace_root : str
        CP-MACE项目的根目录路径
    use_gpu : bool, default=True
        是否使用GPU进行训练和推理
    device : str, optional
        指定设备 ("cuda" 或 "cpu")，默认根据use_gpu自动选择
    default_dtype : str, default="float64"
        默认数据类型
    work_dir : str, optional
        工作目录。如果不指定，将创建临时目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """

    def __init__(
        self,
        cp_mace_root: str,
        use_gpu: bool = True,
        device: Optional[str] = None,
        default_dtype: str = "float64",
        work_dir: Optional[str] = None,
        keep_files: bool = True,
        verbose: bool = True,
    ):
        self.cp_mace_root = os.path.abspath(cp_mace_root)
        self.use_gpu = use_gpu
        self.device = device or ("cuda" if use_gpu else "cpu")
        self.default_dtype = default_dtype
        self.keep_files = keep_files
        self.verbose = verbose

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="cp_mace_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 添加CP-MACE路径到系统路径
        mace_path = os.path.join(self.cp_mace_root, "mace")
        if mace_path not in sys.path:
            sys.path.insert(0, self.cp_mace_root)

        simulation_path = os.path.join(self.cp_mace_root, "simulation")
        if simulation_path not in sys.path:
            sys.path.insert(0, simulation_path)

        # 验证环境
        self._validate_environment()

        # 延迟加载
        self._calculators = None

        if self.verbose:
            print(f"CPMACEPredictor initialized:")
            print(f"  CP-MACE root: {self.cp_mace_root}")
            print(f"  Device: {self.device}")
            print(f"  Work dir: {self.work_dir}")

    def _validate_environment(self):
        """验证运行环境"""
        if not os.path.exists(self.cp_mace_root):
            raise FileNotFoundError(f"CP-MACE root not found: {self.cp_mace_root}")

        mace_init = os.path.join(self.cp_mace_root, "mace", "__init__.py")
        if not os.path.exists(mace_init):
            raise FileNotFoundError(
                f"CP-MACE mace package not found in {self.cp_mace_root}"
            )

        simulation_dir = os.path.join(self.cp_mace_root, "simulation")
        if not os.path.exists(simulation_dir):
            raise FileNotFoundError(
                f"Simulation directory not found: {simulation_dir}"
            )

    def _log(self, message: str):
        """输出日志"""
        if self.verbose:
            print(message)

    def _convert_structure(self, structure) -> "ase.Atoms":
        """将输入结构转换为ASE Atoms"""
        from ase import Atoms
        from ase.io import read

        if isinstance(structure, str):
            if not os.path.exists(structure):
                raise FileNotFoundError(f"Structure file not found: {structure}")
            return read(structure)
        elif isinstance(structure, Atoms):
            return structure.copy()
        else:
            # 尝试从pymatgen Structure转换
            try:
                from pymatgen.io.ase import AseAtomsAdaptor
                from pymatgen.core import Structure
                if isinstance(structure, Structure):
                    return AseAtomsAdaptor.get_atoms(structure)
            except ImportError:
                pass
            raise TypeError(
                f"Unsupported structure type: {type(structure)}. "
                "Expected str (file path), ase.Atoms, or pymatgen.Structure"
            )

    def prepare_training_data(
        self,
        structures: List[Union[str, "ase.Atoms"]],
        electron_numbers: List[float],
        potentials: List[float],
        energies: List[float],
        forces: List[np.ndarray],
        output_file: str,
    ) -> str:
        """
        准备CP-MACE训练数据

        CP-MACE需要xyz格式的训练数据，每个结构需要包含electron和potential标签。

        Parameters
        ----------
        structures : list
            结构列表 (ASE Atoms或文件路径)
        electron_numbers : list of float
            每个结构的电子数或净电荷
        potentials : list of float
            每个结构的费米能级/电极电位
        energies : list of float
            每个结构的参考能量 (eV)
        forces : list of np.ndarray
            每个结构的参考原子力 (eV/Å)
        output_file : str
            输出xyz文件路径

        Returns
        -------
        str
            输出文件路径
        """
        from ase.io import write

        self._log(f"Preparing training data with {len(structures)} structures...")

        atoms_list = []
        for i, struct in enumerate(structures):
            atoms = self._convert_structure(struct)

            # 设置参考数据
            atoms.info['REF_energy'] = energies[i]
            atoms.info['electron'] = electron_numbers[i]
            atoms.info['potential'] = potentials[i]
            atoms.arrays['REF_forces'] = forces[i]

            atoms_list.append(atoms)

        # 写入extended xyz格式
        output_path = os.path.abspath(output_file)
        write(output_path, atoms_list, format='extxyz')

        self._log(f"  Training data saved to: {output_path}")
        return output_path

    def train(
        self,
        train_file: str,
        model_name: str = "CP_MACE_model",
        valid_file: Optional[str] = None,
        valid_fraction: float = 0.05,
        model_type: Literal["FermiMACE", "FermiMACE_2"] = "FermiMACE",
        hidden_irreps: str = "128x0e + 128x1o",
        r_max: float = 5.0,
        batch_size: int = 10,
        max_num_epochs: int = 300,
        energy_weight: float = 1.0,
        forces_weight: float = 100.0,
        potential_weight: float = 10.0,
        E0s: str = "average",
        ema: bool = True,
        ema_decay: float = 0.99,
        amsgrad: bool = True,
        seed: int = 1,
        restart_latest: bool = True,
        output_dir: Optional[str] = None,
        extra_args: Optional[Dict[str, Any]] = None,
    ) -> TrainingResult:
        """
        训练CP-MACE模型

        Parameters
        ----------
        train_file : str
            训练数据文件路径 (extended xyz格式)
        model_name : str, default="CP_MACE_model"
            模型名称
        valid_file : str, optional
            验证数据文件路径。如果不指定，使用valid_fraction从训练数据中划分
        valid_fraction : float, default=0.05
            验证集比例 (当valid_file未指定时使用)
        model_type : str, default="FermiMACE"
            模型类型: "FermiMACE" (推荐，节点增强方法) 或 "FermiMACE_2" (全局状态方法)
        hidden_irreps : str, default="128x0e + 128x1o"
            隐藏层不可约表示
        r_max : float, default=5.0
            截断半径 (Å)
        batch_size : int, default=10
            批次大小
        max_num_epochs : int, default=300
            最大训练轮数
        energy_weight : float, default=1.0
            能量损失权重
        forces_weight : float, default=100.0
            力损失权重
        potential_weight : float, default=10.0
            电位/费米能级损失权重
        E0s : str, default="average"
            E0计算方式
        ema : bool, default=True
            是否使用指数移动平均
        ema_decay : float, default=0.99
            EMA衰减率
        amsgrad : bool, default=True
            是否使用AMSGrad优化器
        seed : int, default=1
            随机种子
        restart_latest : bool, default=True
            是否从最新检查点重启
        output_dir : str, optional
            输出目录
        extra_args : dict, optional
            额外的命令行参数

        Returns
        -------
        TrainingResult
            训练结果
        """
        self._log("\n" + "=" * 60)
        self._log("CP-MACE Model Training")
        self._log("=" * 60)

        # 设置输出目录
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "training", model_name)
        os.makedirs(output_dir, exist_ok=True)

        train_file = os.path.abspath(train_file)
        if not os.path.exists(train_file):
            raise FileNotFoundError(f"Training file not found: {train_file}")

        self._log(f"Training file: {train_file}")
        self._log(f"Model type: {model_type}")
        self._log(f"Output directory: {output_dir}")

        # 构建命令
        cmd = [
            "mace_run_train",
            f"--name={model_name}",
            f"--train_file={train_file}",
            f"--valid_fraction={valid_fraction}",
            "--config_type_weights={'Default':1.0}",
            "--energy_key=REF_energy",
            "--forces_key=REF_forces",
            f"--E0s={E0s}",
            f"--model={model_type}",
            "--loss=fermi_weighted",
            "--error_table=Fermi_PerAtomRMSE",
            f"--forces_weight={forces_weight}",
            f"--energy_weight={energy_weight}",
            f"--potential_weight={potential_weight}",
            f"--hidden_irreps={hidden_irreps}",
            f"--r_max={r_max}",
            f"--batch_size={batch_size}",
            f"--max_num_epochs={max_num_epochs}",
            f"--device={self.device}",
            f"--seed={seed}",
        ]

        if valid_file:
            cmd.append(f"--valid_file={os.path.abspath(valid_file)}")

        if ema:
            cmd.append("--ema")
            cmd.append(f"--ema_decay={ema_decay}")

        if amsgrad:
            cmd.append("--amsgrad")

        if restart_latest:
            cmd.append("--restart_latest")

        # 添加额外参数
        if extra_args:
            for key, value in extra_args.items():
                if isinstance(value, bool):
                    if value:
                        cmd.append(f"--{key}")
                else:
                    cmd.append(f"--{key}={value}")

        self._log(f"\nRunning training command...")
        self._log(" ".join(cmd[:5]) + " ...")

        # 运行训练
        log_file = os.path.join(output_dir, "train.log")
        with open(log_file, 'w') as f:
            result = subprocess.run(
                cmd,
                cwd=output_dir,
                stdout=f if not self.verbose else None,
                stderr=subprocess.STDOUT if not self.verbose else None,
                text=True,
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"Training failed. Check log file: {log_file}"
            )

        # 查找生成的模型文件
        model_path = os.path.join(output_dir, f"{model_name}.model")
        compiled_model_path = os.path.join(output_dir, f"{model_name}_compiled.model")

        if os.path.exists(compiled_model_path):
            final_model_path = compiled_model_path
        elif os.path.exists(model_path):
            final_model_path = model_path
        else:
            # 搜索checkpoints目录
            checkpoints_dir = os.path.join(output_dir, "checkpoints")
            if os.path.exists(checkpoints_dir):
                model_files = list(Path(checkpoints_dir).glob("*.model"))
                if model_files:
                    final_model_path = str(model_files[0])
                else:
                    raise FileNotFoundError(f"No model file found in {output_dir}")
            else:
                raise FileNotFoundError(f"No model file found in {output_dir}")

        self._log(f"\nTraining completed!")
        self._log(f"Model saved to: {final_model_path}")

        return TrainingResult(
            model_name=model_name,
            model_path=final_model_path,
            checkpoint_dir=output_dir,
            train_file=train_file,
            valid_file=valid_file,
            num_epochs=max_num_epochs,
            final_loss=None,  # 可以从日志中解析
            log_file=log_file,
        )

    def _load_calculators(
        self,
        model_paths: List[str],
    ) -> List["Calculator"]:
        """加载MACE计算器"""
        import torch
        from mace.calculators import MACECalculator

        calculators = []
        for path in model_paths:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Model file not found: {path}")
            calc = MACECalculator(
                model_paths=[path],
                device=self.device,
            )
            calculators.append(calc)

        return calculators

    def _create_average_calculator(
        self,
        calculators: List["Calculator"],
    ) -> "Calculator":
        """创建平均力计算器 (用于ensemble)"""
        from ase.calculators.calculator import Calculator
        import copy

        class AverageForceCalculator(Calculator):
            """计算多个模型的平均力、能量和费米能级"""
            implemented_properties = ['forces', 'energy', 'potential']

            def __init__(self, calculators, **kwargs):
                super().__init__(**kwargs)
                self.calculators = calculators

            def calculate(self, atoms=None, properties=['forces', 'energy', 'potential'],
                         system_changes=None):
                super().calculate(atoms, properties, system_changes)

                total_forces = 0
                total_energy = 0
                total_mu = 0
                all_forces = []
                all_mu = []

                for calc in self.calculators:
                    atoms_copy = copy.deepcopy(atoms)
                    atoms_copy.calc = calc
                    calc.calculate(atoms_copy, properties, system_changes)
                    total_forces += calc.results['forces']
                    total_energy += calc.results['energy']
                    total_mu += calc.results['potential']
                    all_forces.append(calc.results['forces'])
                    all_mu.append(calc.results['potential'])

                n = len(self.calculators)
                self.results['forces'] = total_forces / n
                self.results['energy'] = total_energy / n
                self.results['potential'] = total_mu / n
                self.results['mu'] = total_mu / n
                atoms.info['potential'] = total_mu / n

                # 计算不确定性
                all_forces_array = np.array(all_forces)
                force_std = np.std(all_forces_array, axis=0)
                total_std = np.sqrt(np.sum(force_std**2, axis=1))
                self.results['force_std'] = np.max(total_std)
                self.results['mu_std'] = np.std(np.array(all_mu))

            def get_max_std(self):
                return self.results.get('force_std', 0)

            def get_mu_std(self):
                return self.results.get('mu_std', 0)

            def get_mu(self, atoms=None):
                return self.results.get('mu', 0)

        return AverageForceCalculator(calculators)

    def simulate(
        self,
        structure: Union[str, "ase.Atoms", "pymatgen.core.Structure"],
        model_paths: List[str],
        target_potential: float,
        initial_electron_number: Optional[float] = None,
        temperature: float = 300.0,
        timestep: float = 1.0,
        steps: int = 1000,
        integrator: Literal["NoseHoover", "Langevin", "VelocityVerlet"] = "NoseHoover",
        ttime: float = 40.0,
        eta_length: int = 2,
        Mne: Optional[float] = None,
        constraints: Optional[List] = None,
        constraint_increment: float = 0.0,
        read_velocity: bool = False,
        force_threshold: float = 0.15,
        fermi_threshold: float = 0.04,
        save_frequency: int = 1,
        output_dir: Optional[str] = None,
        run_name: Optional[str] = None,
    ) -> SimulationResult:
        """
        运行常电位分子动力学模拟

        Parameters
        ----------
        structure : str or ase.Atoms or pymatgen.Structure
            初始结构（带有吸附物的表面）。需要包含'electron'属性。
        model_paths : list of str
            CP-MACE模型文件路径列表。建议使用2个以上模型进行ensemble计算。
        target_potential : float
            目标电极电位 (V vs. SHE)
        initial_electron_number : float, optional
            初始电子数。如果不指定，从结构的info中读取。
        temperature : float, default=300.0
            模拟温度 (K)
        timestep : float, default=1.0
            时间步长 (fs)
        steps : int, default=1000
            模拟步数
        integrator : str, default="NoseHoover"
            积分器类型: "NoseHoover", "Langevin", "VelocityVerlet"
        ttime : float, default=40.0
            热浴耦合时间 (fs)
        eta_length : int, default=2
            Nose-Hoover链长度
        Mne : float, optional
            电子质量参数。如果不指定，使用初始电子数。
        constraints : list, optional
            约束列表，格式: [[type, atom1, atom2, distance], ...]
            type=0表示固定距离约束
        constraint_increment : float, default=0.0
            每步约束距离增量 (Å)，用于慢增长模拟
        read_velocity : bool, default=False
            是否从结构文件读取初始速度
        force_threshold : float, default=0.15
            力标准差阈值，超过此值收集结构 (eV/Å)
        fermi_threshold : float, default=0.04
            费米能级标准差阈值，超过此值收集结构 (V)
        save_frequency : int, default=1
            保存频率 (步数)
        output_dir : str, optional
            输出目录
        run_name : str, optional
            模拟名称

        Returns
        -------
        SimulationResult
            模拟结果
        """
        import torch
        from ase import units
        from ase.io import Trajectory, write
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

        self._log("\n" + "=" * 60)
        self._log("CP-MACE Constant-Potential MD Simulation")
        self._log("=" * 60)

        # 转换结构
        atoms = self._convert_structure(structure)
        atoms.pbc = [True, True, True]
        formula = atoms.get_chemical_formula()

        self._log(f"Structure: {formula}")
        self._log(f"Number of atoms: {len(atoms)}")
        self._log(f"Models: {len(model_paths)}")
        self._log(f"Target potential: {target_potential:.4f} V")
        self._log(f"Temperature: {temperature:.1f} K")
        self._log(f"Timestep: {timestep:.1f} fs")
        self._log(f"Steps: {steps}")

        # 设置输出目录
        if output_dir is None:
            output_dir = os.path.join(
                self.work_dir,
                run_name or f"cp_mace_sim_{formula}"
            )
        os.makedirs(output_dir, exist_ok=True)

        # 设置电子数
        if initial_electron_number is not None:
            atoms.info['electron'] = initial_electron_number
        elif 'electron' not in atoms.info:
            raise ValueError(
                "Initial electron number not provided and not found in structure. "
                "Please specify initial_electron_number parameter or set "
                "atoms.info['electron'] in the structure."
            )

        if Mne is None:
            Mne = atoms.info['electron']

        self._log(f"Initial electron number: {atoms.info['electron']:.2f}")

        # 设置dtype
        torch.set_default_dtype(torch.float64)

        # 加载计算器
        self._log("\nLoading models...")
        calculators = self._load_calculators(model_paths)

        if len(calculators) > 1:
            avg_calculator = self._create_average_calculator(calculators)
            atoms.calc = avg_calculator
            self._log(f"  Using ensemble of {len(calculators)} models")
        else:
            atoms.calc = calculators[0]
            self._log("  Using single model")

        # 创建积分器配置
        integrator_config = {
            "timestep": timestep * units.fs,
            "temperature": temperature * units.kB,
            "ttime": ttime,
            "Mne": Mne,
            "eta_length": eta_length,
            "targetmu": target_potential,
            "constraints": constraints or [],
            "increm": constraint_increment,
        }

        # 初始化速度
        if not read_velocity:
            MaxwellBoltzmannDistribution(atoms, temperature * units.kB)
            Stationary(atoms)
            self._log(f"  Initialized velocities at {temperature:.1f} K")

        # 创建积分器
        self._log(f"  Integrator: {integrator}")

        # 动态导入积分器
        slow_growth_path = os.path.join(self.cp_mace_root, "simulation", "slow_growth")
        if slow_growth_path not in sys.path:
            sys.path.insert(0, slow_growth_path)

        import integrator as md_integrator

        dyn = getattr(md_integrator, integrator)(atoms, **integrator_config)

        # 设置轨迹和日志
        traj_path = os.path.join(output_dir, "atoms.traj")
        traj = Trajectory(traj_path, 'w', atoms)
        dyn.attach(traj.write, interval=save_frequency)

        # 运行模拟
        self._log(f"\nRunning simulation...")

        energy_history = []
        fermi_history = []
        electron_history = []
        force_std_history = []
        fermi_std_history = []
        high_uncertainty_structures = []
        trajectory = []

        try:
            from tqdm import tqdm
            iterator = tqdm(range(steps))
        except ImportError:
            iterator = range(steps)

        completed_steps = 0
        for step in iterator:
            try:
                dyn.run(1)
                completed_steps += 1

                # 收集数据
                energy = atoms.get_potential_energy()
                energy_history.append(energy)

                calc = atoms.calc
                if hasattr(calc, 'get_mu'):
                    fermi = calc.get_mu()
                elif hasattr(calc, 'results') and 'potential' in calc.results:
                    fermi = calc.results['potential']
                else:
                    fermi = atoms.info.get('potential', 0)
                fermi_history.append(fermi)

                electron_history.append(atoms.info.get('electron', 0))

                # 不确定性
                if hasattr(calc, 'get_max_std'):
                    force_std = calc.get_max_std()
                    fermi_std = calc.get_mu_std()
                else:
                    force_std = 0
                    fermi_std = 0
                force_std_history.append(force_std)
                fermi_std_history.append(fermi_std)

                # 检查阈值，收集高不确定性结构
                if force_std > force_threshold:
                    self._log(f"  Step {step}: Force std {force_std:.4f} > {force_threshold}")
                    high_uncertainty_structures.append(atoms.copy())
                elif fermi_std > fermi_threshold:
                    self._log(f"  Step {step}: Fermi std {fermi_std:.4f} > {fermi_threshold}")
                    high_uncertainty_structures.append(atoms.copy())

                if step % save_frequency == 0:
                    trajectory.append(atoms.copy())

                # 清理GPU内存
                if step % 100 == 0:
                    gc.collect()
                    torch.cuda.empty_cache()

            except Exception as e:
                self._log(f"  Error at step {step}: {e}")
                break

        traj.close()

        # 保存最终结构
        final_path = os.path.join(output_dir, "final.xyz")
        write(final_path, atoms, format='extxyz')

        # 保存高不确定性结构
        if high_uncertainty_structures:
            uncertain_path = os.path.join(output_dir, "high_uncertainty.xyz")
            write(uncertain_path, high_uncertainty_structures, format='extxyz')
            self._log(f"  High uncertainty structures saved: {len(high_uncertainty_structures)}")

        self._log(f"\nSimulation completed!")
        self._log(f"  Completed steps: {completed_steps}/{steps}")
        self._log(f"  Output directory: {output_dir}")

        result = SimulationResult(
            surface_formula=formula,
            model_paths=model_paths,
            target_potential=target_potential,
            temperature=temperature,
            total_steps=steps,
            completed_steps=completed_steps,
            timestep_fs=timestep,
            trajectory=trajectory,
            energy_history=energy_history,
            fermi_history=fermi_history,
            electron_history=electron_history,
            force_std_history=force_std_history,
            fermi_std_history=fermi_std_history,
            final_structure=atoms.copy(),
            output_dir=output_dir,
            high_uncertainty_structures=high_uncertainty_structures,
        )

        self._log(f"\n{result.summary()}")

        return result

    def run_slow_growth(
        self,
        structure: Union[str, "ase.Atoms"],
        model_paths: List[str],
        target_potential: float,
        atom1_index: int,
        atom2_index: int,
        initial_distance: float,
        final_distance: float,
        steps: int = 1000,
        temperature: float = 300.0,
        timestep: float = 1.0,
        output_dir: Optional[str] = None,
        **kwargs,
    ) -> SimulationResult:
        """
        运行慢增长模拟（用于自由能计算）

        Parameters
        ----------
        structure : str or ase.Atoms
            初始结构
        model_paths : list of str
            模型路径列表
        target_potential : float
            目标电极电位 (V)
        atom1_index : int
            约束原子1的索引
        atom2_index : int
            约束原子2的索引
        initial_distance : float
            初始距离 (Å)
        final_distance : float
            最终距离 (Å)
        steps : int, default=1000
            模拟步数
        temperature : float, default=300.0
            温度 (K)
        timestep : float, default=1.0
            时间步长 (fs)
        output_dir : str, optional
            输出目录
        **kwargs
            其他参数传递给simulate()

        Returns
        -------
        SimulationResult
            模拟结果
        """
        increment = (final_distance - initial_distance) / steps

        constraints = [[0, atom1_index, atom2_index, initial_distance]]

        self._log(f"Slow growth simulation:")
        self._log(f"  Reaction coordinate: atoms {atom1_index}-{atom2_index}")
        self._log(f"  Distance: {initial_distance:.4f} -> {final_distance:.4f} Å")
        self._log(f"  Increment per step: {increment:.6f} Å")

        return self.simulate(
            structure=structure,
            model_paths=model_paths,
            target_potential=target_potential,
            temperature=temperature,
            timestep=timestep,
            steps=steps,
            constraints=constraints,
            constraint_increment=increment,
            output_dir=output_dir,
            **kwargs,
        )

    def __del__(self):
        """清理临时目录"""
        if hasattr(self, '_temp_dir') and self._temp_dir and not self.keep_files:
            try:
                shutil.rmtree(self._temp_dir)
            except:
                pass


# 便捷函数
def train_cp_mace(
    train_file: str,
    cp_mace_root: str,
    model_name: str = "CP_MACE_model",
    use_gpu: bool = True,
    output_dir: Optional[str] = None,
    **kwargs,
) -> TrainingResult:
    """
    便捷函数：训练CP-MACE模型

    Parameters
    ----------
    train_file : str
        训练数据文件路径
    cp_mace_root : str
        CP-MACE项目根目录
    model_name : str, default="CP_MACE_model"
        模型名称
    use_gpu : bool, default=True
        是否使用GPU
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数传递给CPMACEPredictor.train()

    Returns
    -------
    TrainingResult
        训练结果
    """
    predictor = CPMACEPredictor(
        cp_mace_root=cp_mace_root,
        use_gpu=use_gpu,
    )
    return predictor.train(
        train_file=train_file,
        model_name=model_name,
        output_dir=output_dir,
        **kwargs,
    )


def run_cp_md(
    structure: Union[str, "ase.Atoms"],
    model_paths: List[str],
    cp_mace_root: str,
    target_potential: float,
    temperature: float = 300.0,
    steps: int = 1000,
    use_gpu: bool = True,
    output_dir: Optional[str] = None,
    **kwargs,
) -> SimulationResult:
    """
    便捷函数：运行常电位分子动力学模拟

    Parameters
    ----------
    structure : str or ase.Atoms
        初始结构
    model_paths : list of str
        模型路径列表
    cp_mace_root : str
        CP-MACE项目根目录
    target_potential : float
        目标电极电位 (V)
    temperature : float, default=300.0
        温度 (K)
    steps : int, default=1000
        模拟步数
    use_gpu : bool, default=True
        是否使用GPU
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数

    Returns
    -------
    SimulationResult
        模拟结果
    """
    predictor = CPMACEPredictor(
        cp_mace_root=cp_mace_root,
        use_gpu=use_gpu,
    )
    return predictor.simulate(
        structure=structure,
        model_paths=model_paths,
        target_potential=target_potential,
        temperature=temperature,
        steps=steps,
        output_dir=output_dir,
        **kwargs,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="CP-MACE Constant-Potential Molecular Dynamics"
    )
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # Train command
    train_parser = subparsers.add_parser("train", help="Train CP-MACE model")
    train_parser.add_argument("train_file", help="Training data file (xyz)")
    train_parser.add_argument("--cp-mace-root",
                              default="deps/CP-MACE",
                              help="CP-MACE root directory")
    train_parser.add_argument("--model-name", default="CP_MACE_model",
                              help="Model name")
    train_parser.add_argument("--max-epochs", type=int, default=300,
                              help="Maximum training epochs")
    train_parser.add_argument("--cpu", action="store_true",
                              help="Use CPU instead of GPU")
    train_parser.add_argument("--output-dir", default=None,
                              help="Output directory")

    # Simulate command
    sim_parser = subparsers.add_parser("simulate", help="Run CP-MD simulation")
    sim_parser.add_argument("structure", help="Initial structure file")
    sim_parser.add_argument("--model-paths", nargs="+", required=True,
                            help="Model file paths")
    sim_parser.add_argument("--cp-mace-root",
                            default="deps/CP-MACE",
                            help="CP-MACE root directory")
    sim_parser.add_argument("--target-potential", type=float, required=True,
                            help="Target electrode potential (V)")
    sim_parser.add_argument("--initial-electron", type=float, default=None,
                            help="Initial electron number")
    sim_parser.add_argument("--temperature", type=float, default=300.0,
                            help="Temperature (K)")
    sim_parser.add_argument("--steps", type=int, default=1000,
                            help="Number of simulation steps")
    sim_parser.add_argument("--timestep", type=float, default=1.0,
                            help="Timestep (fs)")
    sim_parser.add_argument("--cpu", action="store_true",
                            help="Use CPU instead of GPU")
    sim_parser.add_argument("--output-dir", default=None,
                            help="Output directory")

    args = parser.parse_args()

    if args.command == "train":
        predictor = CPMACEPredictor(
            cp_mace_root=args.cp_mace_root,
            use_gpu=not args.cpu,
        )
        result = predictor.train(
            train_file=args.train_file,
            model_name=args.model_name,
            max_num_epochs=args.max_epochs,
            output_dir=args.output_dir,
        )
        print("\n" + "=" * 60)
        print("TRAINING COMPLETE")
        print("=" * 60)
        print(f"Model saved to: {result.model_path}")

    elif args.command == "simulate":
        predictor = CPMACEPredictor(
            cp_mace_root=args.cp_mace_root,
            use_gpu=not args.cpu,
        )
        result = predictor.simulate(
            structure=args.structure,
            model_paths=args.model_paths,
            target_potential=args.target_potential,
            initial_electron_number=args.initial_electron,
            temperature=args.temperature,
            steps=args.steps,
            timestep=args.timestep,
            output_dir=args.output_dir,
        )
        print("\n" + "=" * 60)
        print("SIMULATION COMPLETE")
        print("=" * 60)
        print(result.summary())

    else:
        parser.print_help()
