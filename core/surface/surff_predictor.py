"""
SurFF Predictor - 金属间化合物表面能预测封装类

这个模块提供了一个完整的封装类，用于预测任意bulk结构的表面暴露情况。

使用方式:
    from surff_predictor import SurFFPredictor

    # 初始化预测器
    predictor = SurFFPredictor(
        surff_root="/path/to/SurFF",
        checkpoint_path="/path/to/checkpoint.pt",  # 可选，默认使用SurFF自带的模型
        use_gpu=True,
    )

    # 预测表面暴露
    results = predictor.predict(
        structure="POSCAR",  # 可以是文件路径、ASE Atoms或pymatgen Structure
        top_n=5,  # 返回最可能暴露的前N个表面
    )

作者: Claude
"""

import os
import sys
import shutil
import tempfile
import subprocess
import pickle
from typing import Union, Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class SlabInfo:
    """表面slab信息"""
    slab_id: str
    miller_index: Tuple[int, int, int]
    shift: float
    num_atoms: int
    formula: str


@dataclass
class SurfaceResult:
    """单个表面的预测结果"""
    slab_id: str
    miller_index: Tuple[int, int, int]
    shift: float
    surface_energy: float  # eV/Å²
    area_fraction: float   # 暴露面积比例 (0-1)
    exposure_level: str    # "high", "medium", "low"

    def __repr__(self):
        return (f"Surface({self.miller_index}, shift={self.shift:.3f}, "
                f"E={self.surface_energy:.4f} eV/Å², area={self.area_fraction:.2%}, "
                f"level={self.exposure_level})")


@dataclass
class PredictionResult:
    """完整的预测结果"""
    crystal_id: str
    formula: str
    total_slabs: int
    surfaces: List[SurfaceResult]
    wulff_shape_path: Optional[str] = None
    raw_results: Dict = field(default_factory=dict)

    def get_top_n(self, n: int) -> List[SurfaceResult]:
        """获取暴露面积最大的前N个表面"""
        sorted_surfaces = sorted(self.surfaces, key=lambda x: x.area_fraction, reverse=True)
        return sorted_surfaces[:n]

    def get_exposed_surfaces(self, threshold: float = 0.01) -> List[SurfaceResult]:
        """获取所有暴露面积大于阈值的表面"""
        return [s for s in self.surfaces if s.area_fraction > threshold]

    def to_dataframe(self) -> pd.DataFrame:
        """转换为DataFrame"""
        data = []
        for s in self.surfaces:
            data.append({
                "miller_index": str(s.miller_index),
                "shift": s.shift,
                "surface_energy": s.surface_energy,
                "area_fraction": s.area_fraction,
                "exposure_level": s.exposure_level,
                "slab_id": s.slab_id,
            })
        return pd.DataFrame(data).sort_values("area_fraction", ascending=False)

    def summary(self, top_n: int = 5) -> str:
        """生成预测结果摘要"""
        lines = [
            f"Crystal: {self.formula} (ID: {self.crystal_id})",
            f"Total slabs analyzed: {self.total_slabs}",
            f"",
            f"Top {top_n} most exposed surfaces:",
            "-" * 60,
        ]
        for i, s in enumerate(self.get_top_n(top_n), 1):
            lines.append(f"  {i}. {s.miller_index} (shift={s.shift:.3f})")
            lines.append(f"     Surface energy: {s.surface_energy:.4f} eV/Å²")
            lines.append(f"     Area fraction: {s.area_fraction:.2%}")
            lines.append(f"     Exposure level: {s.exposure_level}")
        return "\n".join(lines)


class SurFFPredictor:
    """
    SurFF表面能预测器

    这个类封装了SurFF项目的所有功能，提供简洁的API用于预测任意bulk结构的表面暴露情况。

    Parameters
    ----------
    surff_root : str
        SurFF项目的根目录路径
    checkpoint_path : str, optional
        模型checkpoint文件路径。如果不指定，将使用SurFF自带的模型
    use_gpu : bool, default=True
        是否使用GPU进行推理
    max_miller_index : int, default=2
        生成slab时的最大Miller指数
    min_slab_size : float, default=10.0
        最小slab厚度 (Å)
    min_vacuum_size : float, default=15.0
        最小真空层厚度 (Å)
    fix_distance : float, default=3.0
        固定原子的距离阈值 (Å)，距离表面超过此距离的原子将被固定
    relaxation_steps : int, default=200
        最大弛豫步数
    relaxation_fmax : float, default=0.05
        弛豫收敛力阈值 (eV/Å)
    max_step : float, default=0.03
        LBFGS优化器的最大步长
    work_dir : str, optional
        工作目录，用于存储中间文件。如果不指定，将创建临时目录
    keep_files : bool, default=False
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """

    def __init__(
        self,
        surff_root: str,
        checkpoint_path: Optional[str] = None,
        use_gpu: bool = True,
        max_miller_index: int = 2,
        min_slab_size: float = 10.0,
        min_vacuum_size: float = 15.0,
        fix_distance: float = 3.0,
        relaxation_steps: int = 200,
        relaxation_fmax: float = 0.05,
        max_step: float = 0.03,
        work_dir: Optional[str] = None,
        keep_files: bool = False,
        verbose: bool = True,
    ):
        self.surff_root = os.path.abspath(surff_root)
        self.ocp_root = os.path.join(self.surff_root, "ocp")
        self.use_gpu = use_gpu
        self.max_miller_index = max_miller_index
        self.min_slab_size = min_slab_size
        self.min_vacuum_size = min_vacuum_size
        self.fix_distance = fix_distance
        self.relaxation_steps = relaxation_steps
        self.relaxation_fmax = relaxation_fmax
        self.max_step = max_step
        self.keep_files = keep_files
        self.verbose = verbose

        # 设置checkpoint路径
        if checkpoint_path is None:
            # 使用默认的checkpoint
            default_ckpt = os.path.join(
                self.surff_root,
                "SurFF_DataFiles/ocp/checkpoints/2024-03-03-23-57-52/best_checkpoint.pt"
            )
            if os.path.exists(default_ckpt):
                self.checkpoint_path = default_ckpt
            else:
                # 尝试ocp目录下的checkpoint
                alt_ckpt = os.path.join(
                    self.ocp_root,
                    "checkpoints/2024-03-03-23-57-52/best_checkpoint.pt"
                )
                if os.path.exists(alt_ckpt):
                    self.checkpoint_path = alt_ckpt
                else:
                    raise FileNotFoundError(
                        f"未找到默认checkpoint文件。请指定checkpoint_path参数。"
                    )
        else:
            self.checkpoint_path = os.path.abspath(checkpoint_path)

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="surff_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 添加SurFF路径到系统路径
        if self.surff_root not in sys.path:
            sys.path.insert(0, self.surff_root)

        # 验证环境
        self._validate_environment()

        if self.verbose:
            print(f"SurFFPredictor initialized:")
            print(f"  SurFF root: {self.surff_root}")
            print(f"  Checkpoint: {self.checkpoint_path}")
            print(f"  Use GPU: {self.use_gpu}")
            print(f"  Work dir: {self.work_dir}")

    def _validate_environment(self):
        """验证运行环境"""
        # 检查必要的目录
        if not os.path.exists(self.surff_root):
            raise FileNotFoundError(f"SurFF root not found: {self.surff_root}")
        if not os.path.exists(self.ocp_root):
            raise FileNotFoundError(f"OCP root not found: {self.ocp_root}")
        if not os.path.exists(self.checkpoint_path):
            raise FileNotFoundError(f"Checkpoint not found: {self.checkpoint_path}")

        # 检查必要的配置文件
        config_path = os.path.join(self.ocp_root, "configs/equiformer_v2_002_relax.yml")
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")

    def _log(self, message: str):
        """输出日志"""
        if self.verbose:
            print(message)

    def _convert_structure(self, structure) -> "pymatgen.core.Structure":
        """将输入结构转换为pymatgen Structure"""
        from pymatgen.core import Structure
        from pymatgen.io.vasp import Poscar
        from pymatgen.io.ase import AseAtomsAdaptor

        if isinstance(structure, str):
            # 文件路径
            if not os.path.exists(structure):
                raise FileNotFoundError(f"Structure file not found: {structure}")
            return Poscar.from_file(structure).structure
        elif isinstance(structure, Structure):
            return structure
        else:
            # 尝试作为ASE Atoms处理
            try:
                from ase import Atoms
                if isinstance(structure, Atoms):
                    return AseAtomsAdaptor.get_structure(structure)
            except ImportError:
                pass
            raise TypeError(
                f"Unsupported structure type: {type(structure)}. "
                "Expected str (file path), pymatgen.Structure, or ase.Atoms"
            )

    def _generate_slabs(
        self,
        structure: "pymatgen.core.Structure",
        crystal_id: str,
        slab_dir: str,
        bulk_dir: str,
    ) -> pd.DataFrame:
        """生成所有可能的slab结构"""
        from pymatgen.core import Structure, Lattice
        from pymatgen.core.surface import generate_all_slabs
        from pymatgen.io.vasp import Poscar
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        self._log("Generating slab structures...")

        # 获取conventional cell
        structure = SpacegroupAnalyzer(structure, symprec=0.01).get_conventional_standard_structure()
        formula = structure.formula

        # 生成所有slabs
        try:
            slabs = generate_all_slabs(
                structure,
                max_index=self.max_miller_index,
                min_slab_size=self.min_slab_size,
                min_vacuum_size=self.min_vacuum_size,
                primitive=False,
                center_slab=True
            )
        except Exception as e:
            raise RuntimeError(f"Failed to generate slabs: {e}")

        if len(slabs) == 0:
            raise RuntimeError("No slabs generated. Check the input structure.")

        self._log(f"  Generated {len(slabs)} unique slabs")

        # 保存slab和bulk结构
        slab_info = []
        for j, slab in enumerate(slabs):
            miller_index = slab.miller_index
            miller_str = ''.join(map(str, miller_index))
            shift = slab.shift
            num_atoms = len(slab.sites)

            slab_id = f"{crystal_id}_{j}"

            # 保存slab (带selective dynamics)
            slab_path = os.path.join(slab_dir, slab_id)
            poscar = self._fix_middle_atoms(slab)
            with open(slab_path, 'w') as f:
                f.write(str(poscar))

            # 保存oriented unit cell
            bulk_path = os.path.join(bulk_dir, slab_id)
            ouc = slab.oriented_unit_cell
            ouc_lattice = Lattice.from_parameters(
                ouc.lattice.a, ouc.lattice.b, ouc.lattice.c,
                ouc.lattice.alpha, ouc.lattice.beta, ouc.lattice.gamma,
            )
            ouc = Structure(ouc_lattice, ouc.species, ouc.frac_coords, coords_are_cartesian=False)
            ouc.to(filename=bulk_path, fmt='poscar')

            slab_info.append({
                'slab_id': slab_id,
                'formula': formula,
                'miller_index': miller_str,
                'shift': shift,
                'num_atom': num_atoms,
            })

        return pd.DataFrame(slab_info)

    def _fix_middle_atoms(self, structure: "pymatgen.core.Structure") -> "pymatgen.io.vasp.Poscar":
        """固定slab中间的原子"""
        from pymatgen.io.vasp import Poscar

        selective_dynamics = []

        A = structure.lattice.matrix[0]
        B = structure.lattice.matrix[1]
        N = np.cross(A, B)
        N = N / np.linalg.norm(N)

        distances = [abs(np.dot(N, site)) for site in structure.cart_coords]
        d_max = max(distances)
        d_min = min(distances)

        for d in distances:
            dist_from_top = d_max - d
            dist_from_bottom = d - d_min

            if dist_from_top > self.fix_distance and dist_from_bottom > self.fix_distance:
                selective_dynamics.append([False, False, False])
            else:
                selective_dynamics.append([True, True, True])

        return Poscar(structure, selective_dynamics=selective_dynamics)

    def _generate_lmdb(self, slab_dir: str, lmdb_path: str, slab_info: pd.DataFrame):
        """生成LMDB数据集"""
        import lmdb
        import torch
        from pymatgen.io.vasp import Poscar
        from torch_geometric.data import Data

        self._log("Generating LMDB dataset...")

        slab_files = [f for f in os.listdir(slab_dir) if not f.endswith('.csv')]
        dataset = [{'slab_id': f, 'POSCAR_pth': os.path.join(slab_dir, f)} for f in slab_files]

        db = lmdb.open(lmdb_path, map_size=10**9, subdir=False, meminit=False, map_async=True)
        txn = db.begin(write=True)

        db_idx = 0
        for item in dataset:
            try:
                poscar_path = item['POSCAR_pth']
                sid = item['slab_id']

                crystal = Poscar.from_file(poscar_path).structure
                natoms = len(crystal)

                pos = torch.Tensor(crystal.cart_coords.tolist())

                try:
                    tags = crystal.site_properties['selective_dynamics']
                    tags = [1 if x[0] else 0 for x in tags]
                    tags = torch.LongTensor(tags)
                except:
                    tags = torch.LongTensor([1] * natoms)

                fixed = (tags == 0).float()
                lattice = torch.Tensor([crystal.lattice.matrix.tolist()])
                atom_fea = torch.Tensor([crystal[i].specie.number for i in range(natoms)])

                data = Data(
                    pos=pos, cell=lattice, atomic_numbers=atom_fea,
                    natoms=natoms, tags=tags, fixed=fixed, sid=sid, fid=0
                )

                txn.put(key=f"{db_idx}".encode("ascii"), value=pickle.dumps(data, protocol=-1))
                db_idx += 1

            except Exception as e:
                self._log(f"  Warning: Error processing {sid}: {e}")

        txn.commit()
        db.sync()

        # 存储长度
        txn = db.begin(write=True)
        txn.put(key='length'.encode("ascii"), value=pickle.dumps(db_idx, protocol=-1))
        txn.commit()
        db.sync()
        db.close()

        self._log(f"  Created LMDB with {db_idx} entries")

    def _run_relaxation(self, lmdb_dir: str, traj_dir: str):
        """运行MLFF弛豫"""
        self._log("Running ML relaxation...")

        # Convert absolute paths to relative paths from ocp root
        lmdb_rel = os.path.relpath(os.path.abspath(lmdb_dir), self.ocp_root)
        traj_rel = os.path.relpath(os.path.abspath(traj_dir), self.ocp_root)

        cmd = [
            "python", "main.py",
            "--mode", "run-relaxations",
            "--config-yml", "configs/equiformer_v2_002_relax.yml",
            "--checkpoint", self.checkpoint_path,
            f"--task.relax_opt.traj_dir={traj_rel}",
            f"--task.relax_opt.maxstep={self.max_step}",
            f"--task.relax_dataset.src={lmdb_rel}",
            f"--task.relaxation_steps={self.relaxation_steps}",
            f"--task.relaxation_fmax={self.relaxation_fmax}",
        ]

        if not self.use_gpu:
            cmd.append("--cpu")

        result = subprocess.run(
            cmd,
            cwd=self.ocp_root,
            capture_output=not self.verbose,
            text=True,
        )

        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else "Unknown error"
            raise RuntimeError(f"Relaxation failed: {error_msg}")

        self._log("  Relaxation completed")

    def _collect_results(self, traj_dir: str) -> pd.DataFrame:
        """收集弛豫结果"""
        from ase.io import Trajectory

        self._log("Collecting relaxation results...")

        traj_files = [f for f in os.listdir(traj_dir) if f.endswith('.traj')]

        results = []
        for traj_file in traj_files:
            slab_id = os.path.splitext(traj_file)[0]
            crystal_id = slab_id.split('_')[0]

            traj = Trajectory(os.path.join(traj_dir, traj_file))
            final_atoms = traj[-1]

            energy = final_atoms.get_potential_energy()
            cell = final_atoms.get_cell()
            area = np.linalg.norm(np.cross(cell[0], cell[1]))
            surface_energy = energy / (2 * area)

            results.append({
                'slab_id': slab_id,
                'crystal_id': crystal_id,
                'surface_energy': surface_energy,
            })

        return pd.DataFrame(results)

    def _calculate_wulff(
        self,
        results_df: pd.DataFrame,
        slab_info: pd.DataFrame,
        crystal_dir: str,
        crystal_id: str,
        save_dir: str,
    ) -> Dict:
        """计算Wulff构型"""
        from pymatgen.analysis.wulff import WulffShape
        from pymatgen.io.vasp import Poscar
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        import re

        self._log("Calculating Wulff shape...")

        # 合并结果
        df = pd.merge(results_df, slab_info, on='slab_id', how='inner')
        df = df.rename(columns={'surface_energy': 'surface_energy_pred'})

        # 读取晶体结构
        crystal_path = os.path.join(crystal_dir, crystal_id)
        crystal_struct = Poscar.from_file(crystal_path).structure
        crystal_struct = SpacegroupAnalyzer(crystal_struct).get_conventional_standard_structure()

        # 解析Miller指数
        def str_to_tuple(s):
            matches = re.findall(r'-?\d', str(s))
            return tuple(int(num) for num in matches)

        miller_indices = [str_to_tuple(m) for m in df['miller_index'].values]
        surface_energies = df['surface_energy_pred'].values.tolist()

        # 计算Wulff形状
        wulff = WulffShape(crystal_struct.lattice, miller_indices, surface_energies)

        # 提取结果
        wulff_results = {
            'miller_indices': miller_indices,
            'surface_energies': surface_energies,
            'area_fractions': [],
            'shifts': df['shift'].values.tolist(),
            'slab_ids': df['slab_id'].values.tolist(),
        }

        # 计算每个slab的暴露面积
        for i, (miller, energy) in enumerate(zip(miller_indices, surface_energies)):
            area = 0.0
            for j, wulff_miller in enumerate(wulff.miller_list):
                if wulff_miller == miller:
                    if abs(wulff.e_surf_list[j] - energy) < 1e-5:
                        area = wulff.color_area[j] / sum(wulff.color_area) if sum(wulff.color_area) > 0 else 0
                        break
            wulff_results['area_fractions'].append(area)

        # 保存Wulff形状图
        wulff_shape_path = None
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            try:
                import matplotlib
                matplotlib.use('Agg')  # 使用非交互式后端
                import matplotlib.pyplot as plt
                ax_3d = wulff.get_plot(
                    color_set='YlGnBu',
                    grid_off=True,
                    axis_off=True,
                    show_area=True,
                )
                wulff_shape_path = os.path.join(save_dir, f"{crystal_id}_wulff.png")
                ax_3d.get_figure().savefig(wulff_shape_path, dpi=150, bbox_inches='tight')
                plt.close('all')
                self._log(f"  Wulff shape saved to: {wulff_shape_path}")
            except Exception as e:
                self._log(f"  Warning: Failed to save Wulff shape: {e}")

        return wulff_results, wulff_shape_path

    def predict(
        self,
        structure: Union[str, "pymatgen.core.Structure", "ase.Atoms"],
        crystal_id: str = "crystal",
        top_n: int = 5,
        save_wulff_shape: bool = True,
        output_dir: Optional[str] = None,
    ) -> PredictionResult:
        """
        预测bulk结构的表面暴露情况

        Parameters
        ----------
        structure : str or pymatgen.Structure or ase.Atoms
            输入的bulk结构，可以是:
            - POSCAR/VASP格式文件路径
            - pymatgen Structure对象
            - ASE Atoms对象
        crystal_id : str, default="crystal"
            晶体ID，用于标识输出文件
        top_n : int, default=5
            返回最可能暴露的前N个表面
        save_wulff_shape : bool, default=True
            是否保存Wulff形状图
        output_dir : str, optional
            输出目录。如果不指定，使用工作目录

        Returns
        -------
        PredictionResult
            预测结果，包含所有表面的信息和Wulff形状
        """
        # 设置输出目录
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, crystal_id)
        os.makedirs(output_dir, exist_ok=True)

        # 创建子目录
        crystal_dir = os.path.join(output_dir, "crystal")
        slab_dir = os.path.join(output_dir, "slabs")
        bulk_dir = os.path.join(output_dir, "bulks")
        lmdb_dir = os.path.join(output_dir, "lmdb")
        traj_dir = os.path.join(output_dir, "traj")
        wulff_dir = os.path.join(output_dir, "wulff") if save_wulff_shape else None

        for d in [crystal_dir, slab_dir, bulk_dir, lmdb_dir, traj_dir]:
            os.makedirs(d, exist_ok=True)

        try:
            # 1. 转换并保存输入结构
            self._log(f"\n{'='*60}")
            self._log(f"SurFF Prediction for: {crystal_id}")
            self._log(f"{'='*60}")

            pmg_structure = self._convert_structure(structure)
            formula = pmg_structure.formula

            # 保存晶体结构
            crystal_path = os.path.join(crystal_dir, crystal_id)
            pmg_structure.to(filename=crystal_path, fmt='poscar')

            self._log(f"Input structure: {formula}")
            self._log(f"Atoms: {len(pmg_structure)}")

            # 2. 生成slab结构
            slab_info = self._generate_slabs(pmg_structure, crystal_id, slab_dir, bulk_dir)

            # 3. 生成LMDB数据集
            lmdb_path = os.path.join(lmdb_dir, "relaxation.lmdb")
            self._generate_lmdb(slab_dir, lmdb_path, slab_info)

            # 4. 运行弛豫
            self._run_relaxation(lmdb_dir, traj_dir)

            # 5. 收集结果
            results_df = self._collect_results(traj_dir)

            # 6. 计算Wulff构型
            wulff_results, wulff_shape_path = self._calculate_wulff(
                results_df, slab_info, crystal_dir, crystal_id, wulff_dir
            )

            # 7. 构建返回结果
            surfaces = []
            for i in range(len(wulff_results['miller_indices'])):
                area = wulff_results['area_fractions'][i]
                if area > 0.1:
                    level = "high"
                elif area > 0.01:
                    level = "medium"
                else:
                    level = "low"

                surfaces.append(SurfaceResult(
                    slab_id=wulff_results['slab_ids'][i],
                    miller_index=wulff_results['miller_indices'][i],
                    shift=wulff_results['shifts'][i],
                    surface_energy=wulff_results['surface_energies'][i],
                    area_fraction=area,
                    exposure_level=level,
                ))

            result = PredictionResult(
                crystal_id=crystal_id,
                formula=formula,
                total_slabs=len(surfaces),
                surfaces=surfaces,
                wulff_shape_path=wulff_shape_path,
                raw_results=wulff_results,
            )

            self._log(f"\n{result.summary(top_n)}")

            return result

        finally:
            # 清理临时文件（如果需要）
            if not self.keep_files and self._temp_dir and output_dir.startswith(self._temp_dir):
                pass  # 保留输出目录的结果

    def predict_batch(
        self,
        structures: List[Union[str, "pymatgen.core.Structure", "ase.Atoms"]],
        crystal_ids: Optional[List[str]] = None,
        top_n: int = 5,
        save_wulff_shape: bool = True,
        output_dir: Optional[str] = None,
    ) -> List[PredictionResult]:
        """
        批量预测多个bulk结构的表面暴露情况

        Parameters
        ----------
        structures : list
            输入的bulk结构列表
        crystal_ids : list of str, optional
            晶体ID列表。如果不指定，将自动生成
        top_n : int, default=5
            返回最可能暴露的前N个表面
        save_wulff_shape : bool, default=True
            是否保存Wulff形状图
        output_dir : str, optional
            输出目录

        Returns
        -------
        list of PredictionResult
            每个结构的预测结果列表
        """
        if crystal_ids is None:
            crystal_ids = [f"crystal_{i}" for i in range(len(structures))]
        elif len(crystal_ids) != len(structures):
            raise ValueError("Length of crystal_ids must match length of structures")

        results = []
        for i, (structure, crystal_id) in enumerate(zip(structures, crystal_ids)):
            self._log(f"\nProcessing {i+1}/{len(structures)}: {crystal_id}")
            try:
                result = self.predict(
                    structure, crystal_id, top_n, save_wulff_shape, output_dir
                )
                results.append(result)
            except Exception as e:
                self._log(f"Error processing {crystal_id}: {e}")
                results.append(None)

        return results

    def __del__(self):
        """清理临时目录"""
        if hasattr(self, '_temp_dir') and self._temp_dir and not self.keep_files:
            try:
                shutil.rmtree(self._temp_dir)
            except:
                pass


# 便捷函数
def predict_surface_exposure(
    structure: Union[str, "pymatgen.core.Structure", "ase.Atoms"],
    surff_root: str,
    checkpoint_path: Optional[str] = None,
    top_n: int = 5,
    use_gpu: bool = True,
    output_dir: Optional[str] = None,
    **kwargs,
) -> PredictionResult:
    """
    便捷函数：预测bulk结构的表面暴露情况

    Parameters
    ----------
    structure : str or pymatgen.Structure or ase.Atoms
        输入的bulk结构
    surff_root : str
        SurFF项目根目录
    checkpoint_path : str, optional
        模型checkpoint路径
    top_n : int, default=5
        返回最可能暴露的前N个表面
    use_gpu : bool, default=True
        是否使用GPU
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数传递给SurFFPredictor

    Returns
    -------
    PredictionResult
        预测结果
    """
    predictor = SurFFPredictor(
        surff_root=surff_root,
        checkpoint_path=checkpoint_path,
        use_gpu=use_gpu,
        keep_files=True if output_dir else False,
        **kwargs,
    )
    return predictor.predict(structure, top_n=top_n, output_dir=output_dir)


if __name__ == "__main__":
    # 示例用法
    import argparse

    parser = argparse.ArgumentParser(description="SurFF Surface Exposure Predictor")
    parser.add_argument("structure", help="Path to structure file (POSCAR format)")
    parser.add_argument("--surff-root", default="deps/SurFF",
                       help="Path to SurFF root directory")
    parser.add_argument("--checkpoint", default=None, help="Path to model checkpoint")
    parser.add_argument("--top-n", type=int, default=5, help="Number of top surfaces to show")
    parser.add_argument("--cpu", action="store_true", help="Use CPU instead of GPU")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--crystal-id", default="crystal", help="Crystal ID")

    args = parser.parse_args()

    predictor = SurFFPredictor(
        surff_root=args.surff_root,
        checkpoint_path=args.checkpoint,
        use_gpu=not args.cpu,
        keep_files=True,
    )

    result = predictor.predict(
        structure=args.structure,
        crystal_id=args.crystal_id,
        top_n=args.top_n,
        output_dir=args.output_dir,
    )

    # 打印结果
    print("\n" + "="*60)
    print("PREDICTION RESULTS")
    print("="*60)
    print(result.summary(args.top_n))

    # 保存CSV
    if args.output_dir:
        csv_path = os.path.join(args.output_dir, f"{args.crystal_id}_results.csv")
        result.to_dataframe().to_csv(csv_path, index=False)
        print(f"\nResults saved to: {csv_path}")
