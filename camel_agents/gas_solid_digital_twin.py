"""
Gas-Solid Catalysis Digital Twin - 气固界面催化数字孪生

这个模块提供了完整的气固界面催化反应模拟流程，从bulk结构到最终的KMC模拟。

完整流程：
1. 从bulk结构生成最可能暴露的表面 (SurFF)
2. 在每个表面上预测初始反应物的吸附位点 (AdsorbDiff)
3. 模拟表面重构 (VSSR-MC, 不考虑电势)
4. 计算反应路径的自由能和能垒 (Pathway/Fairchem)
5. 进行动力学蒙特卡洛模拟 (KMC/CatMAP)

使用方式:
    from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin

    # 初始化
    dt = GasSolidDigitalTwin(
        surff_root="deps/SurFF",
        adsorbdiff_root="deps/AdsorbDiff",
        surface_sampling_root="deps/surface-sampling",
        fairchem_root="deps/fairchem",
    )

    # 运行完整流程
    result = dt.run_complete_workflow(
        bulk_structure="POSCAR",
        reaction_intermediates=["*O", "*OH", "*OOH", "*"],
        initial_reactant="*O",
        temperature=500,  # K
        output_dir="output/gas_solid",
    )

    print(result.summary())

作者: Claude
日期: 2026-01-20
"""

import os
import logging
import tempfile
import shutil
from fractions import Fraction
from math import lcm
from typing import Union, Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
import warnings

import numpy as np
from ase import Atoms
from ase.io import read, write

# 导入各个模块的预测器
from core.surface.surff_predictor import SurFFPredictor, PredictionResult as SurFFResult
from core.reconstruction.adsorbdiff_predictor import AdsorbDiffPredictor, PredictionOutput as AdsorbDiffOutput
from core.reconstruction.adsorbml_predictor import AdsorbMLPredictor
from core.reconstruction.vssr_mc_predictor import VSSRMCPredictor, VSSRMCResult
from core.pathway.pathway_predictor import PathwayPredictor, CompletePathwayResult


# =============================================================================
# 数据结构定义
# =============================================================================

@dataclass
class SurfaceAnalysisResult:
    """单个表面的完整分析结果"""
    miller_index: Tuple[int, int, int]
    surface_energy: float  # eV/Å²
    area_fraction: float  # 表面暴露比例
    pristine_surface: Atoms  # 初始表面
    adsorption_sites: AdsorbDiffOutput  # 吸附位点预测结果
    reconstructed_surface: Optional[VSSRMCResult] = None  # 重构后的表面
    pathway_analysis: Optional[CompletePathwayResult] = None  # 反应路径分析

    def __repr__(self):
        return (f"SurfaceAnalysis({self.miller_index}, "
                f"E_surf={self.surface_energy:.4f} eV/Å², "
                f"area={self.area_fraction:.2%})")


@dataclass
class GasSolidWorkflowResult:
    """气固界面催化数字孪生的完整结果"""
    bulk_formula: str
    reaction_intermediates: List[str]
    temperature: float  # K
    surff_result: SurFFResult  # SurFF表面预测结果
    surface_analyses: List[SurfaceAnalysisResult]  # 每个表面的分析结果
    best_surface: Optional[SurfaceAnalysisResult] = None  # 最优表面
    kmc_result: Optional[Any] = None  # KMC模拟结果
    output_dir: Optional[str] = None

    def summary(self, top_n: int = 3) -> str:
        """生成结果摘要"""
        lines = [
            "=" * 80,
            "Gas-Solid Catalysis Digital Twin Results",
            "=" * 80,
            f"Bulk material: {self.bulk_formula}",
            f"Temperature: {self.temperature} K",
            f"Reaction intermediates: {', '.join(self.reaction_intermediates)}",
            f"",
            f"Surface Analysis:",
            f"-" * 80,
        ]

        for i, surface in enumerate(self.surface_analyses[:top_n], 1):
            lines.append(f"")
            lines.append(f"Surface {i}: {surface.miller_index}")
            lines.append(f"  Surface energy: {surface.surface_energy:.4f} eV/Å²")
            lines.append(f"  Area fraction: {surface.area_fraction:.2%}")

            if surface.adsorption_sites:
                best_site = surface.adsorption_sites.best_result
                lines.append(f"  Best adsorption site: E = {best_site.energy:.3f} eV")

            if surface.reconstructed_surface:
                recon = surface.reconstructed_surface.lowest_energy_structure
                lines.append(f"  Reconstructed surface: E = {recon.energy:.3f} eV")

            if surface.pathway_analysis:
                pathway = surface.pathway_analysis
                lines.append(f"  Pathway analysis:")
                lines.append(f"    Overall ΔE: {pathway.overall_reaction_energy:.3f} eV")
                if pathway.max_barrier:
                    lines.append(f"    Maximum barrier: {pathway.max_barrier:.3f} eV")
                if pathway.rate_determining_step:
                    rds = pathway.rate_determining_step
                    lines.append(f"    RDS: {rds.name} (E_act = {rds.activation_energy:.3f} eV)")

        if self.best_surface:
            lines.append(f"")
            lines.append(f"Best surface for catalysis: {self.best_surface.miller_index}")

        lines.append("=" * 80)
        return "\n".join(lines)


# =============================================================================
# 主类
# =============================================================================

class GasSolidDigitalTwin:
    """
    气固界面催化数字孪生

    整合了从bulk结构到KMC模拟的完整工作流程。
    """

    def __init__(
        self,
        surff_root: str,
        adsorbdiff_root: str,
        surface_sampling_root: str,
        fairchem_root: str,
        surff_checkpoint: Optional[str] = None,
        adsorbdiff_checkpoint: Optional[str] = None,
        fairchem_model: str = "uma-s-1p1",
        fairchem_model_path: Optional[str] = None,
        vssr_mc_model: str = "CHGNetNFF",
        vssr_mc_device: str = "cuda",
        potential_she: Optional[float] = None,
        ph: Optional[float] = None,
        use_gpu: bool = True,
        logger: Optional[logging.Logger] = None,
        use_llm_controller: bool = False,
        llm_model: Optional[str] = None,
        adsorption_backend: str = "adsorbml",
        adsorbml_model: str = "uma-s-1p1",
        adsorbml_num_sites: int = 20,
        adsorbml_placement_mode: str = "random_site_heuristic_placement",
        adsorbml_interstitial_gap: float = 0.1,
        adsorbml_relax_steps: int = 200,
        adsorbml_relax_fmax: float = 0.05,
    ):
        """
        初始化气固界面催化数字孪生

        参数:
            surff_root: SurFF库的路径
            adsorbdiff_root: AdsorbDiff库的路径
            surface_sampling_root: surface-sampling库的路径
            fairchem_root: Fairchem库的路径
            surff_checkpoint: SurFF模型检查点路径（可选）
            adsorbdiff_checkpoint: AdsorbDiff模型检查点路径（可选）
            fairchem_model: Fairchem模型名称 ("uma-s-1p1" 或 "uma-m-1p1")
            fairchem_model_path: Fairchem模型本地路径（可选，如果提供则不从HF下载）
            vssr_mc_model: VSSR-MC使用的力场模型 ("CHGNetNFF" 或 "MACENFF")
            use_gpu: 是否使用GPU
            logger: 日志记录器
            use_llm_controller: 是否使用LLM控制NEB结构准备
            llm_model: LLM模型名称（如"claude-opus-4-5-20251101"）
        """
        self.surff_root = surff_root
        self.adsorbdiff_root = adsorbdiff_root
        self.surface_sampling_root = surface_sampling_root
        self.fairchem_root = fairchem_root
        self.use_gpu = use_gpu
        self.logger = logger or logging.getLogger(__name__)

        # 初始化各个预测器
        self.logger.info("Initializing SurFF predictor...")
        self.surff_predictor = SurFFPredictor(
            surff_root=surff_root,
            checkpoint_path=surff_checkpoint,
            use_gpu=use_gpu,
        )

        self.adsorption_backend = adsorption_backend.lower()
        if self.adsorption_backend not in ("adsorbdiff", "adsorbml"):
            raise ValueError(
                f"adsorption_backend must be 'adsorbdiff' or 'adsorbml', "
                f"got: {adsorption_backend}"
            )

        if self.adsorption_backend == "adsorbml":
            self.logger.info("Initializing AdsorbML predictor (fairchem backend)...")
            self.adsorbdiff_predictor = AdsorbMLPredictor(
                fairchem_model=adsorbml_model,
                fairchem_model_path=fairchem_model_path,
                use_gpu=use_gpu,
                num_sites=adsorbml_num_sites,
                placement_mode=adsorbml_placement_mode,
                interstitial_gap=adsorbml_interstitial_gap,
                relax_steps=adsorbml_relax_steps,
                relax_fmax=adsorbml_relax_fmax,
            )
        else:
            self.logger.info("Initializing AdsorbDiff predictor...")
            self.adsorbdiff_predictor = AdsorbDiffPredictor(
                adsorbdiff_root=adsorbdiff_root,
                checkpoint_path=adsorbdiff_checkpoint,
                use_gpu=use_gpu,
            )

        self.logger.info("Initializing VSSR-MC predictor...")
        self.vssr_mc_predictor = VSSRMCPredictor(
            surface_sampling_root=surface_sampling_root,
            model_type=vssr_mc_model,
            device=vssr_mc_device,
            potential_she=potential_she,
            ph=ph,
        )

        self.logger.info("Initializing Pathway predictor...")
        self.pathway_predictor = PathwayPredictor(
            fairchem_root=fairchem_root,
            model_name=fairchem_model,
            model_path=fairchem_model_path,
            use_gpu=use_gpu,
            use_llm_controller=use_llm_controller,
            llm_model=llm_model,
        )

        if use_llm_controller:
            self.logger.info(f"  LLM-controlled NEB enabled (model: {llm_model or 'default'})")

        self.logger.info("All predictors initialized successfully")

    def generate_surfaces_from_bulk(
        self,
        bulk_structure: Union[str, Atoms],
        top_n: int = 5,
        output_dir: Optional[str] = None,
    ) -> Tuple[SurFFResult, List[Atoms]]:
        """
        步骤1: 从bulk结构生成最可能暴露的表面

        参数:
            bulk_structure: bulk结构文件路径或ASE Atoms对象
            top_n: 生成前N个最可能暴露的表面
            output_dir: 输出目录

        返回:
            (SurFF预测结果, 表面slab列表)
        """
        self.logger.info(f"Step 1: Generating top {top_n} surfaces from bulk structure...")

        # 使用SurFF预测表面
        surff_result = self.surff_predictor.predict(
            structure=bulk_structure,
            top_n=top_n,
            output_dir=output_dir,
        )

        self.logger.info(f"SurFF identified {len(surff_result.surfaces)} surfaces")
        for i, surface in enumerate(surff_result.get_top_n(top_n), 1):
            self.logger.info(f"  {i}. {surface.miller_index}: "
                           f"E_surf = {surface.surface_energy:.4f} eV/Å², "
                           f"area = {surface.area_fraction:.2%}")

        # 生成表面slab结构
        # 注意: SurFF只预测表面能，需要手动生成slab结构
        # 这里我们使用ASE的surface模块生成标准slab
        if isinstance(bulk_structure, str):
            bulk = read(bulk_structure)
        else:
            bulk = bulk_structure

        slabs = []
        for surface in surff_result.get_top_n(top_n):
            miller = surface.miller_index
            slab = self._generate_slab_from_bulk(bulk, miller)
            slabs.append(slab)

            # 保存slab结构
            if output_dir:
                slab_file = os.path.join(output_dir, f"slab_{miller[0]}{miller[1]}{miller[2]}.vasp")
                write(slab_file, slab)
                self.logger.info(f"  Saved slab {miller} to {slab_file}")

        return surff_result, slabs

    def predict_adsorption_sites(
        self,
        surface: Atoms,
        adsorbate: str,
        num_sites: int = 10,
        output_dir: Optional[str] = None,
        min_adsorbate_distance: float = 3.0,
        max_expansion_attempts: int = 2,
    ) -> AdsorbDiffOutput:
        """
        步骤2: 在表面上预测初始反应物的吸附位点

        如果吸附后分子间距离过近（考虑周期性边界），自动扩大表面并重新吸附

        参数:
            surface: 表面slab结构
            adsorbate: 吸附物（SMILES格式或名称，如"*CO", "*O"）
            num_sites: 尝试的吸附位点数
            output_dir: 输出目录
            min_adsorbate_distance: 吸附原子最小允许距离 (Å)，考虑周期性
            max_expansion_attempts: 最大扩胞尝试次数

        返回:
            吸附位点预测结果
        """
        self.logger.info(f"Step 2: Predicting adsorption sites for {adsorbate}...")

        current_surface = surface.copy()
        n_surface_atoms = len(surface)  # 记录原始表面原子数
        expansion_factor = 1  # 追踪扩胞倍数

        for attempt in range(max_expansion_attempts + 1):
            # 进行吸附预测
            result = self.adsorbdiff_predictor.predict(
                surface=current_surface,
                adsorbate=adsorbate,
                num_samples=num_sites,
                output_dir=output_dir,
            )

            self.logger.info(f"  Found {len(result.results)} adsorption configurations")
            self.logger.info(f"  Best site: E = {result.best_result.energy:.3f} eV")

            # 检查最佳吸附构型中吸附原子的周期性距离
            best_structure = result.best_result.final_structure
            is_too_close, min_dist = self._check_adsorbate_periodic_distance(
                surface_before=current_surface,
                surface_after=best_structure,
                threshold=min_adsorbate_distance
            )

            if is_too_close:
                self.logger.warning(f"  ⚠️  Adsorbate atoms too close (periodic): min distance = {min_dist:.2f} Å (threshold = {min_adsorbate_distance:.2f} Å)")

                if attempt < max_expansion_attempts:
                    # 扩大表面并重新吸附
                    expansion_factor += 1
                    self.logger.info(f"  → Expanding surface to {expansion_factor}x{expansion_factor} supercell and re-adsorbing...")
                    from ase.build import make_supercell
                    P = np.diag([expansion_factor, expansion_factor, 1])
                    current_surface = make_supercell(surface, P)
                    n_surface_atoms = len(current_surface)  # 更新表面原子数
                    self.logger.info(f"  → Surface expanded from {len(surface)} to {len(current_surface)} atoms")
                else:
                    self.logger.warning(f"  ⚠️  Max expansion attempts ({max_expansion_attempts}) reached, using current result")
                    break
            else:
                # 距离合适，返回结果
                if expansion_factor > 1:
                    self.logger.info(f"  ✓ Adsorbate distance acceptable after {expansion_factor}x{expansion_factor} expansion")
                else:
                    self.logger.info(f"  ✓ Adsorbate distance acceptable")
                break

        return result

    def _check_adsorbate_periodic_distance(
        self,
        surface_before: Atoms,
        surface_after: Atoms,
        threshold: float = 3.0
    ) -> tuple:
        """
        检查吸附分子与其周期性镜像之间的距离

        通用逻辑：
        1. 对比吸附前后识别吸附原子（新增的原子）
        2. 计算吸附分子质心
        3. 检查cell大小，确保周期性镜像不会太近

        Parameters
        ----------
        surface_before : Atoms
            吸附前的表面
        surface_after : Atoms
            吸附后的结构
        threshold : float
            最小允许距离 (Å)

        Returns
        -------
        (is_too_close, min_distance)
        """
        n_surface = len(surface_before)
        n_total = len(surface_after)

        # 吸附的原子索引（新增的原子）
        ads_indices = list(range(n_surface, n_total))

        if len(ads_indices) == 0:
            return False, float('inf')

        # 获取吸附原子位置
        ads_pos = surface_after.get_positions()[ads_indices]
        cell = surface_after.cell

        # 对于吸附分子（无论单原子还是多原子），检查与周期性镜像的距离
        # 计算吸附分子的质心位置
        ads_center = ads_pos.mean(axis=0)

        # 最小周期性距离约为 min(cell_x, cell_y) / 2
        # 这是分子质心与其周期性镜像之间的最小可能距离
        cell_lengths = cell.lengths()
        min_cell_length = min(cell_lengths[0], cell_lengths[1])  # xy平面
        min_periodic_dist = min_cell_length / 2.0

        # Debug: 记录cell信息
        if min_periodic_dist < threshold:
            self.logger.debug(f"    Cell lengths: x={cell_lengths[0]:.2f} Å, y={cell_lengths[1]:.2f} Å, z={cell_lengths[2]:.2f} Å")
            self.logger.debug(f"    Adsorbate atoms: {len(ads_indices)}, Center: {ads_center}")

        is_too_close = min_periodic_dist < threshold
        return is_too_close, min_periodic_dist

    def simulate_surface_reconstruction(
        self,
        surface: Atoms,
        adsorbates: List[str],
        temperature: float = 500.0,  # K
        total_sweeps: int = 100,
        output_dir: Optional[str] = None,
        surface_indices: Optional[List[int]] = None,
        adsorbate_indices: Optional[List[int]] = None,
        clean_slab: Optional[Atoms] = None,
        num_adsorbates_override: Optional[int] = None,
        chem_pots: Optional[Dict[str, float]] = None,
        use_seed_for_virtual_sites: bool = False,
        adsorbate_exclusion_radius_A: Optional[float] = None,
        existing_atom_exclusion_radius_A: Optional[float] = None,
    ) -> VSSRMCResult:
        """
        步骤3: 使用VSSR-MC模拟表面重构（不考虑电势）

        参数:
            surface: 表面slab结构
            adsorbates: 可能的吸附物种列表（元素符号，如["O", "H"]）
            temperature: 温度 (K)
            total_sweeps: MC总扫描次数
            output_dir: 输出目录

        返回:
            VSSR-MC模拟结果
        """
        self.logger.info(f"Step 3: Simulating surface reconstruction at {temperature} K...")
        self.logger.info(f"  Adsorbates: {', '.join(adsorbates)}")

        # Convert temperature from K to kT units for MC sampling
        k_B = 8.617333262e-5  # eV/K
        mc_temperature = k_B * temperature / 0.025  # Normalize to ~1.0 at 300K
        mc_temperature = max(0.5, min(2.0, mc_temperature))
        self.logger.info(f"  MC temperature: {mc_temperature:.3f} kT")

        # Run VSSR-MC in a subprocess to avoid mpi4py/UCX state corruption
        # from other models (SurFF, AdsorbDiff) loaded in this process
        result = self._run_vssr_mc_subprocess(
            surface=surface,
            adsorbates=adsorbates,
            temperature=mc_temperature,
            total_sweeps=total_sweeps,
            output_dir=output_dir,
            surface_indices=surface_indices,
            adsorbate_indices=adsorbate_indices,
            clean_slab=clean_slab,
            num_adsorbates_override=num_adsorbates_override,
            chem_pots=chem_pots,
            use_seed_for_virtual_sites=use_seed_for_virtual_sites,
            adsorbate_exclusion_radius_A=adsorbate_exclusion_radius_A,
            existing_atom_exclusion_radius_A=existing_atom_exclusion_radius_A,
        )

        self.logger.info(f"  Sampling complete: {len(result.structures)} structures sampled")
        self.logger.info(f"  Lowest energy: {result.lowest_energy_structure.energy:.3f} eV")

        return result

    def _run_vssr_mc_subprocess(
        self,
        surface: Atoms,
        adsorbates: List[str],
        temperature: float,
        total_sweeps: int,
        output_dir: Optional[str],
        surface_indices: Optional[List[int]] = None,
        adsorbate_indices: Optional[List[int]] = None,
        clean_slab: Optional[Atoms] = None,
        chem_pots: Optional[Dict[str, float]] = None,
        use_seed_for_virtual_sites: bool = False,
        num_adsorbates_override: Optional[int] = None,
        adsorbate_exclusion_radius_A: Optional[float] = None,
        existing_atom_exclusion_radius_A: Optional[float] = None,
    ) -> VSSRMCResult:
        """
        Run VSSR-MC in a subprocess to isolate it from GPU/MPI state corruption.
        """
        import json
        import pickle
        import subprocess
        import tempfile

        self.logger.info("  Running VSSR-MC in isolated subprocess...")

        # Create temporary files for communication
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            surface_path = f.name
            pickle.dump(surface, f)

        # Determine which atoms are pre-existing adsorbates:
        #   1. If adsorbate_indices was explicitly provided (including an
        #      empty list), trust it.
        #   2. If surface_indices is provided but adsorbate_indices is None,
        #      assume no adsorbates (LLM may have dropped an empty list
        #      arg; the explicit surface_indices signals an intentional
        #      clean-slab call).
        #   3. Otherwise fall back to the element-composition heuristic,
        #      which works for single-element slabs but misclassifies
        #      minority-element atoms (e.g. Ga in Ni5Ga3) on alloys.
        if adsorbate_indices is not None:
            detected_ads_indices = list(adsorbate_indices)
            self.logger.info(f"  Using provided adsorbate indices: {detected_ads_indices}")
        elif surface_indices is not None:
            detected_ads_indices = []
            self.logger.info(
                "  surface_indices provided without adsorbate_indices; "
                "assuming clean slab (no pre-existing adsorbates)"
            )
        else:
            surface_elements = set(self.catalyst_elements) if hasattr(self, 'catalyst_elements') else set()
            if not surface_elements:
                from collections import Counter
                element_counts = Counter(surface.get_chemical_symbols())
                surface_elements = {element_counts.most_common(1)[0][0]}
            detected_ads_indices = [
                i for i, symbol in enumerate(surface.get_chemical_symbols())
                if symbol not in surface_elements
            ]
            self.logger.info(f"  Detected adsorbate indices (heuristic): {detected_ads_indices}")

        num_adsorbate_atoms = len(detected_ads_indices)
        self.logger.info(f"  Adsorbate atoms: {num_adsorbate_atoms}")

        if clean_slab is not None:
            clean_slab_for_sites = clean_slab.copy()
        else:
            clean_slab_for_sites = surface.copy()
            del clean_slab_for_sites[detected_ads_indices]
        self.logger.info(
            f"  Clean slab: {clean_slab_for_sites.get_chemical_formula()} "
            f"({len(clean_slab_for_sites)} atoms)"
        )

        # 保存clean_slab
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            clean_slab_path = f.name
            pickle.dump(clean_slab_for_sites, f)

        # Determine the VSSR-MC canonical count. This count refers to atoms
        # occupying VSSR virtual adsorption sites (e.g. Ti/O overlayer atoms),
        # not molecular adsorbate atoms such as C3H8 from Agent2. Therefore
        # only an explicit caller override may trigger canonical mode here.
        # adsorbate_indices are used above to remove/protect the molecule when
        # generating clean_slab and constraints; they must not become a Ti/O
        # fixed-count target.
        final_num_ads = (
            int(num_adsorbates_override)
            if num_adsorbates_override is not None
            else 0
        )
        use_canonical = final_num_ads > 0
        # When use_seed_for_virtual_sites=True, skip passing the clean_slab
        # so VSSRMCPredictor enumerates virtual sites on the seed (with
        # overlayer). Otherwise it enumerates on bare metal and overcounts.
        effective_clean_slab_path = None if use_seed_for_virtual_sites else clean_slab_path
        config = {
            "surface_sampling_root": self.vssr_mc_predictor.surface_sampling_root,
            "model_type": getattr(self.vssr_mc_predictor, "model_type", "CHGNetNFF"),
            "device": "cuda",
            "verbose": True,
            "keep_files": True,
            "adsorbates": adsorbates,
            "canonical": use_canonical,
            "num_adsorbates": final_num_ads,
            "total_sweeps": total_sweeps,
            "sweep_size": 20,
            "temperature": temperature,
            "output_dir": output_dir,
            "clean_slab_path": effective_clean_slab_path,
        }
        if adsorbate_exclusion_radius_A is not None:
            config["adsorbate_exclusion_radius_A"] = float(adsorbate_exclusion_radius_A)
        if existing_atom_exclusion_radius_A is not None:
            config["existing_atom_exclusion_radius_A"] = float(existing_atom_exclusion_radius_A)
        if chem_pots:
            # mcmc calculator expects chem_pots for EVERY element in the slab;
            # user config may only supply overrides. Fill defaults to 0.
            all_elements = set(surface.get_chemical_symbols()) | set(adsorbates or [])
            full_chem_pots = {el: 0.0 for el in all_elements}
            full_chem_pots.update({k: float(v) for k, v in chem_pots.items()})
            config["chem_pots"] = full_chem_pots

        # Pass surface/adsorbate indices for proper constraint handling
        if surface_indices is not None:
            config["surface_indices"] = list(surface_indices)
        if adsorbate_indices is not None:
            config["adsorbate_indices"] = list(adsorbate_indices)

        # Electrochemical parameters (liquid-solid Pourbaix mode)
        if getattr(self.vssr_mc_predictor, "potential_she", None) is not None:
            config["potential_she"] = self.vssr_mc_predictor.potential_she
        if getattr(self.vssr_mc_predictor, "ph", None) is not None:
            config["ph"] = self.vssr_mc_predictor.ph

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False, mode="w") as f:
            config_path = f.name
            json.dump(config, f)

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            output_path = f.name

        # Run subprocess
        current_dir = Path(__file__).resolve().parent
        project_root = current_dir.parent
        candidates = [
            current_dir / "reconstruction" / "run_vssr_mc_subprocess.py",
            project_root / "core" / "reconstruction" / "run_vssr_mc_subprocess.py",
        ]
        script_path = next((str(path) for path in candidates if path.exists()), None)
        if script_path is None:
            candidate_text = ", ".join(str(path) for path in candidates)
            raise FileNotFoundError(
                f"Cannot locate run_vssr_mc_subprocess.py. Tried: {candidate_text}"
            )

        try:
            import sys as _sys
            proc = subprocess.run(
                [
                    _sys.executable,
                    script_path,
                    "--surface", surface_path,
                    "--config", config_path,
                    "--output", output_path,
                ],
                capture_output=True,
                text=True,
            )

            # Log subprocess output
            if proc.stdout:
                for line in proc.stdout.strip().split("\n"):
                    self.logger.info(f"  [subprocess] {line}")
            if proc.stderr:
                for line in proc.stderr.strip().split("\n"):
                    if "INFO" in line or "WARNING" in line:
                        self.logger.debug(f"  [subprocess] {line}")
                    elif line.strip():
                        self.logger.warning(f"  [subprocess stderr] {line}")

            if proc.returncode != 0:
                raise RuntimeError(f"VSSR-MC subprocess failed with code {proc.returncode}")

            # Load result
            with open(output_path, "rb") as f:
                output_data = pickle.load(f)

            if not output_data.get("success", False):
                error = output_data.get("error", "Unknown error")
                tb = output_data.get("traceback", "")
                raise RuntimeError(f"VSSR-MC failed: {error}\n{tb}")

            # Reconstruct VSSRMCResult
            from core.reconstruction.vssr_mc_predictor import VSSRMCResult, SampledStructure

            def make_structure(s):
                return SampledStructure(
                    atoms=s["atoms"],
                    energy=s["energy"],
                    sweep_number=s["sweep_number"],
                    num_adsorbates=s["num_adsorbates"],
                    acceptance_rate=s["acceptance_rate"],
                )

            structures = [make_structure(s) for s in output_data["structures"]]
            best_structure = make_structure(output_data["best_structure"])
            lowest_structure = make_structure(output_data["lowest_energy_structure"])

            return VSSRMCResult(
                surface_name=output_data.get("surface_name", "unknown"),
                model_type=output_data.get("model_type", "CHGNetNFF"),
                total_sweeps=output_data.get("total_sweeps", 50),
                sweep_size=output_data.get("sweep_size", 20),
                temperature=output_data.get("temperature", 1.0),
                canonical=output_data.get("canonical", True),
                structures=structures,
                energy_history=output_data.get("energy_history", []),
                acceptance_history=output_data.get("acceptance_history", []),
                adsorption_count_history=output_data.get("adsorption_count_history", []),
                best_structure=best_structure,
                lowest_energy_structure=lowest_structure,
                offset_data=output_data.get("offset_data", {}),
                output_dir=output_data.get("output_dir", output_dir),
            )

        finally:
            # Cleanup temporary files
            for path in [surface_path, config_path, output_path]:
                try:
                    os.unlink(path)
                except Exception:
                    pass

    def analyze_reaction_pathway(
        self,
        surface: Atoms,
        intermediates: List[str],
        calculate_barriers: bool = True,
        num_sites: int = 5,
        n_frames: int = 10,
        fmax: float = 0.05,
        output_dir: Optional[str] = None,
        current_adsorbate_indices: Optional[List[int]] = None,
    ) -> CompletePathwayResult:
        """
        步骤4: 计算反应路径的自由能和能垒

        参数:
            surface: 表面slab结构（可以是重构后的，可能包含吸附物）
            intermediates: 反应中间体列表（按反应顺序，如["*O", "*OH", "*OOH", "*"]）
            calculate_barriers: 是否计算能垒（使用NEB）
            num_sites: 每个中间体尝试的吸附位点数
            n_frames: NEB计算的帧数
            fmax: NEB力收敛标准 (eV/Å)
            output_dir: 输出目录
            current_adsorbate_indices: 当前吸附物的原子索引列表（如果表面已有吸附物）

        返回:
            完整的反应路径分析结果
        """
        self.logger.info(f"Step 4: Analyzing reaction pathway...")
        self.logger.info(f"  Intermediates: {' -> '.join(intermediates)}")

        result = self.pathway_predictor.predict_pathway(
            surface=surface,
            adsorbates=intermediates,
            calculate_barriers=calculate_barriers,
            num_sites=num_sites,
            n_frames=n_frames,
            fmax=fmax,
            output_dir=output_dir,
            current_adsorbate_indices=current_adsorbate_indices,
        )

        self.logger.info(f"  Pathway analysis complete")
        self.logger.info(f"  Overall ΔE: {result.overall_reaction_energy:.3f} eV")
        if result.max_barrier:
            self.logger.info(f"  Maximum barrier: {result.max_barrier:.3f} eV")

        return result

    def run_kmc_simulation(
        self,
        pathway_result: CompletePathwayResult,
        temperature: float,
        pressures: Dict[str, float],
        time_end: float = 1.0,
        output_dir: Optional[str] = None,
    ) -> Any:
        """
        步骤5: 进行动力学蒙特卡洛模拟

        参数:
            pathway_result: 反应路径分析结果
            temperature: 温度 (K)
            pressures: 气相物种分压 (bar)，如{"CO": 1.0, "O2": 0.21}
            time_end: 模拟时间 (s)
            output_dir: 输出目录

        返回:
            KMC模拟结果
        """
        self.logger.info(f"Step 5: Running KMC simulation...")
        self.logger.info(f"  Temperature: {temperature} K")
        self.logger.info(f"  Pressures: {pressures}")

        try:
            # 导入CatMAP预测器
            from core.kmc.catmap_predictor import CatMAPPredictor

            # 创建CatMAP预测器
            predictor = CatMAPPredictor(work_dir=output_dir or tempfile.mkdtemp())

            # 从pathway_result构建反应机理
            self._build_catmap_mechanism(predictor, pathway_result, pressures)

            # Gas mode: prefer ideal_gas if every MOLECULAR (non-monatomic)
            # gas species has frequencies + ideal-gas parameters. Monatomic
            # spectator species (e.g. C_g, O_g used only for element balance)
            # have no vibrational modes by definition — they're checked
            # separately via geometry=='monatomic'. Without entropy, TOF is
            # biased upward by ~10^7 for reactants with significant gas-phase
            # translational/rotational entropy.
            gas_mode_default = 'frozen_gas'
            try:
                gas_iter = list(getattr(predictor, "gases", {}).values())

                def _is_ready(g) -> bool:
                    geom = getattr(g, "geometry", None)
                    sym = getattr(g, "symmetry_number", None)
                    if geom == 'monatomic':
                        # Single atom: only translational entropy; no freqs needed.
                        return sym is not None
                    return (
                        bool(getattr(g, "frequencies", None))
                        and sym is not None
                        and geom is not None
                    )

                if gas_iter and all(_is_ready(g) for g in gas_iter):
                    gas_mode_default = 'ideal_gas'
            except Exception:
                pass
            # Adsorbate mode: harmonic requires per-adsorbate frequencies,
            # which aren't computed by default. Stick with frozen_adsorbate
            # unless the pathway carries them.
            adsorbate_mode_default = 'frozen_adsorbate'
            try:
                has_any_freq = any(
                    bool(getattr(v, "frequencies", None))
                    for v in getattr(predictor, "adsorbates", {}).values()
                )
                if has_any_freq:
                    adsorbate_mode_default = 'harmonic_adsorbate'
            except Exception:
                pass
            self.logger.info(
                "CatMAP thermo modes: gas=%s, adsorbate=%s",
                gas_mode_default, adsorbate_mode_default,
            )
            predictor.set_thermo_modes(
                gas_mode=gas_mode_default,
                adsorbate_mode=adsorbate_mode_default,
            )

            # 设置条件
            predictor.set_conditions(temperature=temperature, pressures=pressures)

            # 运行单点计算（针对当前表面）
            surface_name = pathway_result.surface_formula
            result = predictor.run_single_point(surface=surface_name, output_dir=output_dir)

            self.logger.info(f"  KMC simulation complete")
            self.logger.info(f"  Turnover frequency: {result.get_turnover_frequency():.2e} 1/s")

            return result

        except ImportError as e:
            self.logger.warning(f"Could not import CatMAP: {e}")
            self.logger.warning("KMC simulation skipped. Please ensure CatMAP is installed.")
            return None
        except Exception as e:
            self.logger.error(f"KMC simulation failed: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _build_catmap_mechanism(
        self,
        predictor: Any,
        pathway_result: CompletePathwayResult,
        pressures: Dict[str, float],
    ):
        """
        从pathway结果构建CatMAP反应机理

        参数:
            predictor: CatMAPPredictor实例
            pathway_result: 反应路径分析结果
            pressures: 气相物种分压
        """
        from collections import Counter

        def get_composition(formula: str) -> Counter:
            """从化学式获取元素组成"""
            import re
            # 处理带*号的吸附物种名称
            formula = formula.replace('*', '')
            if not formula:
                return Counter()
            # 解析化学式
            pattern = r'([A-Z][a-z]?)(\d*)'
            composition = Counter()
            for match in re.finditer(pattern, formula):
                element = match.group(1)
                count = int(match.group(2)) if match.group(2) else 1
                composition[element] += count
            return composition

        def get_element_change(reactant_formula: str, product_formula: str) -> dict:
            """计算从反应物到产物的元素变化"""
            reactant_comp = get_composition(reactant_formula)
            product_comp = get_composition(product_formula)
            change = {}
            all_elements = set(reactant_comp.keys()) | set(product_comp.keys())
            for elem in all_elements:
                diff = product_comp[elem] - reactant_comp[elem]
                if diff != 0:
                    change[elem] = diff
            return change

        # Auto-derive gas compositions from pressures dict + pathway species.
        gas_compositions: Dict[str, Dict[str, int]] = {}
        for gas_key in pressures:
            clean = gas_key.replace("_g", "")
            comp = get_composition(clean)
            if comp:
                gas_compositions[clean] = dict(comp)
        # Ensure common diatomic sources exist for element balancing
        for gas in ["H2", "O2", "N2"]:
            if gas not in gas_compositions:
                gas_compositions[gas] = dict(get_composition(gas))
        gas_compositions.setdefault("C", {"C": 1})

        # Normalize pressure keys: strip _g suffix for matching
        pressure_names = set()
        for k in pressures:
            pressure_names.add(k)
            pressure_names.add(k.replace("_g", ""))

        def identify_intact_surface_exchange(reactant_clean: str, product_clean: str) -> Optional[Tuple[str, str]]:
            """Detect pure adsorption/desorption (no chemical change).

            An "intact surface exchange" is a step where the molecule that
            enters/exits the surface is chemically identical on both sides.
            A step like *C3H7 -> *C3H6 shares no atoms of identical formula,
            so it is NOT an exchange — it's a dehydrogenation and must go
            through the element-balance path so that H bookkeeping is correct.
            """
            if not reactant_clean and product_clean and (
                product_clean in gas_compositions or product_clean in pressure_names
            ):
                return ("adsorption", product_clean)
            if not product_clean and reactant_clean and (
                reactant_clean in gas_compositions or reactant_clean in pressure_names
            ):
                return ("desorption", reactant_clean)
            # Adsorbate ↔ gas of the SAME molecule (composition-preserving transition)
            if (
                reactant_clean
                and product_clean
                and reactant_clean == product_clean
                and (product_clean in gas_compositions or product_clean in pressure_names)
            ):
                return ("desorption", product_clean)
            return None

        # 收集路径元素：active 表示该元素在步骤间发生了净增减，
        # conserved 表示路径中存在但总是守恒。
        pathway_elements = set()
        active_elements = set()
        for step in pathway_result.steps:
            reactant_clean = step.reactant_adsorbate.replace('*', '')
            product_clean = step.product_adsorbate.replace('*', '')
            pathway_elements.update(get_composition(step.reactant_adsorbate).keys())
            pathway_elements.update(get_composition(step.product_adsorbate).keys())
            if identify_intact_surface_exchange(reactant_clean, product_clean) is not None:
                continue
            active_elements.update(get_element_change(step.reactant_adsorbate, step.product_adsorbate).keys())
        conserved_elements = set(pathway_elements) - set(active_elements)

        # 由用户输入分压估计元素活度
        element_activity: Dict[str, float] = {}
        for gas_name, pressure in pressures.items():
            comp = gas_compositions.get(gas_name, get_composition(gas_name))
            if not comp:
                continue
            p = max(float(pressure), 1e-20)
            for elem, count in comp.items():
                n = max(int(count), 1)
                activity = p ** (1.0 / float(n))
                element_activity[elem] = max(element_activity.get(elem, 0.0), float(activity))

        for elem in pathway_elements:
            element_activity.setdefault(elem, 1e-12)

        # active 元素用于真实计量修正；conserved 元素仅作为"旁观气体"参考，
        # 避免 CatMAP 在缺失元素参考时报错。
        preferred_active_sources = {
            'H': 'H2',
            'O': 'O2',
            'N': 'N2',
            'C': 'C',
        }
        preferred_conserved_sources = {
            'H': 'H2',
            'O': 'O2',
            'N': 'N2',
            'C': 'CO',
        }

        atomic_sources: Dict[str, str] = {}
        conserved_reference_gases: Dict[str, str] = {}
        effective_pressures: Dict[str, float] = dict(pressures)

        for elem in active_elements:
            source = preferred_active_sources.get(elem, elem)
            atomic_sources[elem] = source

            if source not in gas_compositions:
                gas_compositions[source] = {elem: 1}

            activity = max(float(element_activity.get(elem, 1e-12)), 1e-12)
            atoms_per_mol = int(gas_compositions[source].get(elem, 1) or 1)
            source_pressure = activity ** float(atoms_per_mol)

            prev_pressure = max(float(effective_pressures.get(source, 0.0)), 0.0)
            effective_pressures[source] = max(prev_pressure, source_pressure, 1e-12)

        for elem in conserved_elements:
            source = preferred_conserved_sources.get(elem, elem)
            conserved_reference_gases[elem] = source

            if source not in gas_compositions:
                gas_compositions[source] = {elem: 1}

            activity = max(float(element_activity.get(elem, 1e-12)), 1e-12)
            atoms_per_mol = int(gas_compositions[source].get(elem, 1) or 1)
            source_pressure = activity ** float(atoms_per_mol)

            prev_pressure = max(float(effective_pressures.get(source, 0.0)), 0.0)
            effective_pressures[source] = max(prev_pressure, source_pressure, 1e-12)

        for step in pathway_result.steps:
            reactant_clean = step.reactant_adsorbate.replace('*', '')
            product_clean = step.product_adsorbate.replace('*', '')
            direct_exchange = identify_intact_surface_exchange(reactant_clean, product_clean)
            if direct_exchange is None:
                continue

            _, gas_name = direct_exchange
            if gas_name not in gas_compositions:
                gas_comp = get_composition(gas_name)
                if gas_comp:
                    gas_compositions[gas_name] = dict(gas_comp)
            effective_pressures[gas_name] = max(float(effective_pressures.get(gas_name, 0.0)), 1e-12)

        # 添加气相物种（包含用户输入分压和为质量守恒/元素参考补齐的库气体）
        # Pull vacuum-relaxed frequencies + geometry from the adsorption step
        # so CatMAP's ideal_gas thermo can compute TΔS(T). If an intermediate
        # label (e.g. '*C3H8') matches a gas name, reuse its vacuum data.
        gas_freq_map = dict(getattr(pathway_result, "gas_frequencies", {}) or {})
        gas_best_configs = dict(getattr(pathway_result, "best_configurations", {}) or {})

        # For gas species that appear in gas_pressures but are NOT among the
        # pathway adsorbates (e.g. H2, CH4 in a PDH feed), we still need vacuum
        # frequencies + geometry to activate ideal_gas thermo. Run those
        # on-demand here using the pathway predictor's UMA calculator.
        fp = None
        try:
            if hasattr(self, "pathway_predictor"):
                self.pathway_predictor._init_predictors()
                fp = getattr(self.pathway_predictor, "_fairchem_predictor", None)
        except Exception as exc:
            self.logger.warning("Could not initialise fairchem predictor for gas refs: %s", exc)

        def _resolve_gas_atoms(label: str) -> Optional[Atoms]:
            """Build an ASE Atoms object for a gas-phase molecule label.

            Tries ase.build.molecule first; falls back to the adsorbate
            resolver's SMILES logic. Returns None if unknown.
            """
            try:
                from ase.build import molecule as ase_molecule
                return ase_molecule(label)
            except (KeyError, NotImplementedError):
                pass
            try:
                from core.pathway.adsorbate_resolver import resolve_adsorbate
                atoms, _ = resolve_adsorbate(label)
                return atoms
            except Exception:
                return None

        def _lookup_gas_frequencies(gas_label: str) -> List[float]:
            """Return vibrational frequencies in eV (ASE/CatMAP convention).

            `gas_freq_map` stores values in cm^-1 (from fairchem_predictor);
            convert on lookup. Returns [] when no vacuum vibration is available.
            """
            for key in (gas_label, f"*{gas_label}", f"{gas_label}_g", f"*{gas_label}_g"):
                if key in gas_freq_map:
                    return [float(f) / 8065.54429 for f in gas_freq_map[key]]
            return []

        def _lookup_gas_geometry(gas_label: str):
            for key in (gas_label, f"*{gas_label}", f"{gas_label}_g", f"*{gas_label}_g"):
                cfg = gas_best_configs.get(key)
                if cfg is not None:
                    return cfg
            return None

        def _guess_ideal_gas_params(gas_label: str) -> Tuple[Optional[int], Optional[str], Optional[float]]:
            # Simple heuristics. CatMAP has a built-in table for common molecules
            # (see catmap.data.parameter_data.ideal_gas_params) — only provide
            # overrides here when we'd otherwise be missing from that table.
            builtin = {
                "H2": (2, "linear", 0),
                "N2": (2, "linear", 0),
                "O2": (2, "linear", 1.0),
                "CO": (1, "linear", 0),
                "CO2": (2, "linear", 0),
                "H2O": (2, "nonlinear", 0),
                "CH4": (12, "nonlinear", 0),
                "C3H6": (1, "nonlinear", 0),
                "C3H8": (18, "nonlinear", 0),  # NOT in CatMAP default table
                "C2H4": (4, "nonlinear", 0),
                "C2H6": (6, "nonlinear", 0),
                "NH3": (3, "nonlinear", 0),
            }
            params = builtin.get(gas_label)
            if params is not None:
                return params
            # Monatomic species (C, N, O, ...) are usually virtual element-balance
            # references emitted by _build_catmap_mechanism. They have no
            # vibrational or rotational DOF, only translation.
            comp = get_composition(gas_label)
            total_atoms = sum(comp.values()) if comp else 0
            if total_atoms == 1:
                return (1, "monatomic", 0)
            # Unknown polyatomic — fall back to symmetry=1, nonlinear, spin=0
            # (mild underestimation of rotational entropy, but correct order
            # of magnitude).
            return (1, "nonlinear", 0)

        for gas_name, pressure in effective_pressures.items():
            freq_ev = _lookup_gas_frequencies(gas_name)
            atoms_ref = _lookup_gas_geometry(gas_name)
            # If this gas species wasn't in the pathway adsorbates, launch an
            # on-demand vacuum relax + vibration calc so ideal_gas thermo has
            # everything it needs. Cache is shared across calls via fp's
            # _vacuum_energy_cache keyed by chemical formula.
            if fp is not None and (not freq_ev or atoms_ref is None):
                atoms_for_relax = atoms_ref if atoms_ref is not None else _resolve_gas_atoms(gas_name)
                if atoms_for_relax is not None:
                    try:
                        gas_out = fp._relax_adsorbate_in_vacuum(
                            atoms_for_relax, fmax=0.05, max_steps=200,
                            compute_vibrations=True,
                        )
                        if isinstance(gas_out, tuple):
                            _, freq_cm_new, relaxed_atoms = gas_out
                            if freq_cm_new and not freq_ev:
                                freq_ev = [float(f) / 8065.54429 for f in freq_cm_new]
                            if atoms_ref is None and relaxed_atoms is not None:
                                atoms_ref = relaxed_atoms
                            self.logger.info(
                                "Auto-relaxed sidecar gas %s for ideal_gas thermo "
                                "(%d modes)", gas_name, len(freq_ev),
                            )
                    except Exception as exc:
                        self.logger.warning(
                            "Sidecar vacuum relax for %s failed: %s — "
                            "will fall back to frozen_gas for this species.",
                            gas_name, exc,
                        )
            sym, geom, spin = _guess_ideal_gas_params(gas_name)
            predictor.add_gas(
                name=gas_name,
                formation_energy=0.0,
                pressure=max(float(pressure), 1e-12),
                frequencies=freq_ev,
                symmetry_number=sym,
                geometry=geom,
                spin=spin,
                atoms_geometry=atoms_ref,
            )

        # Build CatMAP reaction expressions.
        # CatMAP format: 'A_s + B_g -> C_s' (adsorbates use _s suffix, gases _g).
        # CatMAP REQUIRES INTEGER stoichiometric coefficients (see
        # catmap.model.expression_string_to_list line ~570 which calls int()).
        # When an elementary step releases or consumes a diatomic (H2, O2, N2, ...)
        # atom-by-atom, the natural coefficient is 1/atoms_per_mol (e.g. 0.5 H2_g).
        # To keep integer coefficients, we split such steps into:
        #   (a) an elementary step with an ATOMIC adsorbate (e.g. H_s)
        #   (b) a recombination/desorption step per element (e.g. 2H_s -> H2_g + 2*_s)
        # This is a general, reaction-agnostic transform: it activates for any gas
        # whose molecular formula has >1 atom of the exchanged element.
        atomic_intermediates_needed: Dict[str, str] = {}  # elem -> source_gas

        for step in pathway_result.steps:
            reactant_clean = step.reactant_adsorbate.replace('*', '')
            product_clean = step.product_adsorbate.replace('*', '')

            # CatMAP species names: adsorbate = X_s, empty site = *_s
            reactant_s = f"{reactant_clean}_s" if reactant_clean else "*_s"
            product_s = f"{product_clean}_s" if product_clean else "*_s"

            # Check for direct adsorption/desorption (adsorbate = gas species)
            direct_exchange = identify_intact_surface_exchange(reactant_clean, product_clean)
            if direct_exchange is not None:
                kind, gas_name = direct_exchange
                if kind == "adsorption":
                    predictor.add_reaction(f"*_s + {gas_name}_g -> {product_s}")
                else:
                    predictor.add_reaction(f"{reactant_s} -> *_s + {gas_name}_g")
                continue

            # Element change between adsorbates
            elem_change = get_element_change(step.reactant_adsorbate, step.product_adsorbate)

            # Build gas-phase / atomic-adsorbate terms for element balance
            reactant_parts = [reactant_s]
            product_parts = [product_s]
            extra_empty_sites = 0  # *_s needed when atoms move to atomic adsorbates

            for elem, change in elem_change.items():
                source_gas = atomic_sources.get(elem)
                if source_gas is None or source_gas not in gas_compositions:
                    continue
                atoms_per_mol = int(gas_compositions[source_gas].get(elem, 1) or 1)
                count = int(abs(change))
                if count == 0:
                    continue

                if atoms_per_mol <= 1 or (count % atoms_per_mol) == 0:
                    # Integer gas coefficient — emit gas directly
                    coeff = count // max(atoms_per_mol, 1)
                    gas_species = f"{source_gas}_g"
                    coeff_str = f"{coeff}{gas_species}" if coeff > 1 else gas_species
                    if change > 0:
                        reactant_parts.append(coeff_str)
                    else:
                        product_parts.append(coeff_str)
                else:
                    # Fractional gas coefficient — use atomic adsorbate intermediate
                    atomic_intermediates_needed[elem] = source_gas
                    atomic_token = f"{count}{elem}_s" if count > 1 else f"{elem}_s"
                    if change > 0:
                        reactant_parts.append(atomic_token)
                    else:
                        product_parts.append(atomic_token)
                        extra_empty_sites += count

            if extra_empty_sites > 0:
                reactant_parts.append(
                    f"{extra_empty_sites}*_s" if extra_empty_sites > 1 else "*_s"
                )

            reaction_expr = " + ".join(reactant_parts) + " -> " + " + ".join(product_parts)
            self.logger.info(f"CatMAP reaction: {reaction_expr}")
            predictor.add_reaction(reaction_expr)

        # Emit one recombination/desorption step per atomic intermediate so the
        # atoms exit to the gas phase. E.g. 2H_s -> H2_g + 2*_s.
        for elem, source_gas in atomic_intermediates_needed.items():
            atoms_per_mol = int(gas_compositions[source_gas].get(elem, 1) or 1)
            if atoms_per_mol < 2:
                continue
            lhs = f"{atoms_per_mol}{elem}_s"
            rhs = f"{source_gas}_g + {atoms_per_mol}*_s"
            recomb_expr = f"{lhs} -> {rhs}"
            self.logger.info(f"CatMAP recombination reaction: {recomb_expr}")
            predictor.add_reaction(recomb_expr)

        # --- Pathway closure: gas adsorption / desorption at chain endpoints ---
        # A microkinetic model is only solvable if every surface intermediate has
        # a gas-phase source OR sink. When the mechanism search produces a chain
        # of surface-surface steps (e.g. *C3H8 -> *C3H7 -> *C3H6), the first
        # adsorbate needs an adsorption step and any intermediate whose molecule
        # is a listed gas species should also be allowed to desorb (for selectivity
        # and mass balance). We add these only if the existing reaction list does
        # not already provide them (so user-authored adsorption steps are kept).
        existing_exprs = set(predictor.reactions)

        def _has_exchange(expr: str, mol: str) -> bool:
            """True if some reaction contains {mol}_g on either side."""
            token_g = f"{mol}_g"
            for e in existing_exprs:
                if token_g in e:
                    return True
            return False

        endpoint_mols = []
        if pathway_result.steps:
            first = pathway_result.steps[0].reactant_adsorbate.replace('*', '').strip()
            if first:
                endpoint_mols.append(first)
        # Collect every intermediate that matches a known gas species — each one
        # is a potential desorption branch (product selectivity channel).
        for step in pathway_result.steps:
            for label in (step.reactant_adsorbate, step.product_adsorbate):
                clean = label.replace('*', '').strip()
                if clean and (clean in gas_compositions or clean in pressure_names):
                    if clean not in endpoint_mols:
                        endpoint_mols.append(clean)

        surface_ads_set = {a.replace('*', '').strip() for a in pathway_result.adsorbate_energies}
        for mol in endpoint_mols:
            if not (mol in gas_compositions or mol in pressure_names):
                # Not declared as a gas species; skip.
                continue
            if mol not in surface_ads_set:
                # No adsorbed form in pathway — can't emit ads/des.
                continue
            if _has_exchange("", mol):
                continue
            # Reversible adsorption/desorption closes the microkinetic network
            # so every surface species has a gas-phase source or sink.
            rev_expr = f"*_s + {mol}_g <-> {mol}_s"
            self.logger.info(f"CatMAP pathway-closure reaction: {rev_expr}")
            predictor.add_reaction(rev_expr)
            existing_exprs.add(rev_expr)

        # Adsorbate formation energies for CatMAP.
        #
        # CatMAP wants one consistent set of formation energies (E_form) such
        # that for any reaction A_s -> B_s + X_g, the CatMAP rxn ΔE matches the
        # physically correct reaction energy. Raw DFT total energies cannot be
        # used (they're ~-700 eV, dwarfing every barrier). Nor can E_ads
        # referenced independently to each adsorbate's parent gas (open-shell
        # radicals like *C3H7 give meaningless E_ads).
        #
        # General approach: build a consistent set by propagating along the
        # reaction chain. Anchor the first adsorbate at its measured E_ads
        # (referenced to a stable closed-shell parent), then use NEB-derived
        # gas-phase-corrected step ΔE to compute the next. Atomic adsorbates
        # (H_s, O_s, ...) are anchored at 0 eV so their recombination steps
        # are thermoneutral with the gas; users can override by supplying
        # explicit E_ads values for "H", "O", etc.
        surface_name = pathway_result.surface_formula
        existing_adsorbate_labels = set()
        ads_binding_energies = getattr(pathway_result, "adsorption_energies", None) or {}

        DEFAULT_ATOMIC_E_S = 0.0  # eV — neutral placeholder
        atomic_formation = {}
        for elem in atomic_intermediates_needed:
            override = ads_binding_energies.get(elem) or ads_binding_energies.get(f"*{elem}")
            atomic_formation[elem] = float(override) if override is not None else DEFAULT_ATOMIC_E_S

        formation = {}

        def _label_key(label: str) -> str:
            return label.replace('*', '').strip()

        # Anchor: first adsorbate in the chain uses its measured binding energy.
        if pathway_result.steps:
            first_label = pathway_result.steps[0].reactant_adsorbate
            first_clean = _label_key(first_label)
            anchor = ads_binding_energies.get(first_label)
            if anchor is None:
                anchor = ads_binding_energies.get(first_clean)
            if anchor is None:
                # Fallback: if we really have no reference, use 0 eV (worst
                # case, puts absolute energies arbitrarily but keeps relative
                # reaction energetics correct within the chain).
                anchor = 0.0
            formation[first_clean] = float(anchor)

        # Propagate along steps using the NEB-derived step.reaction_energy
        # (already element-balance-corrected by
        # _compute_reaction_energy_with_correction in workflow.py).
        for step in pathway_result.steps:
            reactant_key = _label_key(step.reactant_adsorbate)
            product_key = _label_key(step.product_adsorbate)
            if reactant_key not in formation:
                # Branch whose parent hasn't been placed — skip, handled below.
                continue
            # Sum of atomic-adsorbate formations appearing on each side of the
            # CatMAP elementary step (which splits fractional gas into atomic
            # intermediates — see reaction-builder above). Δ(atomic) = (atoms
            # appearing on the product side) × E_form(atom_s) − (atoms on
            # reactant side) × E_form(atom_s).
            elem_change = get_element_change(step.reactant_adsorbate, step.product_adsorbate)
            atomic_delta = 0.0
            for elem, change in elem_change.items():
                e_atom = atomic_formation.get(elem, 0.0)
                source_gas = atomic_sources.get(elem)
                atoms_per_mol = 1
                if source_gas and source_gas in gas_compositions:
                    atoms_per_mol = int(gas_compositions[source_gas].get(elem, 1) or 1)
                count = int(abs(change))
                if count == 0 or atoms_per_mol <= 1 or (count % atoms_per_mol) == 0:
                    # Integer gas coefficient — no atomic adsorbate involved.
                    continue
                # Product released `count` atoms of `elem` to X_s (if change<0),
                # or consumed from X_s (if change>0). Solve:
                #   ΔE_rxn = E_form(P) + count*e_atom*(+1 if released) - E_form(R)
                if change < 0:
                    atomic_delta += count * e_atom
                else:
                    atomic_delta -= count * e_atom
            # Prefer the raw NEB ΔE (same atom count on both endpoints because
            # the staged atoms are integrated in the NEB trajectory) over
            # step.reaction_energy which applies a gas-phase correction that
            # assumes atoms exit to gas — incompatible with the split reactions
            # we emit. Fall back to reaction_energy when NEB data is absent.
            neb_rxn = getattr(step, "neb_reaction_energy", None)
            if neb_rxn is None:
                neb_rxn = step.reaction_energy
            if neb_rxn is not None and product_key not in formation:
                formation[product_key] = (
                    formation[reactant_key]
                    + float(neb_rxn)
                    - atomic_delta
                )

        # Any adsorbate that was not reached by propagation falls back to E_ads.
        for ads_name in pathway_result.adsorbate_energies:
            clean_name = _label_key(ads_name)
            if not clean_name or clean_name in formation:
                continue
            fallback = ads_binding_energies.get(ads_name)
            if fallback is None:
                fallback = ads_binding_energies.get(clean_name)
            if fallback is None:
                self.logger.warning(
                    "Adsorbate %s: no reference available — using 0 eV. Add an "
                    "explicit adsorption calculation for this species for a "
                    "meaningful TOF.",
                    clean_name,
                )
                fallback = 0.0
            formation[clean_name] = float(fallback)

        for clean_name, e_form in formation.items():
            existing_adsorbate_labels.add(clean_name)
            predictor.add_adsorbate(
                name=clean_name,
                formation_energies={surface_name: e_form},
            )

        # Register atomic intermediates that were introduced by the
        # fractional-coefficient split above. `atomic_formation[elem]` is already
        # set (default 0.0 eV, overridden if the user provides H/O/... in
        # adsorbate_energies). Only skip if something else has registered it.
        for elem, e_f in atomic_formation.items():
            if elem in existing_adsorbate_labels:
                continue
            if e_f == DEFAULT_ATOMIC_E_S:
                self.logger.warning(
                    "Atomic adsorbate %s_s: no explicit binding energy — "
                    "using placeholder E_f = %.2f eV. For quantitative TOF, "
                    "add an adsorption calculation of %s* to the pathway.",
                    elem, e_f, elem,
                )
            predictor.add_adsorbate(
                name=elem,
                formation_energies={surface_name: e_f},
            )

        # 添加过渡态（如果有能垒信息）— TS formation energy = E_form(reactant) + Ea_fwd
        # Use the propagated `formation[...]` from the step above (same frame as
        # the adsorbate formation_energies passed to CatMAP); using step.reactant_energy
        # (raw DFT) or ads_binding_energies (radical-referenced) would offset the
        # TS into a different frame and collapse the kinetics.
        for step in pathway_result.steps:
            if step.activation_energy is None or step.activation_energy <= 0.01:
                continue
            ts_name = f"{step.name}-TS"
            reactant_key = _label_key(step.reactant_adsorbate)
            reactant_form = formation.get(reactant_key)
            if reactant_form is None:
                self.logger.warning(
                    "Step %s: no formation energy for reactant %s; skipping TS.",
                    step.name, step.reactant_adsorbate,
                )
                continue
            ts_energy = float(reactant_form) + float(step.activation_energy)
            predictor.add_transition_state(
                name=ts_name,
                formation_energies={surface_name: ts_energy},
                scaling_mode='initial_state',
            )

        # 由 CatMAP 自动从 gas species 推断原子参考。
        # 显式 set_atomic_reservoirs 在不同 CatMAP 版本间存在兼容性差异，
        # 这里避免手动覆盖以提升稳定性。

    def run_complete_workflow(
        self,
        bulk_structure: Union[str, Atoms],
        reaction_intermediates: List[str],
        initial_reactant: str,
        temperature: float = 500.0,  # K
        top_n_surfaces: int = 3,
        num_adsorption_sites: int = 10,
        reconstruction_sweeps: int = 100,
        calculate_barriers: bool = True,
        neb_frames: int = 10,
        neb_fmax: float = 0.05,
        run_kmc: bool = False,
        pressures: Optional[Dict[str, float]] = None,
        output_dir: str = "output/gas_solid",
    ) -> GasSolidWorkflowResult:
        """
        运行完整的气固界面催化数字孪生流程

        参数:
            bulk_structure: bulk结构文件路径或ASE Atoms对象
            reaction_intermediates: 反应中间体列表（如["*O", "*OH", "*OOH", "*"]）
            initial_reactant: 初始反应物（如"*O"）
            temperature: 反应温度 (K)
            top_n_surfaces: 分析前N个最可能暴露的表面
            num_adsorption_sites: 每个表面尝试的吸附位点数
            reconstruction_sweeps: VSSR-MC扫描次数
            calculate_barriers: 是否计算反应能垒
            neb_frames: NEB计算的帧数
            neb_fmax: NEB力收敛标准 (eV/Å)
            run_kmc: 是否运行KMC模拟
            pressures: 气相物种分压（如果运行KMC）
            output_dir: 输出目录

        返回:
            完整的工作流程结果
        """
        self.logger.info("=" * 80)
        self.logger.info("Starting Gas-Solid Catalysis Digital Twin Workflow")
        self.logger.info("=" * 80)

        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)

        # 获取bulk结构的化学式
        if isinstance(bulk_structure, str):
            bulk = read(bulk_structure)
        else:
            bulk = bulk_structure
        bulk_formula = bulk.get_chemical_formula()

        # 步骤1: 生成表面
        surff_result, slabs = self.generate_surfaces_from_bulk(
            bulk_structure=bulk_structure,
            top_n=top_n_surfaces,
            output_dir=output_dir,
        )

        # For VSSR-MC surface reconstruction, we should only add/remove surface atoms
        # NOT the adsorbate molecule atoms (which would break the molecules)
        # Extract surface elements from the bulk structure
        surface_elements = list(set(bulk.get_chemical_symbols()))
        self.logger.info(f"Surface elements for VSSR-MC reconstruction: {surface_elements}")

        # 对每个表面进行完整分析
        surface_analyses = []
        for i, (surface_info, slab) in enumerate(zip(surff_result.get_top_n(top_n_surfaces), slabs)):
            self.logger.info("")
            self.logger.info(f"Analyzing surface {i+1}/{top_n_surfaces}: {surface_info.miller_index}")
            self.logger.info("-" * 80)

            # 创建表面专用输出目录
            surface_dir = os.path.join(output_dir, f"surface_{surface_info.miller_index[0]}{surface_info.miller_index[1]}{surface_info.miller_index[2]}")
            os.makedirs(surface_dir, exist_ok=True)

            # ============================================================
            # Step 3a: Reconstruct clean slab (before adsorption)
            # ============================================================
            if reconstruction_sweeps > 0:
                clean_sweeps = max(1, reconstruction_sweeps // 2)
                self.logger.info(f"Step 3a: Reconstructing clean slab ({clean_sweeps} sweeps)")
                recon_clean_result = self.simulate_surface_reconstruction(
                    surface=slab,
                    adsorbates=surface_elements,
                    temperature=temperature,
                    total_sweeps=clean_sweeps,
                    output_dir=os.path.join(surface_dir, "reconstruction_clean"),
                )
                slab = recon_clean_result.lowest_energy_structure.atoms
                self.logger.info(
                    f"  Clean slab after reconstruction: {slab.get_chemical_formula()} ({len(slab)} atoms)"
                )

            # ============================================================
            # Step 2: Adsorption (on reconstructed surface)
            # ============================================================
            adsorption_result = self.predict_adsorption_sites(
                surface=slab,
                adsorbate=initial_reactant,
                num_sites=num_adsorption_sites,
                output_dir=os.path.join(surface_dir, "adsorption"),
            )
            surface_with_adsorbate = adsorption_result.best_result.final_structure

            # ============================================================
            # Step 3b: Reconstruct with adsorbate present
            # ============================================================
            # Compute surface and adsorbate indices for proper constraint handling
            n_slab = len(slab)
            n_total = len(surface_with_adsorbate)
            pre_mc_adsorbate_indices = list(range(n_slab, n_total))
            slab_tags = slab.get_tags()
            if slab_tags.any():
                pre_mc_surface_indices = [j for j in range(n_slab) if slab_tags[j] == 1]
                if not pre_mc_surface_indices:
                    pre_mc_surface_indices = list(range(n_slab))
            else:
                z_coords = slab.positions[:, 2]
                z_unique = sorted(set(round(z, 1) for z in z_coords), reverse=True)
                top_2_z = set(z_unique[:2]) if len(z_unique) >= 2 else set(z_unique)
                pre_mc_surface_indices = [
                    j for j in range(n_slab) if round(z_coords[j], 1) in top_2_z
                ]
            self.logger.info(
                f"  Pre-MC indices: {len(pre_mc_surface_indices)} surface, "
                f"{len(pre_mc_adsorbate_indices)} adsorbate"
            )

            if reconstruction_sweeps > 0:
                reconstruction_result = self.simulate_surface_reconstruction(
                    surface=surface_with_adsorbate,
                    adsorbates=surface_elements,
                    temperature=temperature,
                    total_sweeps=reconstruction_sweeps,
                    output_dir=os.path.join(surface_dir, "reconstruction"),
                    surface_indices=pre_mc_surface_indices,
                    adsorbate_indices=pre_mc_adsorbate_indices,
                )
                # 使用重构后的最低能量结构（保留吸附物）
                reconstructed_surface = reconstruction_result.lowest_energy_structure.atoms
                # 识别吸附物原子索引，用于确定吸附位点
                adsorbate_indices = [
                    i for i, symbol in enumerate(reconstructed_surface.get_chemical_symbols())
                    if symbol not in surface_elements
                ]
            else:
                self.logger.info("Step 3: Skipping surface reconstruction (reconstruction_sweeps=0)")
                reconstruction_result = None
                # 直接使用吸附后的结构
                reconstructed_surface = surface_with_adsorbate
                # 识别吸附物原子索引
                adsorbate_indices = list(range(len(slab), len(reconstructed_surface)))

            # 步骤4: 反应路径分析（使用重构后的表面，通过替换吸附物进行）
            pathway_result = self.analyze_reaction_pathway(
                surface=reconstructed_surface,
                intermediates=reaction_intermediates,
                calculate_barriers=calculate_barriers,
                num_sites=5,  # 重构后使用较少的位点数以节省时间
                n_frames=neb_frames,
                fmax=neb_fmax,
                output_dir=os.path.join(surface_dir, "pathway"),
                current_adsorbate_indices=adsorbate_indices,  # 传递当前吸附物索引
            )

            # 保存表面分析结果
            analysis = SurfaceAnalysisResult(
                miller_index=surface_info.miller_index,
                surface_energy=surface_info.surface_energy,
                area_fraction=surface_info.area_fraction,
                pristine_surface=slab,
                adsorption_sites=adsorption_result,
                reconstructed_surface=reconstruction_result,
                pathway_analysis=pathway_result,
            )
            surface_analyses.append(analysis)

        # 选择最优表面（基于能垒最低）
        best_surface = None
        if surface_analyses:
            best_surface = min(
                surface_analyses,
                key=lambda x: x.pathway_analysis.max_barrier if x.pathway_analysis and x.pathway_analysis.max_barrier else float('inf')
            )
            self.logger.info("")
            self.logger.info(f"Best surface identified: {best_surface.miller_index}")

        # 步骤5: KMC模拟（可选）
        kmc_result = None
        if run_kmc and best_surface and best_surface.pathway_analysis:
            kmc_result = self.run_kmc_simulation(
                pathway_result=best_surface.pathway_analysis,
                temperature=temperature,
                pressures=pressures or {},
                output_dir=os.path.join(output_dir, "kmc"),
            )

        # 组装最终结果
        result = GasSolidWorkflowResult(
            bulk_formula=bulk_formula,
            reaction_intermediates=reaction_intermediates,
            temperature=temperature,
            surff_result=surff_result,
            surface_analyses=surface_analyses,
            best_surface=best_surface,
            kmc_result=kmc_result,
            output_dir=output_dir,
        )

        self.logger.info("")
        self.logger.info("=" * 80)
        self.logger.info("Workflow Complete!")
        self.logger.info("=" * 80)

        return result

    def _generate_slab_from_bulk(
        self,
        bulk: Atoms,
        miller: Tuple[int, int, int],
        layers: int = 4,
        vacuum: float = 10.0,
        size: Tuple[int, int, int] = (2, 2, 1),  # 默认2x2 supercell
    ) -> Atoms:
        """
        从bulk结构生成slab

        参数:
            bulk: bulk结构
            miller: Miller指数
            layers: 层数
            vacuum: 真空层厚度 (Å)
            size: supercell大小 (a, b, c)，前两个是xy方向的扩展

        返回:
            slab结构
        """
        from ase.build import surface

        try:
            # 使用size参数生成supercell，确保表面足够大
            slab = surface(bulk, miller, layers, vacuum=vacuum)
            # 扩展xy方向
            if size[0] > 1 or size[1] > 1:
                slab = slab.repeat((size[0], size[1], 1))
            return slab
        except Exception as e:
            self.logger.warning(f"Failed to generate slab with ASE: {e}")
            # 备用方案：返回简单的重复结构
            slab = bulk.repeat((size[0] * 2, size[1] * 2, layers))
            slab.center(vacuum=vacuum, axis=2)
            return slab


# =============================================================================
# 便捷函数
# =============================================================================

def run_gas_solid_catalysis(
    bulk_structure: Union[str, Atoms],
    reaction_intermediates: List[str],
    initial_reactant: str,
    temperature: float = 500.0,
    output_dir: str = "output/gas_solid",
    **kwargs
) -> GasSolidWorkflowResult:
    """
    便捷函数：运行气固界面催化数字孪生

    参数:
        bulk_structure: bulk结构文件路径或ASE Atoms对象
        reaction_intermediates: 反应中间体列表（如["*O", "*OH", "*OOH", "*"]）
        initial_reactant: 初始反应物（如"*O"）
        temperature: 反应温度 (K)
        output_dir: 输出目录
        **kwargs: 其他参数传递给GasSolidDigitalTwin.run_complete_workflow

    返回:
        完整的工作流程结果
    """
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # 检查是否有本地模型路径
    fairchem_model_path = kwargs.pop('fairchem_model_path', None)
    fairchem_model = kwargs.pop('fairchem_model', 'uma-s-1p1')

    # 如果没有提供本地路径，尝试使用项目中的模型
    if fairchem_model_path is None:
        project_model_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'deps', 'fairchem_models', f'{fairchem_model}.pt'
        )
        if os.path.exists(project_model_path):
            fairchem_model_path = project_model_path
            print(f"Using local model: {project_model_path}")

    # 初始化数字孪生
    dt = GasSolidDigitalTwin(
        surff_root="deps/SurFF",
        adsorbdiff_root="deps/AdsorbDiff",
        surface_sampling_root="deps/surface-sampling",
        fairchem_root="deps/fairchem",
        fairchem_model=fairchem_model,
        fairchem_model_path=fairchem_model_path,
    )

    # 运行工作流程
    result = dt.run_complete_workflow(
        bulk_structure=bulk_structure,
        reaction_intermediates=reaction_intermediates,
        initial_reactant=initial_reactant,
        temperature=temperature,
        output_dir=output_dir,
        **kwargs
    )

    return result


if __name__ == "__main__":
    # 示例：CO氧化反应
    result = run_gas_solid_catalysis(
        bulk_structure="POSCAR_Pt",
        reaction_intermediates=["*CO", "*O", "*CO2", "*"],
        initial_reactant="*CO",
        temperature=500.0,
        top_n_surfaces=3,
        calculate_barriers=True,
        output_dir="output/gas_solid_co_oxidation",
    )

    print(result.summary())
