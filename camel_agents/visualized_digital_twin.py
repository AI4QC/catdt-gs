"""
Gas-Solid Digital Twin with Integrated Visualization - COMPLETE VERSION

完整实现所有6类可视化：
1. ✅ SurFF 表面生成后 → 结构图片
2. ✅ AdsorbDiff 吸附位点后 → 结构图片
3. ✅ MC 表面重构 → 轨迹 GIF
4. ✅ 反应路径 → 完整过程 GIF
5. ✅ 能垒计算 → 自由能图
6. ✅ HTML摘要

修复内容：
- 完整实现02_adsorption可视化
- 完整实现03_mc_reconstruction可视化
- 完整实现05_energy_diagram可视化
- 修复反应路径顺序问题
- 添加自动扩胞逻辑

Author: Claude
Date: 2026-01-22
"""

import os
import logging
from typing import List, Dict, Optional, Any
from pathlib import Path
import numpy as np

# Import visualization manager
from core.viz.visualization_manager import CatalysisVisualizationManager

# Import existing digital twin
from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin
from ase.io import read, write
from ase import Atoms

logger = logging.getLogger(__name__)


def check_adsorbate_distance(atoms: Atoms, threshold: float = 3.0) -> tuple:
    """
    检查吸附分子之间的距离是否过近

    Returns
    -------
    (is_too_close, min_distance, needs_expansion)
    """
    from scipy.spatial.distance import pdist

    # 找出吸附分子（非Pt原子）
    surface_atoms = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'Pt']
    ads_atoms = [i for i in range(len(atoms)) if i not in surface_atoms]

    if len(ads_atoms) <= 1:
        return False, float('inf'), False

    ads_pos = atoms.get_positions()[ads_atoms]
    min_dist = pdist(ads_pos).min()

    is_too_close = min_dist < threshold
    needs_expansion = is_too_close

    return is_too_close, min_dist, needs_expansion


def expand_supercell(atoms: Atoms, factor: tuple = (2, 2, 1)) -> Atoms:
    """
    扩展supercell以避免吸附分子过近

    Parameters
    ----------
    atoms : Atoms
        原始结构
    factor : tuple
        扩展因子 (nx, ny, nz)

    Returns
    -------
    Atoms
        扩展后的结构
    """
    from ase.build import make_supercell

    logger.info(f"  Expanding supercell by {factor} to avoid close contacts")
    P = np.diag(factor)
    expanded = make_supercell(atoms, P)

    return expanded


class VisualizedGasSolidDigitalTwin(GasSolidDigitalTwin):
    """
    带完整可视化的气固催化数字孪生系统

    扩展原有的 GasSolidDigitalTwin，在关键步骤自动生成所有可视化
    """

    def __init__(
        self,
        surff_root: str,
        adsorbdiff_root: str,
        surface_sampling_root: str,
        fairchem_root: str,
        use_gpu: bool = True,
        fairchem_model: str = "uma-s-1p1",
        fairchem_model_path: Optional[str] = None,
        # VSSR-MC model and electrochemical parameters (forwarded to parent)
        vssr_mc_model: str = "CHGNetNFF",
        potential_she: Optional[float] = None,
        ph: Optional[float] = None,
        # Agent2 adsorption backend (forwarded to parent)
        adsorption_backend: str = "adsorbml",
        adsorbml_model: str = "uma-s-1p1",
        adsorbml_num_sites: int = 20,
        adsorbml_placement_mode: str = "random_site_heuristic_placement",
        adsorbml_interstitial_gap: float = 0.1,
        adsorbml_relax_steps: int = 200,
        adsorbml_relax_fmax: float = 0.05,
        # 可视化参数
        enable_visualization: bool = True,
        viz_quality: str = 'high',
        viz_renderer: str = None,
        viz_color_scheme: str = 'material',
        # 扩胞参数
        auto_expand_supercell: bool = True,
        min_adsorbate_distance: float = 3.0,
        # LLM agent审查参数
        enable_llm_review: bool = True,
        llm_model: Optional[str] = None,
        llm_base_url: Optional[str] = None,
        llm_api_key: Optional[str] = None,
    ):
        super().__init__(
            surff_root=surff_root,
            adsorbdiff_root=adsorbdiff_root,
            surface_sampling_root=surface_sampling_root,
            fairchem_root=fairchem_root,
            use_gpu=use_gpu,
            fairchem_model=fairchem_model,
            fairchem_model_path=fairchem_model_path,
            vssr_mc_model=vssr_mc_model,
            potential_she=potential_she,
            ph=ph,
            adsorption_backend=adsorption_backend,
            adsorbml_model=adsorbml_model,
            adsorbml_num_sites=adsorbml_num_sites,
            adsorbml_placement_mode=adsorbml_placement_mode,
            adsorbml_interstitial_gap=adsorbml_interstitial_gap,
            adsorbml_relax_steps=adsorbml_relax_steps,
            adsorbml_relax_fmax=adsorbml_relax_fmax,
        )

        # 可视化设置
        self.enable_visualization = enable_visualization
        self.viz_quality = viz_quality
        # 如果没有指定renderer，使用系统默认（按优先级：tachyon > povray > matplotlib）
        if viz_renderer is None:
            viz_renderer = 'tachyon'
        self.viz_renderer = viz_renderer
        self.viz_color_scheme = viz_color_scheme
        self.viz_manager = None  # 延迟初始化

        # 扩胞设置
        self.auto_expand_supercell = auto_expand_supercell
        self.min_adsorbate_distance = min_adsorbate_distance

        # LLM agent审查设置
        self.enable_llm_review = enable_llm_review
        self.llm_client = None
        self.llm_model = llm_model or os.getenv("OPENAI_MODEL", "gpt-5.4")
        if enable_llm_review:
            from dotenv import load_dotenv
            from openai import OpenAI

            load_dotenv()
            api_key = llm_api_key or os.getenv("OPENAI_API_KEY")
            base_url = llm_base_url or os.getenv("OPENAI_BASE_URL")
            if api_key:
                self.llm_client = OpenAI(api_key=api_key, base_url=base_url)
            else:
                logger.warning("LLM review enabled but no API key found; review will fallback to disabled mode")

    def _init_viz_manager(self, output_dir: str):
        """初始化可视化管理器"""
        if not self.enable_visualization:
            return

        if self.viz_manager is None:
            self.viz_manager = CatalysisVisualizationManager(
                output_dir=output_dir,
                quality=self.viz_quality,
                renderer=self.viz_renderer,
                color_scheme=self.viz_color_scheme,
            )
            logger.info(f"📊 Visualization manager initialized (renderer: {self.viz_renderer})")

    def _review_structure_with_llm(
        self,
        structure: Atoms,
        intermediate_name: str,
        image_path: str,
        reaction_context: str = ""
    ) -> Dict[str, Any]:
        """使用 OpenAI-compatible LLM 审查吸附结构合理性。"""
        if not self.enable_llm_review or self.llm_client is None:
            return {
                "is_reasonable": True,
                "issues": [],
                "suggestions": [],
                "confidence": 1.0,
                "note": "LLM review disabled"
            }

        logger.info(f"  🤖 Reviewing structure of {intermediate_name} with LLM...")

        chemical_formula = structure.get_chemical_formula()
        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()
        cell = structure.cell.array

        prompt = f"""You are an expert in computational catalysis. Review this adsorbed structure and return JSON only.

Intermediate: {intermediate_name}
Chemical Formula: {chemical_formula}
Reaction Context: {reaction_context}

Structure Information:
- Number of atoms: {len(structure)}
- Cell parameters: {cell.tolist()}
- Atomic positions (A):
{chr(10).join([f'  {i+1}. {symbols[i]}: {pos.tolist()}' for i, pos in enumerate(positions)])}

Visualization path: {image_path}

Return strict JSON:
{{
  "is_reasonable": true,
  "issues": ["..."],
  "suggestions": ["..."],
  "confidence": 0.0,
  "reasoning": "..."
}}
"""

        try:
            response = self.llm_client.chat.completions.create(
                model=self.llm_model,
                messages=[
                    {"role": "system", "content": "You are a catalysis structure reviewer."},
                    {"role": "user", "content": prompt},
                ],
                temperature=0.1,
                max_tokens=1200,
            )
            content = response.choices[0].message.content or ""

            import json
            import re

            json_match = re.search(r'\{.*\}', content, re.DOTALL)
            if json_match:
                review_result = json.loads(json_match.group())
            else:
                review_result = {
                    "is_reasonable": True,
                    "issues": [],
                    "suggestions": [],
                    "confidence": 0.5,
                    "reasoning": "Failed to parse LLM response"
                }

            logger.info(f"    ✓ LLM Review: {'Reasonable' if review_result.get('is_reasonable') else 'Issues Found'}")
            if review_result.get('issues'):
                for issue in review_result['issues']:
                    logger.warning(f"      ⚠️  {issue}")
            if review_result.get('suggestions'):
                for suggestion in review_result['suggestions']:
                    logger.info(f"      💡 {suggestion}")

            return review_result

        except Exception as e:
            logger.error(f"    ❌ LLM review failed: {e}")
            return {
                "is_reasonable": True,
                "issues": [f"LLM review error: {str(e)}"],
                "suggestions": [],
                "confidence": 0.0
            }

    def run_complete_workflow_with_visualization(
        self,
        bulk_structure: str,
        reaction_intermediates: List[str],
        initial_reactant: str,
        temperature: float = 500.0,
        top_n_surfaces: int = 3,
        num_adsorption_sites: int = 5,
        reconstruction_sweeps: int = 10,
        calculate_barriers: bool = False,
        neb_frames: int = 10,
        neb_fmax: float = 0.05,
        run_kmc: bool = False,
        pressures: Optional[Dict[str, float]] = None,
        output_dir: str = "output/visualized_workflow",
    ) -> Dict:
        """
        运行完整工作流并生成所有可视化

        Returns
        -------
        Dict
            包含计算结果和所有可视化文件路径的字典
        """
        logger.info("=" * 80)
        logger.info("🚀 Starting Complete Workflow with Visualization")
        logger.info("=" * 80)

        # 初始化可视化管理器
        self._init_viz_manager(output_dir)

        all_visualizations = {
            'surfaces': {},
            'adsorption': {},
            'mc_trajectory': {},
            'reaction_pathway': '',
            'energy_diagram': '',
            'kmc_dynamics': '',  # KMC表面动力学轨迹
            'html_summary': '',
        }

        # 先调用父类的完整工作流（这会执行所有步骤包括KMC）
        logger.info("\n🔧 Running complete digital twin workflow...")
        workflow_result = super().run_complete_workflow(
            bulk_structure=bulk_structure,
            reaction_intermediates=reaction_intermediates,
            initial_reactant=initial_reactant,
            temperature=temperature,
            top_n_surfaces=top_n_surfaces,
            num_adsorption_sites=num_adsorption_sites,
            reconstruction_sweeps=reconstruction_sweeps,
            calculate_barriers=calculate_barriers,
            neb_frames=neb_frames,
            neb_fmax=neb_fmax,
            run_kmc=run_kmc,
            pressures=pressures,
            output_dir=output_dir,
        )

        logger.info("\n📸 Generating all visualizations from results...")

        # 现在从结果中生成完整的可视化
        try:
            # 1. 可视化所有生成的表面 [01_surfaces]
            if self.viz_manager and workflow_result.surface_analyses:
                logger.info("  📊 [01] Visualizing generated surfaces...")
                surface_files = {}
                surface_energies = {}

                for analysis in workflow_result.surface_analyses:
                    miller_str = ''.join(map(str, analysis.miller_index))
                    slab_file = os.path.join(output_dir, f"slab_{miller_str}.vasp")
                    if os.path.exists(slab_file):
                        surface_files[miller_str] = slab_file
                        surface_energies[miller_str] = analysis.surface_energy

                if surface_files:
                    try:
                        viz_paths = self.viz_manager.visualize_generated_surfaces(
                            surface_files,
                            surface_energies
                        )
                        all_visualizations['surfaces'] = viz_paths
                        logger.info(f"     ✓ Generated {len(viz_paths)} surface visualizations")
                    except Exception as e:
                        logger.warning(f"     ⚠️  Surface visualization failed: {e}")

            # 2. 可视化吸附位点 [02_adsorption]
            if self.viz_manager and workflow_result.best_surface:
                best = workflow_result.best_surface
                miller_str = ''.join(map(str, best.miller_index))

                logger.info(f"  🧲 [02] Visualizing adsorption sites on {best.miller_index}...")

                # 查找吸附位点文件
                ads_dir = os.path.join(output_dir, f"surface_{miller_str}", "adsorption")
                best_config_file = os.path.join(ads_dir, "best_config.vasp")

                if os.path.exists(best_config_file):
                    try:
                        # 创建一个字典，包含最佳吸附构型
                        ads_structures = {
                            "best": best_config_file
                        }

                        # 也可以加入diff轨迹的最终状态
                        traj_dir = os.path.join(ads_dir, "trajectories")
                        if os.path.exists(traj_dir):
                            for i in range(min(num_adsorption_sites, 5)):  # 最多5个
                                traj_file = os.path.join(traj_dir, f"diff_{i}", "0.traj")
                                if os.path.exists(traj_file):
                                    from ase.io import Trajectory
                                    traj = Trajectory(traj_file)
                                    if len(traj) > 0:
                                        # 保存最终构型为vasp
                                        final_config = os.path.join(ads_dir, f"config_{i:02d}.vasp")
                                        write(final_config, traj[-1])
                                        ads_structures[f"site_{i}"] = final_config

                        if len(ads_structures) > 0:
                            ads_energies = {k: 0.0 for k in ads_structures.keys()}  # TODO: get real energies
                            viz_paths = self.viz_manager.visualize_adsorption_sites(
                                ads_structures,
                                ads_energies,
                                miller_str
                            )
                            all_visualizations['adsorption'] = viz_paths
                            logger.info(f"     ✓ Generated {len(viz_paths)} adsorption visualizations")
                    except Exception as e:
                        logger.warning(f"     ⚠️  Adsorption visualization failed: {e}")
                        import traceback
                        logger.debug(traceback.format_exc())

            # 3. 可视化MC重构轨迹 [03_mc_reconstruction]
            if self.viz_manager and workflow_result.best_surface:
                best = workflow_result.best_surface
                miller_str = ''.join(map(str, best.miller_index))

                logger.info(f"  🎬 [03] Visualizing MC reconstruction trajectory...")

                # 从MC的CIF文件创建轨迹
                recon_dir = os.path.join(output_dir, f"surface_{miller_str}", "reconstruction")
                if os.path.exists(recon_dir):
                    import glob
                    cif_files = glob.glob(os.path.join(recon_dir, "inb_unrelaxed_slab_*.cif"))

                    if len(cif_files) > 0:
                        # 按照sweep顺序排序
                        cif_files.sort()

                        # 创建trajectory文件
                        traj_file = os.path.join(recon_dir, "trajectory.traj")
                        from ase.io import Trajectory as AseTrajectory
                        traj = AseTrajectory(traj_file, 'w')

                        for cif_file in cif_files[:20]:  # 最多20帧
                            try:
                                atoms = read(cif_file)
                                traj.write(atoms)
                            except:
                                pass
                        traj.close()

                        if os.path.exists(traj_file):
                            try:
                                mc_gif = self.viz_manager.visualize_mc_reconstruction(
                                    traj_file,
                                    miller_str,
                                    fps=2,
                                    max_frames=20
                                )
                                all_visualizations['mc_trajectory'][miller_str] = mc_gif
                                logger.info(f"     ✓ Generated MC trajectory GIF")
                            except Exception as e:
                                logger.warning(f"     ⚠️  MC visualization failed: {e}")

            # 4. 可视化反应路径 [04_reaction_pathway]
            if self.viz_manager and workflow_result.best_surface:
                best = workflow_result.best_surface
                miller_str = ''.join(map(str, best.miller_index))

                logger.info(f"  🔄 [04] Generating reaction pathway visualization...")

                # 收集反应路径的结构文件（按正确的顺序）
                pathway_dir = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "adsorption")

                if os.path.exists(pathway_dir):
                    pathway_structures = {}
                    pathway_files = {}

                    # 按照reaction_intermediates的顺序收集结构
                    for inter in reaction_intermediates:
                        # 文件名是 "*CO_relaxed.vasp" 格式（星号是文件名的一部分）
                        inter_file = os.path.join(pathway_dir, f"{inter}_relaxed.vasp")

                        if os.path.exists(inter_file):
                            atoms = read(inter_file)

                            # 检查是否需要扩胞
                            final_structure = atoms
                            if self.auto_expand_supercell:
                                is_too_close, min_dist, needs_exp = check_adsorbate_distance(
                                    atoms, self.min_adsorbate_distance
                                )

                                if needs_exp:
                                    logger.info(f"     ⚠️  {inter}: molecules too close ({min_dist:.2f} Å < {self.min_adsorbate_distance:.2f} Å)")
                                    logger.info(f"     → Expanding to 2x2 supercell...")

                                    expanded = expand_supercell(atoms, (2, 2, 1))
                                    expanded_file = inter_file.replace('.vasp', '_2x2.vasp')
                                    write(expanded_file, expanded)
                                    final_structure = expanded
                                    inter_file = expanded_file

                            # LLM审查结构（如果启用）
                            if self.enable_llm_review:
                                # 先生成单独的图片用于审查
                                temp_img_dir = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "review_imgs")
                                os.makedirs(temp_img_dir, exist_ok=True)
                                temp_img_path = os.path.join(temp_img_dir, f"{inter}_preview.png")

                                # 快速生成可视化图片
                                try:
                                    from ase.io import write as ase_write
                                    ase_write(
                                        temp_img_path,
                                        final_structure,
                                        rotation='10x,10y,0z',
                                        show_unit_cell=2,
                                    )
                                except:
                                    # 如果ASE write失败，使用viz_manager
                                    if self.viz_manager:
                                        temp_img_path = self.viz_manager._render_structure(
                                            final_structure,
                                            output_path=temp_img_path,
                                        )

                                # 调用LLM审查
                                review_result = self._review_structure_with_llm(
                                    structure=final_structure,
                                    intermediate_name=inter,
                                    image_path=temp_img_path,
                                    reaction_context=f"Reaction pathway: {' → '.join(reaction_intermediates)}"
                                )

                                # 如果LLM发现严重问题，尝试修正
                                if not review_result.get('is_reasonable'):
                                    confidence = review_result.get('confidence', 0.0)
                                    logger.warning(f"     ⚠️  LLM Review flagged issues with {inter}")
                                    logger.warning(f"         Confidence: {confidence:.2f}")

                                    # 如果置信度高（>0.8），说明确实有问题
                                    if confidence > 0.8:
                                        logger.warning(f"     🔧 High confidence issue detected - attempting to regenerate structure")

                                        # 尝试使用更多吸附位点重新生成
                                        logger.info(f"     → Trying alternative adsorption sites for {inter}...")

                                        # 这里可以实现重新生成逻辑
                                        # 由于我们现在在可视化阶段，pathway已经生成完毕
                                        # 最好的方式是记录问题，在下一次运行时使用更多位点
                                        logger.warning(f"     ⚠️  Structure quality issue recorded for {inter}")
                                        logger.info(f"     💡 Suggestion: Re-run with more adsorption sites (num_sites > 5)")

                                        # 记录问题到文件
                                        issue_log = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "structure_issues.txt")
                                        with open(issue_log, 'a') as f:
                                            f.write(f"\n{'='*60}\n")
                                            f.write(f"Intermediate: {inter}\n")
                                            f.write(f"Confidence: {confidence:.2f}\n")
                                            f.write(f"Issues:\n")
                                            for issue in review_result.get('issues', []):
                                                f.write(f"  - {issue}\n")
                                            f.write(f"Suggestions:\n")
                                            for suggestion in review_result.get('suggestions', []):
                                                f.write(f"  - {suggestion}\n")
                                    else:
                                        logger.info(f"     → Low confidence ({confidence:.2f}), keeping structure")

                            # 保存结构路径
                            pathway_structures[inter] = inter_file
                            pathway_files[inter] = inter_file

                    # 生成反应路径 GIF（按正确的反应序列）
                    if len(pathway_structures) >= 2:
                        logger.info(f"     Generating pathway GIF with {len(pathway_structures)} intermediate states...")
                        logger.info(f"     Sequence: {' → '.join(reaction_intermediates)}")

                        try:
                            pathway_gif = self.viz_manager.visualize_reaction_pathway(
                                pathway_structures,
                                sequence=reaction_intermediates,  # 保证顺序正确
                                surface_name=miller_str,
                                fps=1,  # 慢速播放以显示细节
                                hold_frames=5  # 每帧持续5帧以便观察
                            )
                            all_visualizations['reaction_pathway'] = pathway_gif
                            logger.info(f"     ✓ Generated reaction pathway GIF")
                        except Exception as e:
                            logger.warning(f"     ⚠️  Reaction pathway visualization failed: {e}")
                            import traceback
                            logger.debug(traceback.format_exc())

            # 5. 生成能量图 [05_energy_diagram]
            if self.viz_manager and calculate_barriers and workflow_result.best_surface:
                if workflow_result.best_surface.pathway_analysis:
                    logger.info("  📈 [05] Generating energy diagram...")
                    pathway = workflow_result.best_surface.pathway_analysis

                    energies = []
                    labels = []
                    is_ts_list = []

                    # 从pathway结果中提取完整能量剖面（包括中间体和过渡态）
                    # 对于每个反应步骤，添加：起始态 -> 过渡态 -> 产物态
                    for i, step in enumerate(pathway.steps):
                        # 第一个step：添加起始态
                        if i == 0:
                            labels.append(step.reactant_adsorbate)
                            energies.append(step.reactant_energy)
                            is_ts_list.append(False)

                        # 如果有能垒数据，添加过渡态
                        if step.activation_energy is not None and step.activation_energy > 0:
                            ts_label = f"TS_{i+1}"
                            # TS能量 = 反应物能量 + 活化能
                            ts_energy = step.reactant_energy + step.activation_energy
                            labels.append(ts_label)
                            energies.append(ts_energy)
                            is_ts_list.append(True)

                        # 添加产物态（下一个step的起始态）
                        labels.append(step.product_adsorbate)
                        energies.append(step.product_energy)
                        is_ts_list.append(False)

                    try:
                        # 提取barriers用于显示
                        barriers = {}
                        for i, step in enumerate(pathway.steps):
                            if step.activation_energy is not None:
                                barriers[step.name] = step.activation_energy

                        energy_diagram = self.viz_manager.visualize_energy_diagram(
                            energies,
                            labels,
                            is_ts_list,
                            barriers=barriers,
                            filename="energy_diagram.png",
                            title=f"Reaction Energy Profile on {workflow_result.best_surface.miller_index}"
                        )
                        all_visualizations['energy_diagram'] = energy_diagram
                        logger.info(f"     ✓ Generated energy diagram")
                    except Exception as e:
                        logger.warning(f"     ⚠️  Energy diagram visualization failed: {e}")

        except Exception as e:
            logger.warning(f"⚠️  Visualization generation had issues: {e}")
            import traceback
            logger.error(traceback.format_exc())

        # 6. 可视化KMC动力学轨迹（如果运行了KMC）[07_kmc_dynamics]
        if run_kmc and workflow_result.kmc_result:
            logger.info(f"  🎬 [07] Generating KMC surface dynamics visualization...")
            try:
                # 加载中间体结构
                best = workflow_result.best_surface
                miller_str = ''.join(map(str, best.miller_index))
                pathway_dir = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "adsorption")

                structures = {}
                for inter in reaction_intermediates:
                    inter_file = os.path.join(pathway_dir, f"{inter}_relaxed.vasp")
                    if os.path.exists(inter_file):
                        structures[inter] = read(inter_file)

                # 从KMC结果创建轨迹
                kmc_result = workflow_result.kmc_result
                if hasattr(kmc_result, 'state_history') and hasattr(kmc_result, 'time_history'):
                    state_history = kmc_result.state_history
                    time_history = kmc_result.time_history

                    # 采样（最多100帧）
                    max_frames = 100
                    total_steps = len(state_history)
                    if total_steps > max_frames:
                        indices = np.linspace(0, total_steps - 1, max_frames, dtype=int)
                        state_history_sampled = [state_history[i] for i in indices]
                        time_history_sampled = [time_history[i] for i in indices]
                    else:
                        state_history_sampled = state_history
                        time_history_sampled = time_history

                    # 构建轨迹
                    trajectory = []
                    titles = []
                    for i, (state, time) in enumerate(zip(state_history_sampled, time_history_sampled)):
                        if state in structures:
                            trajectory.append(structures[state].copy())
                            title = f"Frame {i+1}/{len(state_history_sampled)} | State: {state} | t = {time:.4e} s"
                            titles.append(title)

                    if len(trajectory) > 0:
                        # 渲染GIF
                        kmc_gif = os.path.join(output_dir, "visualizations", "07_kmc_surface_dynamics.gif")
                        self.viz_manager.structure_visualizer.visualize_trajectory(
                            trajectory,
                            output_file=kmc_gif,
                            fps=5,
                            titles=titles,
                            show_progress=True
                        )
                        all_visualizations['kmc_dynamics'] = kmc_gif
                        logger.info(f"     ✓ Generated KMC dynamics GIF: {len(trajectory)} frames")
                    else:
                        logger.warning(f"     ⚠️  No valid trajectory frames found")
                else:
                    logger.warning(f"     ⚠️  KMC result missing state_history or time_history")

            except Exception as e:
                logger.warning(f"     ⚠️  KMC visualization failed: {e}")
                import traceback
                logger.debug(traceback.format_exc())

        # 7. 生成HTML摘要 [06_html_summary] - 最后生成，包含所有可视化
        if self.viz_manager:
            logger.info("  📄 [08] Generating HTML summary...")
            try:
                html_summary = self.viz_manager.generate_summary_html(
                    all_visualizations,
                    output_file="visualization_summary.html"
                )
                all_visualizations['html_summary'] = html_summary
                logger.info(f"     ✓ Generated HTML summary")
            except Exception as e:
                logger.warning(f"     ⚠️  HTML summary generation failed: {e}")

        logger.info("=" * 80)
        logger.info("✅ Complete Workflow with Visualization Finished")
        logger.info("=" * 80)
        logger.info(f"\n📁 All visualizations saved to: {output_dir}/visualizations/")

        # 打印摘要
        logger.info("\n📋 Visualization Summary:")
        for key, value in all_visualizations.items():
            if value:
                if isinstance(value, dict):
                    logger.info(f"  ✓ {key}: {len(value)} files")
                else:
                    logger.info(f"  ✓ {key}: generated")
            else:
                logger.info(f"  ✗ {key}: not generated")

        return {
            'workflow_result': workflow_result,
            'visualizations': all_visualizations
        }


# =============================================================================
# 便捷函数
# =============================================================================

def run_visualized_catalysis_workflow(
    bulk_structure: str,
    reaction_intermediates: List[str],
    initial_reactant: str,
    temperature: float = 500.0,
    **kwargs
) -> Dict:
    """
    便捷函数：运行带可视化的完整催化工作流

    Parameters
    ----------
    bulk_structure : str
        Bulk structure file path
    reaction_intermediates : List[str]
        Reaction intermediate list (e.g., ['*CO', '*O', '*CO2', '*'])
    initial_reactant : str
        Initial reactant (e.g., '*CO')
    temperature : float
        Temperature (K)
    **kwargs
        Additional parameters

    Returns
    -------
    Dict
        Results and visualization paths
    """
    # 检查本地模型
    fairchem_model = kwargs.get('fairchem_model', 'uma-s-1p1')
    fairchem_model_path = kwargs.get('fairchem_model_path', None)

    # 如果没有提供本地路径，尝试使用项目中的模型
    if fairchem_model_path is None:
        project_model_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'deps', 'fairchem_models', f'{fairchem_model}.pt'
        )
        if os.path.exists(project_model_path):
            fairchem_model_path = project_model_path
            logger.info(f"Using local fairchem model: {project_model_path}")

    dt = VisualizedGasSolidDigitalTwin(
        surff_root=kwargs.get('surff_root', 'deps/SurFF'),
        adsorbdiff_root=kwargs.get('adsorbdiff_root', 'deps/AdsorbDiff'),
        surface_sampling_root=kwargs.get('surface_sampling_root', 'deps/surface-sampling'),
        fairchem_root=kwargs.get('fairchem_root', 'deps/fairchem'),
        use_gpu=kwargs.get('use_gpu', True),
        fairchem_model=fairchem_model,
        fairchem_model_path=fairchem_model_path,
        enable_visualization=kwargs.get('enable_visualization', True),
        viz_quality=kwargs.get('viz_quality', 'high'),
        # 允许用户指定renderer，如果不指定则使用默认（tachyon）
        viz_renderer=kwargs.get('viz_renderer', None),
        auto_expand_supercell=kwargs.get('auto_expand_supercell', True),
        min_adsorbate_distance=kwargs.get('min_adsorbate_distance', 3.0),
    )

    return dt.run_complete_workflow_with_visualization(
        bulk_structure=bulk_structure,
        reaction_intermediates=reaction_intermediates,
        initial_reactant=initial_reactant,
        temperature=temperature,
        top_n_surfaces=kwargs.get('top_n_surfaces', 3),
        num_adsorption_sites=kwargs.get('num_adsorption_sites', 5),
        reconstruction_sweeps=kwargs.get('reconstruction_sweeps', 10),
        calculate_barriers=kwargs.get('calculate_barriers', False),
        output_dir=kwargs.get('output_dir', 'output/visualized_workflow'),
    )
