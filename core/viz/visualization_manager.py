"""
Unified Visualization Manager for Gas-Solid Catalysis Digital Twin

自动在以下关键步骤生成可视化：
1. SurFF 表面生成后 → 结构图片
2. AdsorbDiff 吸附位点后 → 结构图片
3. MC 表面重构模拟 → 轨迹 GIF
4. 反应路径中间体 → 完整反应 GIF
5. 能垒计算后 → 自由能图

Author: Claude
Date: 2026-01-21
"""

import os
import logging
from typing import List, Dict, Optional, Tuple, Union
from pathlib import Path
import numpy as np

from ase import Atoms
from ase.io import read, write, Trajectory

# Import viz modules
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle, EnergyStep


logger = logging.getLogger(__name__)


class CatalysisVisualizationManager:
    """
    统一的可视化管理器

    在催化工作流的关键步骤自动生成高质量可视化：
    - 表面结构图
    - 吸附构型图
    - MC轨迹动画
    - 反应路径动画
    - 自由能图
    """

    def __init__(
        self,
        output_dir: str,
        quality: str = 'high',
        renderer: str = 'tachyon',
        elevation: float = 30,
        azimuth: float = 45,
        color_scheme: str = 'material',
        enable_all: bool = True,
    ):
        """
        Initialize visualization manager

        Parameters
        ----------
        output_dir : str
            输出目录
        quality : str
            渲染质量 ('low', 'medium', 'high', 'ultra')
        renderer : str
            渲染器 ('tachyon', 'ospray', 'opengl', 'anari')
        elevation : float
            视角仰角
        azimuth : float
            视角方位角
        color_scheme : str
            能量图配色方案 ('material', 'nature', 'science', 'elegant')
        enable_all : bool
            是否启用所有可视化（可以单独禁用某些步骤）
        """
        self.output_dir = Path(output_dir)
        self.quality = quality
        self.renderer = renderer
        self.elevation = elevation
        self.azimuth = azimuth
        self.color_scheme = color_scheme

        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.viz_dir = self.output_dir / "visualizations"
        self.viz_dir.mkdir(exist_ok=True)

        # 可视化开关
        self.enable_surface_viz = enable_all
        self.enable_adsorption_viz = enable_all
        self.enable_mc_trajectory_viz = enable_all
        self.enable_pathway_viz = enable_all
        self.enable_energy_diagram_viz = enable_all

        # 初始化可视化器
        self.structure_visualizer = CatalystSurfaceVisualizer(
            elevation=elevation,
            azimuth=azimuth,
            quality=quality,
            renderer=renderer,
            show_cell=False,
            auto_expand=True,
        )

        self.energy_plotter = None  # 延迟初始化

        logger.info(f"📊 Visualization Manager initialized")
        logger.info(f"   Output: {self.viz_dir}")
        logger.info(f"   Quality: {quality}")
        logger.info(f"   Renderer: {renderer}")

    # ========================================================================
    # 1. SurFF 表面生成后可视化
    # ========================================================================

    def visualize_generated_surfaces(
        self,
        surface_files: Dict[str, str],
        surface_energies: Optional[Dict[str, float]] = None,
    ) -> Dict[str, str]:
        """
        可视化 SurFF 生成的所有表面

        Parameters
        ----------
        surface_files : Dict[str, str]
            表面文件字典 {miller_index: filepath}
            例如: {'111': 'slab_111.vasp', '110': 'slab_110.vasp'}
        surface_energies : Dict[str, float], optional
            表面能 {miller_index: energy_eV_per_A2}

        Returns
        -------
        Dict[str, str]
            生成的图片路径 {miller_index: image_path}
        """
        if not self.enable_surface_viz:
            logger.info("Surface visualization disabled")
            return {}

        logger.info("📸 Visualizing generated surfaces...")

        surface_dir = self.viz_dir / "01_surfaces"
        surface_dir.mkdir(exist_ok=True)

        output_paths = {}

        for miller_idx, filepath in surface_files.items():
            try:
                # 读取结构
                slab = read(filepath)

                # 生成标题
                title = f"Surface ({miller_idx})"
                if surface_energies and miller_idx in surface_energies:
                    energy = surface_energies[miller_idx]
                    title += f"\nγ = {energy:.4f} eV/Å²"

                # 输出路径
                output_path = surface_dir / f"surface_{miller_idx}.png"

                # 可视化
                self.structure_visualizer.visualize_structure(
                    slab,
                    output_file=str(output_path),
                    title=title,
                    show=False
                )

                output_paths[miller_idx] = str(output_path)
                logger.info(f"  ✓ {miller_idx}: {output_path}")

            except Exception as e:
                logger.error(f"  ✗ Failed to visualize surface {miller_idx}: {e}")

        logger.info(f"📸 Generated {len(output_paths)} surface images")
        return output_paths

    # ========================================================================
    # 2. AdsorbDiff 吸附位点后可视化
    # ========================================================================

    def visualize_adsorption_sites(
        self,
        adsorbate_structures: Dict[str, str],
        adsorption_energies: Optional[Dict[str, float]] = None,
        surface_name: str = "surface",
    ) -> Dict[str, str]:
        """
        可视化吸附构型

        Parameters
        ----------
        adsorbate_structures : Dict[str, str]
            吸附构型文件 {adsorbate_name: filepath}
            例如: {'CO': 'CO_adsorbed.vasp', 'O': 'O_adsorbed.vasp'}
        adsorption_energies : Dict[str, float], optional
            吸附能 {adsorbate_name: energy_eV}
        surface_name : str
            表面名称（用于文件命名）

        Returns
        -------
        Dict[str, str]
            生成的图片路径 {adsorbate_name: image_path}
        """
        if not self.enable_adsorption_viz:
            logger.info("Adsorption visualization disabled")
            return {}

        logger.info("📸 Visualizing adsorption sites...")

        ads_dir = self.viz_dir / "02_adsorption"
        ads_dir.mkdir(exist_ok=True)

        output_paths = {}

        for ads_name, filepath in adsorbate_structures.items():
            try:
                # 读取结构
                structure = read(filepath)

                # 生成标题
                title = f"{ads_name} adsorbed"
                if adsorption_energies and ads_name in adsorption_energies:
                    energy = adsorption_energies[ads_name]
                    title += f"\nE_ads = {energy:.3f} eV"

                # 输出路径
                output_path = ads_dir / f"{surface_name}_{ads_name}_adsorbed.png"

                # 可视化
                self.structure_visualizer.visualize_structure(
                    structure,
                    output_file=str(output_path),
                    title=title,
                    show=False
                )

                output_paths[ads_name] = str(output_path)
                logger.info(f"  ✓ {ads_name}: {output_path}")

            except Exception as e:
                logger.error(f"  ✗ Failed to visualize {ads_name}: {e}")

        logger.info(f"📸 Generated {len(output_paths)} adsorption images")
        return output_paths

    # ========================================================================
    # 3. MC 表面重构轨迹可视化
    # ========================================================================

    def visualize_mc_reconstruction(
        self,
        trajectory_file: Union[str, List[Atoms]],
        surface_name: str = "surface",
        fps: int = 5,
        max_frames: Optional[int] = None,
    ) -> str:
        """
        可视化 MC 表面重构轨迹为 GIF

        Parameters
        ----------
        trajectory_file : Union[str, List[Atoms]]
            轨迹文件路径(.traj)或直接传入的 Atoms 帧列表
        surface_name : str
            表面名称
        fps : int
            GIF 帧率
        max_frames : int, optional
            最大帧数（用于长轨迹采样）

        Returns
        -------
        str
            GIF 文件路径
        """
        if not self.enable_mc_trajectory_viz:
            logger.info("MC trajectory visualization disabled")
            return ""

        logger.info("🎬 Generating MC reconstruction GIF...")

        mc_dir = self.viz_dir / "03_mc_reconstruction"
        mc_dir.mkdir(exist_ok=True)

        try:
            if isinstance(trajectory_file, list):
                frames = [atoms.copy() for atoms in trajectory_file if isinstance(atoms, Atoms)]
                n_frames = len(frames)
            else:
                traj = Trajectory(trajectory_file)
                n_frames = len(traj)
                frames = [atoms for atoms in traj]

            logger.info(f"   Trajectory: {n_frames} frames")

            if n_frames == 0:
                logger.warning("   Empty MC trajectory, skip visualization")
                return ""

            if max_frames and n_frames > max_frames:
                indices = np.linspace(0, n_frames - 1, max_frames, dtype=int)
                frames = [frames[i] for i in indices]
                logger.info(f"   Sampled to {len(frames)} frames")

            titles = [f"MC Step {i+1}/{len(frames)}" for i in range(len(frames))]
            output_path = mc_dir / f"{surface_name}_mc_reconstruction.gif"

            self.structure_visualizer.visualize_trajectory(
                frames,
                output_file=str(output_path),
                fps=fps,
                titles=titles,
                show_progress=True,
            )

            logger.info(f"🎬 Generated MC GIF: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"✗ Failed to generate MC GIF: {e}")
            return ""

    # ========================================================================
    # 4. 反应路径完整过程可视化
    # ========================================================================

    def visualize_reaction_pathway(
        self,
        pathway_structures: Dict[str, str],
        sequence: List[str],
        surface_name: str = "surface",
        fps: int = 2,
        hold_frames: int = 3,
    ) -> str:
        """
        可视化完整反应路径为 GIF

        Parameters
        ----------
        pathway_structures : Dict[str, str]
            路径结构 {intermediate_name: filepath}
            例如: {'*CO': 'CO_ads.vasp', '*O': 'O_ads.vasp', '*CO2': 'CO2_ads.vasp'}
        sequence : List[str]
            反应顺序
            例如: ['*CO', '*O', '*CO2', '*']
        surface_name : str
            表面名称
        fps : int
            GIF 帧率
        hold_frames : int
            每个中间体停留帧数

        Returns
        -------
        str
            GIF 文件路径
        """
        if not self.enable_pathway_viz:
            logger.info("Pathway visualization disabled")
            return ""

        logger.info("🎬 Generating reaction pathway GIF...")

        pathway_dir = self.viz_dir / "04_reaction_pathway"
        pathway_dir.mkdir(exist_ok=True)

        try:
            # 读取所有结构
            structures = {}
            for name, filepath in pathway_structures.items():
                if os.path.exists(filepath):
                    structures[name] = read(filepath)

            # 构建帧序列（每个中间体重复几帧）
            frames = []
            titles = []

            for intermediate in sequence:
                if intermediate in structures:
                    for i in range(hold_frames):
                        frames.append(structures[intermediate])
                        titles.append(f"{intermediate} ({i+1}/{hold_frames})")
                else:
                    logger.warning(f"  ⚠ Structure not found: {intermediate}")

            if not frames:
                logger.error("✗ No frames to visualize")
                return ""

            logger.info(f"   Total frames: {len(frames)}")

            # 输出路径
            output_path = pathway_dir / f"{surface_name}_reaction_pathway.gif"

            # 生成 GIF
            self.structure_visualizer.visualize_trajectory(
                frames,
                output_file=str(output_path),
                fps=fps,
                titles=titles,
                show_progress=True
            )

            logger.info(f"🎬 Generated pathway GIF: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"✗ Failed to generate pathway GIF: {e}")
            return ""

    # ========================================================================
    # 5. 自由能图可视化
    # ========================================================================

    def visualize_energy_diagram(
        self,
        energies: List[float],
        labels: List[str],
        is_ts_list: Optional[List[bool]] = None,
        barriers: Optional[Dict[str, float]] = None,
        filename: str = "energy_diagram.png",
        title: Optional[str] = None,
    ) -> str:
        """
        生成自由能图

        Parameters
        ----------
        energies : List[float]
            能量值列表 (eV)
        labels : List[str]
            状态标签
        is_ts_list : List[bool], optional
            过渡态标记
        barriers : Dict[str, float], optional
            能垒信息 {step_name: barrier_eV}
        filename : str
            输出文件名
        title : str, optional
            图表标题

        Returns
        -------
        str
            图片文件路径
        """
        if not self.enable_energy_diagram_viz:
            logger.info("Energy diagram visualization disabled")
            return ""

        logger.info("📈 Generating energy diagram...")

        energy_dir = self.viz_dir / "05_energy_diagram"
        energy_dir.mkdir(exist_ok=True)

        try:
            # 创建样式
            style = DiagramStyle(
                color_scheme=self.color_scheme,
                figsize=(12, 7),
                ylabel='Free Energy',
                ylabel_units='eV',
                show_energies=True,
                energy_format='.2f',
                dpi=300,
            )

            # 创建绘图器
            plotter = EnergyDiagramPlotter(style=style)

            # 绘制
            plotter.plot(
                energies=energies,
                labels=labels,
                is_ts_list=is_ts_list,
                reference_index=0,
                show=False
            )

            # 添加标题（如果提供）
            if title:
                plotter.ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

            # 保存
            output_path = energy_dir / filename
            plotter.save(str(output_path), dpi=300)

            # 同时导出数据
            csv_path = energy_dir / filename.replace('.png', '.csv')
            plotter.export_data(str(csv_path))

            logger.info(f"📈 Generated energy diagram: {output_path}")
            logger.info(f"📊 Exported data: {csv_path}")

            # 打印关键信息
            if barriers:
                logger.info("   Activation barriers:")
                for step, barrier in barriers.items():
                    logger.info(f"     {step}: {barrier:.3f} eV")

            return str(output_path)

        except Exception as e:
            logger.error(f"✗ Failed to generate energy diagram: {e}")
            import traceback
            traceback.print_exc()
            return ""

    # ========================================================================
    # 便捷方法：完整工作流可视化
    # ========================================================================

    def visualize_complete_workflow(
        self,
        workflow_result,
        surface_name: str = "surface",
    ) -> Dict[str, str]:
        """
        对完整的工作流结果进行可视化

        Parameters
        ----------
        workflow_result : WorkflowResult
            工作流的完整结果对象
        surface_name : str
            表面名称

        Returns
        -------
        Dict[str, str]
            所有生成的可视化文件路径
            {
                'surfaces': [path1, path2, ...],
                'adsorption': [path1, path2, ...],
                'mc_trajectory': path,
                'reaction_pathway': path,
                'energy_diagram': path
            }
        """
        logger.info("=" * 80)
        logger.info("🎨 Starting complete workflow visualization")
        logger.info("=" * 80)

        viz_outputs = {
            'surfaces': [],
            'adsorption': [],
            'mc_trajectory': None,
            'reaction_pathway': None,
            'energy_diagram': None,
        }

        # TODO: 根据实际的 workflow_result 结构调用各个可视化方法
        # 这里提供模板

        # 1. 表面可视化
        # if hasattr(workflow_result, 'surface_files'):
        #     viz_outputs['surfaces'] = self.visualize_generated_surfaces(
        #         workflow_result.surface_files,
        #         workflow_result.surface_energies
        #     )

        # 2. 吸附可视化
        # if hasattr(workflow_result, 'adsorbate_structures'):
        #     viz_outputs['adsorption'] = self.visualize_adsorption_sites(
        #         workflow_result.adsorbate_structures,
        #         workflow_result.adsorption_energies,
        #         surface_name
        #     )

        # 3. MC轨迹
        # if hasattr(workflow_result, 'mc_trajectory_file'):
        #     viz_outputs['mc_trajectory'] = self.visualize_mc_reconstruction(
        #         workflow_result.mc_trajectory_file,
        #         surface_name
        #     )

        # 4. 反应路径
        # if hasattr(workflow_result, 'pathway_structures'):
        #     viz_outputs['reaction_pathway'] = self.visualize_reaction_pathway(
        #         workflow_result.pathway_structures,
        #         workflow_result.intermediate_sequence,
        #         surface_name
        #     )

        # 5. 能量图
        # if hasattr(workflow_result, 'pathway_energies'):
        #     viz_outputs['energy_diagram'] = self.visualize_energy_diagram(
        #         workflow_result.pathway_energies,
        #         workflow_result.intermediate_labels,
        #         workflow_result.is_ts_list
        #     )

        logger.info("=" * 80)
        logger.info("🎨 Complete workflow visualization finished")
        logger.info("=" * 80)

        return viz_outputs

    def generate_summary_html(
        self,
        viz_outputs: Dict[str, str],
        output_file: str = "visualization_summary.html"
    ) -> str:
        """
        生成HTML摘要页面，展示所有可视化结果

        Parameters
        ----------
        viz_outputs : Dict[str, str]
            可视化输出路径
        output_file : str
            HTML文件名

        Returns
        -------
        str
            HTML文件路径
        """
        import datetime
        html_path = self.viz_dir / output_file

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Catalysis Workflow Visualization Summary</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #1976D2;
            border-bottom: 3px solid #1976D2;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #388E3C;
            margin-top: 30px;
        }}
        .viz-section {{
            background-color: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            margin: 10px 0;
        }}
        .image-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
        }}
    </style>
</head>
<body>
    <h1>🎨 Catalysis Workflow Visualization Summary</h1>
    <p>Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

    <!-- Add sections for each visualization type -->
    <!-- This is a template - customize based on actual viz_outputs -->

</body>
</html>
"""

        with open(html_path, 'w') as f:
            f.write(html_content)

        logger.info(f"📄 Generated HTML summary: {html_path}")
        return str(html_path)

    # ========================================================================
    # 6. NEB轨迹可视化
    # ========================================================================

    def visualize_neb_trajectory(
        self,
        trajectory: List[Atoms],
        reaction_name: str,
        output_filename: str = None,
        fps: int = 2,
        is_optimized: bool = False,
    ) -> str:
        """
        可视化NEB轨迹（插值路径或优化后的路径）

        Parameters
        ----------
        trajectory : List[Atoms]
            NEB轨迹帧列表
        reaction_name : str
            反应名称（如 "*CO_to_*CHO"）
        output_filename : str, optional
            输出文件名（默认自动生成）
        fps : int
            GIF帧率
        is_optimized : bool
            是否为优化后的路径（True）还是插值路径（False）

        Returns
        -------
        str
            GIF文件路径
        """
        logger.info(f"🎬 Generating NEB trajectory visualization for {reaction_name}...")

        neb_dir = self.viz_dir / "06_neb_trajectories"
        neb_dir.mkdir(exist_ok=True)

        try:
            # 自动生成文件名
            if output_filename is None:
                prefix = "optimized" if is_optimized else "interpolated"
                output_filename = f"{reaction_name}_{prefix}_path.gif"

            output_path = neb_dir / output_filename

            # 生成标题
            n_frames = len(trajectory)
            if is_optimized:
                # 优化后的路径：提取最后n帧（收敛路径）
                frame_titles = [
                    f"{reaction_name} - Reactant (Frame 0/{n_frames-1})"
                ]
                for i in range(1, n_frames-1):
                    frame_titles.append(
                        f"{reaction_name} - Image {i} (Frame {i}/{n_frames-1})"
                    )
                frame_titles.append(
                    f"{reaction_name} - Product (Frame {n_frames-1}/{n_frames-1})"
                )
            else:
                # 插值路径
                frame_titles = [
                    f"{reaction_name} - Initial (Frame 0/{n_frames-1})"
                ]
                for i in range(1, n_frames-1):
                    frame_titles.append(
                        f"{reaction_name} - Interpolated {i} (Frame {i}/{n_frames-1})"
                    )
                frame_titles.append(
                    f"{reaction_name} - Final (Frame {n_frames-1}/{n_frames-1})"
                )

            # 生成GIF
            self.structure_visualizer.visualize_trajectory(
                trajectory=trajectory,
                output_file=str(output_path),
                fps=fps,
                titles=frame_titles,
                show_progress=False  # 不显示进度条，避免干扰日志
            )

            logger.info(f"✓ Saved NEB trajectory: {output_path}")
            logger.info(f"  Frames: {n_frames}, Type: {'Optimized' if is_optimized else 'Interpolated'}")

            return str(output_path)

        except Exception as e:
            logger.error(f"✗ Failed to generate NEB trajectory GIF: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return ""

    def visualize_neb_from_file(
        self,
        traj_file: str,
        reaction_name: str = None,
        n_images: int = 5,
        output_filename: str = None,
        fps: int = 2,
    ) -> str:
        """
        从.traj文件读取并可视化NEB轨迹（提取收敛路径）

        ASE NEB 轨迹每个优化步写入整个 band（n_images 帧），因此收敛路径
        是最后一个完整 band，即最后 n_images 帧。n_images 必须等于 NEB
        band 的图像数，否则切片会混合两个优化迭代的帧。

        Parameters
        ----------
        traj_file : str
            轨迹文件路径
        reaction_name : str, optional
            反应名称（从文件名推断）
        n_images : int
            NEB band 的图像数（每个优化步写入的帧数）；
            提取最后 n_images 帧作为收敛路径
        output_filename : str, optional
            输出文件名
        fps : int
            GIF帧率

        Returns
        -------
        str
            GIF文件路径
        """
        from ase.io import read as ase_read

        logger.info(f"🎬 Visualizing NEB trajectory from {traj_file}...")

        try:
            # 读取轨迹
            full_trajectory = ase_read(traj_file, index=':')
            total_frames = len(full_trajectory)

            logger.info(f"  Total frames in trajectory: {total_frames}")

            # 提取收敛路径（最后一个完整 band = 最后 n_images 帧）
            if total_frames >= n_images:
                if total_frames % n_images != 0:
                    logger.warning(
                        f"  Trajectory length {total_frames} is not divisible by "
                        f"n_images={n_images}; n_images probably does not match "
                        f"the NEB band size, so the last {n_images} frames may mix "
                        f"two optimization iterations. Falling back to slicing the "
                        f"last {n_images} frames anyway."
                    )
                else:
                    logger.info(
                        f"  Trajectory contains {total_frames // n_images} bands "
                        f"of {n_images} images"
                    )
                converged_path = full_trajectory[-n_images:]
                logger.info(f"  Extracting last {n_images} frames as converged path")
            else:
                converged_path = full_trajectory
                logger.warning(f"  Trajectory has only {total_frames} frames, using all")

            # 推断反应名称
            if reaction_name is None:
                reaction_name = Path(traj_file).stem.replace('_neb', '')

            # 可视化
            return self.visualize_neb_trajectory(
                trajectory=converged_path,
                reaction_name=reaction_name,
                output_filename=output_filename,
                fps=fps,
                is_optimized=True,
            )

        except Exception as e:
            logger.error(f"✗ Failed to visualize NEB from file: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return ""


# =============================================================================
# 测试代码
# =============================================================================

def test_visualization_manager():
    """测试可视化管理器"""
    from ase.build import fcc111, molecule, add_adsorbate
    from ase.io import write, Trajectory

    print("=" * 80)
    print("Testing Catalysis Visualization Manager")
    print("=" * 80)

    # 创建测试目录
    test_dir = "/tmp/viz_test"
    os.makedirs(test_dir, exist_ok=True)

    # 创建管理器
    viz_manager = CatalysisVisualizationManager(
        output_dir=test_dir,
        quality='high',
        renderer='tachyon',
        color_scheme='material'
    )

    # 1. 测试表面可视化
    print("\n1. Testing surface visualization...")
    slab_111 = fcc111('Pt', size=(4, 4, 4), vacuum=10.0)
    slab_110 = fcc111('Pt', size=(4, 4, 4), vacuum=10.0)

    write(f"{test_dir}/slab_111.vasp", slab_111)
    write(f"{test_dir}/slab_110.vasp", slab_110)

    surface_files = {
        '111': f"{test_dir}/slab_111.vasp",
        '110': f"{test_dir}/slab_110.vasp"
    }
    surface_energies = {'111': 0.0928, '110': 0.1055}

    viz_manager.visualize_generated_surfaces(surface_files, surface_energies)

    # 2. 测试吸附可视化
    print("\n2. Testing adsorption visualization...")
    slab_with_co = slab_111.copy()
    co = molecule('CO')
    add_adsorbate(slab_with_co, co, height=2.0, position='ontop')
    write(f"{test_dir}/CO_ads.vasp", slab_with_co)

    ads_structures = {'*CO': f"{test_dir}/CO_ads.vasp"}
    ads_energies = {'*CO': -1.5}

    viz_manager.visualize_adsorption_sites(ads_structures, ads_energies, "Pt111")

    # 3. 测试能量图
    print("\n3. Testing energy diagram...")
    energies = [0.0, 1.2, 0.3, 0.9, -1.5]
    labels = ['*CO + *O', 'TS1', '*CO-O', 'TS2', '*CO2 + *']
    is_ts = [False, True, False, True, False]

    viz_manager.visualize_energy_diagram(
        energies, labels, is_ts,
        filename="test_energy_diagram.png",
        title="CO Oxidation on Pt(111)"
    )

    print("\n" + "=" * 80)
    print("✓ All visualization tests completed!")
    print(f"Check output: {test_dir}/visualizations/")
    print("=" * 80)


if __name__ == '__main__':
    test_visualization_manager()
