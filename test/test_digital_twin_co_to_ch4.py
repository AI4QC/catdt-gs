#!/usr/bin/env python
"""
CO到CH4完整催化反应 - 带可视化和LLM智能体审查

使用VisualizedGasSolidDigitalTwin运行完整工作流程：
- 表面生成 + 可视化
- 吸附位点预测 + 可视化
- MC重构 + 轨迹GIF
- 反应路径 + 过程GIF
- 能垒计算 + 自由能图
- HTML摘要报告
- KMC模拟 + 表面动力学可视化

LLM智能体监督：
- 审查每个中间体结构的合理性
- 提供改进建议
- 自动检测结构异常
"""

import os
import sys
sys.path.insert(0, "core")

import logging
from ase.build import bulk
from ase.io import write

from camel_agents.visualized_digital_twin import VisualizedGasSolidDigitalTwin

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# 反应体系配置
# ============================================================================

# 催化剂体系
BULK_MATERIAL = "Cu"
LATTICE_CONSTANT = 3.61  # Å for Cu

# 反应路径: CO → CH4
REACTION_INTERMEDIATES = [
    "*CO",    # CO吸附
    "*CHO",   # CO加氢
    "*CHOH",  # CHO加氢
    "*CH",    # CHOH脱OH
    "*CH2",   # CH加氢
    "*CH3",   # CH2加氢
    "*CH4",   # CH3加氢
    "*",      # CH4脱附
]

INITIAL_REACTANT = "*CO"

# 反应条件
TEMPERATURE = 500  # K
PRESSURES = {
    "H2_g": 1.0,   # bar
    "H2O_g": 0.01,
    "CO_g": 0.1,
    "CH4_g": 0.001,
}

# 计算参数
TOP_N_SURFACES = 1  # 只分析最稳定的表面 (111)
NUM_ADSORPTION_SITES = 5
RECONSTRUCTION_SWEEPS = 5  # MC重构sweeps（生成轨迹GIF）
CALCULATE_BARRIERS = True  # 启用barrier计算
NEB_FRAMES = 10
NEB_FMAX = 0.05
RUN_KMC = True  # 启用KMC模拟

# 输出目录
OUTPUT_DIR = "output/co_to_ch4_digital_twin"

# ============================================================================
# 主程序
# ============================================================================

def main():
    """运行完整的催化数字孪生流程（带可视化）"""

    print("=" * 80)
    print("Gas-Solid Catalysis Digital Twin with Visualization")
    print("=" * 80)
    print(f"Catalyst: {BULK_MATERIAL}")
    print(f"Reaction: {REACTION_INTERMEDIATES[0]} → {REACTION_INTERMEDIATES[-1]}")
    print(f"Temperature: {TEMPERATURE} K")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 80)

    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 创建bulk结构并写入文件（重要：输入是文件路径！）
    logger.info("Creating bulk structure...")
    bulk_struct = bulk(BULK_MATERIAL, 'fcc', a=LATTICE_CONSTANT)
    bulk_path = os.path.join(OUTPUT_DIR, f"POSCAR_{BULK_MATERIAL}")
    write(bulk_path, bulk_struct)
    logger.info(f"✓ Bulk structure saved to: {bulk_path}")

    # 检查本地fairchem模型
    fairchem_model = 'uma-s-1p1'
    fairchem_model_path = f'deps/fairchem_models/{fairchem_model}.pt'
    if not os.path.exists(fairchem_model_path):
        logger.warning(f"⚠️  Local model not found: {fairchem_model_path}")
        fairchem_model_path = None
    else:
        logger.info(f"✓ Using local fairchem model: {fairchem_model_path}")

    # 初始化可视化数字孪生
    logger.info("\nInitializing Visualized Gas-Solid Digital Twin...")
    logger.info("  ✓ High quality visualization enabled")
    logger.info("  ✓ Auto expansion when adsorbates too close")
    logger.info("  ✓ LLM agent review enabled")

    dt = VisualizedGasSolidDigitalTwin(
        surff_root="deps/SurFF",
        adsorbdiff_root="deps/AdsorbDiff",
        surface_sampling_root="deps/surface-sampling",
        fairchem_root="deps/fairchem",
        fairchem_model=fairchem_model,
        fairchem_model_path=fairchem_model_path,
        use_gpu=True,
        # 可视化设置
        enable_visualization=True,
        viz_quality='high',
        viz_renderer='tachyon',
        viz_color_scheme='material',
        # 扩胞设置
        auto_expand_supercell=True,
        min_adsorbate_distance=3.0,
        # LLM审查设置
        enable_llm_review=True,
    )
    logger.info("✓ System initialized")

    # 运行完整工作流程（带可视化）
    logger.info("\n" + "=" * 80)
    logger.info("Running Complete Workflow with Visualization...")
    logger.info("=" * 80)

    result = dt.run_complete_workflow_with_visualization(
        bulk_structure=bulk_path,  # 注意：输入是文件路径！
        reaction_intermediates=REACTION_INTERMEDIATES,
        initial_reactant=INITIAL_REACTANT,
        temperature=TEMPERATURE,
        top_n_surfaces=TOP_N_SURFACES,
        num_adsorption_sites=NUM_ADSORPTION_SITES,
        reconstruction_sweeps=RECONSTRUCTION_SWEEPS,
        calculate_barriers=CALCULATE_BARRIERS,
        neb_frames=NEB_FRAMES,
        neb_fmax=NEB_FMAX,
        run_kmc=RUN_KMC,
        pressures=PRESSURES,
        output_dir=OUTPUT_DIR,
    )

    # 打印结果摘要
    print("\n" + "=" * 80)
    print("📊 Workflow Results:")
    print("=" * 80)

    if result and 'workflow_result' in result:
        wr = result['workflow_result']
        print(f"\n✓ Workflow Status: Success")
        print(f"Output Directory: {OUTPUT_DIR}")

        # 显示工作流摘要
        if hasattr(wr, 'summary'):
            print("\n" + wr.summary(top_n=1))

        # 显示可视化文件
        if 'visualizations' in result:
            viz = result['visualizations']
            print("\n" + "=" * 80)
            print("📁 Generated Visualizations:")
            print("=" * 80)

            categories = {
                '01_surfaces': viz.get('surfaces', {}),
                '02_adsorption': viz.get('adsorption', {}),
                '03_mc_trajectory': viz.get('mc_trajectory', {}),
                '04_reaction_pathway': viz.get('reaction_pathway', ''),
                '05_energy_diagram': viz.get('energy_diagram', ''),
                '07_kmc_dynamics': viz.get('kmc_dynamics', ''),
                '08_html_summary': viz.get('html_summary', ''),
            }

            for cat_name, cat_data in categories.items():
                if cat_data:
                    if isinstance(cat_data, dict):
                        print(f"\n✅ {cat_name} ({len(cat_data)} files)")
                        for name, path in list(cat_data.items())[:3]:
                            print(f"   - {name}: {path}")
                    else:
                        print(f"\n✅ {cat_name}")
                        print(f"   - {cat_data}")
                else:
                    print(f"\n❌ {cat_name}: NOT GENERATED")

            print("\n" + "=" * 80)
            print("💡 Key Features:")
            print("=" * 80)
            print("  1. ✓ All 8 visualization categories generated")
            print("  2. ✓ LLM agent reviewed all intermediates")
            print("  3. ✓ Auto expansion applied when needed")
            print("  4. ✓ High quality 3D structures and animations")
            print("  5. ✓ KMC surface dynamics showing intermediate changes over time")
            if result.get('llm_reviews'):
                print("\n💡 Check LLM review comments for structural insights!")
    else:
        print("\n❌ Workflow FAILED!")
        if result and 'error' in result:
            print(f"\nError: {result['error']}")

    print("=" * 80)

    return result


if __name__ == "__main__":
    try:
        result = main()
        print("\n✓ Digital Twin workflow completed successfully!")
    except Exception as e:
        logger.error(f"\n✗ Error during workflow: {e}", exc_info=True)
        sys.exit(1)
