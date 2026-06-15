#!/usr/bin/env python
"""
CO到CH4催化数字孪生测试（简化版）

只运行surface generation和adsorption prediction，跳过pathway和barrier计算
"""

import os
import sys
sys.path.insert(0, "core")

import logging
from ase.build import bulk

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin

# ============================================================================
# 配置参数
# ============================================================================

# 催化剂体系
BULK_MATERIAL = "Cu"
LATTICE_CONSTANT = 3.61  # Å for Cu

# 输出目录
OUTPUT_DIR = "output/co_to_ch4_simple"

def main():
    """运行简化的催化数字孪生流程"""

    print("=" * 80)
    print("Gas-Solid Catalysis Digital Twin (Simplified)")
    print("=" * 80)
    print(f"Catalyst: {BULK_MATERIAL}")
    print(f"Adsorbate: CO")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 80)

    # 创建bulk结构
    logger.info("Creating bulk structure...")
    bulk_struct = bulk(BULK_MATERIAL, 'fcc', a=LATTICE_CONSTANT)

    # 初始化数字孪生
    logger.info("Initializing Gas-Solid Digital Twin...")
    dt = GasSolidDigitalTwin(
        surff_root="deps/SurFF",
        adsorbdiff_root="deps/AdsorbDiff",
        surface_sampling_root="deps/surface-sampling",
        fairchem_root="deps/fairchem",
        fairchem_model="uma-s-1p1",
        vssr_mc_model="CHGNetNFF",
        use_gpu=True,
        logger=logger,
    )

    # Step 1: 生成表面
    logger.info("\nStep 1: Generating surfaces from bulk...")
    surff_result, slabs = dt.generate_surfaces_from_bulk(
        bulk_structure=bulk_struct,
        top_n=1,  # 只生成Cu(111)表面
        output_dir=f"{OUTPUT_DIR}/surfaces",
    )
    logger.info(f"✓ Generated {len(slabs)} surface(s)")

    # Step 2: 预测吸附位点
    logger.info("\nStep 2: Predicting CO adsorption sites...")
    adsorption_result = dt.predict_adsorption_sites(
        surface=slabs[0],
        adsorbate="*CO",
        num_sites=5,
        output_dir=f"{OUTPUT_DIR}/adsorption",
    )
    logger.info(f"✓ Best adsorption energy: {adsorption_result.best_result.energy:.3f} eV")

    # 打印摘要
    print("\n" + "=" * 80)
    print("Workflow Summary")
    print("=" * 80)
    print(f"Surface: Cu(111)")
    print(f"  Surface energy: {surff_result.surfaces[0].surface_energy:.4f} eV/Å²")
    print(f"  Area fraction: {surff_result.surfaces[0].area_fraction:.2%}")
    print(f"\nCO Adsorption:")
    print(f"  Best site energy: {adsorption_result.best_result.energy:.3f} eV")
    print(f"  Number of sites tested: {adsorption_result.num_sites}")
    print(f"\nOutput saved to: {OUTPUT_DIR}/")
    print("=" * 80)

    return surff_result, adsorption_result


if __name__ == "__main__":
    try:
        surff_result, adsorption_result = main()
        print("\n✓ Simplified digital twin workflow completed successfully!")
    except Exception as e:
        logger.error(f"\n✗ Error during workflow: {e}", exc_info=True)
        sys.exit(1)
