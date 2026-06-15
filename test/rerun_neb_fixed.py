#!/usr/bin/env python
"""
重新运行NEB能垒计算 - 使用修复后的代码
测试改进的H原子放置（0.6 Å偏移 vs 旧的2.5 Å）
"""

import os
import sys
sys.path.insert(0, "core")

import logging
from ase.io import read, write
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 配置
OUTPUT_DIR = "output/co_to_ch4_digital_twin"

# 强制设置LLM模型为opus
os.environ["OPENAI_MODEL"] = "claude-opus-4-5-20251101"

def main():
    """测试单步NEB计算（*CO -> *CHO）"""

    print("=" * 80)
    print("Test NEB Barrier Calculation with Fixed H Positioning")
    print("Using 0.6 Å z-offset (instead of old 2.5 Å)")
    print("=" * 80)

    # 1. 加载结构
    logger.info("Loading structures...")

    pathway_dir = os.path.join(OUTPUT_DIR, "surface_111", "pathway")
    adsorption_dir = os.path.join(pathway_dir, "adsorption")

    reactant = read(os.path.join(adsorption_dir, "*CO_relaxed.vasp"))
    product = read(os.path.join(adsorption_dir, "*CHO_relaxed.vasp"))

    logger.info(f"  Reactant: {len(reactant)} atoms, {reactant.get_chemical_formula()}")
    logger.info(f"  Product: {len(product)} atoms, {product.get_chemical_formula()}")

    # 检查产物中H的位置
    prod_symbols = product.get_chemical_symbols()
    prod_positions = product.get_positions()
    for i, (sym, pos) in enumerate(zip(prod_symbols, prod_positions)):
        if sym == 'H':
            logger.info(f"  Product H at index {i}: z = {pos[2]:.3f} Å")
            logger.info(f"  Expected reactant H z (with fix): {pos[2] + 0.6:.3f} Å")

    # 2. 初始化predictor
    logger.info("\nInitializing barrier predictor...")

    neb_output_dir = os.path.join(OUTPUT_DIR, "neb_fixed_test")
    os.makedirs(neb_output_dir, exist_ok=True)

    from pathway.fairchem_predictor import FairchemPredictor
    from pathway.barrier_predictor import BarrierPredictor

    fairchem_model_path = 'deps/fairchem_models/uma-s-1p1.pt'
    if not os.path.exists(fairchem_model_path):
        logger.error(f"FairChem model not found: {fairchem_model_path}")
        return

    fairchem_predictor = FairchemPredictor(
        fairchem_root="deps/fairchem",
        model_name="uma-s-1p1",
        model_path=fairchem_model_path,
        device="cuda",
    )

    barrier_predictor = BarrierPredictor(
        fairchem_predictor=fairchem_predictor,
        work_dir=neb_output_dir,
    )

    # 3. 运行NEB (只测试第一步)
    logger.info("\n" + "=" * 60)
    logger.info("Running NEB: *CO -> *CHO")
    logger.info("=" * 60)

    # 获取吸附分子索引
    def get_adsorbate_indices(atoms):
        symbols = atoms.get_chemical_symbols()
        return [i for i, s in enumerate(symbols) if s not in {'Cu'}]

    reactant_adsorbate_indices = get_adsorbate_indices(reactant)
    product_adsorbate_indices = get_adsorbate_indices(product)

    logger.info(f"  Reactant adsorbate indices: {reactant_adsorbate_indices}")
    logger.info(f"  Product adsorbate indices: {product_adsorbate_indices}")

    # 使用较少的帧数和步数进行快速测试
    result = barrier_predictor.predict_from_structures(
        reactant=reactant.copy(),
        product=product.copy(),
        n_frames=8,       # 减少帧数加快测试
        fmax=0.1,         # 放宽收敛标准
        max_steps=100,    # 减少步数
        relax_endpoints=False,
        reaction_name="*CO_to_*CHO_fixed",
        reactant_adsorbate_indices=reactant_adsorbate_indices,
        product_adsorbate_indices=product_adsorbate_indices,
    )

    # 4. 打印结果
    print("\n" + "=" * 60)
    print("NEB Result for *CO -> *CHO (with fixed H positioning)")
    print("=" * 60)
    print(f"E_act (forward) = {result.activation_energy_forward:.3f} eV")
    print(f"E_act (reverse) = {result.activation_energy_reverse:.3f} eV")
    print(f"ΔE = {result.reaction_energy:.3f} eV")
    print(f"TS at frame {result.transition_state_index}")
    print(f"Converged: {result.converged}")
    print(f"Final fmax: {result.fmax_final:.4f} eV/Å")
    print(f"Optimization steps: {result.optimization_steps}")

    # 检查能垒是否合理
    if result.activation_energy_forward > 0:
        print(f"\n✓ SUCCESS: Found non-zero forward barrier!")
        print(f"  Expected range for CO hydrogenation: 0.3-1.0 eV")
        print(f"  Calculated: {result.activation_energy_forward:.3f} eV")
    else:
        print(f"\n⚠ WARNING: Forward barrier is still 0 eV")
        print("  This may indicate NEB path issues")

    # 保存NEB轨迹
    if result.neb_frames:
        traj_file = os.path.join(neb_output_dir, "*CO_to_*CHO_fixed_neb.traj")
        from ase.io import Trajectory
        with Trajectory(traj_file, 'w') as traj:
            for frame in result.neb_frames:
                traj.write(frame)
        logger.info(f"\nSaved trajectory: {traj_file}")

        # 打印能量曲线
        if result.energies:
            print("\nEnergy profile:")
            e0 = result.energies[0]
            for i, e in enumerate(result.energies):
                rel_e = e - e0
                bar = "#" * int(max(0, rel_e) * 20)
                print(f"  Frame {i}: {rel_e:+.3f} eV {bar}")

    print("\n" + "=" * 60)
    print("Test Complete!")
    print("=" * 60)

    return result


if __name__ == "__main__":
    result = main()
