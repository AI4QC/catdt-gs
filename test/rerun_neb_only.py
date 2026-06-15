#!/usr/bin/env python
"""
只重新运行NEB能垒计算 - 使用已有的中间体结构
使用改进的LLM提示词（claude-opus-4-5-20251101）

主要修改：
1. LLM现在看到完整的POSCAR结构
2. LLM知道哪些原子是吸附分子
3. H原子放置更接近产物位置（z方向只高0.5-1.0 Å）
"""

import os
import sys
sys.path.insert(0, "core")

import logging
from typing import List, Dict
from ase.io import read, write
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 配置
OUTPUT_DIR = "output/co_to_ch4_digital_twin"
TEMPERATURE = 500  # K
REACTION_INTERMEDIATES = ["*CO", "*CHO", "*CHOH", "*CH", "*CH2", "*CH3", "*CH4", "*"]

# NEB参数
NEB_FRAMES = 8        # 减少帧数加快计算
NEB_FMAX = 0.08       # 稍微放宽收敛标准
NEB_MAX_STEPS = 150   # 适中的步数


def main():
    """只重新运行NEB能垒计算"""

    print("=" * 80)
    print("Rerun NEB Barrier Calculations Only")
    print("Using Improved LLM Prompts (claude-opus-4-5-20251101)")
    print("=" * 80)

    # 1. 加载中间体结构
    logger.info("Loading intermediate structures...")

    pathway_dir = os.path.join(OUTPUT_DIR, "surface_111", "pathway")
    adsorption_dir = os.path.join(pathway_dir, "adsorption")

    if not os.path.exists(adsorption_dir):
        logger.error(f"Pathway results not found: {adsorption_dir}")
        return

    structures = {}
    for inter in REACTION_INTERMEDIATES:
        inter_file = os.path.join(adsorption_dir, f"{inter}_relaxed.vasp")
        if os.path.exists(inter_file):
            structures[inter] = read(inter_file)
            logger.info(f"  ✓ Loaded {inter}: {len(structures[inter])} atoms")
        else:
            logger.warning(f"  ⚠ Missing {inter}")

    if len(structures) < 2:
        logger.error("Not enough structures for NEB calculation")
        return

    # 2. 初始化barrier predictor
    logger.info("\nInitializing barrier predictor...")

    neb_output_dir = os.path.join(OUTPUT_DIR, "neb_recalc")
    os.makedirs(neb_output_dir, exist_ok=True)

    # 强制设置LLM模型为opus
    os.environ["OPENAI_MODEL"] = "claude-opus-4-5-20251101"

    from pathway.fairchem_predictor import FairchemPredictor
    from pathway.barrier_predictor import BarrierPredictor

    # 检查本地fairchem模型
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

    # 3. 运行NEB计算
    logger.info("\n" + "=" * 80)
    logger.info("Running NEB Barrier Calculations")
    logger.info("=" * 80)

    results = {}

    for i in range(len(REACTION_INTERMEDIATES) - 1):
        reactant_name = REACTION_INTERMEDIATES[i]
        product_name = REACTION_INTERMEDIATES[i + 1]

        if reactant_name not in structures or product_name not in structures:
            logger.warning(f"Skipping {reactant_name} -> {product_name}: missing structures")
            continue

        reaction_name = f"{reactant_name}_to_{product_name}"
        logger.info(f"\n{'='*60}")
        logger.info(f"Step {i+1}: {reactant_name} -> {product_name}")
        logger.info(f"{'='*60}")

        reactant = structures[reactant_name].copy()
        product = structures[product_name].copy()

        logger.info(f"  Reactant: {len(reactant)} atoms, {reactant.get_chemical_formula()}")
        logger.info(f"  Product: {len(product)} atoms, {product.get_chemical_formula()}")

        try:
            # 获取吸附分子索引
            # 由于VASP格式不保留ASE tags，我们根据元素类型识别吸附分子
            # Cu是表面原子，其他（C, H, O等）是吸附分子
            def get_adsorbate_indices(atoms):
                """根据元素类型获取吸附分子索引"""
                symbols = atoms.get_chemical_symbols()
                surface_elements = {'Cu'}  # 表面元素
                indices = [i for i, s in enumerate(symbols) if s not in surface_elements]
                return indices

            reactant_adsorbate_indices = get_adsorbate_indices(reactant)
            product_adsorbate_indices = get_adsorbate_indices(product)

            logger.info(f"  Reactant adsorbate indices: {reactant_adsorbate_indices}")
            logger.info(f"  Product adsorbate indices: {product_adsorbate_indices}")

            # 运行NEB
            result = barrier_predictor.predict_from_structures(
                reactant=reactant,
                product=product,
                n_frames=NEB_FRAMES,
                fmax=NEB_FMAX,
                max_steps=NEB_MAX_STEPS,
                relax_endpoints=False,  # 端点已经优化过了
                reaction_name=reaction_name,
                reactant_adsorbate_indices=reactant_adsorbate_indices,
                product_adsorbate_indices=product_adsorbate_indices,
            )

            results[reaction_name] = result

            logger.info(f"\n  ✓ NEB completed for {reaction_name}")
            logger.info(f"    E_act (forward) = {result.activation_energy_forward:.3f} eV")
            logger.info(f"    E_act (reverse) = {result.activation_energy_reverse:.3f} eV")
            logger.info(f"    ΔE = {result.reaction_energy:.3f} eV")
            logger.info(f"    TS at frame {result.transition_state_index}")
            logger.info(f"    Converged: {result.converged}")

            # 保存NEB轨迹
            traj_file = os.path.join(neb_output_dir, f"{reaction_name}_neb.traj")
            from ase.io import Trajectory
            with Trajectory(traj_file, 'w') as traj:
                for frame in result.neb_frames:
                    traj.write(frame)
            logger.info(f"    Saved trajectory: {traj_file}")

        except Exception as e:
            logger.error(f"  ✗ NEB failed for {reaction_name}: {e}")
            import traceback
            logger.debug(traceback.format_exc())

    # 4. 打印结果汇总
    print("\n" + "=" * 80)
    print("NEB Barrier Calculation Results")
    print("=" * 80)

    print(f"\n{'Reaction':<25} {'E_act (eV)':<12} {'ΔE (eV)':<12} {'Converged':<10}")
    print("-" * 60)

    for name, result in results.items():
        print(f"{name:<25} {result.activation_energy_forward:>10.3f} {result.reaction_energy:>10.3f} {'Yes' if result.converged else 'No':>10}")

    # 找出RDS
    if results:
        rds_name = max(results, key=lambda x: results[x].activation_energy_forward)
        rds_barrier = results[rds_name].activation_energy_forward
        print("-" * 60)
        print(f"\nRate Determining Step: {rds_name}")
        print(f"Maximum Barrier: {rds_barrier:.3f} eV")

    # 保存结果到CSV
    csv_file = os.path.join(neb_output_dir, "barrier_results.csv")
    with open(csv_file, 'w') as f:
        f.write("Reaction,E_act_forward,E_act_reverse,Delta_E,TS_frame,Converged\n")
        for name, result in results.items():
            f.write(f"{name},{result.activation_energy_forward:.4f},{result.activation_energy_reverse:.4f},{result.reaction_energy:.4f},{result.transition_state_index},{result.converged}\n")
    print(f"\nResults saved to: {csv_file}")

    print("\n" + "=" * 80)
    print("NEB Recalculation Complete!")
    print("=" * 80)

    return results


if __name__ == "__main__":
    results = main()
