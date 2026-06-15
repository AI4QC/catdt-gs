#!/usr/bin/env python
"""
测试LLM控制的NEB结构准备

这个脚本演示如何使用新的LLM控制器来准备任意反应的NEB计算结构。

LLM完全控制：
1. 分析反应类型（氢化、脱氢、键断裂等）
2. 决定是否需要添加/删除原子
3. 确定添加原子的位置
4. 验证准备好的结构是否合理
5. 如果不合理则迭代改进

关键原则：吸附位点已固定，LLM只负责原子数匹配和初始路径构建
"""

import os
import sys
sys.path.insert(0, "core")

import logging
from ase.io import read

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 配置
OUTPUT_DIR = "output/co_to_ch4_digital_twin"
USE_LLM_CONTROLLER = True  # 启用LLM控制器

# 强制设置LLM模型为opus
os.environ["OPENAI_MODEL"] = "claude-opus-4-5-20251101"


def test_llm_controller_single_step():
    """测试单步NEB计算（*CO -> *CHO）使用LLM控制器"""

    print("=" * 80)
    print("Test LLM-Controlled NEB Structure Preparation")
    print("Reaction: *CO + H → *CHO")
    print("=" * 80)

    # 1. 加载结构
    logger.info("Loading structures...")

    pathway_dir = os.path.join(OUTPUT_DIR, "surface_111", "pathway")
    adsorption_dir = os.path.join(pathway_dir, "adsorption")

    reactant = read(os.path.join(adsorption_dir, "*CO_relaxed.vasp"))
    product = read(os.path.join(adsorption_dir, "*CHO_relaxed.vasp"))

    logger.info(f"  Reactant: {len(reactant)} atoms, {reactant.get_chemical_formula()}")
    logger.info(f"  Product: {len(product)} atoms, {product.get_chemical_formula()}")

    # 获取吸附分子索引
    def get_adsorbate_indices(atoms):
        symbols = atoms.get_chemical_symbols()
        return [i for i, s in enumerate(symbols) if s not in {'Cu'}]

    reactant_adsorbate_indices = get_adsorbate_indices(reactant)
    product_adsorbate_indices = get_adsorbate_indices(product)

    logger.info(f"  Reactant adsorbate indices: {reactant_adsorbate_indices}")
    logger.info(f"  Product adsorbate indices: {product_adsorbate_indices}")

    # 2. 初始化带LLM控制器的barrier predictor
    logger.info("\nInitializing barrier predictor with LLM controller...")

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

    # 关键：启用LLM控制器
    barrier_predictor = BarrierPredictor(
        fairchem_predictor=fairchem_predictor,
        work_dir="output/llm_neb_test",
        use_llm_controller=USE_LLM_CONTROLLER,  # 启用LLM控制器
        llm_model="claude-opus-4-5-20251101",   # 指定LLM模型
    )

    # 3. 运行NEB
    logger.info("\n" + "=" * 60)
    logger.info("Running NEB with LLM-controlled structure preparation")
    logger.info("=" * 60)

    result = barrier_predictor.predict_from_structures(
        reactant=reactant.copy(),
        product=product.copy(),
        n_frames=8,
        fmax=0.1,
        max_steps=100,
        relax_endpoints=False,
        reaction_name="*CO_to_*CHO",
        reactant_adsorbate_indices=reactant_adsorbate_indices,
        product_adsorbate_indices=product_adsorbate_indices,
    )

    # 4. 打印结果
    print("\n" + "=" * 60)
    print("NEB Result (*CO → *CHO)")
    print("=" * 60)
    print(f"E_act (forward) = {result.activation_energy_forward:.3f} eV")
    print(f"E_act (reverse) = {result.activation_energy_reverse:.3f} eV")
    print(f"ΔE = {result.reaction_energy:.3f} eV")
    print(f"TS at frame {result.transition_state_index}")
    print(f"Converged: {result.converged}")
    print(f"Optimization steps: {result.optimization_steps}")

    if result.activation_energy_forward > 0:
        print(f"\n✓ SUCCESS: Found non-zero forward barrier!")
        print(f"  LLM-controlled approach produced realistic barrier")
    else:
        print(f"\n⚠ WARNING: Forward barrier is 0 eV")

    # 5. 保存轨迹
    if result.neb_frames:
        from ase.io import Trajectory
        traj_file = "output/llm_neb_test/*CO_to_*CHO_llm.traj"
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


def test_general_reaction(reactant_file, product_file, reaction_name):
    """
    测试通用反应

    这个函数演示如何对任意反应使用LLM控制器

    参数：
        reactant_file: 反应物VASP文件路径
        product_file: 产物VASP文件路径
        reaction_name: 反应名称（供LLM参考）
    """
    print("=" * 80)
    print(f"Test General Reaction: {reaction_name}")
    print("=" * 80)

    # 加载结构
    reactant = read(reactant_file)
    product = read(product_file)

    logger.info(f"Reactant: {len(reactant)} atoms, {reactant.get_chemical_formula()}")
    logger.info(f"Product: {len(product)} atoms, {product.get_chemical_formula()}")

    # 自动识别吸附分子（假设表面是Cu）
    def get_adsorbate_indices(atoms, surface_element='Cu'):
        symbols = atoms.get_chemical_symbols()
        return [i for i, s in enumerate(symbols) if s != surface_element]

    reactant_ads_idx = get_adsorbate_indices(reactant)
    product_ads_idx = get_adsorbate_indices(product)

    logger.info(f"Reactant adsorbate: {[reactant.get_chemical_symbols()[i] for i in reactant_ads_idx]}")
    logger.info(f"Product adsorbate: {[product.get_chemical_symbols()[i] for i in product_ads_idx]}")

    # 初始化
    from pathway.fairchem_predictor import FairchemPredictor
    from pathway.barrier_predictor import BarrierPredictor

    fairchem_predictor = FairchemPredictor(
        fairchem_root="deps/fairchem",
        model_name="uma-s-1p1",
        model_path='deps/fairchem_models/uma-s-1p1.pt',
        device="cuda",
    )

    barrier_predictor = BarrierPredictor(
        fairchem_predictor=fairchem_predictor,
        work_dir=f"output/llm_neb_{reaction_name.replace(' ', '_')}",
        use_llm_controller=True,  # 启用LLM控制器
        llm_model="claude-opus-4-5-20251101",
    )

    # 运行NEB
    result = barrier_predictor.predict_from_structures(
        reactant=reactant,
        product=product,
        n_frames=8,
        fmax=0.1,
        max_steps=100,
        relax_endpoints=False,
        reaction_name=reaction_name,
        reactant_adsorbate_indices=reactant_ads_idx,
        product_adsorbate_indices=product_ads_idx,
    )

    print("\n" + "=" * 60)
    print(f"NEB Result: {reaction_name}")
    print("=" * 60)
    print(f"E_act (forward) = {result.activation_energy_forward:.3f} eV")
    print(f"ΔE = {result.reaction_energy:.3f} eV")
    print(f"Converged: {result.converged}")

    return result


if __name__ == "__main__":
    # 测试单步反应
    result = test_llm_controller_single_step()

    # 如果要测试其他反应，取消注释以下代码：
    # test_general_reaction(
    #     "path/to/reactant.vasp",
    #     "path/to/product.vasp",
    #     "My Reaction"
    # )
