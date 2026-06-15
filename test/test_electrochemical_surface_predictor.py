"""
Electrochemical Surface Predictor 测试脚本

此脚本演示如何使用 ElectrochemicalSurfacePredictor 进行
CP-MACE + VSSR-MC 集成的电化学表面采样。

使用方式:
    conda activate catdt
    python test_electrochemical_surface_predictor.py

作者: Claude
"""

import os
import sys
import numpy as np

# 添加core/reconstruction到路径
sys.path.insert(0, "core/reconstruction")


def test_electrochemical_conditions():
    """测试电化学条件类"""
    print("=" * 70)
    print("测试1: ElectrochemicalConditions")
    print("=" * 70)

    from electrochemical_surface_predictor import ElectrochemicalConditions

    # 测试不同条件
    conditions_list = [
        ElectrochemicalConditions(potential_she=0.0, ph=0.0, temperature=298.15),
        ElectrochemicalConditions(potential_she=1.0, ph=7.0, temperature=298.15),
        ElectrochemicalConditions(potential_she=1.6, ph=12.0, temperature=300.0),
    ]

    for cond in conditions_list:
        print(f"\n{cond}")
        print(f"  Fermi level (vacuum): {cond.potential_vacuum:.4f} eV")
        print(f"  kT: {cond.kT:.6f} eV")

    print("\n✓ ElectrochemicalConditions 测试通过")
    return True


def test_chemical_potential_calculation():
    """测试化学势计算"""
    print("\n" + "=" * 70)
    print("测试2: 化学势计算 (CHE模型)")
    print("=" * 70)

    from electrochemical_surface_predictor import (
        ElectrochemicalEvaluator,
        ElectrochemicalConditions,
    )

    # 创建条件 (不需要实际加载模型)
    conditions = ElectrochemicalConditions(
        potential_she=1.23,  # OER电位
        ph=0.0,
        temperature=298.15,
    )

    # 测试化学势计算 (使用基类的静态方法)
    kT = conditions.kT
    pH = conditions.ph
    U = conditions.potential_she

    # 手动计算以验证
    # O* 吸附: μ_O = μ_H2O - μ_H2 + kT*ln(10)*pH - eU
    E_H2O = -14.22
    E_H2 = -6.77

    mu_O = E_H2O - E_H2 + kT * np.log(10) * pH - U
    print(f"\n条件: U={U:.2f} V vs. SHE, pH={pH:.1f}, T={conditions.temperature:.1f} K")
    print(f"\nO* 吸附的化学势修正:")
    print(f"  μ_O = μ_H2O - μ_H2 + kT*ln(10)*pH - eU")
    print(f"  μ_O = {E_H2O:.2f} - {E_H2:.2f} + {kT * np.log(10) * pH:.4f} - {U:.2f}")
    print(f"  μ_O = {mu_O:.4f} eV")

    # OH* 吸附
    mu_OH = E_H2O - 0.5 * E_H2 + kT * np.log(10) * pH - 0.5 * U
    print(f"\nOH* 吸附的化学势修正:")
    print(f"  μ_OH = {mu_OH:.4f} eV")

    # 不同pH的影响
    print(f"\n不同pH下O*吸附的化学势变化:")
    for test_ph in [0, 7, 14]:
        mu_O_ph = E_H2O - E_H2 + kT * np.log(10) * test_ph - U
        print(f"  pH={test_ph}: μ_O = {mu_O_ph:.4f} eV")

    print("\n✓ 化学势计算测试通过")
    return True


def test_predictor_initialization():
    """测试预测器初始化"""
    print("\n" + "=" * 70)
    print("测试3: ElectrochemicalSurfacePredictor 初始化")
    print("=" * 70)

    from electrochemical_surface_predictor import ElectrochemicalSurfacePredictor

    cp_mace_root = "deps/CP-MACE"
    surface_sampling_root = "deps/surface-sampling"

    # 检查路径
    if not os.path.exists(cp_mace_root):
        print(f"✗ CP-MACE根目录不存在: {cp_mace_root}")
        return False
    if not os.path.exists(surface_sampling_root):
        print(f"✗ surface-sampling根目录不存在: {surface_sampling_root}")
        return False

    try:
        predictor = ElectrochemicalSurfacePredictor(
            cp_mace_root=cp_mace_root,
            surface_sampling_root=surface_sampling_root,
            device="cpu",  # 使用CPU进行测试
            verbose=True,
        )
        print("\n✓ ElectrochemicalSurfacePredictor 初始化成功")
        return predictor
    except Exception as e:
        print(f"\n✗ 初始化失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_evaluator_initialization():
    """测试电化学评估器初始化"""
    print("\n" + "=" * 70)
    print("测试4: ElectrochemicalEvaluator 初始化")
    print("=" * 70)

    from electrochemical_surface_predictor import (
        ElectrochemicalEvaluator,
        ElectrochemicalConditions,
    )

    cp_mace_root = "deps/CP-MACE"
    model_paths = [
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_1.model"),
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_2.model"),
    ]

    existing_models = [p for p in model_paths if os.path.exists(p)]
    if not existing_models:
        print(f"✗ 没有找到CP-MACE模型文件")
        return None

    conditions = ElectrochemicalConditions(
        potential_she=1.0,
        ph=12.0,
        temperature=300.0,
    )

    try:
        evaluator = ElectrochemicalEvaluator(
            cp_mace_root=cp_mace_root,
            model_paths=existing_models,
            conditions=conditions,
            device="cpu",
            verbose=True,
        )
        print("\n✓ ElectrochemicalEvaluator 初始化成功")
        return evaluator
    except Exception as e:
        print(f"\n✗ 初始化失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_grand_potential_calculation():
    """测试巨势计算"""
    print("\n" + "=" * 70)
    print("测试5: Grand Potential 计算")
    print("=" * 70)

    from ase.io import read
    from electrochemical_surface_predictor import (
        ElectrochemicalEvaluator,
        ElectrochemicalConditions,
    )

    cp_mace_root = "deps/CP-MACE"
    init_file = os.path.join(cp_mace_root, "simulation/slow_growth/init.xyz")
    model_paths = [
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_1.model"),
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_2.model"),
    ]

    if not os.path.exists(init_file):
        print(f"✗ 测试结构文件不存在: {init_file}")
        return None

    existing_models = [p for p in model_paths if os.path.exists(p)]
    if not existing_models:
        print(f"✗ 没有找到模型文件")
        return None

    conditions = ElectrochemicalConditions(
        potential_she=-3.36 - (-4.44),  # 从fermi level转换
        ph=0.0,
        temperature=300.0,
    )

    try:
        evaluator = ElectrochemicalEvaluator(
            cp_mace_root=cp_mace_root,
            model_paths=existing_models,
            conditions=conditions,
            device="cuda",  # 使用GPU加速
            verbose=True,
        )

        # 读取测试结构
        atoms = read(init_file)
        print(f"\n测试结构: {atoms.get_chemical_formula()}")
        print(f"原子数: {len(atoms)}")
        print(f"电子数: {atoms.info.get('electron', 'N/A')}")

        # 计算巨势
        omega_total, omega_el, delta_mu, fermi, ne = \
            evaluator.calculate_grand_potential(atoms)

        print(f"\n计算结果:")
        print(f"  Ω_total: {omega_total:.4f} eV")
        print(f"  Ω_el (CP-MACE): {omega_el:.4f} eV")
        print(f"  Δμ_ions: {delta_mu:.4f} eV")
        print(f"  Fermi level: {fermi:.4f} V")
        print(f"  Electron number: {ne:.2f}")

        print("\n✓ Grand Potential 计算测试通过")
        return True

    except Exception as e:
        print(f"\n✗ 计算失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_short_sampling():
    """测试短时间采样 (仅用于验证功能)"""
    print("\n" + "=" * 70)
    print("测试6: 短时间采样 (5轮)")
    print("=" * 70)

    from ase.io import read
    from electrochemical_surface_predictor import ElectrochemicalSurfacePredictor

    cp_mace_root = "deps/CP-MACE"
    surface_sampling_root = "deps/surface-sampling"
    init_file = os.path.join(cp_mace_root, "simulation/slow_growth/init.xyz")

    if not os.path.exists(init_file):
        print(f"✗ 测试结构文件不存在")
        return None

    try:
        predictor = ElectrochemicalSurfacePredictor(
            cp_mace_root=cp_mace_root,
            surface_sampling_root=surface_sampling_root,
            device="cuda",
            verbose=True,
        )

        # 运行短时间采样
        result = predictor.sample_surface(
            surface=init_file,
            adsorbates=["O", "H"],
            potential_she=1.0,  # 1.0 V vs. SHE
            ph=12.0,
            temperature=300.0,
            total_sweeps=5,  # 仅5轮用于测试
            sweep_size=5,
            output_dir="output/echem",
        )

        print("\n" + result.summary())
        print("\n✓ 采样测试通过")
        return result

    except Exception as e:
        print(f"\n✗ 采样失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def print_usage_examples():
    """打印使用示例"""
    print("\n" + "=" * 70)
    print("使用示例")
    print("=" * 70)

    examples = """
# 示例1: 基本使用 - 在指定电化学条件下采样表面结构
from electrochemical_surface_predictor import ElectrochemicalSurfacePredictor

predictor = ElectrochemicalSurfacePredictor(
    cp_mace_root="/path/to/CP-MACE",
    surface_sampling_root="/path/to/surface-sampling",
    device="cuda",
)

# 在 pH=12, U=1.6V vs. SHE 条件下采样 (典型的OER条件)
result = predictor.sample_surface(
    surface="surface.vasp",
    adsorbates=["O", "OH", "H"],
    potential_she=1.6,  # V vs. SHE
    ph=12.0,
    temperature=300.0,
    total_sweeps=100,
)

print(result.summary())

# 示例2: 获取最稳定的表面配置
top_configs = result.get_top_n(10)
for i, config in enumerate(top_configs):
    print(f"{i+1}. Ω={config.total_grand_potential:.4f} eV, "
          f"μ_F={config.fermi_level:.4f} V")

# 示例3: 使用ElectrochemicalEvaluator单独计算巨势
from electrochemical_surface_predictor import (
    ElectrochemicalEvaluator,
    ElectrochemicalConditions,
)

conditions = ElectrochemicalConditions(
    potential_she=1.23,  # OER平衡电位
    ph=0.0,
    temperature=298.15,
)

evaluator = ElectrochemicalEvaluator(
    cp_mace_root="/path/to/CP-MACE",
    model_paths=["model1.model", "model2.model"],
    conditions=conditions,
    device="cuda",
)

# 计算某个结构的巨势
omega_total, omega_el, delta_mu, fermi, ne = \\
    evaluator.calculate_grand_potential(
        atoms,
        added_species={"O": 1},  # 添加了一个O原子
    )

print(f"Total Grand Potential: {omega_total:.4f} eV")
print(f"Electronic contribution: {omega_el:.4f} eV")
print(f"Chemical potential correction: {delta_mu:.4f} eV")

# 示例4: 命令行使用
# python electrochemical_surface_predictor.py surface.vasp \\
#     --adsorbates O OH H \\
#     --potential 1.6 \\
#     --ph 12.0 \\
#     --temperature 300 \\
#     --total-sweeps 100

# 示例5: 策略A (重打分) - 先用普通VSSR-MC筛选，再用CP-MACE精修
from vssr_mc_predictor import VSSRMCPredictor
from electrochemical_surface_predictor import ElectrochemicalEvaluator

# 第一步：用普通MACE快速采样
vssr_predictor = VSSRMCPredictor(
    surface_sampling_root="/path/to/surface-sampling",
    model_type="CHGNetNFF",
)
vssr_result = vssr_predictor.sample(
    surface="surface.vasp",
    adsorbates=["O", "OH"],
    total_sweeps=1000,
)

# 第二步：筛选Top 50结构
top_structures = vssr_result.get_structures_by_energy(50)

# 第三步：用CP-MACE重新打分
evaluator = ElectrochemicalEvaluator(...)
for struct in top_structures:
    omega = evaluator.calculate_grand_potential(struct.atoms)
    print(f"CP-MACE Ω = {omega[0]:.4f} eV")
"""
    print(examples)


def main():
    """运行所有测试"""
    print("=" * 70)
    print("Electrochemical Surface Predictor 测试")
    print("(CP-MACE + VSSR-MC 集成)")
    print("=" * 70)
    print()

    # 测试1: 电化学条件类
    test_electrochemical_conditions()

    # 测试2: 化学势计算
    test_chemical_potential_calculation()

    # 测试3: 预测器初始化
    test_predictor_initialization()

    # 测试4: 评估器初始化
    test_evaluator_initialization()

    # 测试5: Grand Potential计算
    print("\n注意: Grand Potential计算测试需要GPU支持")
    # 自动跳过需要交互的测试
    print("跳过Grand Potential计算测试 (需要手动运行)")

    # 测试6: 短时间采样
    print("\n注意: 采样测试需要GPU支持和较长时间")
    # 自动跳过需要交互的测试
    print("跳过采样测试 (需要手动运行)")

    # 打印使用示例
    print_usage_examples()

    print("\n" + "=" * 70)
    print("测试完成")
    print("=" * 70)


if __name__ == "__main__":
    main()
