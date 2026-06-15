"""
CP-MACE Predictor 测试脚本

此脚本演示如何使用CPMACEPredictor进行常电位分子动力学模拟。

使用方式:
    # 使用catdt环境
    conda activate catdt
    python test_cp_mace_predictor.py

作者: Claude
"""

import os
import sys
import numpy as np

# 添加core/reconstruction到路径
sys.path.insert(0, "core/reconstruction")

from cp_mace_predictor import CPMACEPredictor


def test_initialization():
    """测试初始化"""
    print("=" * 60)
    print("测试1: 初始化CPMACEPredictor")
    print("=" * 60)

    cp_mace_root = "deps/CP-MACE"

    try:
        predictor = CPMACEPredictor(
            cp_mace_root=cp_mace_root,
            use_gpu=True,
            verbose=True,
        )
        print("✓ CPMACEPredictor初始化成功")
        return predictor
    except Exception as e:
        print(f"✗ 初始化失败: {e}")
        return None


def test_structure_loading():
    """测试结构加载"""
    print("\n" + "=" * 60)
    print("测试2: 结构加载")
    print("=" * 60)

    cp_mace_root = "deps/CP-MACE"
    predictor = CPMACEPredictor(
        cp_mace_root=cp_mace_root,
        use_gpu=False,  # 使用CPU进行简单测试
        verbose=True,
    )

    # 测试从CP-MACE自带的示例文件加载
    test_file = os.path.join(cp_mace_root, "simulation/slow_growth/init.xyz")

    if os.path.exists(test_file):
        try:
            atoms = predictor._convert_structure(test_file)
            print(f"✓ 成功加载结构: {atoms.get_chemical_formula()}")
            print(f"  原子数: {len(atoms)}")
            print(f"  晶胞: {atoms.get_cell()}")
            if 'electron' in atoms.info:
                print(f"  电子数: {atoms.info['electron']}")
            return atoms
        except Exception as e:
            print(f"✗ 结构加载失败: {e}")
            return None
    else:
        print(f"✗ 测试文件不存在: {test_file}")
        return None


def test_model_loading():
    """测试模型加载"""
    print("\n" + "=" * 60)
    print("测试3: 模型加载")
    print("=" * 60)

    cp_mace_root = "deps/CP-MACE"
    predictor = CPMACEPredictor(
        cp_mace_root=cp_mace_root,
        use_gpu=True,
        verbose=True,
    )

    # 测试加载预训练模型
    model_paths = [
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_1.model"),
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_2.model"),
    ]

    existing_models = [p for p in model_paths if os.path.exists(p)]

    if existing_models:
        try:
            calculators = predictor._load_calculators(existing_models)
            print(f"✓ 成功加载 {len(calculators)} 个模型")
            return calculators
        except Exception as e:
            print(f"✗ 模型加载失败: {e}")
            import traceback
            traceback.print_exc()
            return None
    else:
        print(f"✗ 没有找到可用的模型文件")
        print(f"  期望路径: {model_paths}")
        return None


def test_short_simulation():
    """测试短时间模拟"""
    print("\n" + "=" * 60)
    print("测试4: 短时间模拟 (仅5步)")
    print("=" * 60)

    cp_mace_root = "deps/CP-MACE"

    # 检查必要文件
    init_file = os.path.join(cp_mace_root, "simulation/slow_growth/init.xyz")
    model_files = [
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_1.model"),
        os.path.join(cp_mace_root, "simulation/slow_growth/MACE_model_compiled_2.model"),
    ]

    if not os.path.exists(init_file):
        print(f"✗ 初始结构文件不存在: {init_file}")
        return None

    existing_models = [m for m in model_files if os.path.exists(m)]
    if not existing_models:
        print(f"✗ 没有找到模型文件")
        return None

    print(f"  初始结构: {init_file}")
    print(f"  模型数量: {len(existing_models)}")

    try:
        predictor = CPMACEPredictor(
            cp_mace_root=cp_mace_root,
            use_gpu=True,
            verbose=True,
        )

        # 运行非常短的模拟
        result = predictor.simulate(
            structure=init_file,
            model_paths=existing_models,
            target_potential=-3.36,  # 从inputs.yml中获取
            temperature=300.0,
            timestep=1.0,
            steps=5,  # 仅5步用于测试
            save_frequency=1,
            output_dir="output/cp_mace",
        )

        print("\n✓ 模拟完成!")
        print(f"  完成步数: {result.completed_steps}")
        print(f"  最终能量: {result.energy_history[-1]:.4f} eV" if result.energy_history else "  N/A")
        print(f"  最终费米能级: {result.fermi_history[-1]:.4f} V" if result.fermi_history else "  N/A")
        print(f"  输出目录: {result.output_dir}")

        return result

    except Exception as e:
        print(f"✗ 模拟失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_training_data_preparation():
    """测试训练数据准备"""
    print("\n" + "=" * 60)
    print("测试5: 训练数据准备")
    print("=" * 60)

    from ase import Atoms
    from ase.build import fcc111

    cp_mace_root = "deps/CP-MACE"
    predictor = CPMACEPredictor(
        cp_mace_root=cp_mace_root,
        use_gpu=False,
        verbose=True,
    )

    # 创建简单的测试结构
    slab = fcc111('Au', size=(2, 2, 3), vacuum=10.0)

    # 创建不同电子数的配置
    structures = [slab.copy() for _ in range(3)]
    electron_numbers = [100.0, 100.5, 101.0]
    potentials = [-3.0, -3.2, -3.4]
    energies = [-50.0, -50.1, -50.2]  # 假设的参考能量
    forces = [np.zeros((len(slab), 3)) for _ in range(3)]  # 假设的参考力

    try:
        output_file = "output/cp_mace/test_train.xyz"
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        result = predictor.prepare_training_data(
            structures=structures,
            electron_numbers=electron_numbers,
            potentials=potentials,
            energies=energies,
            forces=forces,
            output_file=output_file,
        )

        print(f"✓ 训练数据准备成功")
        print(f"  输出文件: {result}")

        # 验证文件内容
        from ase.io import read
        loaded = read(output_file, index=':')
        print(f"  加载的结构数: {len(loaded)}")
        print(f"  第一个结构的电子数: {loaded[0].info.get('electron', 'N/A')}")
        print(f"  第一个结构的电位: {loaded[0].info.get('potential', 'N/A')}")

        return result

    except Exception as e:
        print(f"✗ 训练数据准备失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def print_usage_examples():
    """打印使用示例"""
    print("\n" + "=" * 60)
    print("使用示例")
    print("=" * 60)

    examples = """
# 示例1: 使用预训练模型进行模拟
from cp_mace_predictor import CPMACEPredictor

predictor = CPMACEPredictor(
    cp_mace_root="deps/CP-MACE",
    use_gpu=True,
)

result = predictor.simulate(
    structure="init.xyz",  # 带有electron信息的结构
    model_paths=["model1.model", "model2.model"],
    target_potential=-3.36,  # 目标电极电位 (V)
    temperature=300.0,
    steps=1000,
)

print(result.summary())

# 示例2: 训练新模型
train_result = predictor.train(
    train_file="train.xyz",  # 包含electron和potential的训练数据
    model_name="my_model",
    max_num_epochs=300,
    potential_weight=10.0,
)

# 示例3: 慢增长模拟 (用于计算自由能)
result = predictor.run_slow_growth(
    structure="init.xyz",
    model_paths=["model1.model", "model2.model"],
    target_potential=-3.36,
    atom1_index=137,  # 约束原子1
    atom2_index=204,  # 约束原子2
    initial_distance=1.36,
    final_distance=1.50,
    steps=1000,
)

# 示例4: 命令行使用
# 训练:
#   python cp_mace_predictor.py train train.xyz --model-name my_model

# 模拟:
#   python cp_mace_predictor.py simulate init.xyz \\
#       --model-paths model1.model model2.model \\
#       --target-potential -3.36 \\
#       --temperature 300 \\
#       --steps 1000
"""
    print(examples)


def main():
    """运行所有测试"""
    print("=" * 60)
    print("CP-MACE Predictor 测试")
    print("=" * 60)
    print()

    # 测试1: 初始化
    test_initialization()

    # 测试2: 结构加载
    test_structure_loading()

    # 测试3: 模型加载
    test_model_loading()

    # 测试4: 短时间模拟
    # 注意: 这需要GPU和有效的模型文件
    print("\n注意: 短时间模拟测试需要GPU支持和有效的模型文件")
    # 自动跳过需要交互的测试
    print("跳过模拟测试 (需要手动运行)")

    # 测试5: 训练数据准备
    test_training_data_preparation()

    # 打印使用示例
    print_usage_examples()

    print("\n" + "=" * 60)
    print("测试完成")
    print("=" * 60)


if __name__ == "__main__":
    main()
