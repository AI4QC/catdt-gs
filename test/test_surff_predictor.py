"""
测试SurFFPredictor封装类

使用Cu和CuAl示例结构测试SurFF封装类的功能
"""

import os
import sys

# 添加路径
sys.path.insert(0, "core/surface")

from surff_predictor import SurFFPredictor, predict_surface_exposure


def test_with_poscar():
    """使用POSCAR文件测试"""
    print("="*60)
    print("Test 1: Predict from POSCAR file")
    print("="*60)

    # 使用SurFF自带的示例结构
    structure_path = "deps/SurFF/SurFF_DataFiles/Data/example/crystal_structures/104"  # TiAl

    predictor = SurFFPredictor(
        surff_root="deps/SurFF",
        use_gpu=True,
        relaxation_steps=50,  # 减少步数加快测试
        verbose=True,
        keep_files=True,
    )

    result = predictor.predict(
        structure=structure_path,
        crystal_id="TiAl_test",
        top_n=5,
        save_wulff_shape=True,
        output_dir="output/surff/TiAl",
    )

    print("\n" + "="*60)
    print("Test 1 Results:")
    print("="*60)
    print(result.summary(5))

    # 打印DataFrame
    print("\nFull results DataFrame:")
    print(result.to_dataframe())

    return result


def test_with_pymatgen():
    """使用pymatgen Structure测试"""
    print("\n" + "="*60)
    print("Test 2: Predict from pymatgen Structure (Cu bulk)")
    print("="*60)

    from pymatgen.core import Structure, Lattice

    # 创建一个简单的Cu FCC结构
    a = 3.615  # Cu晶格常数 (Å)
    lattice = Lattice.cubic(a)
    cu_structure = Structure(
        lattice,
        ["Cu", "Cu", "Cu", "Cu"],
        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    )

    predictor = SurFFPredictor(
        surff_root="deps/SurFF",
        use_gpu=True,
        relaxation_steps=50,
        verbose=True,
        keep_files=True,
    )

    result = predictor.predict(
        structure=cu_structure,
        crystal_id="Cu_fcc",
        top_n=5,
        save_wulff_shape=True,
        output_dir="output/surff/Cu",
    )

    print("\n" + "="*60)
    print("Test 2 Results:")
    print("="*60)
    print(result.summary(5))

    return result


def test_with_ase():
    """使用ASE Atoms测试"""
    print("\n" + "="*60)
    print("Test 3: Predict from ASE Atoms (Al bulk)")
    print("="*60)

    from ase.build import bulk

    # 创建Al FCC结构
    al_atoms = bulk('Al', 'fcc', a=4.05)

    predictor = SurFFPredictor(
        surff_root="deps/SurFF",
        use_gpu=True,
        relaxation_steps=50,
        verbose=True,
        keep_files=True,
    )

    result = predictor.predict(
        structure=al_atoms,
        crystal_id="Al_fcc",
        top_n=5,
        save_wulff_shape=True,
        output_dir="output/surff/Al",
    )

    print("\n" + "="*60)
    print("Test 3 Results:")
    print("="*60)
    print(result.summary(5))

    return result


def test_convenience_function():
    """测试便捷函数"""
    print("\n" + "="*60)
    print("Test 4: Using convenience function")
    print("="*60)

    structure_path = "deps/SurFF/SurFF_DataFiles/Data/example/crystal_structures/101"  # VRu

    result = predict_surface_exposure(
        structure=structure_path,
        surff_root="deps/SurFF",
        top_n=3,
        use_gpu=True,
        output_dir="output/surff/VRu",
        relaxation_steps=50,
        verbose=True,
    )

    print("\n" + "="*60)
    print("Test 4 Results:")
    print("="*60)
    print(result.summary(3))

    return result


def main():
    """运行所有测试"""
    print("="*60)
    print("SurFFPredictor Test Suite")
    print("="*60)

    # 创建输出目录
    os.makedirs("output/surff", exist_ok=True)

    # 测试1: POSCAR文件
    result1 = test_with_poscar()

    # 测试2: pymatgen Structure
    # result2 = test_with_pymatgen()

    # 测试3: ASE Atoms
    # result3 = test_with_ase()

    # 测试4: 便捷函数
    # result4 = test_convenience_function()

    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)
    print(f"\nOutput directory: output/surff/")

    # 显示输出文件
    for root, dirs, files in os.walk("output/surff"):
        level = root.replace("output/surff", "").count(os.sep)
        indent = " " * 2 * level
        print(f"{indent}{os.path.basename(root)}/")
        subindent = " " * 2 * (level + 1)
        for file in files[:5]:  # 只显示前5个文件
            print(f"{subindent}{file}")
        if len(files) > 5:
            print(f"{subindent}... and {len(files) - 5} more files")


if __name__ == "__main__":
    main()
