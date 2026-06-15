"""
测试AdsorbDiffPredictor封装类

使用Cu和Al4Cu2示例结构测试AdsorbDiff封装类的功能
"""

import os
import sys

# 添加路径
sys.path.insert(0, "core/reconstruction")

from adsorbdiff_predictor import AdsorbDiffPredictor, predict_adsorption_site


def test_with_slab_and_adsorbate():
    """测试方式1: 提供表面结构和吸附物（分开提供）"""
    print("="*60)
    print("Test 1: Surface + Adsorbate (separate)")
    print("="*60)

    # 使用pymatgen创建一个简单的Cu(111)表面
    from pymatgen.core import Structure, Lattice
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.ase import AseAtomsAdaptor

    # 创建Cu FCC结构
    a = 3.615
    lattice = Lattice.cubic(a)
    cu_bulk = Structure(
        lattice,
        ["Cu", "Cu", "Cu", "Cu"],
        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    )

    # 生成111表面
    slabgen = SlabGenerator(
        cu_bulk,
        miller_index=(1, 1, 1),
        min_slab_size=8.0,
        min_vacuum_size=15.0,
        center_slab=True,
    )
    slabs = slabgen.get_slabs()
    slab = slabs[0]  # 取第一个slab

    # 转换为ASE Atoms
    slab_atoms = AseAtomsAdaptor.get_atoms(slab)

    # 设置tags (表面原子tag=1)
    positions = slab_atoms.get_positions()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    tags = [1 if z > z_max - 2.5 else 0 for z in z_coords]
    slab_atoms.set_tags(tags)

    print(f"Surface: Cu(111), {len(slab_atoms)} atoms")

    # 初始化预测器
    predictor = AdsorbDiffPredictor(
        adsorbdiff_root="deps/AdsorbDiff",
        relax_checkpoint_path="deps/AdsorbDiff/eq2_153M_ec4_allmd.pt",  # 使用fairchem弛豫模型
        use_gpu=True,
        num_sites=3,  # 采样3个位点加快测试
        diffusion_steps=50,  # 减少步数加快测试
        relax_steps=30,
        verbose=True,
        keep_files=True,
    )

    # 预测CO吸附
    result = predictor.predict(
        surface=slab_atoms,
        adsorbate="*CO",  # 使用SMILES格式
        has_adsorbate=False,
        output_dir="output/adsorbdiff/Cu111_CO",
    )

    print("\n" + "="*60)
    print("Test 1 Results:")
    print("="*60)
    print(result.summary(3))

    return result


def test_with_existing_adsorbate():
    """测试方式2: 表面已经包含吸附物"""
    print("\n" + "="*60)
    print("Test 2: Surface with existing adsorbate")
    print("="*60)

    from ase.build import fcc111, add_adsorbate
    from ase import Atoms

    # 创建Cu(111)表面
    slab = fcc111('Cu', size=(3, 3, 4), vacuum=15.0)

    # 设置tags
    positions = slab.get_positions()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    tags = [1 if z > z_max - 2.5 else 0 for z in z_coords]
    slab.set_tags(tags)

    # 添加CO吸附物
    co = Atoms('CO', positions=[[0, 0, 0], [0, 0, 1.128]])
    add_adsorbate(slab, co, height=1.8, position='ontop')

    # 将吸附物tag设为2
    new_tags = list(slab.get_tags())
    new_tags[-2] = 2  # C
    new_tags[-1] = 2  # O
    slab.set_tags(new_tags)

    print(f"Surface with adsorbate: Cu(111)+CO, {len(slab)} atoms")
    print(f"Adsorbate atoms (tag=2): {sum(1 for t in slab.get_tags() if t == 2)}")

    # 初始化预测器
    predictor = AdsorbDiffPredictor(
        adsorbdiff_root="deps/AdsorbDiff",
        use_gpu=True,
        diffusion_steps=50,
        relax_steps=30,
        verbose=True,
        keep_files=True,
    )

    # 预测（吸附物已在结构中）
    result = predictor.predict(
        surface=slab,
        has_adsorbate=True,
        adsorbate_tag=2,
        output_dir="output/adsorbdiff/Cu111_CO_existing",
    )

    print("\n" + "="*60)
    print("Test 2 Results:")
    print("="*60)
    print(result.summary())

    return result


def test_with_file_input():
    """测试方式3: 从文件读取表面结构"""
    print("\n" + "="*60)
    print("Test 3: Surface from file")
    print("="*60)

    # 使用AdsorbDiff数据库中的结构
    from ase.build import fcc111
    from ase.io import write

    # 创建一个Al(100)表面
    slab = fcc111('Al', size=(3, 3, 4), vacuum=15.0, a=4.05)

    # 设置tags
    positions = slab.get_positions()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    tags = [1 if z > z_max - 2.5 else 0 for z in z_coords]
    slab.set_tags(tags)

    # 保存为VASP格式
    slab_path = "output/adsorbdiff/Al111_slab.vasp"
    os.makedirs(os.path.dirname(slab_path), exist_ok=True)
    write(slab_path, slab, format='vasp')

    print(f"Surface file: {slab_path}")

    # 初始化预测器
    predictor = AdsorbDiffPredictor(
        adsorbdiff_root="deps/AdsorbDiff",
        use_gpu=True,
        num_sites=3,
        diffusion_steps=50,
        relax_steps=30,
        verbose=True,
        keep_files=True,
    )

    # 预测O吸附
    result = predictor.predict(
        surface=slab_path,
        adsorbate="*O",
        has_adsorbate=False,
        output_dir="output/adsorbdiff/Al111_O",
    )

    print("\n" + "="*60)
    print("Test 3 Results:")
    print("="*60)
    print(result.summary(3))

    return result


def test_list_adsorbates():
    """测试列出可用吸附物"""
    print("\n" + "="*60)
    print("Test 4: List available adsorbates")
    print("="*60)

    predictor = AdsorbDiffPredictor(
        adsorbdiff_root="deps/AdsorbDiff",
        use_gpu=False,  # 不需要GPU
        verbose=False,
    )

    adsorbates = predictor.list_available_adsorbates()

    print(f"Found {len(adsorbates)} adsorbates in database:")
    for ads in adsorbates[:20]:
        print(f"  ID={ads['id']}: {ads['smiles']} ({ads['formula']}, {ads['num_atoms']} atoms)")

    if len(adsorbates) > 20:
        print(f"  ... and {len(adsorbates) - 20} more")

    return adsorbates


def test_convenience_function():
    """测试便捷函数"""
    print("\n" + "="*60)
    print("Test 5: Using convenience function")
    print("="*60)

    from ase.build import fcc100

    # 创建Cu(100)表面
    slab = fcc100('Cu', size=(3, 3, 4), vacuum=15.0)

    # 设置tags
    positions = slab.get_positions()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    tags = [1 if z > z_max - 2.5 else 0 for z in z_coords]
    slab.set_tags(tags)

    result = predict_adsorption_site(
        surface=slab,
        adsorbate="*H",
        adsorbdiff_root="deps/AdsorbDiff",
        use_gpu=True,
        num_samples=3,
        output_dir="output/adsorbdiff/Cu100_H",
        diffusion_steps=50,
        relax_steps=30,
        verbose=True,
    )

    print("\n" + "="*60)
    print("Test 5 Results:")
    print("="*60)
    print(result.summary(3))

    return result


def main():
    """运行所有测试"""
    print("="*60)
    print("AdsorbDiffPredictor Test Suite")
    print("="*60)

    # 创建输出目录
    os.makedirs("output/adsorbdiff", exist_ok=True)

    # 测试4: 列出可用吸附物（快速测试，不需要GPU）
    test_list_adsorbates()

    # 测试1: 表面+吸附物分开提供
    result1 = test_with_slab_and_adsorbate()

    # 测试2: 已有吸附物的表面（跳过，因为使用已有配置）
    # result2 = test_with_existing_adsorbate()

    # 测试3: 从文件读取
    # result3 = test_with_file_input()

    # 测试5: 便捷函数
    # result5 = test_convenience_function()

    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)
    print(f"\nOutput directory: output/adsorbdiff/")

    # 显示输出文件
    for root, dirs, files in os.walk("output/adsorbdiff"):
        level = root.replace("output/adsorbdiff", "").count(os.sep)
        indent = " " * 2 * level
        print(f"{indent}{os.path.basename(root)}/")
        subindent = " " * 2 * (level + 1)
        for file in files[:5]:
            print(f"{subindent}{file}")
        if len(files) > 5:
            print(f"{subindent}... and {len(files) - 5} more files")


if __name__ == "__main__":
    main()
