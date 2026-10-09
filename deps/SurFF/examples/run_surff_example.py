"""
SurFF 示例脚本
表面能预测和Wulff构型计算

运行方式：
    conda activate catdt
    cd $CATDT_ROOT/deps/SurFF
    python examples/run_surff_example.py

功能:
1. 从晶体结构生成表面slab
2. 使用MLFF模型进行结构弛豫
3. 计算表面能和Wulff形状
"""

import os
import sys
import pandas as pd
import subprocess

# 添加项目路径
SURFF_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SURFF_ROOT)

from ase.io import read
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt


def setup_directories():
    """设置必要的目录结构"""
    dirs = [
        os.path.join(SURFF_ROOT, "Data/example/slabs"),
        os.path.join(SURFF_ROOT, "Data/example/bulks"),
        os.path.join(SURFF_ROOT, "ocp/data/example"),
        os.path.join(SURFF_ROOT, "ocp/traj/example"),
        os.path.join(SURFF_ROOT, "results/example"),
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

    # 创建checkpoint符号链接
    ckpt_src = os.path.join(SURFF_ROOT, "SurFF_DataFiles/ocp/checkpoints/2024-03-03-23-57-52")
    ckpt_dst = os.path.join(SURFF_ROOT, "ocp/checkpoints/2024-03-03-23-57-52")
    if os.path.exists(ckpt_src) and not os.path.exists(ckpt_dst):
        os.makedirs(os.path.dirname(ckpt_dst), exist_ok=True)
        os.symlink(ckpt_src, ckpt_dst)


def step1_generate_slabs(crystal_dir, slab_dir, bulk_dir):
    """
    步骤1: 从晶体结构生成表面slab
    """
    from functionals.gen_bulk_slab import process_structure_sep

    print("="*60)
    print("步骤1: 从晶体结构生成表面slab")
    print("="*60)

    files = [f for f in os.listdir(crystal_dir) if not f.startswith('.')]
    print(f"找到 {len(files)} 个晶体结构: {files}")

    info = []
    for poscar in files:
        try:
            slab_info = process_structure_sep(
                crystal_id=int(poscar),
                slab_dir=slab_dir,
                bulk_dir=bulk_dir,
                source_dir=crystal_dir
            )
            info.extend(slab_info)
            print(f"  处理 {poscar}: 生成 {len(slab_info)} 个slab")
        except Exception as e:
            print(f"  处理 {poscar} 时出错: {e}")

    info_df = pd.DataFrame(info, columns=['slab_id', 'formula', 'miller_index', 'shift', 'num_atom'])
    print(f"\n总共生成 {len(info_df)} 个slab结构")
    print(info_df.head(10))

    # 保存info
    info_df.to_csv(os.path.join(slab_dir, "info.csv"), index=False)

    return info_df


def step2_generate_lmdb(slab_dir, lmdb_path):
    """
    步骤2: 将slab结构编译成LMDB数据集
    """
    from functionals.gen_lmdb_relaxation import gen_lmdb

    print("\n" + "="*60)
    print("步骤2: 生成LMDB数据集")
    print("="*60)

    slab_files = [f for f in os.listdir(slab_dir) if not f.endswith('.csv') and not f.startswith('.')]
    dataset = [{'slab_id': f, "POSCAR_pth": os.path.join(slab_dir, f)} for f in slab_files]

    gen_lmdb(dataset=dataset, DB_path=lmdb_path)
    print(f"LMDB数据集保存到: {lmdb_path}")


def step3_run_relaxation(use_gpu=False):
    """
    步骤3: 使用MLFF模型进行结构弛豫
    """
    print("\n" + "="*60)
    print("步骤3: 运行MLFF结构弛豫")
    print("="*60)

    ocp_dir = os.path.join(SURFF_ROOT, "ocp")

    cmd = [
        "python", "main.py",
        "--mode", "run-relaxations",
        "--config-yml", "configs/equiformer_v2_002_relax.yml",
        "--checkpoint", "checkpoints/2024-03-03-23-57-52/best_checkpoint.pt",
        "--task.relax_opt.traj_dir=traj/example",
        "--task.relax_dataset.src=data/example",
        "--task.relaxation_steps=50",
        "--task.relaxation_fmax=0.05",
    ]

    if not use_gpu:
        cmd.append("--cpu")

    print(f"运行命令: {' '.join(cmd)}")
    print(f"工作目录: {ocp_dir}")

    result = subprocess.run(cmd, cwd=ocp_dir, capture_output=False)

    if result.returncode == 0:
        print("弛豫完成!")
    else:
        print(f"弛豫失败，返回码: {result.returncode}")

    return result.returncode == 0


def step4_calculate_wulff(traj_dir, info_df, crystal_dir, save_dir):
    """
    步骤4: 计算表面能和Wulff形状
    """
    print("\n" + "="*60)
    print("步骤4: 计算表面能和Wulff形状")
    print("="*60)

    # 复制run_wulff.py到results目录
    wulff_src = os.path.join(SURFF_ROOT, "SurFF_DataFiles/results/run_wulff.py")
    wulff_dst = os.path.join(SURFF_ROOT, "results/run_wulff.py")
    if os.path.exists(wulff_src) and not os.path.exists(wulff_dst):
        import shutil
        shutil.copy(wulff_src, wulff_dst)

    from results.run_wulff import run_wulff

    run_wulff(traj_dir, info=info_df, crystal_dir=crystal_dir, save_dir=save_dir)

    # 读取结果
    results_csv = os.path.join(save_dir, "wulff_results.csv")
    if os.path.exists(results_csv):
        results = pd.read_csv(results_csv)
        print("\n表面能预测结果:")
        print(results)
        return results
    return None


def visualize_results(save_dir):
    """
    可视化Wulff形状结果
    """
    print("\n" + "="*60)
    print("可视化结果")
    print("="*60)

    wulff_shape_dir = os.path.join(save_dir, "wulff_shape")
    if os.path.exists(wulff_shape_dir):
        images = [f for f in os.listdir(wulff_shape_dir) if f.endswith('.png')]
        print(f"生成的Wulff形状图: {images}")
    else:
        print("未找到Wulff形状图目录")


def main():
    """主函数"""
    print("="*60)
    print("SurFF 示例运行")
    print("金属间化合物表面能预测")
    print("="*60)

    # 设置目录
    setup_directories()

    # 路径配置
    crystal_dir = os.path.join(SURFF_ROOT, "SurFF_DataFiles/Data/example/crystal_structures")
    slab_dir = os.path.join(SURFF_ROOT, "Data/example/slabs")
    bulk_dir = os.path.join(SURFF_ROOT, "Data/example/bulks")
    lmdb_path = os.path.join(SURFF_ROOT, "ocp/data/example/surface_relaxation.lmdb")
    traj_dir = os.path.join(SURFF_ROOT, "ocp/traj/example")
    save_dir = os.path.join(SURFF_ROOT, "results/example")

    # 检查是否已有预处理的数据
    use_existing_lmdb = os.path.exists(os.path.join(SURFF_ROOT, "SurFF_DataFiles/ocp/data/example/surface_relaxation.lmdb"))
    use_existing_traj = len([f for f in os.listdir(traj_dir) if f.endswith('.traj')]) > 0 if os.path.exists(traj_dir) else False

    print(f"\n配置:")
    print(f"  晶体结构目录: {crystal_dir}")
    print(f"  使用已有LMDB: {use_existing_lmdb}")
    print(f"  使用已有轨迹: {use_existing_traj}")

    # 步骤1: 生成slab (如果需要)
    info_csv = os.path.join(slab_dir, "info.csv")
    if os.path.exists(info_csv):
        info_df = pd.read_csv(info_csv)
        print(f"\n使用已有的slab信息: {len(info_df)} 个slab")
    else:
        info_df = step1_generate_slabs(crystal_dir, slab_dir, bulk_dir)

    # 步骤2: 生成LMDB (如果需要)
    if use_existing_lmdb:
        # 使用预处理的LMDB
        src_lmdb = os.path.join(SURFF_ROOT, "SurFF_DataFiles/ocp/data/example/surface_relaxation.lmdb")
        if not os.path.exists(lmdb_path):
            os.symlink(src_lmdb, lmdb_path)
            os.symlink(src_lmdb + "-lock", lmdb_path + "-lock")
        print(f"\n使用已有LMDB: {src_lmdb}")
    else:
        step2_generate_lmdb(slab_dir, lmdb_path)

    # 步骤3: 运行弛豫 (如果需要)
    if not use_existing_traj:
        success = step3_run_relaxation(use_gpu=True)
        if not success:
            print("弛豫失败，请检查错误信息")
            return
    else:
        print(f"\n使用已有轨迹文件: {traj_dir}")

    # 步骤4: 计算Wulff形状
    results = step4_calculate_wulff(traj_dir, info_df, crystal_dir, save_dir)

    # 可视化
    visualize_results(save_dir)

    print("\n" + "="*60)
    print("SurFF 示例运行完成!")
    print(f"结果保存在: {save_dir}")
    print("="*60)


if __name__ == "__main__":
    main()
