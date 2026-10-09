"""
AdsorbDiff 示例脚本
使用Cu和CuAl结构运行扩散采样和弛豫优化

运行方式：
    conda activate catdt
    cd $CATDT_ROOT/deps/AdsorbDiff
    python examples/run_example.py
"""

import os
import sys

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from ase.optimize import BFGS
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

from adsorbdiff import AdsorbDiffCalculator
from adsorbdiff.placement import (
    Adsorbate,
    AdsorbateSlabConfig,
    Bulk,
    Slab,
)


def run_adsorbdiff_example(
    bulk_id: int,
    adsorbate_id: int,
    diff_ckpt_path: str,
    relax_ckpt_path: str = None,
    miller_indices: tuple = (1, 1, 1),
    output_dir: str = "output",
    use_gpu: bool = True,
):
    """
    运行AdsorbDiff模型的示例函数

    Args:
        bulk_id: bulk结构在数据库中的索引 (12=Cu, 292=Al4Cu2)
        adsorbate_id: 吸附物在数据库中的索引 (0=O, 1=H, 5=CO)
        diff_ckpt_path: 扩散模型checkpoint路径
        relax_ckpt_path: 弛豫模型checkpoint路径 (可选)
        miller_indices: Miller指数，默认(1,1,1)
        output_dir: 输出目录
        use_gpu: 是否使用GPU
    """
    os.makedirs(output_dir, exist_ok=True)

    # 设置数据库路径
    main_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ads_db_path = os.path.join(main_path, "adsorbdiff/placement/pkls/adsorbates.pkl")
    bulk_db_path = os.path.join(main_path, "adsorbdiff/placement/pkls/bulks.pkl")

    print(f"="*60)
    print(f"AdsorbDiff 示例运行")
    print(f"="*60)

    # 1. 创建Bulk对象
    print(f"\n[1] 加载Bulk结构 (id={bulk_id})...")
    bulk = Bulk(bulk_id_from_db=bulk_id, bulk_db_path=bulk_db_path)
    print(f"    Bulk: {bulk}")
    print(f"    化学式: {bulk.atoms.get_chemical_formula()}")
    print(f"    原子数: {len(bulk.atoms)}")

    # 2. 创建Slab对象
    print(f"\n[2] 创建Slab结构 (Miller指数={miller_indices})...")
    slabs = Slab.from_bulk_get_specific_millers(
        specific_millers=miller_indices,
        bulk=bulk,
        min_ab=8.0,
    )

    if len(slabs) == 0:
        print(f"    警告: 未找到Miller指数为{miller_indices}的slab，尝试随机选择...")
        slab = Slab.from_bulk_get_random_slab(bulk=bulk, max_miller=2, min_ab=8.0)
    else:
        slab = slabs[0]

    print(f"    Slab: {slab}")
    print(f"    原子数: {len(slab.atoms)}")

    # 3. 创建Adsorbate对象
    print(f"\n[3] 加载吸附物 (id={adsorbate_id})...")
    adsorbate = Adsorbate(adsorbate_id_from_db=adsorbate_id, adsorbate_db_path=ads_db_path)
    print(f"    Adsorbate: {adsorbate}")

    # 4. 创建吸附物-表面配置
    print(f"\n[4] 创建吸附物-表面配置...")
    adslab_config = AdsorbateSlabConfig(
        slab=slab,
        adsorbate=adsorbate,
        num_sites=1,
        num_augmentations_per_site=1,
        interstitial_gap=2.0,
        mode="heuristic",
    )

    adslab = adslab_config.atoms_list[0]
    print(f"    总原子数: {len(adslab)}")
    print(f"    表面原子 (tag=1): {sum(adslab.get_tags() == 1)}")
    print(f"    吸附物原子 (tag=2): {sum(adslab.get_tags() == 2)}")

    # 保存初始结构图
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_atoms(adslab, ax=ax, rotation="-75x, 45y, 10z")
    ax.set_title(f"Initial: {bulk.atoms.get_chemical_formula()} + {adsorbate.atoms.get_chemical_formula()}")
    fig.savefig(os.path.join(output_dir, "initial_structure.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    已保存初始结构图: {output_dir}/initial_structure.png")

    # 5. 加载扩散模型并运行扩散
    print(f"\n[5] 加载扩散模型: {diff_ckpt_path}")
    calc_diff = AdsorbDiffCalculator(
        checkpoint_path=diff_ckpt_path,
        cpu=not use_gpu,
    )

    print(f"\n[6] 运行扩散采样...")
    diff_traj_dir = os.path.join(output_dir, "diff_traj")
    os.makedirs(diff_traj_dir, exist_ok=True)

    diffused_adslab = calc_diff.run_diffusion(adslab, trajectory=diff_traj_dir)
    print(f"    扩散完成！轨迹保存在: {diff_traj_dir}")

    # 获取扩散后的吸附位点
    diffused_adsorption_site = diffused_adslab.get_positions()[diffused_adslab.get_tags() == 2]

    # 重新放置吸附物确保正确的几何结构
    diffused_adslab = AdsorbateSlabConfig(
        slab=slab,
        adsorbate=adsorbate,
        sites=diffused_adsorption_site,
        interstitial_gap=0.1
    ).atoms_list[0]

    # 保存扩散后结构图
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_atoms(diffused_adslab, ax=ax, rotation="-75x, 45y, 10z")
    ax.set_title(f"After Diffusion")
    fig.savefig(os.path.join(output_dir, "diffused_structure.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    已保存扩散后结构图: {output_dir}/diffused_structure.png")

    # 6. 可选：使用MLFF进行弛豫优化
    if relax_ckpt_path and os.path.exists(relax_ckpt_path):
        print(f"\n[7] 加载弛豫模型: {relax_ckpt_path}")
        calc_opt = AdsorbDiffCalculator(
            checkpoint_path=relax_ckpt_path,
            cpu=not use_gpu,
        )

        print(f"\n[8] 运行结构优化...")
        diffused_adslab.calc = calc_opt
        opt_traj_path = os.path.join(output_dir, "opt.traj")
        opt_log_path = os.path.join(output_dir, "opt.log")

        opt = BFGS(diffused_adslab, trajectory=opt_traj_path, logfile=opt_log_path)
        converged = opt.run(fmax=0.05, steps=100)

        print(f"    优化完成！收敛: {converged}")
        print(f"    轨迹文件: {opt_traj_path}")

        # 保存优化后结构图
        fig, ax = plt.subplots(figsize=(8, 8))
        plot_atoms(diffused_adslab, ax=ax, rotation="-75x, 45y, 10z")
        ax.set_title(f"After Optimization")
        fig.savefig(os.path.join(output_dir, "optimized_structure.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"    已保存优化后结构图: {output_dir}/optimized_structure.png")

    print(f"\n{'='*60}")
    print(f"运行完成！所有输出保存在: {output_dir}")
    print(f"{'='*60}")

    return diffused_adslab


if __name__ == "__main__":
    # Checkpoint路径
    main_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # 使用PaiNN模型 (推荐，推理速度快)
    diff_ckpt_painn = os.path.join(main_path, "PT_zeroshot_painn.pt")

    # 使用EquiformerV2模型 (精度更高)
    diff_ckpt_eqv2 = os.path.join(main_path, "PT_fewshot_eqv2_cond.pt")

    # 选择使用的checkpoint
    if os.path.exists(diff_ckpt_painn):
        diff_ckpt = diff_ckpt_painn
        model_name = "PaiNN"
    elif os.path.exists(diff_ckpt_eqv2):
        diff_ckpt = diff_ckpt_eqv2
        model_name = "EquiformerV2"
    else:
        raise FileNotFoundError("未找到checkpoint文件！请确保PT_zeroshot_painn.pt或PT_fewshot_eqv2_cond.pt存在")

    print(f"使用模型: {model_name}")
    print(f"Checkpoint: {diff_ckpt}")

    # ============================================================
    # 示例1: Cu(111) + CO吸附
    # ============================================================
    print("\n" + "="*60)
    print("示例1: Cu(111) + CO吸附")
    print("="*60)

    run_adsorbdiff_example(
        bulk_id=12,           # Cu
        adsorbate_id=5,       # CO
        diff_ckpt_path=diff_ckpt,
        miller_indices=(1, 1, 1),
        output_dir="output/Cu111_CO",
        use_gpu=True,
    )

    # ============================================================
    # 示例2: Al4Cu2(111) + O吸附
    # ============================================================
    print("\n" + "="*60)
    print("示例2: Al4Cu2(111) + O吸附")
    print("="*60)

    run_adsorbdiff_example(
        bulk_id=292,          # Al4Cu2
        adsorbate_id=0,       # O
        diff_ckpt_path=diff_ckpt,
        miller_indices=(1, 1, 1),
        output_dir="output/Al4Cu2_111_O",
        use_gpu=True,
    )

    # ============================================================
    # 示例3: Cu(100) + H吸附
    # ============================================================
    print("\n" + "="*60)
    print("示例3: Cu(100) + H吸附")
    print("="*60)

    run_adsorbdiff_example(
        bulk_id=12,           # Cu
        adsorbate_id=1,       # H
        diff_ckpt_path=diff_ckpt,
        miller_indices=(1, 0, 0),
        output_dir="output/Cu100_H",
        use_gpu=True,
    )
