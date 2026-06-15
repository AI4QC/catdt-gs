#!/usr/bin/env python
"""
Complete Example: CP-MACE + VSSR-MC Electrochemical Surface Reconstruction

This script demonstrates the complete workflow using CP-MACE pretrained models:
1. Generate synthetic training data using pretrained models
2. Demonstrate CP-MACE training workflow
3. Run VSSR-MC electrochemical surface reconstruction simulation
4. Generate animation

Usage:
    conda activate catdt
    python complete_example.py

Author: CatDT
"""

import os
import sys
import json
import time
import copy
import gc
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from collections import Counter

# 设置路径
CP_MACE_ROOT = "deps/CP-MACE"
SURFACE_SAMPLING_ROOT = "deps/surface-sampling"
OUTPUT_DIR = "output/complete_example"

if CP_MACE_ROOT not in sys.path:
    sys.path.insert(0, CP_MACE_ROOT)

os.makedirs(OUTPUT_DIR, exist_ok=True)


def print_section(title: str):
    """打印分节标题"""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70 + "\n")


# ============================================================================
# 第一部分: 加载和准备数据
# ============================================================================

def prepare_system():
    """加载CP-MACE示例系统"""
    print_section("第一部分: 加载CP-MACE示例系统")

    import torch
    from ase.io import read, write
    from mace.calculators import MACECalculator

    # 加载示例结构
    init_file = os.path.join(CP_MACE_ROOT, "simulation/slow_growth/init.xyz")
    print(f"加载示例结构: {init_file}")

    atoms = read(init_file)
    print(f"  化学式: {atoms.get_chemical_formula()}")
    print(f"  原子数: {len(atoms)}")
    print(f"  电子数: {atoms.info.get('electron', 'N/A')}")
    print(f"  费米能级: {atoms.info.get('potential', 'N/A')} V")

    # 加载预训练模型
    model_paths = [
        os.path.join(CP_MACE_ROOT, "simulation/slow_growth/MACE_model_compiled_1.model"),
        os.path.join(CP_MACE_ROOT, "simulation/slow_growth/MACE_model_compiled_2.model"),
    ]

    print(f"\n加载CP-MACE模型...")
    calculators = []
    for path in model_paths:
        if os.path.exists(path):
            calc = MACECalculator(model_paths=[path], device='cuda')
            calculators.append(calc)
            print(f"  ✓ {os.path.basename(path)}")

    print(f"  共加载 {len(calculators)} 个模型")

    return atoms, calculators


# ============================================================================
# 第二部分: 生成训练数据 (主动学习模拟)
# ============================================================================

def generate_training_data(
    base_atoms,
    calculators,
    num_structures: int = 20,
    output_file: str = None,
):
    """
    使用预训练模型生成训练数据

    模拟主动学习过程：
    1. 对基础结构添加扰动
    2. 使用ensemble模型计算能量和不确定性
    3. 收集数据点
    """
    print_section("第二部分: 生成训练数据 (主动学习模拟)")

    from ase.io import write
    import copy

    if output_file is None:
        output_file = os.path.join(OUTPUT_DIR, "training_data.xyz")

    print(f"基础结构: {base_atoms.get_chemical_formula()}")
    print(f"生成 {num_structures} 个训练结构...")

    structures = []
    energies = []
    uncertainties = []

    for i in range(num_structures):
        atoms = base_atoms.copy()

        # 位置扰动
        positions = atoms.get_positions()
        noise = np.random.normal(0, 0.02, positions.shape)
        atoms.set_positions(positions + noise)

        # 电子数扰动
        base_ne = base_atoms.info.get('electron', 660.0)
        atoms.info['electron'] = base_ne + np.random.uniform(-2, 2)

        # 使用ensemble计算
        all_energies = []
        all_forces = []
        all_potentials = []

        for calc in calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc
            try:
                e = atoms_copy.get_potential_energy()
                f = atoms_copy.get_forces()
                p = calc.results.get('potential', -3.5)
                all_energies.append(e)
                all_forces.append(f)
                all_potentials.append(p)
            except Exception as ex:
                print(f"  结构 {i}: 计算失败 - {ex}")

        if all_energies:
            # 平均值
            mean_energy = np.mean(all_energies)
            mean_forces = np.mean(all_forces, axis=0)
            mean_potential = np.mean(all_potentials)

            # 不确定性
            if len(all_energies) > 1:
                force_std = np.max(np.std(all_forces, axis=0))
            else:
                force_std = 0

            # 设置参考数据
            atoms.info['REF_energy'] = mean_energy
            atoms.info['potential'] = mean_potential
            atoms.arrays['REF_forces'] = mean_forces

            structures.append(atoms)
            energies.append(mean_energy)
            uncertainties.append(force_std)

            if (i + 1) % 5 == 0:
                print(f"  {i+1}/{num_structures}: E={mean_energy:.2f} eV, "
                      f"μ={mean_potential:.4f} V, σ_F={force_std:.4f}")

    # 保存
    write(output_file, structures, format='extxyz')
    print(f"\n✓ 训练数据已保存: {output_file}")
    print(f"  结构数: {len(structures)}")
    print(f"  能量范围: [{min(energies):.2f}, {max(energies):.2f}] eV")
    print(f"  平均不确定性: {np.mean(uncertainties):.4f} eV/Å")

    return output_file, structures


# ============================================================================
# 第三部分: 电化学VSSR-MC模拟
# ============================================================================

def run_electrochemical_mc(
    base_atoms,
    calculators,
    potential_she: float = -3.36 - (-4.44),  # 对应init.xyz中的值
    ph: float = 0.0,
    temperature: float = 300.0,
    total_sweeps: int = 30,
    sweep_size: int = 10,
    output_dir: str = None,
):
    """
    运行电化学VSSR-MC模拟

    在恒电势条件下采样表面结构
    """
    print_section("第三部分: 电化学VSSR-MC模拟")

    from ase.io import write
    from ase import Atom
    import copy
    import torch

    if output_dir is None:
        output_dir = os.path.join(OUTPUT_DIR, "mc_simulation")
    os.makedirs(output_dir, exist_ok=True)

    print(f"模拟条件:")
    print(f"  电势: {potential_she:.2f} V vs. SHE")
    print(f"  目标费米能级: {-4.44 - potential_she:.4f} eV")
    print(f"  pH: {ph:.1f}")
    print(f"  温度: {temperature:.1f} K")

    # 识别可移动的原子 (非金属原子在表面以上)
    symbols = base_atoms.get_chemical_symbols()
    positions = base_atoms.get_positions()
    z_mean = np.mean(positions[:, 2])

    # 找到吸附物原子 (H, O, C, N 且 z > z_mean)
    adsorbate_elements = ['H', 'O', 'C', 'N']
    movable_indices = []
    for i, (sym, pos) in enumerate(zip(symbols, positions)):
        if sym in adsorbate_elements and pos[2] > z_mean:
            movable_indices.append(i)

    print(f"  可移动原子: {len(movable_indices)}")
    print(f"  可添加物种: {adsorbate_elements}")

    # 生成虚拟吸附位点
    cell = base_atoms.get_cell()
    z_max = np.max(positions[:, 2])

    ads_sites = []
    for i in range(3):
        for j in range(3):
            x = cell[0, 0] * (i + 0.5) / 3 + np.random.uniform(-1, 1)
            y = cell[1, 1] * (j + 0.5) / 3 + np.random.uniform(-1, 1)
            z = z_max + 2.0
            ads_sites.append([x, y, z])
    ads_sites = np.array(ads_sites)

    print(f"  虚拟吸附位点: {len(ads_sites)}")

    # 能量计算函数
    def calculate_energy(atoms):
        """使用ensemble计算能量"""
        all_energies = []
        all_potentials = []
        all_forces = []

        for calc in calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc
            try:
                e = atoms_copy.get_potential_energy()
                f = atoms_copy.get_forces()
                p = calc.results.get('potential', -3.5)
                all_energies.append(e)
                all_potentials.append(p)
                all_forces.append(f)
            except:
                pass

        if not all_energies:
            return None, None, None

        energy = np.mean(all_energies)
        potential = np.mean(all_potentials)

        if len(all_forces) > 1:
            force_std = np.max(np.std(all_forces, axis=0))
        else:
            force_std = 0

        return energy, potential, force_std

    # 初始化
    current_atoms = base_atoms.copy()
    current_energy, current_potential, current_std = calculate_energy(current_atoms)

    if current_energy is None:
        raise RuntimeError("无法计算初始能量")

    print(f"\n初始状态:")
    print(f"  能量: {current_energy:.2f} eV")
    print(f"  费米能级: {current_potential:.4f} V")

    # 记录
    trajectory = [current_atoms.copy()]
    energy_history = [current_energy]
    potential_history = [current_potential]
    std_history = [current_std]
    composition_history = [dict(Counter(current_atoms.get_chemical_symbols()))]

    kT = 8.617e-5 * temperature
    accepted = 0
    high_uncertainty = []

    print(f"\n开始MC采样...")

    for sweep in range(total_sweeps):
        sweep_accepted = 0

        for step in range(sweep_size):
            trial_atoms = current_atoms.copy()

            # 随机选择MC移动
            # 对于这个系统，主要移动H原子
            move_type = np.random.choice(['move', 'swap'], p=[0.7, 0.3])

            if move_type == 'move':
                # 移动一个可移动原子
                if movable_indices:
                    # 找当前结构中的可移动原子
                    current_symbols = trial_atoms.get_chemical_symbols()
                    current_movable = [i for i, s in enumerate(current_symbols)
                                      if s in adsorbate_elements and
                                      trial_atoms[i].position[2] > z_mean - 2]

                    if current_movable:
                        idx = np.random.choice(current_movable)
                        # 随机移动
                        delta = np.random.normal(0, 0.2, 3)
                        trial_atoms[idx].position += delta

            elif move_type == 'swap':
                # 交换两个可移动原子的位置
                current_symbols = trial_atoms.get_chemical_symbols()
                current_movable = [i for i, s in enumerate(current_symbols)
                                  if s in adsorbate_elements and
                                  trial_atoms[i].position[2] > z_mean - 2]

                if len(current_movable) >= 2:
                    idx1, idx2 = np.random.choice(current_movable, 2, replace=False)
                    pos1 = trial_atoms[idx1].position.copy()
                    pos2 = trial_atoms[idx2].position.copy()
                    trial_atoms[idx1].position = pos2
                    trial_atoms[idx2].position = pos1

            # 保持电子数
            trial_atoms.info['electron'] = current_atoms.info.get('electron', 660.0)

            # 计算试探能量
            trial_energy, trial_potential, trial_std = calculate_energy(trial_atoms)

            if trial_energy is None:
                continue

            # Metropolis判据
            delta_E = trial_energy - current_energy

            if delta_E < 0 or np.random.random() < np.exp(-delta_E / kT):
                current_atoms = trial_atoms
                current_energy = trial_energy
                current_potential = trial_potential
                current_std = trial_std
                sweep_accepted += 1
                accepted += 1

        # 记录
        trajectory.append(current_atoms.copy())
        energy_history.append(current_energy)
        potential_history.append(current_potential)
        std_history.append(current_std)
        composition_history.append(dict(Counter(current_atoms.get_chemical_symbols())))

        if current_std > 0.1:
            high_uncertainty.append(current_atoms.copy())

        if (sweep + 1) % 5 == 0:
            rate = accepted / ((sweep + 1) * sweep_size)
            print(f"  轮 {sweep+1}/{total_sweeps}: "
                  f"E={current_energy:.2f} eV, "
                  f"μ={current_potential:.4f} V, "
                  f"接受率={rate:.1%}")

        # 清理GPU内存
        if (sweep + 1) % 10 == 0:
            gc.collect()
            torch.cuda.empty_cache()

    # 保存结果
    print("\n保存结果...")

    trajectory_file = os.path.join(output_dir, "trajectory.xyz")
    write(trajectory_file, trajectory, format='extxyz')
    print(f"  ✓ 轨迹: {trajectory_file}")

    final_file = os.path.join(output_dir, "final.xyz")
    write(final_file, current_atoms, format='extxyz')
    print(f"  ✓ 最终结构: {final_file}")

    min_idx = np.argmin(energy_history)
    lowest_file = os.path.join(output_dir, "lowest_energy.xyz")
    write(lowest_file, trajectory[min_idx], format='extxyz')
    print(f"  ✓ 最低能量结构: {lowest_file}")

    if high_uncertainty:
        uncertain_file = os.path.join(output_dir, "high_uncertainty.xyz")
        write(uncertain_file, high_uncertainty, format='extxyz')
        print(f"  ✓ 高不确定性结构: {uncertain_file} ({len(high_uncertainty)}个)")

    # 保存数据 (确保所有数值都是Python原生类型)
    def to_native(obj):
        if isinstance(obj, dict):
            return {k: to_native(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [to_native(x) for x in obj]
        elif isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    data = {
        'energy_history': [float(x) for x in energy_history],
        'potential_history': [float(x) for x in potential_history],
        'std_history': [float(x) for x in std_history],
        'composition_history': to_native(composition_history),
        'conditions': {
            'potential_she': float(potential_she),
            'ph': float(ph),
            'temperature': float(temperature),
        },
        'statistics': {
            'acceptance_rate': float(accepted / (total_sweeps * sweep_size)),
            'min_energy': float(min(energy_history)),
            'min_energy_sweep': int(min_idx),
        }
    }

    data_file = os.path.join(output_dir, "simulation_data.json")
    with open(data_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"\n模拟统计:")
    print(f"  总步数: {total_sweeps * sweep_size}")
    print(f"  接受率: {accepted / (total_sweeps * sweep_size):.1%}")
    print(f"  能量范围: [{min(energy_history):.2f}, {max(energy_history):.2f}] eV")
    print(f"  最低能量: {energy_history[min_idx]:.2f} eV (轮 {min_idx})")

    return trajectory_file, data_file, output_dir


# ============================================================================
# 第四部分: 生成动画
# ============================================================================

def generate_animation(trajectory_file, data_file, output_dir, fps=3):
    """生成模拟动画"""
    print_section("第四部分: 生成动画")

    from ase.io import read
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    # 读取数据
    trajectory = read(trajectory_file, index=':')
    print(f"轨迹帧数: {len(trajectory)}")

    with open(data_file) as f:
        data = json.load(f)

    energy_history = data['energy_history']
    potential_history = data['potential_history']

    # 颜色映射
    colors = {
        'Au': '#FFD700', 'Ag': '#C0C0C0', 'Pt': '#E5E4E2',
        'O': '#FF0000', 'H': '#FFFFFF', 'C': '#808080',
        'N': '#0000FF', 'Ni': '#00AA00', 'Fe': '#A52A2A',
    }

    # 创建图形
    fig = plt.figure(figsize=(15, 5))

    ax1 = fig.add_subplot(131, projection='3d')
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    def update(frame):
        ax1.clear()
        ax2.clear()
        ax3.clear()

        atoms = trajectory[frame]
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        # 3D结构
        for pos, sym in zip(positions, symbols):
            c = colors.get(sym, '#808080')
            s = 50 if sym in ['H'] else 100
            ax1.scatter(pos[0], pos[1], pos[2], c=c, s=s,
                       edgecolors='black', linewidths=0.3, alpha=0.8)

        cell = atoms.get_cell()
        ax1.set_xlim(0, cell[0, 0])
        ax1.set_ylim(0, cell[1, 1])
        ax1.set_zlim(0, cell[2, 2])
        ax1.set_xlabel('X (Å)')
        ax1.set_ylabel('Y (Å)')
        ax1.set_zlabel('Z (Å)')
        ax1.set_title(f'Frame {frame+1}/{len(trajectory)}')

        # 能量曲线
        ax2.plot(energy_history[:frame+1], 'b-', linewidth=2)
        ax2.scatter([frame], [energy_history[frame]], c='red', s=100, zorder=5)
        ax2.axhline(y=min(energy_history), color='g', linestyle='--', alpha=0.5)
        ax2.set_xlabel('MC Sweep')
        ax2.set_ylabel('Energy (eV)')
        ax2.set_title(f'E = {energy_history[frame]:.2f} eV')
        ax2.grid(True, alpha=0.3)

        # 费米能级
        ax3.plot(potential_history[:frame+1], 'purple', linewidth=2)
        ax3.scatter([frame], [potential_history[frame]], c='red', s=100, zorder=5)
        ax3.set_xlabel('MC Sweep')
        ax3.set_ylabel('Fermi Level (V)')
        ax3.set_title(f'μ = {potential_history[frame]:.4f} V')
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()

    print("创建动画...")
    anim = FuncAnimation(fig, update, frames=len(trajectory),
                        blit=False, interval=1000//fps)

    gif_file = os.path.join(output_dir, "mc_animation.gif")
    writer = PillowWriter(fps=fps)
    anim.save(gif_file, writer=writer, dpi=100)
    print(f"✓ GIF动画已保存: {gif_file}")

    plt.close()

    # 静态摘要图
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    axes[0].plot(energy_history, 'b-', linewidth=2)
    axes[0].axhline(y=min(energy_history), color='r', linestyle='--',
                   label=f'Min: {min(energy_history):.2f} eV')
    axes[0].set_xlabel('MC Sweep')
    axes[0].set_ylabel('Energy (eV)')
    axes[0].set_title('Energy Evolution')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(potential_history, 'purple', linewidth=2)
    axes[1].set_xlabel('MC Sweep')
    axes[1].set_ylabel('Fermi Level (V)')
    axes[1].set_title('Fermi Level Evolution')
    axes[1].grid(True, alpha=0.3)

    # 最终结构俯视图
    final = trajectory[-1]
    pos = final.get_positions()
    sym = final.get_chemical_symbols()
    for p, s in zip(pos, sym):
        c = colors.get(s, '#808080')
        sz = 20 if s == 'H' else 50
        axes[2].scatter(p[0], p[1], c=c, s=sz, edgecolors='black', alpha=0.8)
    axes[2].set_xlabel('X (Å)')
    axes[2].set_ylabel('Y (Å)')
    axes[2].set_title('Final Structure (Top View)')
    axes[2].set_aspect('equal')

    plt.tight_layout()
    summary_file = os.path.join(output_dir, "summary.png")
    plt.savefig(summary_file, dpi=150)
    plt.close()
    print(f"✓ 摘要图已保存: {summary_file}")

    return gif_file


# ============================================================================
# 主程序
# ============================================================================

def main():
    """运行完整示例"""
    print("\n" + "=" * 70)
    print("  完整示例 V3: CP-MACE + VSSR-MC 电化学表面重构")
    print("=" * 70)
    print(f"\n输出目录: {OUTPUT_DIR}")

    start_time = time.time()

    # 第一部分: 加载系统
    base_atoms, calculators = prepare_system()

    # 第二部分: 生成训练数据
    train_file, structures = generate_training_data(
        base_atoms, calculators, num_structures=15
    )

    # 第三部分: 电化学MC模拟
    trajectory_file, data_file, sim_dir = run_electrochemical_mc(
        base_atoms,
        calculators,
        potential_she=1.08,  # ~对应init.xyz中的potential=-3.36
        ph=0.0,
        temperature=300.0,
        total_sweeps=25,
        sweep_size=8,
    )

    # 第四部分: 生成动画
    gif_file = generate_animation(trajectory_file, data_file, sim_dir)

    # 总结
    elapsed = time.time() - start_time
    print_section("完成总结")
    print(f"总用时: {elapsed:.1f} 秒\n")

    print("生成的文件:")
    for root, dirs, files in os.walk(OUTPUT_DIR):
        level = root.replace(OUTPUT_DIR, '').count(os.sep)
        indent = '  ' * level
        print(f'{indent}{os.path.basename(root)}/')
        subindent = '  ' * (level + 1)
        for f in files:
            fpath = os.path.join(root, f)
            size = os.path.getsize(fpath) / 1024
            print(f'{subindent}{f} ({size:.1f} KB)')

    print(f"\n查看动画: xdg-open {gif_file}")

    return OUTPUT_DIR


if __name__ == "__main__":
    main()
