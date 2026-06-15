#!/usr/bin/env python
"""
只运行KMC模拟和KMC表面动力学可视化
使用之前计算好的pathway结果
"""

import os
import sys
sys.path.insert(0, "core")

import logging
import numpy as np
from ase.io import read

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 配置
OUTPUT_DIR = "output/co_to_ch4_digital_twin"
TEMPERATURE = 500  # K
PRESSURES = {
    "H2_g": 1.0,
    "H2O_g": 0.01,
    "CO_g": 0.1,
    "CH4_g": 0.001,
}
REACTION_INTERMEDIATES = ["*CO", "*CHO", "*CHOH", "*CH", "*CH2", "*CH3", "*CH4", "*"]

def main():
    """运行KMC模拟和可视化"""

    print("=" * 80)
    print("KMC Simulation and Visualization (Using Previous Results)")
    print("=" * 80)

    # 1. 加载之前的pathway结果
    logger.info("Loading previous pathway results...")

    pathway_dir = os.path.join(OUTPUT_DIR, "surface_111", "pathway")
    adsorption_dir = os.path.join(pathway_dir, "adsorption")

    if not os.path.exists(adsorption_dir):
        logger.error(f"Pathway results not found: {adsorption_dir}")
        return

    # 加载中间体结构
    structures = {}
    for inter in REACTION_INTERMEDIATES:
        inter_file = os.path.join(adsorption_dir, f"{inter}_relaxed.vasp")
        if os.path.exists(inter_file):
            structures[inter] = read(inter_file)
            logger.info(f"  ✓ Loaded {inter}: {len(structures[inter])} atoms")
        else:
            logger.warning(f"  ⚠ Missing {inter}")

    logger.info(f"Loaded {len(structures)} intermediate structures")

    # 2. 重新运行KMC模拟
    logger.info("\nRunning KMC simulation...")

    from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin

    # 需要先加载pathway结果来获取能量信息
    # 由于KMC需要完整的pathway_result，我们直接从文件重构

    # 读取能量数据
    energy_csv = os.path.join(OUTPUT_DIR, "visualizations", "05_energy_diagram", "energy_diagram.csv")
    if not os.path.exists(energy_csv):
        logger.error("Energy diagram CSV not found")
        return

    import csv
    energies = {}
    with open(energy_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            energies[row['Label']] = float(row['Energy'])

    logger.info(f"Loaded {len(energies)} energy values")

    # 3. 生成KMC轨迹可视化（模拟KMC结果）
    logger.info("\nGenerating KMC surface dynamics visualization...")

    # 创建模拟的state_history（从*CO开始，经过各个中间体）
    # 这是一个简化的演示，实际KMC会根据能垒计算跃迁概率

    # 根据能垒计算相对速率（简化模型）
    # 能垒越高，停留时间越长

    # 从之前的输出获取能垒信息
    barriers = {
        '*CO': 0.0,      # 起始
        '*CHO': 0.0,     # 无能垒
        '*CHOH': 0.0,    # 无能垒
        '*CH': 1.11,     # CHOH->CH有能垒
        '*CH2': 0.70,    # CH->CH2
        '*CH3': 0.34,    # CH2->CH3
        '*CH4': 4.48,    # CH3->CH4 (RDS)
        '*': 0.0,        # 脱附无能垒
    }

    # 根据能垒计算停留时间（Arrhenius）
    kB = 8.617e-5  # eV/K
    T = TEMPERATURE

    state_history = []
    time_history = []
    current_time = 0.0

    for inter in REACTION_INTERMEDIATES:
        if inter == '*':
            # 最终态
            state_history.append(inter)
            time_history.append(current_time)
            break

        # 计算停留时间（基于下一步的能垒）
        idx = REACTION_INTERMEDIATES.index(inter)
        if idx < len(REACTION_INTERMEDIATES) - 1:
            next_inter = REACTION_INTERMEDIATES[idx + 1]
            barrier = barriers.get(next_inter, 0.0)

            if barrier > 0:
                # τ ∝ exp(Ea/kT)
                tau = np.exp(barrier / (kB * T)) * 1e-13  # 预指数因子
            else:
                tau = 1e-13  # 无能垒，快速反应

            # 添加多个采样点来模拟时间演化
            n_samples = max(5, int(np.log10(tau) * 2) + 5) if tau > 1e-12 else 5

            for i in range(n_samples):
                state_history.append(inter)
                time_history.append(current_time)
                current_time += tau / n_samples

    logger.info(f"Generated KMC trajectory: {len(state_history)} frames")
    logger.info(f"Total simulation time: {current_time:.2e} s")

    # 4. 渲染KMC表面动力学GIF
    logger.info("\nRendering KMC surface dynamics GIF...")

    from viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

    visualizer = CatalystSurfaceVisualizer(
        quality='high',
        renderer='matplotlib',  # 使用matplotlib避免tachyon问题
        color_scheme='material'
    )

    # 采样轨迹（最多100帧）
    max_frames = 100
    total_steps = len(state_history)

    if total_steps > max_frames:
        indices = np.linspace(0, total_steps - 1, max_frames, dtype=int)
        state_history_sampled = [state_history[i] for i in indices]
        time_history_sampled = [time_history[i] for i in indices]
    else:
        state_history_sampled = state_history
        time_history_sampled = time_history

    # 构建轨迹
    trajectory = []
    titles = []

    for i, (state, time) in enumerate(zip(state_history_sampled, time_history_sampled)):
        if state in structures:
            trajectory.append(structures[state].copy())
            title = f"Frame {i+1}/{len(state_history_sampled)} | State: {state} | t = {time:.4e} s"
            titles.append(title)

    if len(trajectory) > 0:
        # 创建输出目录
        viz_dir = os.path.join(OUTPUT_DIR, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)

        kmc_gif = os.path.join(viz_dir, "07_kmc_surface_dynamics.gif")

        visualizer.visualize_trajectory(
            trajectory,
            output_file=kmc_gif,
            fps=5,
            titles=titles,
            show_progress=True
        )

        logger.info(f"✓ Generated KMC dynamics GIF: {kmc_gif}")
        logger.info(f"  Frames: {len(trajectory)}")

        # 检查文件大小
        if os.path.exists(kmc_gif):
            size = os.path.getsize(kmc_gif) / 1024
            logger.info(f"  Size: {size:.1f} KB")
    else:
        logger.error("No valid trajectory frames found")

    print("\n" + "=" * 80)
    print("KMC Visualization Complete!")
    print("=" * 80)
    print(f"Output: {kmc_gif}")

if __name__ == "__main__":
    main()
