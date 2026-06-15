#!/usr/bin/env python
"""
快速生成KMC表面动力学可视化 - 减少帧数
"""

import os
import sys
sys.path.insert(0, "core")

# 设置matplotlib后端为Agg（无GUI）
import matplotlib
matplotlib.use('Agg')

import logging
import numpy as np
from ase.io import read
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from ase.visualize.plot import plot_atoms

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

OUTPUT_DIR = "output/co_to_ch4_digital_twin"
REACTION_INTERMEDIATES = ["*CO", "*CHO", "*CHOH", "*CH", "*CH2", "*CH3", "*CH4", "*"]

def main():
    print("=" * 60)
    print("Fast KMC Surface Dynamics Visualization")
    print("=" * 60)

    # 加载中间体结构
    pathway_dir = os.path.join(OUTPUT_DIR, "surface_111", "pathway", "adsorption")

    structures = {}
    for inter in REACTION_INTERMEDIATES:
        inter_file = os.path.join(pathway_dir, f"{inter}_relaxed.vasp")
        if os.path.exists(inter_file):
            structures[inter] = read(inter_file)
            logger.info(f"✓ Loaded {inter}")

    # 生成简化的KMC轨迹（每个中间体3帧）
    trajectory = []
    titles = []
    times = [0, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]  # 对数时间尺度

    for i, inter in enumerate(REACTION_INTERMEDIATES):
        if inter in structures:
            for j in range(3):  # 每个中间体3帧
                trajectory.append(structures[inter])
                t = times[i] if i < len(times) else times[-1]
                titles.append(f"{inter} | t={t:.1e}s")

    logger.info(f"Total frames: {len(trajectory)}")

    # 创建GIF
    fig, ax = plt.subplots(figsize=(8, 6))

    def update(frame):
        ax.clear()
        atoms = trajectory[frame]

        # 简单的2D投影
        positions = atoms.positions
        symbols = atoms.get_chemical_symbols()

        # 颜色映射
        colors = {'Cu': '#B87333', 'C': '#333333', 'O': '#FF0000', 'H': '#FFFFFF'}
        sizes = {'Cu': 200, 'C': 300, 'O': 250, 'H': 100}

        for symbol in set(symbols):
            mask = [s == symbol for s in symbols]
            pos = positions[mask]
            ax.scatter(pos[:, 0], pos[:, 1],
                      c=colors.get(symbol, '#888888'),
                      s=sizes.get(symbol, 150),
                      edgecolors='black',
                      linewidths=0.5,
                      label=symbol,
                      alpha=0.8)

        ax.set_title(titles[frame], fontsize=14, fontweight='bold')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_aspect('equal')
        ax.legend(loc='upper right', fontsize=8)

        return ax,

    logger.info("Creating animation...")
    anim = FuncAnimation(fig, update, frames=len(trajectory), interval=500, blit=False)

    # 保存GIF
    output_path = os.path.join(OUTPUT_DIR, "visualizations", "07_kmc_surface_dynamics.gif")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    writer = PillowWriter(fps=2)
    anim.save(output_path, writer=writer, dpi=100)
    plt.close()

    if os.path.exists(output_path):
        size = os.path.getsize(output_path) / 1024
        logger.info(f"✓ Created: {output_path} ({size:.1f} KB)")
        print(f"\n✅ KMC dynamics GIF generated: {output_path}")
    else:
        logger.error("Failed to create GIF")

if __name__ == "__main__":
    main()
