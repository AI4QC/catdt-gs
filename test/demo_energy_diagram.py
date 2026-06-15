#!/usr/bin/env python
"""
Demonstration of improved energy diagram features

Shows the improvements:
1. No shadow beneath lines (cleaner look)
2. No "Reaction Coordinate" label on x-axis
3. State labels centered with platforms
4. TS energy labels at peak positions

Author: Claude
Date: 2026-01-18
"""
import sys
sys.path.insert(0, "core/viz")

from energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle
import matplotlib.pyplot as plt


def demo_improved_features():
    """Demonstrate all the improvements."""
    print("=" * 70)
    print("Improved Energy Diagram Features Demo")
    print("=" * 70)

    # Sample CO oxidation pathway
    energies = [0.0, 0.95, 0.15, 0.85, -2.35]
    labels = ['CO* + O*', 'TS₁', 'CO-O*', 'TS₂', 'CO₂ + 2*']
    is_ts = [False, True, False, True, False]

    # Create improved diagram
    style = DiagramStyle(
        color_scheme='nature',
        ylabel='Free Energy',
        ylabel_units='eV',
        figsize=(12, 7),
        use_shadow=False,       # ✓ No shadow
        show_energies=True,
        line_width=2.5
    )

    plotter = EnergyDiagramPlotter(style=style)
    plotter.plot(energies, labels, is_ts)

    # Add annotations to highlight improvements
    ax = plotter.ax

    # Annotation 1: No shadow
    ax.annotate(
        '✓ No shadow beneath lines',
        xy=(0.5, energies[0] - 0.3),
        xytext=(0.5, energies[0] - 1.0),
        fontsize=10,
        color='green',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.7),
        arrowprops=dict(arrowstyle='->', color='green', lw=1.5)
    )

    # Annotation 2: Labels centered
    ax.annotate(
        '✓ Labels centered with platforms',
        xy=(0.5, ax.get_ylim()[0] - 0.3),
        xytext=(1.5, ax.get_ylim()[0] - 1.5),
        fontsize=10,
        color='blue',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.7),
        arrowprops=dict(arrowstyle='->', color='blue', lw=1.5)
    )

    # Annotation 3: TS labels at peak
    if hasattr(plotter, 'ts_positions') and len(plotter.ts_positions) > 0:
        ts_x, ts_y, _ = plotter.ts_positions[0]
        ax.annotate(
            '✓ TS energy at peak',
            xy=(ts_x, ts_y + 0.15),
            xytext=(ts_x + 0.8, ts_y + 0.5),
            fontsize=10,
            color='red',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6', alpha=0.7),
            arrowprops=dict(arrowstyle='->', color='red', lw=1.5)
        )

    # Add title
    ax.set_title(
        'Improved Energy Diagram Features',
        fontsize=14,
        fontweight='bold',
        pad=15
    )

    plotter.save('output/demo_improved_features.png', dpi=300)
    print("\n✓ Created: demo_improved_features.png")
    print("\nKey improvements:")
    print("  1. ✓ No shadow beneath lines (cleaner, more professional)")
    print("  2. ✓ No 'Reaction Coordinate' label (states speak for themselves)")
    print("  3. ✓ State labels centered with platforms")
    print("  4. ✓ TS energy labels positioned at peak tops")
    print("=" * 70)


def demo_publication_ready():
    """Create a publication-ready example."""
    print("\n" + "=" * 70)
    print("Publication-Ready Example")
    print("=" * 70)

    # More realistic data with units
    steps = [
        {'label': 'CO* + O*', 'energy': 0.00, 'is_ts': False},
        {'label': 'TS-COOH', 'energy': 0.95, 'is_ts': True},
        {'label': 'CO-O*', 'energy': 0.15, 'is_ts': False},
        {'label': 'TS-CO₂', 'energy': 0.85, 'is_ts': True},
        {'label': 'CO₂ + 2*', 'energy': -2.35, 'is_ts': False},
    ]

    style = DiagramStyle(
        color_scheme='nature',
        figsize=(10, 6),
        line_width=2.5,
        use_shadow=False,
        show_energies=True,
        energy_format='.2f',
        ylabel='Free Energy',
        ylabel_units='eV',
        label_fontsize=13,
        energy_fontsize=11
    )

    plotter = EnergyDiagramPlotter(style=style)
    plotter.plot(steps)

    # Add professional title
    plotter.ax.set_title(
        'CO Oxidation Mechanism on Pt(111)',
        fontsize=14,
        fontweight='bold',
        pad=15
    )

    # Add reaction equation
    plotter.ax.text(
        0.02, 0.98,
        'CO* + O* → CO₂ + 2*',
        transform=plotter.ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                 edgecolor='#0C4B8E', linewidth=1.5, alpha=0.9)
    )

    plotter.save('demo_publication_ready.png', dpi=300)
    plotter.save('demo_publication_ready.pdf')  # Vector format

    print("\n✓ Created: demo_publication_ready.png")
    print("✓ Created: demo_publication_ready.pdf")
    print("\nThis figure is ready for:")
    print("  - Journal submissions (Nature, Science, etc.)")
    print("  - Presentations")
    print("  - Posters")
    print("  - Thesis/dissertations")
    print("=" * 70)


def demo_quick_comparison():
    """Show different styles side by side."""
    print("\n" + "=" * 70)
    print("Style Comparison")
    print("=" * 70)

    energies = [0.0, 1.2, -0.5, 0.8, -1.8]
    labels = ['A', 'TS₁', 'B', 'TS₂', 'C']
    is_ts = [False, True, False, True, False]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')
    axes = axes.flatten()

    schemes = ['material', 'nature', 'science', 'elegant']
    titles = ['Material Design', 'Nature', 'Science', 'Elegant']

    for i, (scheme, title) in enumerate(zip(schemes, titles)):
        style = DiagramStyle(
            color_scheme=scheme,
            ylabel='ΔG',
            ylabel_units='eV',
            use_shadow=False,
            show_energies=True
        )

        plotter = EnergyDiagramPlotter(style=style)
        plotter.fig = fig
        plotter.ax = axes[i]
        plotter.plot(energies, labels, is_ts)

        axes[i].set_title(title, fontsize=13, fontweight='bold', pad=10)

    plt.tight_layout()
    plt.savefig('demo_style_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("\n✓ Created: demo_style_comparison.png")
    print("\nAvailable color schemes:")
    print("  - material: Modern, colorful (presentations)")
    print("  - nature: Professional blue/red (Nature journals)")
    print("  - science: Deep blue/red (Science journals)")
    print("  - elegant: Sophisticated gray-blue (general use)")
    print("=" * 70)


if __name__ == '__main__':
    print("\n" + "=" * 70)
    print("ENERGY DIAGRAM PLOTTER - DEMONSTRATION")
    print("=" * 70)

    demo_improved_features()
    demo_publication_ready()
    demo_quick_comparison()

    print("\n" + "=" * 70)
    print("All demonstrations completed!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  1. demo_improved_features.png    - Annotated feature highlights")
    print("  2. demo_publication_ready.png    - Publication-quality example")
    print("  3. demo_publication_ready.pdf    - Vector format")
    print("  4. demo_style_comparison.png     - All 4 color schemes")
    print("\nYou can now use these as templates for your own diagrams!")
    print("=" * 70)
