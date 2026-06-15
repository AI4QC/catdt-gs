#!/usr/bin/env python
"""
Test with 10 reaction steps to verify platform consistency
"""
import sys
sys.path.insert(0, "core/viz")

from energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

# Create 10-step reaction pathway
# Pattern: State -> TS -> State -> TS -> ... (5 barriers)
energies = [
    0.0,    # State 0
    0.8,    # TS1
    -0.3,   # State 1
    1.2,    # TS2
    0.1,    # State 2
    0.9,    # TS3
    -0.5,   # State 3
    1.1,    # TS4
    -0.2,   # State 4
    0.7,    # TS5
    -1.5    # State 5 (final)
]

labels = [
    'S0', 'TS1', 'S1', 'TS2', 'S2',
    'TS3', 'S3', 'TS4', 'S4', 'TS5', 'S5'
]

is_ts = [
    False, True, False, True, False,
    True, False, True, False, True, False
]

# Create diagram
style = DiagramStyle(
    color_scheme='nature',
    ylabel='Free Energy',
    ylabel_units='eV',
    figsize=(16, 7),
    use_shadow=False,
    show_energies=True,
    platform_length=1.0,  # Explicitly set to 1.0
    ts_peak_width=0.8
)

plotter = EnergyDiagramPlotter(style=style)
plotter.plot(energies, labels, is_ts)

# Add title
plotter.ax.set_title(
    '10-Step Reaction Pathway - Platform Length Test',
    fontsize=14,
    fontweight='bold',
    pad=15
)

plotter.save('output/test_10_steps.png', dpi=300)
print("✓ Created: test_10_steps.png")
print(f"Total steps: {len(energies)}")
print(f"Number of barriers: {sum(is_ts)}")
print(f"Platform length setting: {style.platform_length}")
