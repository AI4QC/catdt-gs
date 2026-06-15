#!/usr/bin/env python
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'deps', 'catplot'))

"""
Energy Diagram Plotter - High-quality free energy diagram visualization
"""

import numpy as np
from typing import List, Dict, Optional, Union, Tuple
from dataclasses import dataclass, field
import re

plt = None
rcParams = None


def _ensure_matplotlib():
    global plt, rcParams
    if plt is not None and rcParams is not None:
        return plt, rcParams
    try:
        import matplotlib.pyplot as _plt
        from matplotlib import rcParams as _rc_params
    except Exception as exc:
        raise RuntimeError(f"matplotlib import failed in energy_diagram_plotter: {exc}") from exc
    plt = _plt
    rcParams = _rc_params
    return plt, rcParams

# Try to import catplot for smooth interpolation
try:
    from catplot.ep_components.ep_lines import ElementaryLine
    from catplot.interpolate import get_potential_energy_points
    HAS_CATPLOT = True
except ImportError:
    HAS_CATPLOT = False
    print("Warning: catplot not found. Using fallback interpolation.")


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class EnergyStep:
    """Single step in energy diagram."""
    label: str                          # Step label (e.g., 'H₂ + O₂', 'TS1')
    energy: float                       # Energy value (kcal/mol, kJ/mol, eV, etc.)
    is_ts: bool = False                 # Whether this is a transition state
    color: Optional[str] = None         # Custom color for this step
    line_style: Optional[str] = None    # Custom line style

@dataclass
class DiagramStyle:
    """Styling configuration for energy diagram."""
    # Colors
    color_scheme: str = 'material'      # 'material', 'nature', 'science', 'custom'
    primary_color: str = '#1976D2'      # Default blue
    ts_color: str = '#D32F2F'           # Red for transition states
    background_color: str = '#FFFFFF'   # White background

    # Line properties
    line_width: float = 2.5
    ts_line_width: float = 2.0
    platform_length: float = 1.0        # Length of horizontal platforms
    ts_peak_width: float = 0.8          # Width of TS peak

    # Shadow effects
    use_shadow: bool = False            # Disabled by default
    shadow_color: str = '#CCCCCC'
    shadow_offset: float = 0.05         # Offset as fraction of y-range

    # Interpolation
    interp_method: str = 'spline'       # 'spline' or 'quadratic'
    interp_points: int = 100            # Points per segment

    # Labels
    label_fontsize: int = 12
    energy_fontsize: int = 10
    show_energies: bool = True
    energy_format: str = '.2f'
    state_label_rotation: float = 0.0

    # Axes
    xlabel: str = ''                    # No xlabel by default
    ylabel: str = 'Free energy'
    ylabel_units: str = 'kcal/mol'
    grid: bool = False

    # Figure
    figsize: Tuple[float, float] = (10, 6)
    dpi: int = 300


# =============================================================================
# Color Schemes
# =============================================================================

COLOR_SCHEMES = {
    'material': {
        'primary': '#1976D2',      # Blue 700
        'ts': '#D32F2F',           # Red 700
        'secondary': '#388E3C',     # Green 700
        'accent': '#F57C00',       # Orange 700
        'background': '#FAFAFA',
        'text': '#212121',
    },
    'nature': {
        'primary': '#0C4B8E',      # Nature blue
        'ts': '#C3423F',           # Nature red
        'secondary': '#7CB342',     # Nature green
        'accent': '#FFB300',       # Nature yellow
        'background': '#FFFFFF',
        'text': '#000000',
    },
    'science': {
        'primary': '#003C71',      # Science blue
        'ts': '#B31B1B',           # Science red
        'secondary': '#00A79D',     # Science teal
        'accent': '#FF6F00',       # Science orange
        'background': '#FFFFFF',
        'text': '#000000',
    },
    'elegant': {
        'primary': '#34495E',      # Dark gray-blue
        'ts': '#E74C3C',           # Soft red
        'secondary': '#16A085',     # Soft teal
        'accent': '#F39C12',       # Soft orange
        'background': '#ECF0F1',    # Light gray
        'text': '#2C3E50',         # Dark blue-gray
    },
}


# =============================================================================
# Utility Functions
# =============================================================================

def format_chemical_formula(formula: str) -> str:
    """
    Convert chemical formula to LaTeX format with subscripts and superscripts.

    Examples:
        'H2O' -> 'H$_2$O'
        'CO2' -> 'CO$_2$'
        'H+' -> 'H$^+$'
        'O2-' -> 'O$_2^-$'
        '1/2 O2' -> '$\\frac{1}{2}$ O$_2$'

    Args:
        formula: Chemical formula string

    Returns:
        LaTeX formatted string
    """
    # Handle fractions
    formula = re.sub(r'(\d+)/(\d+)', r'$\\frac{\1}{\2}$ ', formula)

    # Handle charges (+ or -)
    formula = re.sub(r'(\d*)([+-])', r'$^{\1\2}$', formula)

    # Handle subscripts (numbers after letters)
    formula = re.sub(r'([A-Za-z])(\d+)', r'\1$_{\2}$', formula)

    # Clean up multiple $ signs
    formula = re.sub(r'\$+', '$', formula)

    return formula


def parse_energy_input(
    energies: Union[List[float], List[Dict], List[EnergyStep]],
    labels: Optional[List[str]] = None,
    is_ts_list: Optional[List[bool]] = None
) -> List[EnergyStep]:
    """
    Parse various input formats into list of EnergyStep objects.

    Args:
        energies: Energy values or dict/EnergyStep list
        labels: Optional labels for each step
        is_ts_list: Optional list indicating which steps are TS

    Returns:
        List of EnergyStep objects
    """
    # Case 1: Already EnergyStep objects
    if all(isinstance(e, EnergyStep) for e in energies):
        return energies

    # Case 2: List of dictionaries
    if all(isinstance(e, dict) for e in energies):
        return [
            EnergyStep(
                label=e.get('label', e.get('name', f'State {i}')),
                energy=e['energy'],
                is_ts=e.get('is_ts', False),
                color=e.get('color'),
                line_style=e.get('line_style')
            )
            for i, e in enumerate(energies)
        ]

    # Case 3: List of floats
    if all(isinstance(e, (int, float)) for e in energies):
        steps = []
        for i, energy in enumerate(energies):
            # Auto-detect labels
            if labels and i < len(labels):
                label = labels[i]
            else:
                label = f'State {i}'

            # Auto-detect TS
            if is_ts_list and i < len(is_ts_list):
                is_ts = is_ts_list[i]
            else:
                # Heuristic: odd indices might be TS if energy is higher than neighbors
                is_ts = False
                if i > 0 and i < len(energies) - 1:
                    if energy > energies[i-1] and energy > energies[i+1]:
                        is_ts = True

            steps.append(EnergyStep(label=label, energy=energy, is_ts=is_ts))

        return steps

    raise ValueError("Invalid energy input format")


def smooth_barrier_curve(
    E_IS: float,
    E_TS: float,
    E_FS: float,
    n_points: int = 100,
    platform_length: float = 1.0,
    peak_width: float = 0.8,
    method: str = 'spline'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate smooth energy barrier curve using catplot's interpolation.

    Args:
        E_IS: Initial state energy
        E_TS: Transition state energy
        E_FS: Final state energy
        n_points: Number of interpolation points
        platform_length: Length of IS and FS platforms
        peak_width: Width of TS peak
        method: 'spline' or 'quadratic'

    Returns:
        (x, y) arrays for the barrier curve
    """
    if HAS_CATPLOT:
        # Use catplot's high-quality interpolation
        try:
            x, y = get_potential_energy_points(
                energies=np.array([E_IS, E_TS, E_FS]),
                n=n_points,
                hline_length=platform_length,
                peak_width=peak_width,
                kind=method
            )
            return x, y
        except ValueError as e:
            # Catplot fails for cases where TS is not highest (e.g., desorption)
            # Fall back to scipy interpolation
            pass

    # Fallback: simple cubic spline (works for all energy profiles)
    from scipy.interpolate import CubicSpline

    # Define key points: IS platform -> barrier peak -> FS platform
    x_key = np.array([
        0,                                          # Start of IS platform
        platform_length,                            # End of IS platform / start of barrier
        platform_length + peak_width / 2,          # TS peak position
        platform_length + peak_width,              # End of barrier / start of FS platform
        2 * platform_length + peak_width           # End of FS platform
    ])
    y_key = np.array([E_IS, E_IS, E_TS, E_FS, E_FS])

    # Interpolate
    cs = CubicSpline(x_key, y_key, bc_type='clamped')
    x = np.linspace(0, x_key[-1], n_points * 3)
    y = cs(x)

    return x, y


# =============================================================================
# Main Plotter Class
# =============================================================================

class EnergyDiagramPlotter:
    """
    High-quality energy diagram plotter.

    Combines the best features of pMuTT and catplot to create publication-ready
    energy diagrams with minimal user input.

    Example:
        >>> plotter = EnergyDiagramPlotter()
        >>>
        >>> # Simple input: just energies and labels
        >>> energies = [0.0, 1.2, -0.5, 0.8, -2.5]
        >>> labels = ['H₂ + O₂', 'TS1', 'H₂O₂', 'TS2', '2 H₂O']
        >>>
        >>> plotter.plot(energies, labels)
        >>> plotter.save('energy_diagram.png', dpi=300)
    """

    def __init__(self, style: Optional[DiagramStyle] = None, apply_matplotlib_style: bool = True):
        """
        Initialize plotter with optional custom style.

        Args:
            style: DiagramStyle object with custom styling
            apply_matplotlib_style: Whether to override global matplotlib rcParams
        """
        self.style = style or DiagramStyle()
        self.fig = None
        self.ax = None
        self.steps = []
        self.apply_matplotlib_style = apply_matplotlib_style

        # Apply color scheme
        self._apply_color_scheme()

        # Set matplotlib style
        if self.apply_matplotlib_style:
            self._set_matplotlib_style()

    def _apply_color_scheme(self):
        """Apply predefined color scheme to style."""
        if self.style.color_scheme in COLOR_SCHEMES:
            scheme = COLOR_SCHEMES[self.style.color_scheme]
            self.style.primary_color = scheme['primary']
            self.style.ts_color = scheme['ts']
            self.style.background_color = scheme['background']

    def _set_matplotlib_style(self):
        """Set matplotlib rcParams for publication quality."""
        _, rc_params = _ensure_matplotlib()
        rc_params['font.family'] = 'sans-serif'
        rc_params['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
        rc_params['font.size'] = 11
        rc_params['axes.labelsize'] = 12
        rc_params['axes.titlesize'] = 14
        rc_params['xtick.labelsize'] = 10
        rc_params['ytick.labelsize'] = 10
        rc_params['legend.fontsize'] = 10
        rc_params['figure.dpi'] = 100
        rc_params['savefig.dpi'] = self.style.dpi
        rc_params['savefig.bbox'] = 'tight'
        rc_params['axes.linewidth'] = 1.2
        rc_params['axes.edgecolor'] = '#333333'
        rc_params['axes.labelcolor'] = '#333333'
        rc_params['text.usetex'] = False  # Can be set to True if LaTeX is available

    def plot(
        self,
        energies: Union[List[float], List[Dict], List[EnergyStep]],
        labels: Optional[List[str]] = None,
        is_ts_list: Optional[List[bool]] = None,
        reference_index: int = 0,
        show: bool = False,
        ax=None,
    ):
        """
        Create energy diagram.

        Args:
            energies: Energy values, dict list, or EnergyStep list
            labels: Optional labels for each step
            is_ts_list: Optional list indicating which steps are TS
            reference_index: Index of reference state (set to 0 energy)
            show: Whether to display the plot immediately
            ax: Optional matplotlib axis to draw into
        """
        # Parse input
        self.steps = parse_energy_input(energies, labels, is_ts_list)

        # Set reference energy to zero
        ref_energy = self.steps[reference_index].energy
        for step in self.steps:
            step.energy -= ref_energy

        # Create figure
        plt_mod, _ = _ensure_matplotlib()
        if ax is None:
            self.fig, self.ax = plt.subplots(
                figsize=self.style.figsize,
                facecolor=self.style.background_color
            )
        else:
            self.ax = ax
            self.fig = ax.figure
        self.ax.set_facecolor(self.style.background_color)

        # Plot energy profile
        self._plot_energy_profile()

        # Format axes
        self._format_axes()

        # Add labels
        self._add_labels()

        if show:
            plt_mod.show()

    def _plot_energy_profile(self):
        """Plot the main energy profile with barriers using catplot's approach."""
        x_offset = 0.0
        x_positions = [None] * len(self.steps)
        ts_positions = []

        # Color for barrierless connections
        barrierless_color = '#999999'
        barrierless_style = '--'

        # First pass: identify all barriers and handle the first state if it has no barrier before it
        i = 0

        # Handle first state if it doesn't start with a barrier sequence
        if len(self.steps) > 0 and not self.steps[0].is_ts:
            first_step = self.steps[0]
            # Check if first state is followed by a barrier
            if not (len(self.steps) > 1 and self.steps[1].is_ts):
                # First state with no barrier after it - draw it as simple platform
                x_positions[0] = x_offset + self.style.platform_length / 2
                color = first_step.color if first_step.color else self.style.primary_color

                x_platform = np.array([x_offset, x_offset + self.style.platform_length])
                y_platform = np.array([first_step.energy, first_step.energy])

                self.ax.plot(
                    x_platform,
                    y_platform,
                    color=color,
                    linewidth=self.style.line_width,
                    solid_capstyle='round',
                    zorder=2
                )
                x_offset += self.style.platform_length
                i = 1

        # Process all barriers
        while i < len(self.steps):
            step = self.steps[i]

            # Skip TS steps - they're handled as part of barriers
            if step.is_ts:
                i += 1
                continue

            # Check if next step is a TS (forming a barrier)
            has_barrier = (i + 1 < len(self.steps) and
                          self.steps[i + 1].is_ts and
                          i + 2 < len(self.steps))

            if has_barrier:
                # Draw barrier: IS (i) -> TS (i+1) -> FS (i+2)
                ts_step = self.steps[i + 1]
                fs_step = self.steps[i + 2]

                # Generate full reaction curve (IS platform + barrier + FS platform)
                x_full, y_full = smooth_barrier_curve(
                    E_IS=step.energy,
                    E_TS=ts_step.energy,
                    E_FS=fs_step.energy,
                    n_points=self.style.interp_points,
                    platform_length=self.style.platform_length,
                    peak_width=self.style.ts_peak_width,
                    method=self.style.interp_method
                )

                # Check if this IS was already drawn as the FS of previous barrier
                is_already_drawn = (x_positions[i] is not None)

                if is_already_drawn:
                    # Skip the IS platform part, only draw from barrier start
                    platform_end_idx = np.argmax(x_full >= self.style.platform_length)
                    x_plot = x_full[platform_end_idx:] - self.style.platform_length + x_offset
                    y_plot = y_full[platform_end_idx:]
                else:
                    # Draw the full curve including IS platform
                    x_plot = x_full + x_offset
                    y_plot = y_full
                    # Store IS position
                    x_positions[i] = x_offset + self.style.platform_length / 2

                # Find TS peak position
                peak_idx = np.argmax(y_plot)
                ts_x = x_plot[peak_idx]
                ts_y = y_plot[peak_idx]
                x_positions[i + 1] = ts_x
                ts_positions.append((ts_x, ts_y, ts_step))

                # Store FS position (center of FS platform)
                x_positions[i + 2] = x_plot[-1] - self.style.platform_length / 2

                # Determine color
                color = step.color if step.color else self.style.primary_color
                line_width = self.style.line_width

                # Plot the curve
                self.ax.plot(
                    x_plot,
                    y_plot,
                    color=color,
                    linewidth=line_width,
                    solid_capstyle='round',
                    zorder=2
                )

                # Update x_offset to end of this curve
                x_offset = x_plot[-1]

                # Move to FS (i+2) for next iteration
                i += 2

            else:
                # State without barrier after it
                # Check if there's a next state to connect to
                if i + 1 < len(self.steps) and not self.steps[i + 1].is_ts:
                    # Draw current platform if not drawn yet
                    if x_positions[i] is None:
                        x_positions[i] = x_offset + self.style.platform_length / 2
                        color = step.color if step.color else self.style.primary_color
                        x_platform = np.array([x_offset, x_offset + self.style.platform_length])
                        y_platform = np.array([step.energy, step.energy])

                        self.ax.plot(
                            x_platform,
                            y_platform,
                            color=color,
                            linewidth=self.style.line_width,
                            solid_capstyle='round',
                            zorder=2
                        )
                        x_offset += self.style.platform_length

                    # Draw barrierless connection with dashed line
                    next_step = self.steps[i + 1]
                    gap = 0.3
                    x_next_start = x_offset + gap

                    self.ax.plot(
                        [x_offset, x_next_start],
                        [step.energy, next_step.energy],
                        color=barrierless_color,
                        linestyle=barrierless_style,
                        linewidth=self.style.line_width * 0.8,
                        zorder=1
                    )

                    x_offset = x_next_start
                    i += 1

                else:
                    # Last state or followed by something else - just draw platform
                    if x_positions[i] is None:
                        x_positions[i] = x_offset + self.style.platform_length / 2
                        color = step.color if step.color else self.style.primary_color
                        x_platform = np.array([x_offset, x_offset + self.style.platform_length])
                        y_platform = np.array([step.energy, step.energy])

                        self.ax.plot(
                            x_platform,
                            y_platform,
                            color=color,
                            linewidth=self.style.line_width,
                            solid_capstyle='round',
                            zorder=2
                        )
                        x_offset += self.style.platform_length
                    i += 1

        self.x_positions = x_positions
        self.ts_positions = ts_positions

    def _format_axes(self):
        """Format axis labels, ticks, and grid."""
        # X-axis
        if self.style.xlabel:  # Only set if not empty
            self.ax.set_xlabel(
                self.style.xlabel,
                fontsize=self.style.label_fontsize,
                color='black'
            )
        self.ax.set_xticks([])  # Remove x-ticks for cleaner look

        # Y-axis
        ylabel = self.style.ylabel
        if self.style.ylabel_units:
            ylabel += f' ({self.style.ylabel_units})'
        self.ax.set_ylabel(
            ylabel,
            fontsize=self.style.label_fontsize,
            color='black'
        )

        # Grid
        if self.style.grid:
            self.ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

        # Spine styling
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['left'].set_linewidth(1.5)
        self.ax.spines['bottom'].set_linewidth(1.5)

    def _add_labels(self):
        """Add state labels and energy values."""
        if not hasattr(self, 'x_positions'):
            return

        # Add labels for non-TS states
        for i, (step, x_pos) in enumerate(zip(self.steps, self.x_positions)):
            # Skip TS in label positions (they're part of barriers)
            if step.is_ts:
                continue

            # Format label (chemical formula with subscripts)
            formatted_label = format_chemical_formula(step.label)

            # Add label below x-axis with appropriate spacing
            y_range = self.ax.get_ylim()[1] - self.ax.get_ylim()[0]
            y_label_pos = self.ax.get_ylim()[0] - y_range * 0.06
            self.ax.text(
                x_pos,
                y_label_pos,
                formatted_label,
                ha='center',
                va='top',
                rotation=self.style.state_label_rotation,
                rotation_mode='anchor',
                fontsize=self.style.label_fontsize,
                color='black',
                transform=self.ax.transData,
                clip_on=False
            )

            # Add energy value above platform (closer to curve)
            if self.style.show_energies:
                energy_text = f'{step.energy:{self.style.energy_format}}'
                self.ax.text(
                    x_pos,
                    step.energy + y_range * 0.015,
                    energy_text,
                    ha='center',
                    va='bottom',
                    fontsize=self.style.energy_fontsize,
                    color='black'
                )

        # Add labels for TS states (at peak positions)
        if hasattr(self, 'ts_positions') and self.style.show_energies:
            y_range = self.ax.get_ylim()[1] - self.ax.get_ylim()[0]
            for ts_x, ts_y, ts_step in self.ts_positions:
                # Find the index of this TS in self.steps
                ts_idx = None
                for idx, step in enumerate(self.steps):
                    if step is ts_step:
                        ts_idx = idx
                        break

                # Calculate activation barrier (energy difference from previous intermediate)
                if ts_idx is not None and ts_idx > 0:
                    # Find the previous non-TS state
                    prev_energy = None
                    for j in range(ts_idx - 1, -1, -1):
                        if not self.steps[j].is_ts:
                            prev_energy = self.steps[j].energy
                            break

                    if prev_energy is not None:
                        # Display activation barrier instead of absolute energy
                        barrier = ts_step.energy - prev_energy
                        energy_text = f'{barrier:{self.style.energy_format}}'
                    else:
                        # Fallback to absolute energy if can't find previous state
                        energy_text = f'{ts_step.energy:{self.style.energy_format}}'
                else:
                    # Fallback to absolute energy
                    energy_text = f'{ts_step.energy:{self.style.energy_format}}'

                self.ax.text(
                    ts_x,
                    ts_y + y_range * 0.015,
                    energy_text,
                    ha='center',
                    va='bottom',
                    fontsize=self.style.energy_fontsize,
                    color='black'
                )

    def save(self, filename: str, dpi: Optional[int] = None, **kwargs):
        """
        Save figure to file.

        Args:
            filename: Output filename
            dpi: Resolution (overrides style.dpi)
            **kwargs: Additional arguments passed to plt.savefig
        """
        if self.fig is None:
            raise ValueError("No figure to save. Call plot() first.")

        save_dpi = dpi or self.style.dpi
        self.fig.savefig(filename, dpi=save_dpi, bbox_inches='tight', **kwargs)
        print(f"Saved energy diagram to: {filename}")

    def export_data(self, filename: str):
        """
        Export energy data to CSV file.

        Args:
            filename: Output CSV filename
        """
        import csv

        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Label', 'Energy', 'Is_TS'])
            for step in self.steps:
                writer.writerow([step.label, step.energy, step.is_ts])

        print(f"Exported energy data to: {filename}")


# =============================================================================
# Convenience Functions
# =============================================================================

def quick_plot(
    energies: List[float],
    labels: List[str],
    is_ts_list: Optional[List[bool]] = None,
    output: str = 'energy_diagram.png',
    color_scheme: str = 'material',
    **style_kwargs
):
    """
    Quick one-line energy diagram plotting.

    Args:
        energies: List of energy values
        labels: List of state labels
        is_ts_list: Optional list indicating TS states
        output: Output filename
        color_scheme: Color scheme name
        **style_kwargs: Additional style parameters

    Example:
        >>> quick_plot(
        ...     energies=[0, 1.2, -0.5, 0.8, -2.5],
        ...     labels=['H2 + O2', 'TS1', 'H2O2', 'TS2', '2 H2O'],
        ...     output='co_oxidation.png'
        ... )
    """
    style = DiagramStyle(color_scheme=color_scheme, **style_kwargs)
    plotter = EnergyDiagramPlotter(style=style)
    plotter.plot(energies, labels, is_ts_list)
    plotter.save(output)


# =============================================================================
# Main entry point for testing
# =============================================================================

if __name__ == '__main__':
    print("Energy Diagram Plotter - Testing Examples")
    print("=" * 70)

    # ==========================================================================
    # Test 1: Simple CO oxidation pathway
    # ==========================================================================
    print("\nTest 1: CO Oxidation on Pt(111)")
    print("-" * 70)

    energies_1 = [0.0, 0.95, 0.15, 0.85, -2.35]
    labels_1 = ['CO* + O*', 'TS1', 'CO-O*', 'TS2', 'CO₂ + 2*']
    is_ts_1 = [False, True, False, True, False]

    plotter1 = EnergyDiagramPlotter(style=DiagramStyle(color_scheme='material'))
    plotter1.plot(energies_1, labels_1, is_ts_1)
    plotter1.save('test_co_oxidation.png', dpi=300)
    plotter1.export_data('test_co_oxidation.csv')

    print("✓ Created: test_co_oxidation.png")

    # ==========================================================================
    # Test 2: H2O formation (using dict input)
    # ==========================================================================
    print("\nTest 2: H₂O Formation (dict input)")
    print("-" * 70)

    steps_2 = [
        {'label': 'H₂ + ½ O₂', 'energy': 0.0, 'is_ts': False},
        {'label': 'TS-H₂O₂', 'energy': 1.2, 'is_ts': True},
        {'label': 'H₂O₂', 'energy': -0.3, 'is_ts': False},
        {'label': 'TS-2H₂O', 'energy': 0.5, 'is_ts': True},
        {'label': '2 H₂O', 'energy': -2.8, 'is_ts': False},
    ]

    style2 = DiagramStyle(
        color_scheme='nature',
        figsize=(12, 7),
        use_shadow=True,
        show_energies=True
    )

    plotter2 = EnergyDiagramPlotter(style=style2)
    plotter2.plot(steps_2)
    plotter2.save('test_h2o_formation.png')

    print("✓ Created: test_h2o_formation.png")

    # ==========================================================================
    # Test 3: Multiple color schemes comparison
    # ==========================================================================
    print("\nTest 3: Color Scheme Comparison")
    print("-" * 70)

    energies_3 = [0.0, 0.8, -0.5, 1.2, -1.8]
    labels_3 = ['A', 'TS1', 'B', 'TS2', 'C']

    for scheme in ['material', 'nature', 'science', 'elegant']:
        style = DiagramStyle(color_scheme=scheme, figsize=(8, 5))
        plotter = EnergyDiagramPlotter(style=style)
        plotter.plot(energies_3, labels_3, [False, True, False, True, False])
        plotter.save(f'test_scheme_{scheme}.png')
        print(f"✓ Created: test_scheme_{scheme}.png")

    # ==========================================================================
    # Test 4: Quick plot function
    # ==========================================================================
    print("\nTest 4: Quick Plot Function")
    print("-" * 70)

    quick_plot(
        energies=[0, 1.5, -1.0, 0.5, -3.0],
        labels=['Reactants', 'TS1', 'Intermediate', 'TS2', 'Products'],
        output='test_quick_plot.png',
        color_scheme='elegant',
        ylabel_units='kJ/mol'
    )

    print("✓ Created: test_quick_plot.png")

    # ==========================================================================
    # Summary
    # ==========================================================================
    print("\n" + "=" * 70)
    print("All tests completed successfully!")
    print("=" * 70)
    print("\nGenerated files:")
    print("  - test_co_oxidation.png")
    print("  - test_co_oxidation.csv")
    print("  - test_h2o_formation.png")
    print("  - test_scheme_material.png")
    print("  - test_scheme_nature.png")
    print("  - test_scheme_science.png")
    print("  - test_scheme_elegant.png")
    print("  - test_quick_plot.png")
    print("=" * 70)
