#!/usr/bin/env python
import sys

sys.path.insert(0, "core/viz")

from energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle


def test_state_label_rotation_option_applies_to_state_labels():
    energies = [0.0, 0.6, -0.2, 0.4, -0.5]
    labels = ["A", "TS1", "B", "TS2", "C"]
    is_ts = [False, True, False, True, False]

    style = DiagramStyle(
        color_scheme="nature",
        ylabel="Free Energy",
        ylabel_units="eV",
        state_label_rotation=45,
    )

    plotter = EnergyDiagramPlotter(style=style)
    plotter.plot(energies, labels, is_ts)

    state_texts = [
        t for t in plotter.ax.texts
        if t.get_text() in {"A", "B", "C"}
    ]

    assert len(state_texts) == 3
    assert all(abs(t.get_rotation() - 45) < 1e-6 for t in state_texts)

