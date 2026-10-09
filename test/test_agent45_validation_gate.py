from __future__ import annotations

import sys
from pathlib import Path

from ase import Atoms

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from camel_agents.tools import CatDTTools
from camel_agents.prompts import TaskPrompts


def _make_step_dict(
    tools: CatDTTools,
    name: str,
    reactant: Atoms,
    product: Atoms,
    reactant_formula: str = "*CO",
    product_formula: str = "*CO",
) -> dict:
    workflow = tools._workflow_bridge
    reactant_ads = list(range(len(reactant)))
    product_ads = list(range(len(product)))
    return {
        "name": name,
        "reactant": reactant,
        "product": product,
        "reactant_formula": reactant_formula,
        "product_formula": product_formula,
        "reactant_adsorbate_indices": reactant_ads,
        "product_adsorbate_indices": product_ads,
        "reactant_fixed_atom_count": len(reactant),
        "product_fixed_atom_count": len(product),
        "reactant_fixed_reference_signature": workflow._fixed_coordinate_signature(reactant, len(reactant)),
        "product_fixed_reference_signature": workflow._fixed_coordinate_signature(product, len(product)),
        "atoms_to_add_count": 0,
        "atoms_to_remove_count": 0,
    }


def test_programmatic_validation_fails_severe_interpolation_overlap_in_strict_mode(monkeypatch):
    monkeypatch.setenv("CATDT_STRICT_INTERPOLATION_FATAL", "1")
    tools = CatDTTools(output_base_dir="output/_pytest_agent45_gate")
    workflow = tools._workflow_bridge

    reactant = Atoms(
        symbols=["C", "O"],
        positions=[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    product = Atoms(
        symbols=["C", "O"],
        positions=[[1.2, 0.0, 0.0], [0.0, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    step = _make_step_dict(tools, "swap_crossing", reactant, product)

    report = tools.programmatic_validation(workflow=workflow, step_structures=[step])

    assert report["status"] == "FAIL"
    assert any("interpolated path has severe overlap" in item for item in report["fatal_issues"])


def test_programmatic_validation_cross_step_continuity_is_order_invariant():
    tools = CatDTTools(output_base_dir="output/_pytest_agent45_gate")
    workflow = tools._workflow_bridge

    step1_reactant = Atoms(
        symbols=["C", "O"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    step1_product = Atoms(
        symbols=["C", "O"],
        positions=[[0.1, 0.0, 0.0], [0.0, 0.0, 2.1]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    step2_reactant = Atoms(
        symbols=["O", "C"],
        positions=[[0.0, 0.0, 2.1], [0.1, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    step2_product = Atoms(
        symbols=["O", "C"],
        positions=[[0.0, 0.0, 2.2], [0.2, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )

    step_structures = [
        _make_step_dict(tools, "step_1", step1_reactant, step1_product),
        _make_step_dict(tools, "step_2", step2_reactant, step2_product),
    ]

    report = tools.programmatic_validation(workflow=workflow, step_structures=step_structures)

    assert not any("previous step product core does not match current step reactant core" in item for item in report["warning_issues"])


def test_agent5_prompt_contains_precheck_issue_list():
    prompt = TaskPrompts.pathway_validation(
        steps_details="step=1:swap_crossing",
        structures_payload="[reactant_adsorbate] ...",
        precheck_summary="step_1: endpoint_min=1.1",
        precheck_issues="- swap_crossing: interpolated path has close approaches (min 0.67 Å)",
    )
    assert "Pre-check Issues" in prompt
    assert "interpolated path has close approaches" in prompt


def test_programmatic_validation_does_not_fail_on_migrating_oh_fragment():
    tools = CatDTTools(output_base_dir="output/_pytest_agent45_gate")
    workflow = tools._workflow_bridge

    reactant = Atoms(
        symbols=["C", "H", "H", "O"],
        positions=[
            [3.4035, 3.4035, 17.5100],
            [4.1451, 3.8643, 18.2309],
            [1.7867, 2.4782, 17.8048],
            [2.4918, 2.9076, 18.3386],
        ],
        cell=[10.21, 10.21, 27.61],
        pbc=[True, True, True],
    )
    product = Atoms(
        symbols=["C", "H", "H", "O"],
        positions=[
            [3.4035, 3.4035, 17.5100],
            [3.1727, 3.5898, 18.5979],
            [5.1100, 1.2800, 18.5900],
            [5.1100, 1.2800, 17.6100],
        ],
        cell=[10.21, 10.21, 27.61],
        pbc=[True, True, True],
    )
    step = {
        "name": "*CHOH_to_*CH",
        "reactant": reactant,
        "product": product,
        "reactant_formula": "*CHOH",
        "product_formula": "*CH",
        "reactant_adsorbate_indices": [0, 1, 2, 3],
        "product_adsorbate_indices": [0, 1, 2, 3],
        "reactant_fixed_atom_count": 4,
        "product_fixed_atom_count": 2,
        "reactant_fixed_reference_signature": workflow._fixed_coordinate_signature(reactant, 4),
        "product_fixed_reference_signature": workflow._fixed_coordinate_signature(product, 2),
        "reactant_staged_indices": [],
        "product_staged_indices": [2, 3],
        "atoms_to_add_count": 2,
        "atoms_to_remove_count": 0,
        "staging_hints": [
            {
                "side": "product",
                "fragment_elements": ["O", "H"],
                "fragment_label": "OH",
                "anchor_symbol": "C",
                "candidate_sites": [
                    {"rank": 1, "position": [1.28, 5.11, 17.61], "site_label": "bridge", "dist_from_anchor": 2.70},
                    {"rank": 2, "position": [5.11, 1.28, 17.61], "site_label": "bridge", "dist_from_anchor": 2.70},
                    {"rank": 3, "position": [1.28, 1.28, 17.61], "site_label": "bridge", "dist_from_anchor": 3.00},
                ],
            }
        ],
    }
    tools._agent45_site_retry_state = {
        ("*CHOH_to_*CH", "product", "OH"): 2,
    }

    report = tools.programmatic_validation(workflow=workflow, step_structures=[step])

    assert report["status"] == "PASS"
    assert not any("severe overlap" in item for item in report["fatal_issues"])


def test_energy_gate_rejects_high_force_endpoint_before_neb(monkeypatch):
    tools = CatDTTools(output_base_dir="output/_pytest_agent45_gate")
    workflow = tools._workflow_bridge
    geometry = tools._agent45_geometry(workflow)

    reactant = Atoms(
        symbols=["C", "O"],
        positions=[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    product = reactant.copy()
    step = _make_step_dict(tools, "high_force_endpoint", reactant, product)

    class FakePredictor:
        def _load_model(self):
            return None

    calls = []

    def fake_evaluate(predictor, atoms, adsorbate_indices, focus_indices=None):
        calls.append(len(calls))
        fmax = 0.04 if len(calls) == 1 else 1.60
        return {
            "energy": -10.0,
            "max_force_all": fmax,
            "max_force_adsorbate": fmax,
            "mean_force_adsorbate": fmax / 2.0,
            "max_force_focus": fmax,
            "mean_force_focus": fmax / 2.0,
            "min_pair_distance": 1.20,
            "min_adsorbate_surface_distance": 2.00,
        }

    monkeypatch.setattr(geometry, "_get_uma_predictor_for_energy_gate", lambda output_dir="": FakePredictor())
    monkeypatch.setattr(geometry, "_evaluate_endpoint_with_uma", fake_evaluate)

    report = geometry.run_agent45_energy_gate([step])

    assert report["status"] == "FAIL"
    assert any("max adsorbate force 1.60" in item for item in report["fatal_issues"])
