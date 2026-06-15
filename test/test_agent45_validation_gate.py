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


def test_programmatic_validation_blocks_severe_interpolated_collision():
    tools = CatDTTools(output_base_dir="output/_pytest_agent45_gate")
    workflow = tools._workflow_bridge

    reactant = Atoms(
        symbols=["C", "O"],
        positions=[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    product = Atoms(
        symbols=["C", "O"],
        positions=[[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=[True, True, True],
    )
    step = _make_step_dict(tools, "swap_crossing", reactant, product)

    report = tools.programmatic_validation(workflow=workflow, step_structures=[step])

    assert report["status"] == "FAIL"
    assert any("interpolated path has close approaches" in item for item in report["fatal_issues"])


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
        staging_plan_text="(none)",
    )
    assert "程序预检问题明细" in prompt
    assert "interpolated path has close approaches" in prompt
