import sys
from pathlib import Path
_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from camel_agents.workflow import CatDTCamelWorkflow
from camel_agents.schemas import PathwayStepSpec


def test_resolve_atom_index_rejects_surface_global_index():
    wf = CatDTCamelWorkflow(initialize_agents=False)
    ads = [64, 65]
    assert wf._resolve_atom_index(0, ads, 70) == 64
    assert wf._resolve_atom_index(65, ads, 70) == 65
    assert wf._resolve_atom_index(10, ads, 70) == -1


def test_resolve_removal_indices_only_from_adsorbates():
    wf = CatDTCamelWorkflow(initialize_agents=False)
    step = PathwayStepSpec(
        step_name="x",
        reactant_formula="*A",
        product_formula="*B",
        reaction_type="elementary_transition",
        atoms_to_add=[],
        atoms_to_modify=[],
        atoms_to_remove=[0, 1, 10, {"index": 64}, {"index": 2}, {"indices": [65, 11]}],
        confidence=0.5,
    )
    resolved = wf._resolve_removal_indices(step, [64, 65])
    assert resolved == [64, 65]
