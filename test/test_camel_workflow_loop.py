import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from camel_agents.workflow import CatDTCamelWorkflow, WorkflowState, CatDTConfig
from ase import Atoms


def test_iteration_cap(monkeypatch):
    wf = CatDTCamelWorkflow(initialize_agents=False)
    state = WorkflowState()
    wf._max_pathway_iters = 2

    def fake_validate(*args, **kwargs):
        return {"status": "FAIL", "feedback": "retry"}

    wf._validate_pathway_with_agent5 = fake_validate
    base = Atoms('H', positions=[[0, 0, 0]])
    cfg = CatDTConfig(reaction_description="test")
    wf._run_pathway_iteration(state, base, cfg)
    assert wf._last_iteration_count == 2
