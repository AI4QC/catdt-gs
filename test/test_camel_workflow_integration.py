import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from ase import Atoms
from camel_agents.workflow import CatDTCamelWorkflow, WorkflowState, CatDTConfig


class DummyTools:
    pass


def test_pipeline_wires_state():
    tools = DummyTools()
    wf = CatDTCamelWorkflow(tools=tools, initialize_agents=False)
    state = WorkflowState(run_id="test", output_base_dir="/tmp")
    wf._max_pathway_iters = 1
    wf._validate_pathway_with_agent5 = lambda s: {"status": "PASS"}
    base = Atoms('H', positions=[[0, 0, 0]])
    cfg = CatDTConfig(reaction_description="test")
    wf._run_pathway_iteration(state, base, cfg)
    assert wf._last_iteration_count == 1
    assert wf.tools is tools
