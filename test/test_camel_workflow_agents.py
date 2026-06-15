import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from camel_agents.workflow import CatDTCamelWorkflow


def test_workflow_has_agents():
    wf = CatDTCamelWorkflow(initialize_agents=False)
    assert hasattr(wf, "agent1")
    assert hasattr(wf, "agent4")
    assert hasattr(wf, "agent5")
