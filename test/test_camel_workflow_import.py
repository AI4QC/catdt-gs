import importlib
import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

def test_camel_workflow_module_import():
    mod = importlib.import_module("camel_agents.workflow")
    assert hasattr(mod, "CatDTCamelWorkflow")
