import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from camel_agents.tools import CatDTTools


def test_tools_has_new_methods():
    tools = CatDTTools(output_base_dir="/tmp")
    assert hasattr(tools, "compute_adsorption_energies")
    assert hasattr(tools, "run_neb_for_steps")
