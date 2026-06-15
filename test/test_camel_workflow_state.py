import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from camel_agents.workflow import WorkflowState, ReactionContext


def test_workflow_state_defaults():
    state = WorkflowState()
    assert state.reaction_context is None


def test_reaction_context_fields():
    ctx = ReactionContext(
        reaction_type="hydrogenation",
        initial_adsorbate="*CO",
        intermediates=["*CO", "*CHO", "*CH2"],
        constraints="use user constraints",
    )
    assert ctx.initial_adsorbate == "*CO"
