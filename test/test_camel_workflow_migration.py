import importlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_workflow_module_exposes_camel_entrypoints():
    mod = importlib.import_module("camel_agents.workflow")
    assert hasattr(mod, "CatDTCamelWorkflow")
    assert not hasattr(mod, "CatDTCrewWorkflow")


def test_workflow_module_uses_migrated_model_backend_location():
    source = Path("camel_agents/workflow.py").read_text(encoding="utf-8")
    assert "from camel_agents.camel_model_backend import" in source
    assert "core.camel_model_backend" not in source


def test_evolvability_policy_updates_and_persists(tmp_path):
    mod = importlib.import_module("camel_agents.workflow")
    policy_path = tmp_path / "policy.json"

    policy = mod.EvolvablePolicy(policy_path=policy_path)
    strategy = policy.choose_strategy()
    assert strategy in policy.available_strategies

    policy.update(strategy=strategy, reward=0.75)
    assert policy_path.exists()

    payload = json.loads(policy_path.read_text(encoding="utf-8"))
    assert "q_values" in payload
    assert strategy in payload["q_values"]


def test_agent45_prompts_are_generic_and_multi_round_capable():
    from camel_agents.prompts import TaskPrompts

    prompt4 = TaskPrompts.universal_agent4_pathway_design(
        reaction_description="A catalytic surface reaction from initial state to final state",
        surface_structure_info="surface: generic",
        current_adsorbate_info="adsorbate atoms: generic",
        feedback="",
        history="",
    )
    prompt5 = TaskPrompts.universal_agent5_validation(
        pathway_design="{}",
        structures_info="none",
        atomic_distances="none",
    )

    forbidden = ["CO", "CH4", '"element": "H"', "C-O", "Pt-C", "Pt-O"]
    for token in forbidden:
        assert token not in prompt4
        assert token not in prompt5

    assert "10" in prompt4
    assert "历史" in prompt4
    assert "反馈" in prompt4
