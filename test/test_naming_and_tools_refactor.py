import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_core_root_has_no_standalone_files():
    core_root = Path('core')
    files = sorted(p.name for p in core_root.iterdir() if p.is_file())
    assert files == []


def test_core_root_keeps_expected_domain_directories():
    core_root = Path('core')
    dirs = {p.name for p in core_root.iterdir() if p.is_dir()}
    expected = {'kmc', 'pathway', 'reconstruction', 'surface', 'viz'}
    assert expected.issubset(dirs)


def test_camel_agents_workflow_entry_exists_without_legacy_aliases():
    mod = importlib.import_module('camel_agents.workflow')
    assert hasattr(mod, 'CatDTCamelWorkflow')
    assert not hasattr(mod, 'CatDTCrewWorkflow')
    assert not hasattr(mod, 'LegacyCatDTWorkflow')


def test_tools_exposes_camel_tool_registry():
    mod = importlib.import_module('camel_agents.tools')
    tools = mod.CatDTTools(output_base_dir='output/test_camel_tools_refactor')
    registry = tools.to_camel_tools()
    assert isinstance(registry, dict)
    assert 'generate_surfaces' in registry
    assert 'run_neb_for_steps' in registry


def test_refactored_modules_no_longer_depend_on_core_root_files():
    targets = [
        'camel_agents/workflow.py',
        'camel_agents/runtime.py',
        'camel_agents/policy.py',
        'camel_agents/tools.py',
        'camel_agents/resumable.py',
    ]
    forbidden = [
        'core.camel_model_backend',
        'core.workflow_checkpoint',
        'from workflow_checkpoint import',
    ]

    for rel in targets:
        source = Path(rel).read_text(encoding='utf-8')
        for pattern in forbidden:
            assert pattern not in source
