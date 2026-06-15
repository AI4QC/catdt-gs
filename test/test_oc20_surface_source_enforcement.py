import argparse
import importlib.util
from pathlib import Path
import sys

import pytest


def _load_memento_module():
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    module_path = repo_root / "test_memento_multireaction.py"
    spec = importlib.util.spec_from_file_location("test_memento_multireaction", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mmr = _load_memento_module()


def test_validate_surface_source_requires_oc20():
    args = argparse.Namespace(surface_source="bulk")
    with pytest.raises(ValueError, match="OC20"):
        mmr._validate_surface_source_policy(args)


def test_assert_surface_from_oc20_accepts_in_dir(tmp_path: Path):
    oc20_dir = tmp_path / "slabs"
    oc20_dir.mkdir(parents=True, exist_ok=True)
    good = oc20_dir / "101_0"
    good.write_text("dummy")

    normalized = mmr._assert_surface_from_oc20(good, oc20_dir)
    assert Path(normalized) == good.resolve()


def test_assert_surface_from_oc20_rejects_outside_dir(tmp_path: Path):
    oc20_dir = tmp_path / "slabs"
    oc20_dir.mkdir(parents=True, exist_ok=True)
    bad = tmp_path / "outside_slab"
    bad.write_text("dummy")

    with pytest.raises(ValueError, match="outside"):
        mmr._assert_surface_from_oc20(bad, oc20_dir)
