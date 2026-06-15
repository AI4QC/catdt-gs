from __future__ import annotations

import pathlib
import stat
import sys

import ase.io
import numpy as np
import yaml
from ase import Atoms

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from core.active_learning.framework import (
    initial_train_and_recommend,
    retrain_with_feedback,
)


def _make_atoms(al_id: str, energy: float, electron: float, potential: float) -> Atoms:
    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    atoms.info["al_id"] = al_id
    atoms.info["REF_energy"] = energy
    atoms.info["electron"] = electron
    atoms.info["potential"] = potential
    atoms.arrays["REF_forces"] = np.zeros((2, 3), dtype=float)
    return atoms


def _write_exec_script(path: pathlib.Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR)


def test_cp_mace_format_extxyz_runs_init_and_retrain(tmp_path: pathlib.Path) -> None:
    data_dir = tmp_path / "data"
    scripts_dir = tmp_path / "mock"
    work_dir = tmp_path / "workspace"
    data_dir.mkdir(parents=True)
    scripts_dir.mkdir(parents=True)

    train_xyz = data_dir / "train.xyz"
    valid_xyz = data_dir / "valid.xyz"
    pool_xyz = data_dir / "pool.xyz"

    ase.io.write(
        str(train_xyz), [_make_atoms("train_0", -1.0, 0.0, 0.1)], format="extxyz"
    )
    ase.io.write(
        str(valid_xyz), [_make_atoms("valid_0", -0.9, 0.0, 0.1)], format="extxyz"
    )
    ase.io.write(
        str(pool_xyz),
        [
            _make_atoms("pool_0", -0.8, 0.1, 0.2),
            _make_atoms("pool_1", -0.7, 0.2, 0.3),
            _make_atoms("pool_2", -0.6, 0.3, 0.4),
        ],
        format="extxyz",
    )

    train_script = scripts_dir / "mock_train.py"
    score_script = scripts_dir / "mock_score.py"

    _write_exec_script(
        train_script,
        """
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--run_dir', required=True)
parser.add_argument('--seed', required=True)
args = parser.parse_args()
run_dir = Path(args.run_dir)
run_dir.mkdir(parents=True, exist_ok=True)
(run_dir / f'mock_seed_{args.seed}.model').write_text('ok', encoding='utf-8')
""".strip(),
    )

    _write_exec_script(
        score_script,
        """
import argparse
import json
import ase.io

parser = argparse.ArgumentParser()
parser.add_argument('--pool_file', required=True)
parser.add_argument('--score_file', required=True)
_ = parser.add_argument('--checkpoints', nargs='*')
args = parser.parse_args()

atoms_list = ase.io.read(args.pool_file, index=':')
rows = []
for i, a in enumerate(atoms_list):
    rows.append({
        'structure_id': str(a.info.get('al_id', f'pool_{i}')),
        'uncertainty_force': float(i + 1),
        'uncertainty_energy': float(i + 1) * 0.1,
        'uncertainty_potential': float(a.info.get('potential', 0.0)),
        'descriptor': [float(i), float(i)],
        'metadata': {'index': i}
    })

with open(args.score_file, 'w', encoding='utf-8') as f:
    json.dump(rows, f, ensure_ascii=False, indent=2)
""".strip(),
    )

    cfg_path = tmp_path / "al.yaml"
    cfg = {
        "workspace_dir": str(work_dir),
        "train_file": str(train_xyz),
        "valid_file": str(valid_xyz),
        "pool_file": str(pool_xyz),
        "committee_seeds": [1, 2],
        "select_top_k": 2,
        "preselect_k": 3,
        "score_weights": {"force": 0.6, "energy": 0.3, "potential": 0.1},
        "diversity_lambda": 0.3,
        "adapter": {
            "type": "command",
            "params": {
                "train_command_template": f"python {train_script} --run_dir {{run_dir}} --seed {{seed}}",
                "score_command_template": f"python {score_script} --pool_file {{pool_file}} --score_file {{score_file}} --checkpoints {{checkpoints}}",
                "checkpoint_glob": "*.model",
            },
        },
    }
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    state0 = initial_train_and_recommend(str(cfg_path))
    rec0 = pathlib.Path(state0["last_recommendation"])
    assert rec0.exists()

    rec_atoms = ase.io.read(str(rec0), index=":")
    assert len(rec_atoms) == 2
    assert all("al_id" in a.info for a in rec_atoms)

    feedback = data_dir / "feedback_iter0.xyz"
    fb0 = rec_atoms[0]
    fb0.info["REF_energy"] = -1.23
    fb0.arrays["REF_forces"] = np.zeros((len(fb0), 3), dtype=float)
    ase.io.write(str(feedback), [fb0], format="extxyz")

    state1 = retrain_with_feedback(str(cfg_path), str(feedback))
    rec1 = pathlib.Path(state1["last_recommendation"])
    assert rec1.exists()

    train1_atoms = ase.io.read(state1["train_file"], index=":")
    assert any("electron" in a.info and "potential" in a.info for a in train1_atoms)
    assert any(
        "REF_energy" in a.info and "REF_forces" in a.arrays for a in train1_atoms
    )


def test_cp_mace_parser_accepts_same_extxyz_format(tmp_path: pathlib.Path) -> None:
    cp_root = PROJECT_ROOT / "deps" / "CP-MACE"
    if str(cp_root) not in sys.path:
        sys.path.insert(0, str(cp_root))

    from mace.data.utils import load_from_xyz

    xyz_path = tmp_path / "cp_format.xyz"
    atoms = _make_atoms("cpfmt_0", energy=-1.1, electron=0.25, potential=0.33)
    ase.io.write(str(xyz_path), [atoms], format="extxyz")

    _, cfgs = load_from_xyz(
        file_path=str(xyz_path),
        config_type_weights={"Default": 1.0},
        energy_key="REF_energy",
        forces_key="REF_forces",
        stress_key="REF_stress",
        virials_key="REF_virials",
        dipole_key="REF_dipole",
        charges_key="REF_charges",
        head_name="Default",
    )

    assert len(cfgs) == 1
    cfg = cfgs[0]
    assert float(cfg.energy) == -1.1
    assert float(cfg.electron) == 0.25
    assert float(cfg.potential) == 0.33
