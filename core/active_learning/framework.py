from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import ase.io
import numpy as np
import yaml

from .acquisition import rank_candidates
from .types import CandidateScore


@dataclass
class AdapterConfig:
    type: str
    params: dict[str, Any]


@dataclass
class ALConfig:
    workspace_dir: Path
    train_file: Path
    valid_file: Path
    pool_file: Path
    adapter: AdapterConfig
    committee_seeds: list[int]
    select_top_k: int
    preselect_k: int
    score_weights: dict[str, float]
    diversity_lambda: float


class BaseAdapter:
    def train_committee(
        self, cfg: ALConfig, train_file: Path, valid_file: Path, iteration: int
    ) -> list[Path]:
        raise NotImplementedError

    def score_pool(
        self, cfg: ALConfig, checkpoints: list[Path], pool_file: Path
    ) -> list[CandidateScore]:
        raise NotImplementedError


class CommandAdapter(BaseAdapter):
    def _run_shell(self, command: str, cwd: Path) -> None:
        proc = subprocess.run(command, shell=True, cwd=str(cwd), check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"Command failed with code {proc.returncode}: {command}")

    def train_committee(
        self, cfg: ALConfig, train_file: Path, valid_file: Path, iteration: int
    ) -> list[Path]:
        train_tpl = cfg.adapter.params["train_command_template"]
        ckpt_glob = cfg.adapter.params.get("checkpoint_glob", "*.model")
        run_root = cfg.workspace_dir / "runs" / f"iter_{iteration}"
        run_root.mkdir(parents=True, exist_ok=True)

        checkpoints: list[Path] = []
        for seed in cfg.committee_seeds:
            run_dir = run_root / f"seed_{seed}"
            run_dir.mkdir(parents=True, exist_ok=True)
            command = train_tpl.format(
                train_file=str(train_file),
                valid_file=str(valid_file),
                seed=seed,
                run_dir=str(run_dir),
                iteration=iteration,
            )
            self._run_shell(command, cwd=run_dir)
            models = sorted(
                run_dir.rglob(ckpt_glob), key=lambda p: p.stat().st_mtime, reverse=True
            )
            if not models:
                raise FileNotFoundError(
                    f"No checkpoints matched '{ckpt_glob}' in {run_dir}"
                )
            checkpoints.append(models[0])
        return checkpoints

    def score_pool(
        self, cfg: ALConfig, checkpoints: list[Path], pool_file: Path
    ) -> list[CandidateScore]:
        score_tpl = cfg.adapter.params["score_command_template"]
        score_file = cfg.workspace_dir / "artifacts" / "scores.json"
        score_file.parent.mkdir(parents=True, exist_ok=True)
        command = score_tpl.format(
            pool_file=str(pool_file),
            checkpoints=" ".join(str(x) for x in checkpoints),
            score_file=str(score_file),
        )
        self._run_shell(command, cwd=cfg.workspace_dir)
        payload = json.loads(score_file.read_text(encoding="utf-8"))
        scores: list[CandidateScore] = []
        for row in payload:
            descriptor = (
                np.array(row["descriptor"], dtype=float)
                if row.get("descriptor") is not None
                else None
            )
            scores.append(
                CandidateScore(
                    structure_id=str(row["structure_id"]),
                    uncertainty_force=float(row.get("uncertainty_force", 0.0)),
                    uncertainty_energy=float(row.get("uncertainty_energy", 0.0)),
                    uncertainty_potential=float(row.get("uncertainty_potential", 0.0)),
                    descriptor=descriptor,
                    metadata=row.get("metadata", {}),
                )
            )
        return scores


class CPMACEAdapter(BaseAdapter):
    def train_committee(
        self, cfg: ALConfig, train_file: Path, valid_file: Path, iteration: int
    ) -> list[Path]:
        command_adapter = CommandAdapter()
        return command_adapter.train_committee(cfg, train_file, valid_file, iteration)

    def score_pool(
        self, cfg: ALConfig, checkpoints: list[Path], pool_file: Path
    ) -> list[CandidateScore]:
        cp_mace_root = Path(cfg.adapter.params["cp_mace_root"])
        if str(cp_mace_root) not in sys.path:
            sys.path.insert(0, str(cp_mace_root))

        from mace.calculators import MACECalculator

        device = cfg.adapter.params.get("device", "cuda")
        dtype = cfg.adapter.params.get("default_dtype", "float64")
        calc = MACECalculator(
            model_paths=[str(x) for x in checkpoints],
            device=device,
            default_dtype=dtype,
        )

        structures = ase.io.read(str(pool_file), index=":")
        scores: list[CandidateScore] = []
        for idx, atoms in enumerate(structures):
            atoms.calc = calc
            _ = atoms.get_potential_energy()

            forces_comm = calc.results.get("forces_comm", None)
            energy_var = calc.results.get("energy_var", 0.0)

            if forces_comm is not None:
                force_var = np.var(forces_comm, axis=0)
                force_std = np.sqrt(np.sum(force_var, axis=1))
                uncertainty_force = float(np.max(force_std))
            else:
                uncertainty_force = 0.0
            uncertainty_energy = float(np.sqrt(max(float(energy_var), 0.0)))

            descriptor = None
            try:
                desc = calc.get_descriptors(atoms.copy(), invariants_only=True)
                if isinstance(desc, np.ndarray) and desc.size > 0:
                    descriptor = (
                        np.mean(desc, axis=(0, 1))
                        if desc.ndim == 3
                        else np.mean(desc, axis=0)
                    )
            except Exception:
                descriptor = None

            structure_id = str(atoms.info.get("al_id", f"pool_{idx:08d}"))
            scores.append(
                CandidateScore(
                    structure_id=structure_id,
                    uncertainty_force=uncertainty_force,
                    uncertainty_energy=uncertainty_energy,
                    uncertainty_potential=0.0,
                    descriptor=descriptor,
                    metadata={"index": idx},
                )
            )
        return scores


def _load_config(config_path: Path) -> ALConfig:
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    workspace = Path(raw["workspace_dir"]).resolve()
    return ALConfig(
        workspace_dir=workspace,
        train_file=Path(raw["train_file"]).resolve(),
        valid_file=Path(raw["valid_file"]).resolve(),
        pool_file=Path(raw["pool_file"]).resolve(),
        adapter=AdapterConfig(
            type=raw["adapter"]["type"], params=dict(raw["adapter"].get("params", {}))
        ),
        committee_seeds=list(raw.get("committee_seeds", [1, 2, 3])),
        select_top_k=int(raw.get("select_top_k", 30)),
        preselect_k=int(raw.get("preselect_k", 200)),
        score_weights=dict(
            raw.get("score_weights", {"force": 0.6, "energy": 0.3, "potential": 0.1})
        ),
        diversity_lambda=float(raw.get("diversity_lambda", 0.4)),
    )


def _adapter_from_config(cfg: ALConfig) -> BaseAdapter:
    if cfg.adapter.type == "cp_mace":
        return CPMACEAdapter()
    if cfg.adapter.type == "command":
        return CommandAdapter()
    raise ValueError(f"Unsupported adapter type: {cfg.adapter.type}")


def _read_pool(pool_file: Path) -> list[Any]:
    return ase.io.read(str(pool_file), index=":")


def _write_recommendations(
    cfg: ALConfig, iteration: int, selected: list[CandidateScore], pool_file: Path
) -> Path:
    rec_dir = cfg.workspace_dir / "recommendations"
    rec_dir.mkdir(parents=True, exist_ok=True)
    pool = _read_pool(pool_file)

    selected_id_set = {s.structure_id for s in selected}
    by_index: dict[int, CandidateScore] = {}
    for s in selected:
        idx = int(s.metadata.get("index", -1))
        if idx >= 0:
            by_index[idx] = s

    output_atoms = []
    for idx, atoms in enumerate(pool):
        sid = str(atoms.info.get("al_id", f"pool_{idx:08d}"))
        if sid in selected_id_set or idx in by_index:
            score = by_index.get(idx)
            atoms.info["al_id"] = sid
            if score is not None:
                atoms.info["al_unc_force"] = score.uncertainty_force
                atoms.info["al_unc_energy"] = score.uncertainty_energy
                atoms.info["al_unc_potential"] = score.uncertainty_potential
            output_atoms.append(atoms)

    rec_xyz = rec_dir / f"iter_{iteration:03d}_recommend.xyz"
    ase.io.write(str(rec_xyz), output_atoms, format="extxyz")
    report = [
        {
            "structure_id": x.structure_id,
            "uncertainty_force": x.uncertainty_force,
            "uncertainty_energy": x.uncertainty_energy,
            "uncertainty_potential": x.uncertainty_potential,
            "metadata": x.metadata,
        }
        for x in selected
    ]
    (rec_dir / f"iter_{iteration:03d}_recommend.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return rec_xyz


def _save_state(cfg: ALConfig, state: dict[str, Any]) -> None:
    cfg.workspace_dir.mkdir(parents=True, exist_ok=True)
    (cfg.workspace_dir / "state.json").write_text(
        json.dumps(state, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _load_state(cfg: ALConfig) -> dict[str, Any]:
    state_file = cfg.workspace_dir / "state.json"
    if not state_file.exists():
        raise FileNotFoundError(f"State file not found: {state_file}")
    return json.loads(state_file.read_text(encoding="utf-8"))


def _append_xyz(base_xyz: Path, extra_xyz: Path, out_xyz: Path) -> None:
    base_atoms = ase.io.read(str(base_xyz), index=":")
    extra_atoms = ase.io.read(str(extra_xyz), index=":")
    merged = list(base_atoms) + list(extra_atoms)
    out_xyz.parent.mkdir(parents=True, exist_ok=True)
    ase.io.write(str(out_xyz), merged, format="extxyz")


def _prune_pool_by_feedback(
    pool_file: Path, feedback_xyz: Path, out_pool: Path
) -> None:
    pool_atoms = _read_pool(pool_file)
    fb_atoms = ase.io.read(str(feedback_xyz), index=":")
    fb_ids = {str(a.info["al_id"]) for a in fb_atoms if "al_id" in a.info}
    if not fb_ids:
        raise ValueError(
            "Feedback xyz must contain al_id in Atoms.info for each labeled structure"
        )

    kept = []
    for idx, atoms in enumerate(pool_atoms):
        sid = str(atoms.info.get("al_id", f"pool_{idx:08d}"))
        if sid not in fb_ids:
            kept.append(atoms)
    out_pool.parent.mkdir(parents=True, exist_ok=True)
    ase.io.write(str(out_pool), kept, format="extxyz")


def initial_train_and_recommend(config_path: str) -> dict[str, Any]:
    cfg = _load_config(Path(config_path).resolve())
    adapter = _adapter_from_config(cfg)

    iteration = 0
    checkpoints = adapter.train_committee(
        cfg, cfg.train_file, cfg.valid_file, iteration
    )
    scores = adapter.score_pool(cfg, checkpoints, cfg.pool_file)
    selected = rank_candidates(
        scores,
        top_k=cfg.select_top_k,
        preselect_k=cfg.preselect_k,
        weights=cfg.score_weights,
        diversity_lambda=cfg.diversity_lambda,
    )
    rec_xyz = _write_recommendations(cfg, iteration, selected, cfg.pool_file)

    state = {
        "iteration": iteration,
        "train_file": str(cfg.train_file),
        "valid_file": str(cfg.valid_file),
        "pool_file": str(cfg.pool_file),
        "checkpoints": [str(x) for x in checkpoints],
        "last_recommendation": str(rec_xyz),
    }
    _save_state(cfg, state)
    return state


def retrain_with_feedback(config_path: str, feedback_xyz: str) -> dict[str, Any]:
    cfg = _load_config(Path(config_path).resolve())
    adapter = _adapter_from_config(cfg)
    state = _load_state(cfg)

    prev_iter = int(state["iteration"])
    next_iter = prev_iter + 1

    train_prev = Path(state["train_file"]).resolve()
    pool_prev = Path(state["pool_file"]).resolve()
    feedback_file = Path(feedback_xyz).resolve()

    train_next = cfg.workspace_dir / "data" / f"train_iter_{next_iter:03d}.xyz"
    pool_next = cfg.workspace_dir / "data" / f"pool_iter_{next_iter:03d}.xyz"

    _append_xyz(train_prev, feedback_file, train_next)
    _prune_pool_by_feedback(pool_prev, feedback_file, pool_next)

    checkpoints = adapter.train_committee(cfg, train_next, cfg.valid_file, next_iter)
    scores = adapter.score_pool(cfg, checkpoints, pool_next)
    selected = rank_candidates(
        scores,
        top_k=cfg.select_top_k,
        preselect_k=cfg.preselect_k,
        weights=cfg.score_weights,
        diversity_lambda=cfg.diversity_lambda,
    )
    rec_xyz = _write_recommendations(cfg, next_iter, selected, pool_next)

    new_state = {
        "iteration": next_iter,
        "train_file": str(train_next),
        "valid_file": str(cfg.valid_file),
        "pool_file": str(pool_next),
        "checkpoints": [str(x) for x in checkpoints],
        "last_recommendation": str(rec_xyz),
    }
    _save_state(cfg, new_state)
    return new_state
