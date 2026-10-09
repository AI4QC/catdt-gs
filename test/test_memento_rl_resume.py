import importlib
import importlib.util
import json
import random
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _load_rl_monitor_module():
    module_path = ROOT / "scripts" / "monitor_agent45_rl_progress.py"
    spec = importlib.util.spec_from_file_location("monitor_agent45_rl_progress", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rng_state_round_trip_is_lossless():
    mod = importlib.import_module("test_memento_multireaction")
    rng = random.Random(20260309)
    prefix = [rng.random() for _ in range(5)]
    encoded = mod._serialize_rng_state(rng)
    restored = random.Random()
    restored.setstate(mod._deserialize_rng_state(encoded))

    assert prefix
    assert [restored.random() for _ in range(5)] == [rng.random() for _ in range(5)]


def test_resume_fallback_uses_casebank_not_partial_episode_dirs(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")

    (tmp_path / "memento_her_i001_e0001").mkdir()
    (tmp_path / "memento_oer_i001_e0002").mkdir()
    (tmp_path / "memento_co2rr_i001_e0003").mkdir()

    casebank = tmp_path / "agent45_memento_casebank.jsonl"
    with casebank.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps({"run_id": "memento_her_i001_e0001"}) + "\n")
        fh.write(json.dumps({"run_id": "memento_oer_i001_e0002"}) + "\n")

    state = mod._resolve_rl_resume_state(
        output_base_dir=tmp_path,
        resume=True,
        train_batch_size=10,
    )

    assert state["completed_episode_count"] == 2
    assert state["next_global_episode_idx"] == 3
    assert state["next_train_iteration"] == 1
    assert state["next_batch_index"] == 3


def test_resume_state_prefers_checkpoint_when_present(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")
    rng = random.Random(7)
    checkpoint = {
        "completed_episode_count": 17,
        "next_global_episode_idx": 18,
        "next_train_iteration": 2,
        "next_batch_index": 8,
        "post_iteration_state": {
            "iteration": 2,
            "retriever_done": False,
            "bank_extract_done": False,
            "eval_done": False,
        },
        "rng_state": mod._serialize_rng_state(rng),
    }

    mod._save_rl_resume_checkpoint(tmp_path, checkpoint)
    state = mod._resolve_rl_resume_state(
        output_base_dir=tmp_path,
        resume=True,
        train_batch_size=10,
    )

    assert state["completed_episode_count"] == 17
    assert state["next_global_episode_idx"] == 18
    assert state["next_train_iteration"] == 2
    assert state["next_batch_index"] == 8
    assert state["post_iteration_state"]["iteration"] == 2


def test_resume_replays_pending_post_iteration_hooks_for_completed_iteration(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")

    casebank = tmp_path / "agent45_memento_casebank.jsonl"
    with casebank.open("w", encoding="utf-8") as fh:
        for episode_idx in range(1, 51):
            case = ["her", "oer", "co2rr", "nrr"][(episode_idx - 1) % 4]
            train_iter = ((episode_idx - 1) // 10) + 1
            fh.write(json.dumps({"run_id": f"memento_{case}_i{train_iter:03d}_e{episode_idx:04d}"}) + "\n")

    state = mod._resolve_rl_resume_state(
        output_base_dir=tmp_path,
        resume=True,
        train_batch_size=10,
    )

    assert state["completed_episode_count"] == 50
    assert state["next_global_episode_idx"] == 51
    assert state["next_train_iteration"] == 5
    assert state["next_batch_index"] == 11
    assert state["post_iteration_state"]["iteration"] == 5
    assert state["post_iteration_state"]["eval_done"] is False


def test_checkpoint_resume_replays_pending_post_iteration_hooks_for_completed_iteration(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")

    checkpoint = {
        "completed_episode_count": 50,
        "next_global_episode_idx": 51,
        "next_train_iteration": 6,
        "next_batch_index": 1,
        "post_iteration_state": {
            "iteration": 6,
            "retriever_done": False,
            "bank_extract_done": False,
            "eval_done": False,
        },
    }
    mod._save_rl_resume_checkpoint(tmp_path, checkpoint)

    state = mod._resolve_rl_resume_state(
        output_base_dir=tmp_path,
        resume=True,
        train_batch_size=10,
    )

    assert state["completed_episode_count"] == 50
    assert state["next_global_episode_idx"] == 51
    assert state["next_train_iteration"] == 5
    assert state["next_batch_index"] == 11
    assert state["post_iteration_state"]["iteration"] == 5
    assert state["post_iteration_state"]["retriever_done"] is False
    assert state["post_iteration_state"]["bank_extract_done"] is False
    assert state["post_iteration_state"]["eval_done"] is False


def test_resolve_candidate_split_uses_cached_split_without_scanning(tmp_path, monkeypatch):
    mod = importlib.import_module("test_memento_multireaction")

    cache_path = tmp_path / "train_eval_split.json"
    payload = {
        "surface_source": "oc20_is2re_val",
        "oc20_lmdb_path": "/tmp/fake/data.lmdb",
        "seed": 42,
        "eval_set_size": 50,
        "oc20_min_atoms": 36,
        "oc20_max_atoms": 0,
        "oc20_max_scan_entries": 0,
        "train_candidates": [
            {"key": "1", "sid": 11, "natoms": 64, "composition": "Cu64"},
            {"key": "2", "sid": 22, "natoms": 72, "composition": "Cu72"},
        ],
        "eval_candidates": [
            {"key": "9", "sid": 99, "natoms": 80, "composition": "Cu80"},
        ],
    }
    cache_path.write_text(json.dumps(payload), encoding="utf-8")

    monkeypatch.setattr(
        mod,
        "_list_oc20_candidates",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("LMDB scan should not run when cache exists")),
    )

    train_candidates, eval_candidates = mod._resolve_candidate_split(
        output_base_dir=tmp_path,
        lmdb_path=Path("/tmp/fake/data.lmdb"),
        min_atoms=36,
        max_atoms=0,
        max_scan_entries=0,
        eval_size=50,
        seed=42,
    )

    assert train_candidates == payload["train_candidates"]
    assert eval_candidates == payload["eval_candidates"]


def test_resolve_candidate_split_builds_and_persists_cache_when_missing(tmp_path, monkeypatch):
    mod = importlib.import_module("test_memento_multireaction")

    listed_candidates = [
        {"key": "1", "sid": 11, "natoms": 64, "composition": "Cu64"},
        {"key": "2", "sid": 22, "natoms": 72, "composition": "Cu72"},
        {"key": "3", "sid": 33, "natoms": 80, "composition": "Cu80"},
    ]
    split_train = listed_candidates[:2]
    split_eval = listed_candidates[2:]

    monkeypatch.setattr(mod, "_list_oc20_candidates", lambda **kwargs: list(listed_candidates))
    monkeypatch.setattr(mod, "_split_train_eval_candidates", lambda **kwargs: (list(split_train), list(split_eval)))

    train_candidates, eval_candidates = mod._resolve_candidate_split(
        output_base_dir=tmp_path,
        lmdb_path=Path("/tmp/fake/data.lmdb"),
        min_atoms=36,
        max_atoms=0,
        max_scan_entries=0,
        eval_size=50,
        seed=42,
    )

    assert train_candidates == split_train
    assert eval_candidates == split_eval

    written = json.loads((tmp_path / "train_eval_split.json").read_text(encoding="utf-8"))
    assert written["train_candidates"] == split_train
    assert written["eval_candidates"] == split_eval
    assert written["seed"] == 42
    assert written["eval_set_size"] == 50


def test_run_eval_block_reads_memory_but_disables_write(tmp_path, monkeypatch):
    mod = importlib.import_module("test_memento_multireaction")
    seen = {}

    @contextmanager
    def _fake_memento_mode(disable_read: bool, disable_write: bool):
        seen["disable_read"] = bool(disable_read)
        seen["disable_write"] = bool(disable_write)
        yield

    monkeypatch.setattr(mod, "_temporary_memento_mode", _fake_memento_mode)
    monkeypatch.setattr(
        mod,
        "_materialize_clean_slab_from_candidate",
        lambda **kwargs: {
            "seed_surface_path": str(tmp_path / "seed.vasp"),
            "seed_surface_key": "123",
            "seed_surface_sid": 456,
            "surface_atom_count": 12,
        },
    )
    monkeypatch.setattr(
        mod,
        "_prepare_case_base_structure",
        lambda **kwargs: {
            "prepared_base_path": str(tmp_path / "prepared.vasp"),
            "surface_atom_count": 12,
            "random_fixed_site": None,
        },
    )
    monkeypatch.setattr(mod, "_visualize_structure_png", lambda **kwargs: str(tmp_path / "eval.png"))
    monkeypatch.setattr(
        mod,
        "_run_agent45_only_episode",
        lambda **kwargs: {
            "validation_status": "PASS",
            "episode_reward": 1.0,
            "selected_strategy": "balanced",
            "pathway_iterations": 1,
            "metrics": {"neb_converged_ratio": 1.0, "nonzero_barrier_steps": 1, "max_barrier": 0.5},
            "neb_results_summary": {},
            "policy_update": {"enabled": False},
        },
    )
    monkeypatch.setattr(
        mod,
        "_collect_llm_structure_visualizations",
        lambda **kwargs: {
            "llm_step_visualization_dir": str(tmp_path / "viz"),
            "llm_step_structure_files": [],
            "llm_step_png_files": [],
        },
    )

    result = mod._run_eval_block(
        workflow=object(),
        output_base_dir=tmp_path,
        oc20_lmdb_path=tmp_path / "dummy.lmdb",
        eval_candidates=[{"key": "123", "sid": 456}],
        eval_cases=["her"],
        train_iteration_idx=5,
        rng=random.Random(42),
        neb_max_steps=100,
        max_feedback_rounds=1,
        use_random_site=False,
    )

    assert seen == {"disable_read": False, "disable_write": True}
    assert result["records"][0]["memory_mode"]["read_disabled"] is False
    assert result["records"][0]["memory_mode"]["write_disabled"] is True


def test_backend_model_config_preserves_seed():
    mod = importlib.import_module("camel_agents.camel_model_backend")
    cfg = mod._build_model_config(0.2, max_tokens=1024, top_p=0.9, seed=20260310)
    assert cfg["seed"] == 20260310
    assert cfg["top_p"] == 0.9


def test_build_agent45_bank_output_models_validate_when_loaded_from_script_path():
    script_path = ROOT / "scripts" / "build_agent45_banks.py"
    spec = importlib.util.spec_from_file_location("build_agent45_banks_test_mod", script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)

    knowledge = module.KnowledgeOutput.model_validate(
        {
            "items": [
                {
                    "title": "k1",
                    "reaction_type": "HER",
                    "transition_signature": ["*h->h2(g)"],
                    "knowledge": "k",
                    "guidance": "g",
                    "pattern_type": "mixed",
                    "failure_mode": "",
                    "confidence": 0.8,
                    "source_case_ids": ["case1"],
                }
            ]
        }
    )
    skill = module.SkillOutput.model_validate(
        {
            "items": [
                {
                    "title": "s1",
                    "reaction_type": "HER",
                    "transition_signature": ["*h->h2(g)"],
                    "skill": "s",
                    "action_template": "1. do x",
                    "applicability": "general",
                    "pattern_type": "mixed",
                    "failure_mode": "",
                    "confidence": 0.8,
                    "source_case_ids": ["case1"],
                }
            ]
        }
    )

    assert len(knowledge.items) == 1
    assert len(skill.items) == 1


def test_evolvable_policy_persists_rng_state(tmp_path):
    mod = importlib.import_module("camel_agents.policy")
    policy_path = tmp_path / "policy.json"
    policy = mod.EvolvablePolicy(
        policy_path=policy_path,
        available_strategies=["balanced", "conservative", "exploratory"],
        epsilon=1.0,
        seed=20260310,
    )
    first = policy.choose_strategy()
    second = policy.choose_strategy()
    policy._save()

    reloaded = mod.EvolvablePolicy(
        policy_path=policy_path,
        available_strategies=["balanced", "conservative", "exploratory"],
        epsilon=1.0,
        seed=999,
    )
    expected_next = policy.choose_strategy()
    actual_next = reloaded.choose_strategy()

    assert first in {"balanced", "conservative", "exploratory"}
    assert second in {"balanced", "conservative", "exploratory"}
    assert actual_next == expected_next


def test_write_generation_casebank_snapshot_truncates_at_target_episode(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")

    source_dir = tmp_path / "source"
    source_dir.mkdir()
    casebank_path = source_dir / "agent45_memento_casebank.jsonl"
    rows = [
        {"run_id": "memento_her_i001_e0001", "feedback": "a"},
        {"run_id": "memento_oer_i001_e0002", "feedback": "b"},
        {"run_id": "memento_co2rr_i005_e0050", "feedback": "c"},
        {"run_id": "memento_co2rr_i005_e0050", "feedback": "c2"},
        {"run_id": "memento_nrr_i006_e0051", "feedback": "d"},
    ]
    with casebank_path.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")

    snapshot_dir = tmp_path / "snapshot"
    info = mod._write_generation_casebank_snapshot(
        source_output_base_dir=source_dir,
        snapshot_output_base_dir=snapshot_dir,
        max_episode=50,
    )

    snapshot_casebank = snapshot_dir / "agent45_memento_casebank.jsonl"
    assert snapshot_casebank.exists()
    written = [json.loads(line) for line in snapshot_casebank.read_text().splitlines() if line.strip()]
    assert len(written) == 4
    assert all(mod._parse_run_id(item["run_id"]).get("episode_idx", 0) <= 50 for item in written)
    assert info["kept_case_count"] == 4
    assert info["max_episode"] == 50


def test_write_generation_casebank_window_snapshot_keeps_lower_bound_episode(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")

    source_dir = tmp_path / "source"
    source_dir.mkdir()
    casebank_path = source_dir / "agent45_memento_casebank.jsonl"
    rows = [
        {"run_id": "memento_her_i001_e0009", "feedback": "ep9"},
        {"run_id": "memento_her_i001_e0010", "feedback": "ep10"},
        {"run_id": "memento_her_i002_e0015", "feedback": "ep15"},
        {"run_id": "memento_her_i002_e0020", "feedback": "ep20"},
        {"run_id": "memento_her_i003_e0021", "feedback": "ep21"},
    ]
    with casebank_path.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")

    snapshot_dir = tmp_path / "snapshot"
    info = mod._write_generation_casebank_window_snapshot(
        source_output_base_dir=source_dir,
        snapshot_output_base_dir=snapshot_dir,
        from_episode=10,
        to_episode=20,
    )

    snapshot_casebank = snapshot_dir / "agent45_memento_casebank.jsonl"
    written = [json.loads(line) for line in snapshot_casebank.read_text().splitlines() if line.strip()]
    written_episodes = [mod._parse_run_id(item["run_id"]).get("episode_idx", 0) for item in written]

    assert written_episodes == [10, 15, 20]
    assert info["from_episode"] == 10
    assert info["to_episode"] == 20


def test_save_generation_checkpoint_copies_current_retriever_and_banks(tmp_path):
    mod = importlib.import_module("test_memento_multireaction")

    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "agent45_memento_casebank.jsonl").write_text(
        json.dumps({"run_id": "memento_oer_i005_e0050"}) + "\n",
        encoding="utf-8",
    )
    (source_dir / "agent45_knowledge_bank.jsonl").write_text('{"k":1}\n', encoding="utf-8")
    (source_dir / "agent45_skill_bank.jsonl").write_text('{"s":1}\n', encoding="utf-8")
    (source_dir / "agent45_bank_build_manifest.json").write_text(json.dumps({"knowledge_item_count": 1}), encoding="utf-8")
    (source_dir / "evolvability_policy.json").write_text(json.dumps({"q_values": {"balanced": 1.0}}), encoding="utf-8")
    (source_dir / "rl_resume_checkpoint.json").write_text(json.dumps({"completed_episode_count": 71}), encoding="utf-8")

    retriever_dir = source_dir / "agent45_memento_retriever"
    retriever_dir.mkdir()
    (retriever_dir / "best.pt").write_bytes(b"checkpoint")
    (retriever_dir / "train_meta.json").write_text(
        json.dumps({"model_path": str(retriever_dir / "best.pt"), "best_val_acc": 0.9}),
        encoding="utf-8",
    )
    (retriever_dir / "train_history.jsonl").write_text(
        json.dumps({"model_path": str(retriever_dir / "best.pt")}) + "\n",
        encoding="utf-8",
    )
    (source_dir / "parametric_retriever_training.jsonl").write_text(
        json.dumps({"train_iteration": 7, "train_episode_count": 70, "train_result": {"model_path": str(retriever_dir / "best.pt")}}) + "\n",
        encoding="utf-8",
    )

    info = mod._save_generation_checkpoint(
        output_base_dir=source_dir,
        completed_episode_count=50,
        train_iteration=5,
    )

    checkpoint_dir = source_dir / "generation_checkpoints" / "ep_0050"
    assert checkpoint_dir.exists()
    assert (checkpoint_dir / "agent45_memento_casebank.jsonl").exists()
    assert (checkpoint_dir / "agent45_knowledge_bank.jsonl").exists()
    assert (checkpoint_dir / "agent45_skill_bank.jsonl").exists()
    assert (checkpoint_dir / "evolvability_policy.json").exists()
    assert (checkpoint_dir / "rl_resume_checkpoint.json").exists()
    assert (checkpoint_dir / "agent45_memento_retriever" / "best.pt").exists()

    copied_meta = json.loads((checkpoint_dir / "agent45_memento_retriever" / "train_meta.json").read_text())
    assert str(checkpoint_dir / "agent45_memento_retriever" / "best.pt") == copied_meta["model_path"]
    assert info["completed_episode_count"] == 50
    assert info["train_iteration"] == 5


def test_rebuild_generation_full_bank_state_replays_bank_updates_every_10_episodes(tmp_path, monkeypatch):
    mod = importlib.import_module("test_memento_multireaction")

    source_dir = tmp_path / "source"
    source_dir.mkdir()
    casebank_path = source_dir / "agent45_memento_casebank.jsonl"
    with casebank_path.open("w", encoding="utf-8") as fh:
        for episode_idx in range(1, 51):
            case = ["her", "oer", "co2rr", "nrr"][(episode_idx - 1) % 4]
            train_iter = ((episode_idx - 1) // 10) + 1
            fh.write(json.dumps({"run_id": f"memento_{case}_i{train_iter:03d}_e{episode_idx:04d}"}) + "\n")

    calls = []

    monkeypatch.setattr(
        mod,
        "_copy_current_retriever_checkpoint",
        lambda **kwargs: {"source": "current_checkpoint_copy"},
    )

    def _fake_bank_extract(**kwargs):
        calls.append(
            {
                "completed_episode_count": kwargs["completed_episode_count"],
                "train_iteration": kwargs["train_iteration"],
            }
        )
        return {
            "status": "OK",
            "train_iteration": kwargs["train_iteration"],
            "completed_episode_count": kwargs["completed_episode_count"],
            "manifest": {"knowledge_item_count": 1, "skill_item_count": 1},
        }

    monkeypatch.setattr(mod, "_run_periodic_bank_extraction", _fake_bank_extract)

    snapshot_dir = tmp_path / "snapshot"
    manifest = mod._rebuild_generation_full_bank_state(
        source_output_base_dir=source_dir,
        snapshot_output_base_dir=snapshot_dir,
        max_episode=50,
        train_batch_size=10,
        seed=42,
        bank_max_cases=120,
        bank_batch_size=5,
        bank_knowledge_per_group=6,
        bank_skill_per_group=6,
        bank_temperature=0.2,
        bank_max_tokens=4096,
        bank_top_p=0.95,
        parametric_min_cases=30,
        parametric_max_cases=1500,
        parametric_max_pairs=20000,
        parametric_max_pos_per_query=4,
        parametric_max_neg_per_query=8,
        parametric_epochs=1,
        parametric_batch_size=16,
        parametric_learning_rate=2e-5,
        parametric_max_len=256,
        parametric_val_ratio=0.15,
        parametric_model_name="sentence-transformers/all-MiniLM-L6-v2",
        parametric_device="cuda",
        parametric_resume_from_checkpoint=False,
    )

    assert calls == [
        {"completed_episode_count": 10, "train_iteration": 1},
        {"completed_episode_count": 20, "train_iteration": 2},
        {"completed_episode_count": 30, "train_iteration": 3},
        {"completed_episode_count": 40, "train_iteration": 4},
        {"completed_episode_count": 50, "train_iteration": 5},
    ]
    assert manifest["staged_bank_updates"] == [
        {"from_episode": 1, "to_episode": 10, "train_iteration": 1},
        {"from_episode": 10, "to_episode": 20, "train_iteration": 2},
        {"from_episode": 20, "to_episode": 30, "train_iteration": 3},
        {"from_episode": 30, "to_episode": 40, "train_iteration": 4},
        {"from_episode": 40, "to_episode": 50, "train_iteration": 5},
    ]


def test_monitor_phase_is_bank_building_when_bank_extract_pending():
    mon = _load_rl_monitor_module()

    phase = mon.infer_progress_phase(
        checkpoint={
            "completed_episode_count": 100,
            "post_iteration_state": {
                "iteration": 10,
                "retriever_done": True,
                "bank_extract_done": False,
                "eval_done": False,
            },
        },
        has_rl_process=True,
        has_bank_builder_process=False,
    )

    assert phase == "bank_building"


def test_monitor_phase_is_evaluating_when_bank_done_but_eval_pending():
    mon = _load_rl_monitor_module()

    phase = mon.infer_progress_phase(
        checkpoint={
            "completed_episode_count": 100,
            "post_iteration_state": {
                "iteration": 10,
                "retriever_done": True,
                "bank_extract_done": True,
                "eval_done": False,
            },
        },
        has_rl_process=True,
        has_bank_builder_process=False,
    )

    assert phase == "evaluating"


def test_monitor_phase_is_training_when_post_iteration_hooks_are_complete():
    mon = _load_rl_monitor_module()

    phase = mon.infer_progress_phase(
        checkpoint={
            "completed_episode_count": 100,
            "post_iteration_state": {
                "iteration": 11,
                "retriever_done": False,
                "bank_extract_done": False,
                "eval_done": False,
            },
        },
        has_rl_process=True,
        has_bank_builder_process=False,
    )

    assert phase == "training"


def test_monitor_process_summary_filters_wrapper_processes():
    mon = _load_rl_monitor_module()

    summary = mon.collect_process_summary(
        process_rows=[
            {
                "pid": 1,
                "ppid": 0,
                "stat": "Ss",
                "etime": "00:10",
                "cpu": 0.0,
                "mem": 0.0,
                "cmd": "tmux new-session -d -s rl_train_resume bash -lc 'python -u test_memento_multireaction.py ...'",
            },
            {
                "pid": 2,
                "ppid": 1,
                "stat": "Ss",
                "etime": "00:09",
                "cpu": 0.0,
                "mem": 0.0,
                "cmd": "bash -lc cd /repo && python -u test_memento_multireaction.py ...",
            },
            {
                "pid": 3,
                "ppid": 2,
                "stat": "Sl",
                "etime": "00:08",
                "cpu": 10.0,
                "mem": 1.0,
                "cmd": "python -u test_memento_multireaction.py --output-base-dir output/run",
            },
            {
                "pid": 4,
                "ppid": 3,
                "stat": "Sl",
                "etime": "00:03",
                "cpu": 0.2,
                "mem": 0.3,
                "cmd": "/opt/conda/envs/catdt/bin/python /srv/catdt/scripts/build_agent45_banks.py --casebank x",
            },
        ]
    )

    assert [row["pid"] for row in summary["rl_processes"]] == [3]
    assert [row["pid"] for row in summary["bank_builders"]] == [4]
