"""
Workflow persistence: checkpoint system + reporting tools.

Provides:
1. Atomic checkpoint save/load for workflow state
2. Resumable workflow mixin for checkpoint-based restart
3. Step structure export and HTML report generation
"""

import hashlib
import importlib.util
import json
import logging
import os
import pickle
import random
import re
import shutil
import tempfile
import time
from collections import deque
from contextlib import contextmanager
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms
from ase.io import read, write

from .common import logger


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class StepArtifact:
    """步骤产物"""
    name: str
    file_path: str
    file_type: str  # "pickle", "vasp", "json", etc.
    description: str = ""

    def load(self) -> Any:
        """加载产物"""
        path = Path(self.file_path)
        if not path.exists():
            raise FileNotFoundError(f"Artifact not found: {self.file_path}")

        if self.file_type == "pickle":
            with open(path, 'rb') as f:
                return pickle.load(f)
        elif self.file_type == "json":
            with open(path, 'r') as f:
                return json.load(f)
        elif self.file_type == "vasp":
            return read(path)
        else:
            raise ValueError(f"Unknown file type: {self.file_type}")


@dataclass
class WorkflowCheckpoint:
    """工作流检查点"""
    run_id: str
    step_name: str  # 已完成的步骤名称
    step_number: int  # 步骤序号
    timestamp: datetime
    state_data: Dict[str, Any]  # 序列化的状态数据
    artifacts: List[StepArtifact]  # 步骤产物
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        return {
            "run_id": self.run_id,
            "step_name": self.step_name,
            "step_number": self.step_number,
            "timestamp": self.timestamp.isoformat(),
            "state_data": self.state_data,
            "artifacts": [
                {"name": a.name, "path": a.file_path, "type": a.file_type, "desc": a.description}
                for a in self.artifacts
            ],
            "metadata": self.metadata
        }

    @classmethod
    def from_dict(cls, data: Dict) -> "WorkflowCheckpoint":
        """从字典恢复"""
        artifacts = [
            StepArtifact(
                name=a["name"],
                file_path=a["path"],
                file_type=a["type"],
                description=a.get("desc", "")
            )
            for a in data.get("artifacts", [])
        ]

        return cls(
            run_id=data["run_id"],
            step_name=data["step_name"],
            step_number=data["step_number"],
            timestamp=datetime.fromisoformat(data["timestamp"]),
            state_data=data["state_data"],
            artifacts=artifacts,
            metadata=data.get("metadata", {})
        )


# ---------------------------------------------------------------------------
# Checkpoint manager
# ---------------------------------------------------------------------------


class CheckpointManager:
    """
    检查点管理器

    管理检查点的保存、加载和清理
    """

    def __init__(self, checkpoint_dir: str = ".catdt_checkpoints"):
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self._current_checkpoint: Optional[WorkflowCheckpoint] = None

    def save_checkpoint(
        self,
        run_id: str,
        step_name: str,
        step_number: int,
        state: Any,  # WorkflowState 对象
        artifacts: Dict[str, Any] = None,
        metadata: Dict[str, Any] = None
    ) -> Path:
        """
        原子性保存检查点

        先写入临时文件，再重命名，确保不会出现损坏的检查点
        """
        timestamp = datetime.now()

        # 序列化状态
        state_data = self._serialize_state(state)

        # 保存产物
        artifact_list = []
        run_checkpoint_dir = self.checkpoint_dir / run_id
        run_checkpoint_dir.mkdir(parents=True, exist_ok=True)

        if artifacts:
            for name, data in artifacts.items():
                artifact_path = self._save_artifact(
                    run_checkpoint_dir, name, data, step_number
                )
                artifact_list.append(artifact_path)

        # 创建检查点
        checkpoint = WorkflowCheckpoint(
            run_id=run_id,
            step_name=step_name,
            step_number=step_number,
            timestamp=timestamp,
            state_data=state_data,
            artifacts=artifact_list,
            metadata=metadata or {}
        )

        # 原子性保存
        checkpoint_file = run_checkpoint_dir / f"checkpoint_{step_number:02d}_{step_name}.pkl"
        temp_file = checkpoint_file.with_suffix(".tmp")

        try:
            with open(temp_file, 'wb') as f:
                pickle.dump(checkpoint, f)

            # 原子性重命名
            temp_file.rename(checkpoint_file)

            # 同时保存 JSON 版本便于查看
            json_file = checkpoint_file.with_suffix(".json")
            with open(json_file, 'w') as f:
                json.dump(checkpoint.to_dict(), f, indent=2, default=str)

            self._current_checkpoint = checkpoint
            logger.info(f"Checkpoint saved: {checkpoint_file}")

            return checkpoint_file

        except Exception as e:
            # 清理临时文件
            if temp_file.exists():
                temp_file.unlink()
            raise RuntimeError(f"Failed to save checkpoint: {e}")

    def load_latest_checkpoint(self, run_id: str) -> Optional[WorkflowCheckpoint]:
        """加载最新的检查点"""
        run_dir = self.checkpoint_dir / run_id
        if not run_dir.exists():
            return None

        # 找到最新的检查点文件
        checkpoint_files = sorted(
            run_dir.glob("checkpoint_*.pkl"),
            key=lambda p: p.stat().st_mtime,
            reverse=True
        )

        if not checkpoint_files:
            return None

        latest_file = checkpoint_files[0]

        try:
            with open(latest_file, 'rb') as f:
                checkpoint = pickle.load(f)

            logger.info(f"Loaded checkpoint: {latest_file} (step {checkpoint.step_name})")
            self._current_checkpoint = checkpoint
            return checkpoint

        except Exception as e:
            logger.error(f"Failed to load checkpoint {latest_file}: {e}")
            return None

    def load_checkpoint_by_step(self, run_id: str, step_name: str) -> Optional[WorkflowCheckpoint]:
        """按步骤名称加载检查点"""
        run_dir = self.checkpoint_dir / run_id
        if not run_dir.exists():
            return None

        matches = list(run_dir.glob(f"checkpoint_*_{step_name}.pkl"))

        if not matches:
            return None

        with open(matches[0], 'rb') as f:
            return pickle.load(f)

    def list_checkpoints(self, run_id: str) -> List[Dict]:
        """列出所有检查点"""
        run_dir = self.checkpoint_dir / run_id
        if not run_dir.exists():
            return []

        checkpoints = []
        for json_file in sorted(run_dir.glob("checkpoint_*.json")):
            with open(json_file, 'r') as f:
                data = json.load(f)
                checkpoints.append({
                    "step": data["step_name"],
                    "number": data["step_number"],
                    "timestamp": data["timestamp"],
                    "artifacts": len(data.get("artifacts", []))
                })

        return checkpoints

    def clean_old_checkpoints(self, run_id: str, keep_last: int = 3):
        """清理旧检查点，只保留最近的几个"""
        run_dir = self.checkpoint_dir / run_id
        if not run_dir.exists():
            return

        checkpoint_files = sorted(
            run_dir.glob("checkpoint_*.pkl"),
            key=lambda p: p.stat().st_mtime,
            reverse=True
        )

        for old_file in checkpoint_files[keep_last:]:
            old_file.unlink()
            # 同时删除对应的 JSON 文件
            json_file = old_file.with_suffix(".json")
            if json_file.exists():
                json_file.unlink()
            logger.debug(f"Cleaned old checkpoint: {old_file}")

    def _serialize_state(self, state: Any) -> Dict[str, Any]:
        """序列化状态对象"""
        if hasattr(state, '__dict__'):
            return {
                "_type": type(state).__name__,
                "_module": type(state).__module__,
                "data": {
                    k: self._serialize_value(v)
                    for k, v in state.__dict__.items()
                }
            }
        return {"_raw": self._serialize_value(state)}

    def _serialize_value(self, value: Any) -> Any:
        """序列化单个值"""
        if isinstance(value, Atoms):
            # ASE Atoms 特殊处理
            return {
                "_type": "Atoms",
                "symbols": value.get_chemical_symbols(),
                "positions": value.get_positions().tolist(),
                "cell": value.cell.tolist() if value.cell else None,
                "pbc": value.pbc.tolist() if hasattr(value.pbc, 'tolist') else value.pbc,
                "info": {k: v for k, v in value.info.items() if self._is_serializable(v)}
            }
        elif isinstance(value, (list, tuple)):
            return [self._serialize_value(v) for v in value]
        elif isinstance(value, dict):
            return {k: self._serialize_value(v) for k, v in value.items()}
        elif isinstance(value, (int, float, str, bool, type(None))):
            return value
        else:
            # 尝试 pickle
            try:
                pickle.dumps(value)
                return {"_pickled": pickle.dumps(value).hex()}
            except:
                return {"_str": str(value)}

    def _is_serializable(self, value: Any) -> bool:
        """检查值是否可序列化"""
        try:
            json.dumps(value)
            return True
        except:
            return False

    def _save_artifact(
        self,
        run_dir: Path,
        name: str,
        data: Any,
        step_number: int
    ) -> StepArtifact:
        """保存产物"""
        artifact_dir = run_dir / f"step_{step_number:02d}_artifacts"
        artifact_dir.mkdir(exist_ok=True)

        # 根据数据类型选择保存方式
        if isinstance(data, Atoms):
            file_path = artifact_dir / f"{name}.vasp"
            write(file_path, data)
            return StepArtifact(name, str(file_path), "vasp", f"ASE Atoms structure")

        elif isinstance(data, (dict, list)):
            file_path = artifact_dir / f"{name}.json"
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2, default=str)
            return StepArtifact(name, str(file_path), "json", f"JSON data")

        else:
            # 默认 pickle
            file_path = artifact_dir / f"{name}.pkl"
            with open(file_path, 'wb') as f:
                pickle.dump(data, f)
            return StepArtifact(name, str(file_path), "pickle", f"Pickled object")


# ---------------------------------------------------------------------------
# Resumable workflow mixin
# ---------------------------------------------------------------------------


class ResumableWorkflowMixin:
    """
    可恢复工作流的 Mixin 类

    添加到工作流类中以支持断点续传
    """

    def __init__(self, *args, enable_checkpointing: bool = True, **kwargs):
        super().__init__(*args, **kwargs)
        self.enable_checkpointing = enable_checkpointing
        self.checkpoint_manager = CheckpointManager() if enable_checkpointing else None
        self._current_step = 0
        self._step_names = []

    def run_with_resume(self, config: Any, start_from: Optional[str] = None) -> Any:
        """
        支持断点续传的运行

        Parameters
        ----------
        config : Any
            工作流配置
        start_from : str, optional
            从指定步骤开始（如果提供，则忽略检查点）
        """
        run_id = getattr(config, 'run_id', 'default_run')

        # 确定起始步骤
        if start_from:
            start_step = start_from
            state = None
        else:
            checkpoint = self.checkpoint_manager.load_latest_checkpoint(run_id)
            if checkpoint:
                start_step = self._get_next_step(checkpoint.step_name)
                state = self._restore_state_from_checkpoint(checkpoint)
                logger.info(f"Resuming from step: {start_step}")
            else:
                start_step = self._step_names[0] if self._step_names else "start"
                state = None

        # 执行剩余步骤
        remaining_steps = self._get_steps_from(start_step)

        for step_name in remaining_steps:
            logger.info(f"Executing step: {step_name}")

            try:
                # 执行步骤
                state = self._execute_step(step_name, state, config)

                # 保存检查点
                if self.enable_checkpointing:
                    self._save_step_checkpoint(run_id, step_name, state)

            except Exception as e:
                logger.error(f"Step {step_name} failed: {e}")
                raise

        return state

    def _save_step_checkpoint(self, run_id: str, step_name: str, state: Any):
        """保存步骤检查点"""
        self._current_step += 1

        artifacts = self._collect_step_artifacts(step_name, state)
        metadata = {
            "step_duration": None,  # 可以添加计时
            "step_status": "completed"
        }

        self.checkpoint_manager.save_checkpoint(
            run_id=run_id,
            step_name=step_name,
            step_number=self._current_step,
            state=state,
            artifacts=artifacts,
            metadata=metadata
        )

    def _get_next_step(self, completed_step: str) -> str:
        """获取下一步"""
        if not self._step_names:
            return "end"

        try:
            idx = self._step_names.index(completed_step)
            if idx + 1 < len(self._step_names):
                return self._step_names[idx + 1]
            return "end"
        except ValueError:
            return self._step_names[0]

    def _get_steps_from(self, start_step: str) -> List[str]:
        """获取从某步骤开始的所有步骤"""
        if start_step == "end":
            return []

        try:
            idx = self._step_names.index(start_step)
            return self._step_names[idx:]
        except ValueError:
            return self._step_names

    def _restore_state_from_checkpoint(self, checkpoint: WorkflowCheckpoint) -> Any:
        """从检查点恢复状态 - 子类应重写此方法"""
        raise NotImplementedError("Subclasses must implement _restore_state_from_checkpoint")

    def _execute_step(self, step_name: str, state: Any, config: Any) -> Any:
        """执行单个步骤 - 子类应重写此方法"""
        raise NotImplementedError("Subclasses must implement _execute_step")

    def _collect_step_artifacts(self, step_name: str, state: Any) -> Dict[str, Any]:
        """收集步骤产物 - 子类可重写此方法"""
        return {}


def with_checkpoint(step_name: str, step_number: int):
    """
    装饰器：自动为函数保存检查点

    Example:
        >>> @with_checkpoint("surff", 1)
        ... def run_surff(state, config):
        ...     # 执行步骤
        ...     return new_state
    """
    def decorator(func: Callable):
        def wrapper(self, *args, **kwargs):
            # 执行函数
            result = func(self, *args, **kwargs)

            # 保存检查点
            if hasattr(self, 'checkpoint_manager') and self.checkpoint_manager:
                run_id = getattr(kwargs.get('config'), 'run_id', 'default_run')
                self.checkpoint_manager.save_checkpoint(
                    run_id=run_id,
                    step_name=step_name,
                    step_number=step_number,
                    state=result,
                    metadata={"function": func.__name__}
                )

            return result
        return wrapper
    return decorator


# ---------------------------------------------------------------------------
# Agent4/5 Memento-style case memory (zero-parameter, memory-only tuning)
# ---------------------------------------------------------------------------


class Agent45MementoToolsMixin:
    """
    Memento-style memory utilities for Agent4/5.

    Design goals:
    1) Zero-parameter adaptation (no model weight updates).
    2) Persist successful/failed trajectories as retrievable cases.
    3) Provide positive/negative in-context examples to improve next iterations.
    """

    def _agent45_memento_casebank_path(self) -> Path:
        path = self.output_base_dir / "agent45_memento_casebank.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _agent45_knowledge_bank_path(self) -> Path:
        override = str(os.getenv("CATDT_AGENT45_KNOWLEDGE_BANK_PATH", "")).strip()
        if override:
            return Path(override).expanduser()
        path = self.output_base_dir / "agent45_knowledge_bank.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _agent45_skill_bank_path(self) -> Path:
        override = str(os.getenv("CATDT_AGENT45_SKILL_BANK_PATH", "")).strip()
        if override:
            return Path(override).expanduser()
        path = self.output_base_dir / "agent45_skill_bank.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def _agent45_memento_env_flag(name: str, default: str = "0") -> bool:
        value = str(os.getenv(name, default)).strip().lower()
        return value in {"1", "true", "yes", "on"}

    def _load_agent45_memento_cases(self) -> List[Dict[str, Any]]:
        path = self._agent45_memento_casebank_path()
        return self._load_agent45_jsonl_items(path=path)

    def _load_agent45_jsonl_items(self, path: Path) -> List[Dict[str, Any]]:
        if not path.exists():
            return []

        items: List[Dict[str, Any]] = []
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                raw = line.strip()
                if not raw:
                    continue
                try:
                    obj = json.loads(raw)
                except Exception:
                    continue
                if isinstance(obj, dict):
                    items.append(obj)
        return items

    def _load_agent45_knowledge_items(self) -> List[Dict[str, Any]]:
        return self._load_agent45_jsonl_items(path=self._agent45_knowledge_bank_path())

    def _load_agent45_skill_items(self) -> List[Dict[str, Any]]:
        return self._load_agent45_jsonl_items(path=self._agent45_skill_bank_path())

    def _memento_extract_pairs(
        self,
        items: List[Dict[str, Any]],
        key_field: str = "question",
        value_field: str = "plan",
    ) -> List[Tuple[str, Any, int]]:
        """
        Reuse Memento's pair extraction if available; otherwise use compatible fallback.
        """
        np_memory_path = Path(__file__).resolve().parents[2] / "deps" / "Memento" / "memory" / "np_memory.py"
        if np_memory_path.exists():
            try:
                spec = importlib.util.spec_from_file_location("memento_np_memory", str(np_memory_path))
                if spec and spec.loader:
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                    extract_pairs = getattr(module, "extract_pairs", None)
                    if callable(extract_pairs):
                        return list(extract_pairs(items, key_field, value_field))
            except Exception:
                pass

        pairs: List[Tuple[str, Any, int]] = []
        for idx, obj in enumerate(items):
            if not isinstance(obj, dict):
                continue
            if key_field in obj and value_field in obj:
                pairs.append((str(obj.get(key_field, "")), obj.get(value_field), idx))
        return pairs

    @staticmethod
    def _tokenize_case_text(text: str) -> List[str]:
        cleaned = str(text or "").strip().lower()
        if not cleaned:
            return []
        return re.findall(r"[a-z0-9\*\+\-\>\(\)_\/]+", cleaned)

    @staticmethod
    def _extract_transition_tokens(text: str) -> List[str]:
        raw = str(text or "")
        transitions = re.findall(r"(\*?[A-Za-z0-9\(\)\+\-]+)\s*(?:->|→)\s*(\*?[A-Za-z0-9\(\)\+\-]+)", raw)
        return [f"{lhs}->{rhs}".lower() for lhs, rhs in transitions]

    @staticmethod
    def _agent45_memento_mode() -> str:
        mode = str(os.getenv("CATDT_AGENT45_MEMENTO_MODE", "parametric")).strip().lower()
        if mode in {"parametric", "non_parametric", "non-parametric"}:
            return "parametric" if mode == "parametric" else "non_parametric"
        return "parametric"

    def _agent45_parametric_model_name(self) -> str:
        return str(
            os.getenv(
                "CATDT_AGENT45_PARAMETRIC_MODEL_NAME",
                "sentence-transformers/all-MiniLM-L6-v2",
            )
        ).strip()

    def _agent45_parametric_model_path(self) -> Path:
        override = str(os.getenv("CATDT_AGENT45_PARAMETRIC_MODEL_PATH", "")).strip()
        if override:
            return Path(override).expanduser()
        return self.output_base_dir / "agent45_memento_retriever" / "best.pt"

    @staticmethod
    def _agent45_parametric_device() -> str:
        return str(os.getenv("CATDT_AGENT45_PARAMETRIC_DEVICE", "cuda")).strip().lower()

    @staticmethod
    def _agent45_case_query_text(case: Dict[str, Any]) -> str:
        text = str(case.get("question", "")).strip()
        if text:
            return text
        reaction_type = str(case.get("reaction_type", "")).strip()
        transitions = list(case.get("transition_signature", []) or [])
        intermediates = list(case.get("intermediates", []) or [])
        return (
            f"reaction_type={reaction_type}; "
            f"intermediates={intermediates}; transitions={transitions}"
        )

    @staticmethod
    def _agent45_case_plan_text(case: Dict[str, Any]) -> str:
        raw = case.get("plan")
        if isinstance(raw, (dict, list)):
            try:
                return json.dumps(raw, ensure_ascii=False)
            except Exception:
                return str(raw)
        if isinstance(raw, str) and raw.strip():
            return raw.strip()
        payload = {
            "status": str(case.get("validation_status", "")),
            "issues": list(case.get("issues", []) or [])[:8],
            "feedback": str(case.get("feedback", "")),
            "design_outline": list(case.get("design_outline", []) or [])[:12],
        }
        try:
            return json.dumps(payload, ensure_ascii=False)
        except Exception:
            return str(payload)

    def _agent45_case_to_icl_text(self, case: Dict[str, Any]) -> str:
        case_text = self._agent45_case_query_text(case)
        plan_text = self._agent45_case_plan_text(case)
        parts = [f"[CASE]\n{case_text}\n[PLAN]\n{plan_text}"]

        # Add surface context for parametric retriever training
        surf = case.get("surface_snapshot", {})
        if isinstance(surf, dict) and surf.get("cell"):
            top2 = surf.get("top2_layers", [])
            from collections import Counter
            elem_counts = dict(Counter(a.get("elem", "?") for a in top2)) if top2 else {}
            parts.append(
                f"[SURFACE] cell={surf['cell']} top_z={surf.get('surface_top_z')} "
                f"composition={elem_counts}"
            )

        # Add per-step geometry (compact) for retriever training
        geo_parts = []
        for snap in list(case.get("step_geometry_snapshot", []) or []):
            name = snap.get("name", "")
            for ep in ("reactant", "product"):
                ep_data = snap.get(ep, {})
                if not isinstance(ep_data, dict):
                    continue
                staged = ep_data.get("staged", [])
                if staged:
                    staged_str = " ".join(
                        f"{a.get('elem','')}_d={a.get('dist_to_nearest_ads','?')}"
                        for a in staged
                    )
                    geo_parts.append(f"{name}:{ep} staged={staged_str}")
        if geo_parts:
            parts.append(f"[GEOMETRY] {'; '.join(geo_parts[:8])}")

        return "\n".join(parts).strip()

    def _agent45_parametric_cache(self) -> Dict[str, Any]:
        cache = getattr(self, "_agent45_parametric_retriever_cache", None)
        if not isinstance(cache, dict):
            cache = {}
            setattr(self, "_agent45_parametric_retriever_cache", cache)
        return cache

    def _load_agent45_parametric_retriever(
        self,
        model_path: Optional[Path] = None,
        model_name: Optional[str] = None,
        device: Optional[str] = None,
    ) -> Dict[str, Any]:
        model_path = Path(model_path or self._agent45_parametric_model_path()).expanduser()
        model_name = str(model_name or self._agent45_parametric_model_name()).strip()
        device = str(device or self._agent45_parametric_device()).strip().lower()

        if device == "cuda":
            import torch

            if not torch.cuda.is_available():
                raise RuntimeError("Parametric CBR requires CUDA, but torch.cuda.is_available() is False.")

        if not model_path.exists():
            raise FileNotFoundError(
                f"Parametric retriever checkpoint not found: {model_path}. "
                "Train with `train_agent45_parametric_retriever` first."
            )

        mtime = float(model_path.stat().st_mtime)
        cache = self._agent45_parametric_cache()
        cache_key = (
            str(model_path.resolve()),
            model_name,
            device,
            mtime,
        )
        if cache.get("cache_key") == cache_key and cache.get("model") is not None:
            return cache

        import torch
        from transformers import AutoModel, AutoTokenizer

        module_path = Path(__file__).resolve().parents[2] / "deps" / "Memento" / "memory" / "train_memory_retriever.py"
        if not module_path.exists():
            raise FileNotFoundError(f"Memento train module missing: {module_path}")

        spec = importlib.util.spec_from_file_location("memento_train_memory_retriever", str(module_path))
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Failed to load module spec: {module_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        model_cls = getattr(module, "MemoryRetrieverClassifier", None)
        if model_cls is None:
            raise RuntimeError("MemoryRetrieverClassifier not found in Memento train module.")

        tokenizer = AutoTokenizer.from_pretrained(pretrained_model_name_or_path=model_name, use_fast=True)
        backbone = AutoModel.from_pretrained(pretrained_model_name_or_path=model_name, use_safetensors=True)
        model = model_cls(backbone).to(device)
        state_dict = torch.load(str(model_path), map_location=device)
        model.load_state_dict(state_dict)
        model.eval()

        cache.clear()
        cache.update(
            {
                "cache_key": cache_key,
                "tokenizer": tokenizer,
                "model": model,
                "device": device,
                "model_path": str(model_path),
                "model_name": model_name,
            }
        )
        return cache

    def _score_agent45_cases_parametric(
        self,
        query_text: str,
        cases: List[Dict[str, Any]],
        model_path: Optional[Path] = None,
        model_name: Optional[str] = None,
        device: Optional[str] = None,
    ) -> List[float]:
        if not cases:
            return []
        cache = self._load_agent45_parametric_retriever(
            model_path=model_path,
            model_name=model_name,
            device=device,
        )
        tokenizer = cache["tokenizer"]
        model = cache["model"]
        target_device = str(cache["device"])

        import torch

        case_texts = [self._agent45_case_to_icl_text(case) for case in cases]
        natural_texts = [str(query_text)] * len(case_texts)

        batch_size = max(1, int(os.getenv("CATDT_AGENT45_PARAMETRIC_BATCH_SIZE", "32")))
        max_len = max(64, int(os.getenv("CATDT_AGENT45_PARAMETRIC_MAX_LEN", "256")))

        probs: List[float] = []
        with torch.inference_mode():
            for start in range(0, len(case_texts), batch_size):
                end = min(start + batch_size, len(case_texts))
                t1 = tokenizer(
                    case_texts[start:end],
                    padding=True,
                    truncation=True,
                    max_length=max_len,
                    return_tensors="pt",
                )
                t2 = tokenizer(
                    natural_texts[start:end],
                    padding=True,
                    truncation=True,
                    max_length=max_len,
                    return_tensors="pt",
                )
                ids1 = t1["input_ids"].to(target_device)
                mask1 = t1["attention_mask"].to(target_device)
                ids2 = t2["input_ids"].to(target_device)
                mask2 = t2["attention_mask"].to(target_device)
                logits = model(ids1, mask1, ids2, mask2)
                batch_probs = torch.softmax(logits, dim=1)[:, 1].detach().cpu().tolist()
                probs.extend(float(x) for x in batch_probs)
        return probs

    @staticmethod
    def _extract_reaction_type_from_query(query_text: str) -> str:
        match = re.search(r"reaction_type=(\w+)", str(query_text or ""))
        return match.group(1).upper() if match else ""

    def _score_agent45_case(self, query_text: str, case: Dict[str, Any]) -> float:
        query_tokens = set(self._tokenize_case_text(query_text))
        case_tokens = set(
            self._tokenize_case_text(case.get("question", ""))
            + self._tokenize_case_text(case.get("reaction_description", ""))
            + self._tokenize_case_text(case.get("feedback", ""))
        )

        token_score = 0.0
        if query_tokens and case_tokens:
            token_score = len(query_tokens & case_tokens) / max(len(query_tokens | case_tokens), 1)

        query_transitions = set(self._extract_transition_tokens(query_text))
        case_transitions = set(
            list(case.get("transition_signature", []) or [])
            + self._extract_transition_tokens(case.get("question", ""))
        )
        trans_score = 0.0
        if query_transitions and case_transitions:
            trans_score = len(query_transitions & case_transitions) / max(len(query_transitions), 1)

        reward = float(case.get("reward", 0.0) or 0.0)
        reward_bonus = 0.05 if reward > 0 else 0.0

        # Reaction type matching: boost same-type cases, penalize cross-type
        query_rxn = self._extract_reaction_type_from_query(query_text)
        case_rxn = str(case.get("reaction_type", "")).strip().upper()
        rxn_bonus = 0.0
        if query_rxn and case_rxn:
            rxn_bonus = 0.15 if query_rxn == case_rxn else -0.10

        score = 0.50 * token_score + 0.30 * trans_score + reward_bonus + rxn_bonus
        return float(max(0.0, score))

    def _score_agent45_knowledge_or_skill_item(self, query_text: str, item: Dict[str, Any]) -> float:
        query_tokens = set(self._tokenize_case_text(query_text))
        item_tokens = set(
            self._tokenize_case_text(item.get("title", ""))
            + self._tokenize_case_text(item.get("knowledge", ""))
            + self._tokenize_case_text(item.get("skill", ""))
            + self._tokenize_case_text(item.get("guidance", ""))
            + self._tokenize_case_text(item.get("action_template", ""))
            + self._tokenize_case_text(item.get("applicability", ""))
            + self._tokenize_case_text(item.get("reaction_type", ""))
        )
        token_score = 0.0
        if query_tokens and item_tokens:
            token_score = len(query_tokens & item_tokens) / max(len(query_tokens | item_tokens), 1)

        query_transitions = set(self._extract_transition_tokens(query_text))
        item_transitions = set(
            list(item.get("transition_signature", []) or [])
            + list(item.get("transitions", []) or [])
            + self._extract_transition_tokens(item.get("knowledge", ""))
            + self._extract_transition_tokens(item.get("skill", ""))
        )
        trans_score = 0.0
        if query_transitions and item_transitions:
            trans_score = len(query_transitions & item_transitions) / max(len(query_transitions), 1)

        confidence = float(item.get("confidence", item.get("score", 0.0)) or 0.0)
        confidence_bonus = 0.05 * max(0.0, min(confidence, 1.0))

        # Reaction type matching: boost same-type items, penalize cross-type
        query_rxn = self._extract_reaction_type_from_query(query_text)
        item_rxn = str(item.get("reaction_type", "")).strip().upper()
        rxn_bonus = 0.0
        if query_rxn and item_rxn:
            rxn_bonus = 0.12 if query_rxn == item_rxn else -0.08

        score = 0.55 * token_score + 0.25 * trans_score + confidence_bonus + rxn_bonus
        return float(max(0.0, score))

    @staticmethod
    def _agent45_item_short_text(item: Dict[str, Any], keys: List[str], max_len: int = 220) -> str:
        for key in keys:
            value = str(item.get(key, "") or "").strip().replace("\n", " ")
            if value:
                if len(value) > max_len:
                    return value[:max_len] + "..."
                return value
        return ""

    def _build_agent45_generic_bank_prompt_block(
        self,
        bank_name: str,
        query_text: str,
        items: List[Dict[str, Any]],
        max_items: int = 4,
    ) -> str:
        lines: List[str] = []
        lines.append(f"【{bank_name} Top-K】")
        lines.append(f"query={str(query_text or '').strip()}")
        if not items:
            lines.append("(none)")
            return "\n".join(lines)

        for idx, item in enumerate(items[: max(1, int(max_items))], start=1):
            title = self._agent45_item_short_text(item, ["title", "name", "id"], max_len=80) or f"{bank_name}_item_{idx}"
            evidence = self._agent45_item_short_text(
                item,
                ["guidance", "knowledge", "skill", "action_template", "applicability", "notes"],
                max_len=220,
            )
            transitions = list(item.get("transition_signature", []) or [])
            if not transitions:
                transitions = list(item.get("transitions", []) or [])
            lines.append(
                f"- #{idx} score={float(item.get('retrieval_score', 0.0)):.3f}; title={title}; "
                f"reaction_type={item.get('reaction_type', '')}; transitions={transitions[:4]}"
            )
            if evidence:
                lines.append(f"  guidance={evidence}")
        return "\n".join(lines)

    def retrieve_agent45_knowledge_items(
        self,
        query_text: str,
        top_k: int = 4,
        min_score: float = 0.03,
    ) -> Dict[str, Any]:
        if self._agent45_memento_env_flag("CATDT_AGENT45_MEMENTO_DISABLE_READ", "0"):
            return {
                "query": query_text,
                "items": [],
                "prompt_block": "(knowledge read disabled by CATDT_AGENT45_MEMENTO_DISABLE_READ)",
                "bank_path": str(self._agent45_knowledge_bank_path()),
            }

        items = self._load_agent45_knowledge_items()
        if not items:
            return {
                "query": query_text,
                "items": [],
                "prompt_block": "(no knowledge items)",
                "bank_path": str(self._agent45_knowledge_bank_path()),
            }

        scored: List[Dict[str, Any]] = []
        for item in items:
            score = self._score_agent45_knowledge_or_skill_item(query_text=query_text, item=item)
            if score < float(min_score):
                continue
            obj = dict(item)
            obj["retrieval_score"] = round(float(score), 6)
            scored.append(obj)
        scored.sort(key=lambda x: float(x.get("retrieval_score", 0.0)), reverse=True)
        selected = scored[: max(1, int(top_k))]
        return {
            "query": query_text,
            "items": selected,
            "prompt_block": self._build_agent45_generic_bank_prompt_block(
                bank_name="KnowledgeBank",
                query_text=query_text,
                items=selected,
                max_items=max(1, int(top_k)),
            ),
            "bank_path": str(self._agent45_knowledge_bank_path()),
        }

    def retrieve_agent45_skill_items(
        self,
        query_text: str,
        top_k: int = 4,
        min_score: float = 0.03,
    ) -> Dict[str, Any]:
        if self._agent45_memento_env_flag("CATDT_AGENT45_MEMENTO_DISABLE_READ", "0"):
            return {
                "query": query_text,
                "items": [],
                "prompt_block": "(skill read disabled by CATDT_AGENT45_MEMENTO_DISABLE_READ)",
                "bank_path": str(self._agent45_skill_bank_path()),
            }

        items = self._load_agent45_skill_items()
        if not items:
            return {
                "query": query_text,
                "items": [],
                "prompt_block": "(no skill items)",
                "bank_path": str(self._agent45_skill_bank_path()),
            }

        scored: List[Dict[str, Any]] = []
        for item in items:
            score = self._score_agent45_knowledge_or_skill_item(query_text=query_text, item=item)
            if score < float(min_score):
                continue
            obj = dict(item)
            obj["retrieval_score"] = round(float(score), 6)
            scored.append(obj)
        scored.sort(key=lambda x: float(x.get("retrieval_score", 0.0)), reverse=True)
        selected = scored[: max(1, int(top_k))]
        return {
            "query": query_text,
            "items": selected,
            "prompt_block": self._build_agent45_generic_bank_prompt_block(
                bank_name="SkillBank",
                query_text=query_text,
                items=selected,
                max_items=max(1, int(top_k)),
            ),
            "bank_path": str(self._agent45_skill_bank_path()),
        }

    def _retrieve_agent45_memento_cases_non_parametric(
        self,
        query_text: str,
        top_k: int = 6,
        include_negative: bool = True,
        min_score: float = 0.05,
    ) -> Dict[str, Any]:
        cases = self._load_agent45_memento_cases()
        if not cases:
            return {
                "query": query_text,
                "cases": [],
                "positive_cases": [],
                "negative_cases": [],
                "prompt_block": "(no memento cases)",
                "retrieval_mode": "non_parametric",
            }

        _ = self._memento_extract_pairs(cases, key_field="question", value_field="plan")

        scored: List[Dict[str, Any]] = []
        for case in cases:
            score = self._score_agent45_case(query_text=query_text, case=case)
            if score < float(min_score):
                continue
            label = str(case.get("case_label", "")).strip().lower()
            if not label:
                label = "positive" if float(case.get("reward", 0.0) or 0.0) > 0 else "negative"
            if (not include_negative) and label != "positive":
                continue
            item = dict(case)
            item["retrieval_score"] = round(float(score), 6)
            item["case_label"] = label
            scored.append(item)

        scored.sort(key=lambda x: float(x.get("retrieval_score", 0.0)), reverse=True)
        selected = scored[: max(int(top_k), 0)]

        positive = [c for c in selected if str(c.get("case_label", "")).lower() == "positive"]
        negative = [c for c in selected if str(c.get("case_label", "")).lower() != "positive"]

        prompt_block = self.build_agent45_memento_prompt_block(
            query_text=query_text,
            positive_cases=positive,
            negative_cases=negative,
            max_cases_per_group=max(1, min(3, int(top_k))),
        )

        return {
            "query": query_text,
            "cases": selected,
            "positive_cases": positive,
            "negative_cases": negative,
            "prompt_block": prompt_block,
            "retrieval_mode": "non_parametric",
        }

    def _retrieve_agent45_memento_cases_parametric(
        self,
        query_text: str,
        top_k: int = 6,
        include_negative: bool = True,
        min_score: float = 0.05,
    ) -> Dict[str, Any]:
        cases = self._load_agent45_memento_cases()
        if not cases:
            return {
                "query": query_text,
                "cases": [],
                "positive_cases": [],
                "negative_cases": [],
                "prompt_block": "(no memento cases)",
                "retrieval_mode": "parametric",
            }

        model_path = self._agent45_parametric_model_path().expanduser()
        if not model_path.exists():
            return {
                "query": query_text,
                "cases": [],
                "positive_cases": [],
                "negative_cases": [],
                "prompt_block": "(parametric retriever not trained yet)",
                "retrieval_mode": "parametric",
                "model_path": str(model_path),
                "model_name": self._agent45_parametric_model_name(),
            }

        probs = self._score_agent45_cases_parametric(query_text=query_text, cases=cases)
        scored: List[Dict[str, Any]] = []
        for case, prob in zip(cases, probs):
            score = float(prob)
            if score < float(min_score):
                continue
            label = str(case.get("case_label", "")).strip().lower()
            if not label:
                label = "positive" if float(case.get("reward", 0.0) or 0.0) > 0 else "negative"
            if (not include_negative) and label != "positive":
                continue
            item = dict(case)
            item["retrieval_score"] = round(score, 6)
            item["case_label"] = label
            scored.append(item)

        scored.sort(key=lambda x: float(x.get("retrieval_score", 0.0)), reverse=True)
        selected = scored[: max(int(top_k), 0)]
        positive = [c for c in selected if str(c.get("case_label", "")).lower() == "positive"]
        negative = [c for c in selected if str(c.get("case_label", "")).lower() != "positive"]
        prompt_block = self.build_agent45_memento_prompt_block(
            query_text=query_text,
            positive_cases=positive,
            negative_cases=negative,
            max_cases_per_group=max(1, min(3, int(top_k))),
        )
        return {
            "query": query_text,
            "cases": selected,
            "positive_cases": positive,
            "negative_cases": negative,
            "prompt_block": prompt_block,
            "retrieval_mode": "parametric",
            "model_path": str(model_path),
            "model_name": self._agent45_parametric_model_name(),
        }

    def retrieve_agent45_memento_cases(
        self,
        query_text: str,
        top_k: int = 6,
        include_negative: bool = True,
        min_score: float = 0.05,
    ) -> Dict[str, Any]:
        """
        Retrieve Memento-like cases for Agent4/5 prompt conditioning.
        """
        if self._agent45_memento_env_flag("CATDT_AGENT45_MEMENTO_DISABLE_READ", "0"):
            return {
                "query": query_text,
                "cases": [],
                "positive_cases": [],
                "negative_cases": [],
                "prompt_block": "(memento read disabled by CATDT_AGENT45_MEMENTO_DISABLE_READ)",
                "retrieval_mode": self._agent45_memento_mode(),
            }

        mode = self._agent45_memento_mode()
        if mode == "parametric":
            return self._retrieve_agent45_memento_cases_parametric(
                query_text=query_text,
                top_k=top_k,
                include_negative=include_negative,
                min_score=min_score,
            )
        return self._retrieve_agent45_memento_cases_non_parametric(
            query_text=query_text,
            top_k=top_k,
            include_negative=include_negative,
            min_score=min_score,
        )

    def build_agent45_memento_prompt_block(
        self,
        query_text: str,
        positive_cases: List[Dict[str, Any]],
        negative_cases: List[Dict[str, Any]],
        max_cases_per_group: int = 3,
    ) -> str:
        def _short(text: Any, max_len: int = 220) -> str:
            s = str(text or "").strip().replace("\n", " ")
            if len(s) > max_len:
                return s[:max_len] + "..."
            return s

        lines: List[str] = []
        lines.append("【Memento Case Memory】")
        lines.append(f"query={_short(query_text, 320)}")

        if positive_cases:
            lines.append("Positive cases (reuse these patterns):")
            for idx, case in enumerate(positive_cases[:max_cases_per_group], start=1):
                transitions = ", ".join(list(case.get("transition_signature", []) or [])[:4])
                feedback = _short(case.get("feedback", ""), 160)
                lines.append(
                    f"  + #{idx} score={case.get('retrieval_score', 0):.3f}; "
                    f"status={case.get('validation_status', '')}; transitions=[{transitions}]"
                )
                if feedback:
                    lines.append(f"    feedback={feedback}")
                # Surface info (cell + top-2-layer atom count)
                surf = case.get("surface_snapshot", {})
                if surf and surf.get("cell"):
                    n_top2 = len(surf.get("top2_layers", []))
                    lines.append(
                        f"    surface: cell={surf['cell']}, top_z={surf.get('surface_top_z')}, "
                        f"top2_layer_atoms={n_top2}"
                    )
                # Include adsorbate geometry from successful cases
                for snap in list(case.get("step_geometry_snapshot", []) or []):
                    for ep in ("reactant", "product"):
                        ep_data = snap.get(ep, {})
                        if not isinstance(ep_data, dict):
                            continue
                        core = ep_data.get("core", [])
                        staged = ep_data.get("staged", [])
                        if not core and not staged:
                            continue
                        parts = [f"{a['elem']}{a['xyz']}" for a in core]
                        for sa in staged:
                            d = sa.get("dist_to_nearest_ads", "?")
                            parts.append(f"+{sa['elem']}{sa['xyz']}(d={d}Å)")
                        lines.append(
                            f"    [{snap.get('name','')}:{ep}] {' '.join(parts)}"
                        )
        else:
            lines.append("Positive cases: (none)")

        if negative_cases:
            lines.append("Negative cases (avoid these failure patterns):")
            for idx, case in enumerate(negative_cases[:max_cases_per_group], start=1):
                issues = list(case.get("issues", []) or [])
                issue_preview = "; ".join(_short(x, 120) for x in issues[:2])
                transitions = ", ".join(list(case.get("transition_signature", []) or [])[:4])
                lines.append(
                    f"  - #{idx} score={case.get('retrieval_score', 0):.3f}; "
                    f"status={case.get('validation_status', '')}; transitions=[{transitions}]"
                )
                if issue_preview:
                    lines.append(f"    issues={issue_preview}")
                # Surface info
                surf = case.get("surface_snapshot", {})
                if surf and surf.get("cell"):
                    n_top2 = len(surf.get("top2_layers", []))
                    lines.append(
                        f"    surface: cell={surf['cell']}, top_z={surf.get('surface_top_z')}, "
                        f"top2_layer_atoms={n_top2}"
                    )
                # Show adsorbate geometry from failed cases
                for snap in list(case.get("step_geometry_snapshot", []) or []):
                    for ep in ("reactant", "product"):
                        ep_data = snap.get(ep, {})
                        if not isinstance(ep_data, dict):
                            continue
                        staged = ep_data.get("staged", [])
                        if not staged:
                            continue
                        core = ep_data.get("core", [])
                        parts = [f"{a['elem']}{a['xyz']}" for a in core]
                        for sa in staged:
                            d = sa.get("dist_to_nearest_ads", "?")
                            parts.append(f"+BAD {sa['elem']}{sa['xyz']}(d={d}Å)")
                        lines.append(
                            f"    [{snap.get('name','')}:{ep}] {' '.join(parts)}"
                        )
        else:
            lines.append("Negative cases: (none)")

        return "\n".join(lines)

    def record_agent45_memento_case(
        self,
        run_id: str,
        reaction_description: str,
        reaction_type: str,
        intermediates: List[str],
        transition_signature: List[str],
        iteration: int,
        validation_status: str,
        issues: List[str],
        feedback: str,
        design_outline: List[Dict[str, Any]],
        energy_gate_report: Optional[Dict[str, Any]] = None,
        neb_summary: Optional[Dict[str, Any]] = None,
        reward: Optional[float] = None,
        extra_query: str = "",
        step_geometry_snapshot: Optional[List[Dict[str, Any]]] = None,
        surface_snapshot: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Append one memory case in Memento-compatible JSONL format.
        """
        if self._agent45_memento_env_flag("CATDT_AGENT45_MEMENTO_DISABLE_WRITE", "0"):
            return "(memento write disabled by CATDT_AGENT45_MEMENTO_DISABLE_WRITE)"

        status = str(validation_status or "").upper()
        if reward is None:
            reward = 1.0 if status == "PASS" else 0.0

        transition_signature = [str(x).strip().lower() for x in (transition_signature or []) if str(x).strip()]
        query_text = (
            f"reaction_type={reaction_type}; "
            f"intermediates={list(intermediates or [])}; "
            f"transitions={transition_signature}; "
            f"extra={extra_query}"
        )

        case_label = "positive" if float(reward) > 0 else "negative"
        case = {
            "timestamp": datetime.now().isoformat(),
            "run_id": str(run_id),
            "iteration": int(iteration),
            "reaction_description": str(reaction_description or ""),
            "reaction_type": str(reaction_type or ""),
            "intermediates": list(intermediates or []),
            "transition_signature": transition_signature,
            "validation_status": status,
            "issues": list(issues or []),
            "feedback": str(feedback or ""),
            "energy_gate_report": dict(energy_gate_report or {}),
            "neb_summary": dict(neb_summary or {}),
            "reward": float(reward),
            "case_label": case_label,
            # Surface top-2-layer atoms + cell (shared across all steps)
            "surface_snapshot": dict(surface_snapshot or {}),
            # Per-step adsorbate geometry: core + staged atom coords
            "step_geometry_snapshot": list(step_geometry_snapshot or []),
            # Memento-compatible fields
            "question": query_text,
            "plan": json.dumps(
                {
                    "status": status,
                    "issues": list(issues or [])[:8],
                    "feedback": str(feedback or ""),
                    "design_outline": list(design_outline or []),
                },
                ensure_ascii=False,
            ),
        }

        path = self._agent45_memento_casebank_path()
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(case, ensure_ascii=False) + "\n")
        return str(path)

    @staticmethod
    def _agent45_transition_overlap(case_a: Dict[str, Any], case_b: Dict[str, Any]) -> float:
        a = {str(x).strip().lower() for x in list(case_a.get("transition_signature", []) or []) if str(x).strip()}
        b = {str(x).strip().lower() for x in list(case_b.get("transition_signature", []) or []) if str(x).strip()}
        if not a or not b:
            return 0.0
        return float(len(a & b) / max(len(a), 1))

    def _build_agent45_parametric_training_pairs(
        self,
        cases: List[Dict[str, Any]],
        seed: int = 42,
        max_pos_per_query: int = 4,
        max_neg_per_query: int = 8,
        max_pairs: int = 20000,
    ) -> List[Dict[str, Any]]:
        rng = random.Random(int(seed))
        normalized: List[Dict[str, Any]] = []
        for case in cases:
            item = dict(case)
            label = str(item.get("case_label", "")).strip().lower()
            if not label:
                label = "positive" if float(item.get("reward", 0.0) or 0.0) > 0 else "negative"
            item["case_label"] = label
            normalized.append(item)

        positive_cases = [c for c in normalized if str(c.get("case_label", "")) == "positive"]
        if not positive_cases:
            return []

        pairs: List[Dict[str, Any]] = []
        for anchor in normalized:
            query_text = self._agent45_case_query_text(anchor)
            if not query_text:
                continue

            pos_pool: List[Dict[str, Any]] = []
            neg_pool: List[Dict[str, Any]] = []
            for cand in normalized:
                overlap = self._agent45_transition_overlap(anchor, cand)
                is_positive_case = str(cand.get("case_label", "")) == "positive"
                if is_positive_case and (overlap > 0.0 or cand is anchor):
                    pos_pool.append(cand)
                else:
                    neg_pool.append(cand)

            if not pos_pool:
                fallback_pos = [c for c in positive_cases if c is anchor]
                if not fallback_pos and positive_cases:
                    fallback_pos = [rng.choice(positive_cases)]
                pos_pool = fallback_pos

            rng.shuffle(pos_pool)
            rng.shuffle(neg_pool)
            pos_selected = pos_pool[: max(1, int(max_pos_per_query))]
            neg_selected = neg_pool[: max(1, int(max_neg_per_query))]

            for cand in pos_selected:
                pairs.append(
                    {
                        "query": query_text,
                        "icl": self._agent45_case_to_icl_text(cand),
                        "label": 1,
                    }
                )
            for cand in neg_selected:
                pairs.append(
                    {
                        "query": query_text,
                        "icl": self._agent45_case_to_icl_text(cand),
                        "label": 0,
                    }
                )
            if len(pairs) >= int(max_pairs):
                break

        rng.shuffle(pairs)
        return pairs[: int(max_pairs)]

    def train_agent45_parametric_retriever(
        self,
        min_cases: int = 30,
        max_cases: int = 1500,
        max_pairs: int = 20000,
        max_pos_per_query: int = 4,
        max_neg_per_query: int = 8,
        epochs: int = 1,
        batch_size: int = 16,
        learning_rate: float = 2e-5,
        max_len: int = 256,
        val_ratio: float = 0.15,
        seed: int = 42,
        model_name: str = "",
        device: str = "",
        resume_from_checkpoint: bool = False,
    ) -> Dict[str, Any]:
        cases = self._load_agent45_memento_cases()
        if len(cases) < int(min_cases):
            return {
                "status": "SKIP",
                "reason": f"insufficient_cases<{int(min_cases)}",
                "case_count": len(cases),
            }

        subset = list(cases[-max(1, int(max_cases)):])
        pairs = self._build_agent45_parametric_training_pairs(
            cases=subset,
            seed=int(seed),
            max_pos_per_query=int(max_pos_per_query),
            max_neg_per_query=int(max_neg_per_query),
            max_pairs=int(max_pairs),
        )
        if not pairs:
            return {"status": "SKIP", "reason": "no_training_pairs", "case_count": len(subset)}

        n_pos = sum(int(x.get("label", 0)) for x in pairs)
        n_neg = len(pairs) - n_pos
        if n_pos <= 0 or n_neg <= 0:
            return {
                "status": "SKIP",
                "reason": "single_class_pairs",
                "pair_count": len(pairs),
                "positives": n_pos,
                "negatives": n_neg,
            }

        model_name = str(model_name or self._agent45_parametric_model_name()).strip()
        device = str(device or self._agent45_parametric_device()).strip().lower()
        if device == "cuda":
            import torch

            if not torch.cuda.is_available():
                raise RuntimeError("Parametric CBR training requires CUDA, but CUDA is not available.")

        retriever_dir = self._agent45_parametric_model_path().expanduser().parent
        retriever_dir.mkdir(parents=True, exist_ok=True)
        model_path = self._agent45_parametric_model_path().expanduser()
        history_path = retriever_dir / "train_history.jsonl"

        import torch
        from torch import nn
        from torch.utils.data import DataLoader, Dataset
        from transformers import AutoModel, AutoTokenizer

        np.random.seed(int(seed))
        torch.manual_seed(int(seed))
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(int(seed))

        module_path = Path(__file__).resolve().parents[2] / "deps" / "Memento" / "memory" / "train_memory_retriever.py"
        if not module_path.exists():
            raise FileNotFoundError(f"Memento train module missing: {module_path}")
        spec = importlib.util.spec_from_file_location("memento_train_memory_retriever", str(module_path))
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Failed to load module spec: {module_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        model_cls = getattr(module, "MemoryRetrieverClassifier", None)
        if model_cls is None:
            raise RuntimeError("MemoryRetrieverClassifier not found in Memento train module.")

        rng = random.Random(int(seed))
        rng.shuffle(pairs)
        split_idx = int(round(len(pairs) * (1.0 - float(val_ratio))))
        split_idx = max(1, min(split_idx, len(pairs) - 1))
        train_pairs = pairs[:split_idx]
        val_pairs = pairs[split_idx:]

        class _PairDataset(Dataset):
            def __init__(self, rows: List[Dict[str, Any]]):
                self.rows = rows

            def __len__(self) -> int:
                return len(self.rows)

            def __getitem__(self, idx: int) -> Dict[str, Any]:
                return self.rows[idx]

        tokenizer = AutoTokenizer.from_pretrained(pretrained_model_name_or_path=model_name, use_fast=True)
        loader_generator = torch.Generator()
        loader_generator.manual_seed(int(seed))

        def _collate(rows: List[Dict[str, Any]]) -> Tuple[Any, Any, Any, Any, Any]:
            case_texts = [str(r["icl"]) for r in rows]
            query_texts = [str(r["query"]) for r in rows]
            labels = torch.tensor([int(r["label"]) for r in rows], dtype=torch.long)
            t1 = tokenizer(case_texts, padding=True, truncation=True, max_length=max_len, return_tensors="pt")
            t2 = tokenizer(query_texts, padding=True, truncation=True, max_length=max_len, return_tensors="pt")
            return t1["input_ids"], t1["attention_mask"], t2["input_ids"], t2["attention_mask"], labels

        train_loader = DataLoader(
            _PairDataset(train_pairs),
            batch_size=max(1, int(batch_size)),
            shuffle=True,
            num_workers=0,
            collate_fn=_collate,
            generator=loader_generator,
        )
        val_loader = DataLoader(
            _PairDataset(val_pairs),
            batch_size=max(1, int(batch_size)),
            shuffle=False,
            num_workers=0,
            collate_fn=_collate,
        )

        backbone = AutoModel.from_pretrained(pretrained_model_name_or_path=model_name, use_safetensors=True)
        model = model_cls(backbone).to(device)

        if bool(resume_from_checkpoint) and model_path.exists():
            state_dict = torch.load(str(model_path), map_location=device)
            model.load_state_dict(state_dict)

        optimizer = torch.optim.AdamW(model.parameters(), lr=float(learning_rate))
        class_weight = torch.tensor(
            [1.0, max(1.0, float(n_neg) / max(float(n_pos), 1.0))],
            device=device,
        )
        criterion = nn.CrossEntropyLoss(weight=class_weight)

        best_val_acc = -1.0
        best_state: Optional[Dict[str, Any]] = None
        train_loss_curve: List[float] = []
        val_acc_curve: List[float] = []
        val_loss_curve: List[float] = []

        for _epoch in range(max(1, int(epochs))):
            model.train()
            epoch_losses: List[float] = []
            for ids1, mask1, ids2, mask2, labels in train_loader:
                ids1 = ids1.to(device)
                mask1 = mask1.to(device)
                ids2 = ids2.to(device)
                mask2 = mask2.to(device)
                labels = labels.to(device)

                optimizer.zero_grad(set_to_none=True)
                logits = model(ids1, mask1, ids2, mask2)
                loss = criterion(logits, labels)
                loss.backward()
                optimizer.step()
                epoch_losses.append(float(loss.detach().cpu().item()))

            train_loss_curve.append(float(sum(epoch_losses) / max(len(epoch_losses), 1)))

            model.eval()
            correct = 0
            total = 0
            val_loss_items: List[float] = []
            with torch.inference_mode():
                for ids1, mask1, ids2, mask2, labels in val_loader:
                    ids1 = ids1.to(device)
                    mask1 = mask1.to(device)
                    ids2 = ids2.to(device)
                    mask2 = mask2.to(device)
                    labels = labels.to(device)
                    logits = model(ids1, mask1, ids2, mask2)
                    vloss = criterion(logits, labels)
                    val_loss_items.append(float(vloss.detach().cpu().item()))
                    pred = torch.argmax(logits, dim=1)
                    correct += int((pred == labels).sum().detach().cpu().item())
                    total += int(labels.shape[0])
            val_acc = float(correct / max(total, 1))
            val_loss = float(sum(val_loss_items) / max(len(val_loss_items), 1))
            val_acc_curve.append(val_acc)
            val_loss_curve.append(val_loss)
            if val_acc >= best_val_acc:
                best_val_acc = val_acc
                best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}

        if best_state is None:
            best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}

        torch.save(best_state, str(model_path))

        meta = {
            "timestamp": datetime.now().isoformat(),
            "mode": "parametric",
            "model_name": model_name,
            "model_path": str(model_path),
            "device": device,
            "case_count": len(subset),
            "pair_count": len(pairs),
            "train_pairs": len(train_pairs),
            "val_pairs": len(val_pairs),
            "positives": n_pos,
            "negatives": n_neg,
            "epochs": int(epochs),
            "batch_size": int(batch_size),
            "learning_rate": float(learning_rate),
            "max_len": int(max_len),
            "best_val_acc": float(best_val_acc),
            "train_loss_curve": train_loss_curve,
            "val_loss_curve": val_loss_curve,
            "val_acc_curve": val_acc_curve,
        }
        with open(retriever_dir / "train_meta.json", "w", encoding="utf-8") as fh:
            json.dump(meta, fh, ensure_ascii=False, indent=2)
        with open(history_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(meta, ensure_ascii=False) + "\n")

        # Invalidate cached retriever to ensure new checkpoint is loaded on next read.
        cache = self._agent45_parametric_cache()
        cache.clear()

        return {
            "status": "OK",
            "mode": "parametric",
            "model_path": str(model_path),
            "model_name": model_name,
            "best_val_acc": float(best_val_acc),
            "pair_count": len(pairs),
            "case_count": len(subset),
            "train_meta_path": str(retriever_dir / "train_meta.json"),
            "history_path": str(history_path),
            "train_loss_curve": train_loss_curve,
            "val_loss_curve": val_loss_curve,
            "val_acc_curve": val_acc_curve,
        }

    def get_agent45_parametric_retriever_status(self) -> Dict[str, Any]:
        model_path = self._agent45_parametric_model_path().expanduser()
        meta_path = model_path.parent / "train_meta.json"
        payload: Dict[str, Any] = {
            "mode": self._agent45_memento_mode(),
            "model_path": str(model_path),
            "model_exists": model_path.exists(),
            "meta_path": str(meta_path),
            "meta_exists": meta_path.exists(),
        }
        if meta_path.exists():
            try:
                payload["meta"] = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception as exc:
                payload["meta_error"] = str(exc)
        return payload


# ---------------------------------------------------------------------------
# Reporting tools mixin
# ---------------------------------------------------------------------------


class ReportingToolsMixin:
    def export_step_structures(
        self,
        steps: List[Dict[str, Any]],
        output_dir: str,
        prefix: str = "tool_baseline",
        with_images: bool = True,
    ) -> Dict[str, Any]:
        """Export per-step reactant/product structures for audit before/after LLM edits."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        visualizer = None
        if with_images:
            try:
                from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

                visualizer = CatalystSurfaceVisualizer(renderer="tachyon", quality="low", auto_expand=False)
            except Exception as exc:
                logger.warning("Step-structure visualizer unavailable: %s", exc)
                visualizer = None

        manifest: Dict[str, Any] = {"prefix": prefix, "steps": []}
        for idx, step in enumerate(steps, start=1):
            step_name = str(step.get("name", f"step_{idx}"))
            safe_name = ''.join(ch if (ch.isalnum() or ch in ('_', '-')) else '_' for ch in step_name)

            reactant = step.get("reactant")
            product = step.get("product")
            step_entry: Dict[str, Any] = {
                "index": idx,
                "name": step_name,
                "reactant_path": "",
                "product_path": "",
                "reactant_image": "",
                "product_image": "",
            }

            if isinstance(reactant, Atoms):
                reactant_path = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_reactant.vasp"
                write(reactant_path, reactant)
                step_entry["reactant_path"] = str(reactant_path)
                if visualizer is not None:
                    reactant_png = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_reactant.png"
                    visualizer.visualize_structure(
                        reactant,
                        output_file=str(reactant_png),
                        title=f"{prefix} Step {idx} Reactant",
                        show=False,
                    )
                    step_entry["reactant_image"] = str(reactant_png)

            if isinstance(product, Atoms):
                product_path = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_product.vasp"
                write(product_path, product)
                step_entry["product_path"] = str(product_path)
                if visualizer is not None:
                    product_png = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_product.png"
                    visualizer.visualize_structure(
                        product,
                        output_file=str(product_png),
                        title=f"{prefix} Step {idx} Product",
                        show=False,
                    )
                    step_entry["product_image"] = str(product_png)

            manifest["steps"].append(step_entry)

        manifest_path = out_dir / f"{prefix}_manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, ensure_ascii=False, indent=2)
        manifest["manifest_path"] = str(manifest_path)
        return manifest

    def generate_final_report(
        self,
        run_id: str = "default_run",
        workflow_step: str = "06_report",
        surff_result_path: Optional[str] = None,
        adsorbdiff_result_path: Optional[str] = None,
        vssr_mc_result_path: Optional[str] = None,
        pathway_result_path: Optional[str] = None,
        kmc_result_path: Optional[str] = None,
        llm_review_results: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Generates a comprehensive HTML summary report of the entire CatDT workflow.

        Args:
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "06_report").
            surff_result_path (str, optional): Path to SurFFPredictionResult pickle.
            adsorbdiff_result_path (str, optional): Path to AdsorbDiffPredictionOutput pickle.
            vssr_mc_result_path (str, optional): Path to VSSRMCResult pickle.
            pathway_result_path (str, optional): Path to CompletePathwayResult pickle.
            kmc_result_path (str, optional): Path to CatMAPResult pickle.
            llm_review_results (Dict[str, Any], optional): Dictionary of LLM review outputs for various steps.

        Returns:
            str: Path to the generated HTML summary file.
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Generating final report for run {run_id} into {output_dir}")

        viz_manager = self._get_viz_manager_instance(output_dir)
        all_visualizations = {
            'surfaces': {},
            'adsorption': {},
            'mc_trajectory': {},
            'reaction_pathway': '',
            'energy_diagram': '',
            'kmc_dynamics': '',
        }

        # Collect existing visualization paths from previous steps
        surff_viz_dir = self.output_base_dir / run_id / "01_surfaces" / "visualizations"
        if surff_viz_dir.exists():
            for viz_file in surff_viz_dir.glob("*.png"):
                all_visualizations['surfaces'][viz_file.stem] = str(viz_file)
            for viz_file in surff_viz_dir.glob("*.html"):
                 all_visualizations['surfaces'][viz_file.stem] = str(viz_file)

        ads_viz_dir = self.output_base_dir / run_id / "02_adsorption" / "visualizations"
        if ads_viz_dir.exists():
            for viz_file in ads_viz_dir.glob("*.png"):
                all_visualizations['adsorption'][viz_file.stem] = str(viz_file)

        mc_viz_dir = self.output_base_dir / run_id / "03_reconstruction" / "visualizations"
        if mc_viz_dir.exists():
            for viz_file in mc_viz_dir.glob("*.gif"):
                all_visualizations['mc_trajectory'][viz_file.stem] = str(viz_file)

        pathway_viz_dir = self.output_base_dir / run_id / "04_pathway_analysis" / "visualizations"
        if pathway_viz_dir.exists():
            for viz_file in pathway_viz_dir.glob("reaction_pathway.gif"):
                all_visualizations['reaction_pathway'] = str(viz_file)
            for viz_file in pathway_viz_dir.glob("energy_diagram.png"):
                all_visualizations['energy_diagram'] = str(viz_file)

        kmc_viz_dir = self.output_base_dir / run_id / "05_kmc_simulation" / "visualizations"
        if kmc_viz_dir.exists():
            for viz_file in kmc_viz_dir.glob("kmc_surface_dynamics.gif"):
                all_visualizations['kmc_dynamics'] = str(viz_file)


        # Load actual result objects if paths are provided
        full_results = {}
        if surff_result_path:
            full_results['surff'] = self._load_result_from_pickle(surff_result_path)
        if adsorbdiff_result_path:
            full_results['adsorbdiff'] = self._load_result_from_pickle(adsorbdiff_result_path)
        if vssr_mc_result_path:
            full_results['vssr_mc'] = self._load_result_from_pickle(vssr_mc_result_path)
        if pathway_result_path:
            full_results['pathway'] = self._load_result_from_pickle(pathway_result_path)
        if kmc_result_path:
            full_results['kmc'] = self._load_result_from_pickle(kmc_result_path)

        full_results['llm_review'] = llm_review_results


        html_path = viz_manager.generate_summary_html(
            viz_outputs=all_visualizations,
            output_file="final_report.html",
        )
        logger.info(f"Generated final HTML report: {html_path}")
        return html_path

    def check_checkpoint_exists(self, run_id: str) -> Optional[str]:
        """
        检查是否存在可恢复的检查点

        Returns
        -------
        str or None
            最新的检查点步骤名称，如果没有则返回 None
        """
        checkpoint_manager = CheckpointManager()
        checkpoint = checkpoint_manager.load_latest_checkpoint(run_id)

        if checkpoint:
            logger.info(f"Found checkpoint at step: {checkpoint.step_name}")
            return checkpoint.step_name
        return None

    def list_available_checkpoints(self, run_id: str) -> List[Dict]:
        """
        列出所有可用的检查点

        Returns
        -------
        List[Dict]
            检查点列表，每个包含 step, timestamp 等信息
        """
        checkpoint_manager = CheckpointManager()
        return checkpoint_manager.list_checkpoints(run_id)
