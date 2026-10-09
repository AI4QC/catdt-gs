import sys
from pathlib import Path
from types import SimpleNamespace

_repo_root = Path(__file__).resolve().parents[1]
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import pytest
from camel.utils.token_counting import BaseTokenCounter

import camel_agents.camel_model_backend as backend_mod
from camel_agents.workflow import CatDTCamelWorkflow
from camel_agents.camel_model_backend import get_camel_model_backend


def test_workflow_has_agents():
    wf = CatDTCamelWorkflow(initialize_agents=False)
    assert hasattr(wf, "agent1")
    assert hasattr(wf, "agent4")
    assert hasattr(wf, "agent5")


class DummyTokenCounter(BaseTokenCounter):
    def count_tokens_from_messages(self, messages):
        return len(messages or [])

    def encode(self, text):
        return [1] * len(str(text))

    def decode(self, token_ids):
        return "x" * len(token_ids or [])


def test_workflow_memory_uses_backend_token_counter(monkeypatch):
    wf = CatDTCamelWorkflow(initialize_agents=False)
    backend_token_counter = DummyTokenCounter()
    wf.llm = SimpleNamespace(token_counter=backend_token_counter)

    memory = wf._make_memory(token_limit=1024, window_size=8)
    context_creator = memory.get_context_creator()

    assert memory is not None
    assert context_creator.token_counter is backend_token_counter


def test_get_camel_model_backend_uses_claude_token_counter(monkeypatch):
    captured = {}

    def fake_create(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(
            token_counter=kwargs.get("token_counter"),
            token_limit=200000,
            model_type=kwargs.get("model_type"),
        )

    monkeypatch.setenv("OPENAI_MODEL", "claude-opus-4-6")
    monkeypatch.setenv("OPENAI_API_KEY", "dummy-key")
    monkeypatch.setenv("OPENAI_BASE_URL", "https://example.invalid/v1")
    monkeypatch.setattr(backend_mod.ModelFactory, "create", fake_create)

    backend = get_camel_model_backend()

    assert backend is not None
    assert "token_counter" in captured
    assert captured["token_counter"] is backend.token_counter
    assert type(captured["token_counter"]).__name__ == "LiteLLMTokenCounter"


def test_gpt5_model_config_uses_max_completion_tokens():
    cfg = backend_mod._build_model_config(0.2, model_name="gpt-5.4", max_tokens=1024, top_p=0.9)

    assert "max_tokens" not in cfg
    assert cfg["max_completion_tokens"] == 1024
    assert cfg["top_p"] == 0.9
