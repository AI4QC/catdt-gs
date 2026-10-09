#!/usr/bin/env python
"""CAMEL model backend utilities."""

from __future__ import annotations

import os
from typing import Any, Dict, Optional

from camel.models import ModelFactory
from camel.types import ModelPlatformType, ModelType
from camel.types.unified_model_type import UnifiedModelType
from camel.utils.token_counting import LiteLLMTokenCounter


def _parse_custom_token_limits() -> Dict[str, int]:
    raw = str(os.getenv("CATDT_CUSTOM_MODEL_TOKEN_LIMITS", "")).strip()
    if not raw:
        return {}
    try:
        import json

        obj = json.loads(raw)
        if isinstance(obj, dict):
            out: Dict[str, int] = {}
            for key, value in obj.items():
                try:
                    out[str(key).strip()] = int(value)
                except Exception:
                    continue
            return out
    except Exception:
        pass
    return {}


def _patch_unified_model_token_limit() -> None:
    """
    Patch CAMEL UnifiedModelType.token_limit for custom model aliases.

    Why:
    - CAMEL treats unknown model strings as UnifiedModelType and logs:
      "Unknown model ... context window size not defined".
    - We use OPENAI-compatible model names (e.g. `gpt-5.4`, `claude-opus-4-6`).
    - Add a deterministic local mapping to avoid noisy warning and undefined limit.
    """
    try:
        from camel.types.unified_model_type import UnifiedModelType
    except Exception:
        return

    if getattr(UnifiedModelType, "_catdt_token_limit_patched", False):
        return

    original_prop = getattr(UnifiedModelType, "token_limit", None)
    original_getter = getattr(original_prop, "fget", None)
    if original_getter is None:
        return

    default_map: Dict[str, int] = {
        "claude-opus-4-6": int(os.getenv("CATDT_MODEL_TOKEN_LIMIT_CLAUDE_OPUS_46", "200000")),
        "gpt-5.4": int(os.getenv("CATDT_MODEL_TOKEN_LIMIT_GPT_54", "128000")),
    }
    default_map.update(_parse_custom_token_limits())

    def _patched_token_limit(self: Any) -> int:
        model_name = str(self)
        if model_name in default_map:
            return int(default_map[model_name])
        return int(original_getter(self))

    UnifiedModelType._catdt_token_limit_patched = True  # type: ignore[attr-defined]
    UnifiedModelType._catdt_original_token_limit_getter = original_getter  # type: ignore[attr-defined]
    UnifiedModelType.token_limit = property(_patched_token_limit)  # type: ignore[assignment]


def _uses_max_completion_tokens(model_name: Optional[str]) -> bool:
    name = str(model_name or "").strip().lower()
    return name.startswith("gpt-5")


def _build_model_config(temperature: float, **kwargs: Any) -> Dict[str, Any]:
    model_name = str(kwargs.pop("model_name", "") or "")
    model_config: Dict[str, Any] = {"temperature": temperature}
    max_tokens = kwargs.get("max_tokens")
    if max_tokens is not None:
        if _uses_max_completion_tokens(model_name):
            model_config["max_completion_tokens"] = max_tokens
        else:
            model_config["max_tokens"] = max_tokens
    for key in ("top_p", "frequency_penalty", "presence_penalty", "seed"):
        if key in kwargs and kwargs[key] is not None:
            model_config[key] = kwargs[key]
    return model_config


def _build_token_counter(model_name: str):
    """Build a token counter that tracks the configured Claude model itself.

    CAMEL's `OpenAICompatibleModel.token_counter` currently hardcodes
    `OpenAITokenCounter(ModelType.GPT_4O_MINI)` when no counter is injected.
    For our strict `claude-opus-4-6` backend this is architecturally wrong and
    causes the GPT tokenizer path to appear in runtime diagnostics. Inject a
    model-aware counter explicitly at backend creation time.
    """
    return LiteLLMTokenCounter(UnifiedModelType(model_name))


def get_camel_model_backend(
    model: Optional[str] = None,
    api_key: Optional[str] = None,
    base_url: Optional[str] = None,
    temperature: float = 0.2,
    **kwargs: Any,
):
    """Create a strict CAMEL model backend (no fallback)."""
    if os.getenv("CATDT_LLM_BACKEND", "").strip().lower() == "codex":
        # local codex-api (codex login) instead of an OpenAI-compatible HTTP API
        from camel_agents.codex_backend import codex_backend_from_env
        return codex_backend_from_env()
    _patch_unified_model_token_limit()
    model_name = str(model or os.getenv("OPENAI_MODEL", "gpt-5.4") or "").strip()
    if not model_name:
        raise RuntimeError("OPENAI_MODEL is required (e.g. 'gpt-5.4').")

    key = api_key or os.getenv("OPENAI_API_KEY")
    url = base_url or os.getenv("OPENAI_BASE_URL")

    if not key:
        raise RuntimeError("OPENAI_API_KEY is required. Fallback model is disabled.")

    platform = ModelPlatformType.OPENAI_COMPATIBLE_MODEL if url else ModelPlatformType.DEFAULT
    request_timeout = kwargs.pop("timeout", None)
    if request_timeout is None:
        request_timeout = float(os.getenv("OPENAI_REQUEST_TIMEOUT_SEC", "1800"))

    max_retries = kwargs.pop("max_retries", None)
    if max_retries is None:
        max_retries = int(os.getenv("OPENAI_MAX_RETRIES", "3"))
    if "seed" not in kwargs or kwargs.get("seed") is None:
        env_seed = str(os.getenv("OPENAI_SEED", "")).strip()
        if env_seed:
            kwargs["seed"] = int(env_seed)

    backend = ModelFactory.create(
        model_platform=platform,
        model_type=model_name,
        api_key=key,
        url=url,
        token_counter=_build_token_counter(model_name),
        timeout=float(request_timeout),
        max_retries=int(max_retries),
        model_config_dict=_build_model_config(temperature, model_name=model_name, **kwargs),
    )

    # Some API mirrors (e.g. AICodeMirror) return content=null in non-stream
    # mode but work correctly in stream mode. Monkey-patch the backend's run()
    # to always use streaming and reassemble into a standard ChatCompletion,
    # so the rest of the CAMEL framework sees a normal synchronous response.
    if url and "openai.com" not in url:
        _patch_always_stream(backend)

    return backend


def _patch_always_stream(backend) -> None:
    """Monkey-patch *backend.run* to always stream and reassemble."""
    original_run = backend.run

    def _patched_run(messages, *args, **kwargs):
        saved = backend.model_config_dict.get("stream")
        backend.model_config_dict["stream"] = True
        try:
            stream = original_run(messages, *args, **kwargs)
            return _consume_stream(stream)
        finally:
            if saved is None:
                backend.model_config_dict.pop("stream", None)
            else:
                backend.model_config_dict["stream"] = saved

    backend.run = _patched_run


def _consume_stream(stream):
    """Consume an OpenAI stream → ChatCompletion (supports tool_calls)."""
    from openai.types.chat import ChatCompletion, ChatCompletionMessage
    from openai.types.chat.chat_completion import Choice
    from openai.types.chat.chat_completion_message_tool_call import (
        ChatCompletionMessageToolCall,
        Function,
    )
    from openai.types.completion_usage import CompletionUsage

    parts = []
    finish_reason = "stop"
    usage_data = None
    resp_id = ""
    model = ""
    # Accumulate tool calls: index → {id, type, function.name, function.arguments}
    tool_calls_acc: dict = {}

    for chunk in stream:
        if not chunk.choices:
            if hasattr(chunk, "usage") and chunk.usage:
                usage_data = chunk.usage
            continue
        delta = chunk.choices[0].delta
        if hasattr(delta, "content") and delta.content:
            parts.append(delta.content)
        # Accumulate tool_calls from deltas
        if hasattr(delta, "tool_calls") and delta.tool_calls:
            for tc_delta in delta.tool_calls:
                idx = tc_delta.index
                if idx not in tool_calls_acc:
                    tool_calls_acc[idx] = {
                        "id": "", "type": "function",
                        "name": "", "arguments": "",
                    }
                acc = tool_calls_acc[idx]
                if tc_delta.id:
                    acc["id"] = tc_delta.id
                if tc_delta.type:
                    acc["type"] = tc_delta.type
                if hasattr(tc_delta, "function") and tc_delta.function:
                    if tc_delta.function.name:
                        acc["name"] = tc_delta.function.name
                    if tc_delta.function.arguments:
                        acc["arguments"] += tc_delta.function.arguments
        if chunk.choices[0].finish_reason:
            finish_reason = chunk.choices[0].finish_reason
        if not resp_id and hasattr(chunk, "id"):
            resp_id = chunk.id
        if not model and hasattr(chunk, "model"):
            model = chunk.model

    content = "".join(parts) or None

    # Build tool_calls list if any were accumulated
    tool_calls_list = None
    if tool_calls_acc:
        tool_calls_list = []
        for idx in sorted(tool_calls_acc):
            acc = tool_calls_acc[idx]
            tool_calls_list.append(
                ChatCompletionMessageToolCall(
                    id=acc["id"],
                    type=acc["type"],
                    function=Function(
                        name=acc["name"],
                        arguments=acc["arguments"],
                    ),
                )
            )

    message = ChatCompletionMessage(
        role="assistant",
        content=content,
        tool_calls=tool_calls_list,
    )
    choice = Choice(index=0, message=message, finish_reason=finish_reason)
    usage = None
    if usage_data:
        usage = CompletionUsage(
            prompt_tokens=getattr(usage_data, "prompt_tokens", 0),
            completion_tokens=getattr(usage_data, "completion_tokens", 0),
            total_tokens=getattr(usage_data, "total_tokens", 0),
        )
    return ChatCompletion(
        id=resp_id or "stream-fallback",
        choices=[choice],
        created=0,
        model=model,
        object="chat.completion",
        usage=usage,
    )
