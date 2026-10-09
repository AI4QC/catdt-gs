"""CAMEL runtime primitives for agent task execution."""

from __future__ import annotations

import copy
import json
import logging
import os
import re

logger = logging.getLogger(__name__)
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from camel.agents import ChatAgent
from camel.memories import ChatHistoryMemory
from camel.toolkits import FunctionTool


@dataclass
class TaskOutput:
    raw: str = ""
    json_dict: Dict[str, Any] = field(default_factory=dict)
    pydantic: Optional[Any] = None
    tool_call_record: Optional[Any] = None


@dataclass
class Task:
    description: str
    expected_output: str
    agent: Any
    output_pydantic: Optional[Any] = None
    result_handler: Optional[Any] = None
    allow_tool_calls: Optional[bool] = None


class CamelWorkflowAgent:
    @staticmethod
    def _prepare_tools(tools: Optional[List[FunctionTool]]) -> List[FunctionTool]:
        """Clone tools and freeze validated schemas to avoid runtime schema drift."""
        prepared: List[FunctionTool] = []
        for tool in tools or []:
            if not isinstance(tool, FunctionTool):
                prepared.append(tool)
                continue

            tool_name = tool.get_function_name()
            try:
                schema = copy.deepcopy(tool.get_openai_tool_schema())
                cloned_tool = FunctionTool(tool.func, openai_tool_schema=schema)
                cloned_tool.get_openai_tool_schema = (  # type: ignore[method-assign]
                    lambda _schema=schema: _schema
                )
                cloned_tool.get_openai_function_schema = (  # type: ignore[method-assign]
                    lambda _schema=schema: _schema["function"]
                )
                prepared.append(cloned_tool)
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to prepare CAMEL tool schema for '{tool_name}': {exc}"
                ) from exc
        return prepared

    def __init__(
        self,
        role_name: str,
        goal: str,
        backstory: str,
        model_backend: Any,
        memory: Optional[ChatHistoryMemory] = None,
        tools: Optional[List[FunctionTool]] = None,
        step_timeout: float = 600.0,
        message_window_size: int = 40,
        summarize_threshold: int = 80,
        token_limit: Optional[int] = None,
    ):
        system_message = (
            f"Role: {role_name}\n"
            f"Goal: {goal}\n"
            f"Backstory: {backstory}\n"
            "When the task requires JSON, output exactly one machine-parseable JSON object.\n"
            "Do not output analysis, chain-of-thought, markdown, code fences, or any extra text."
        )
        prepared_tools = self._prepare_tools(tools or [])
        # No max_iteration — let the loop run until LLM stops calling tools.
        # Timeout is controlled by step_timeout (set per-agent).
        self.chat_agent = ChatAgent(
            system_message=system_message,
            model=model_backend,
            memory=memory,
            tools=prepared_tools,
            message_window_size=int(message_window_size),
            summarize_threshold=int(summarize_threshold),
            token_limit=(int(token_limit) if token_limit is not None else None),
        )
        self.chat_agent.step_timeout = step_timeout

    @staticmethod
    def _extract_json(text: str) -> Dict[str, Any]:
        if not isinstance(text, str) or not text.strip():
            return {}

        stripped = text.strip()
        block = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, re.DOTALL | re.IGNORECASE)
        candidate = block.group(1) if block else stripped

        if block:
            try:
                parsed = json.loads(candidate)
                return parsed if isinstance(parsed, dict) else {}
            except Exception:
                pass

        # Try full-text direct parse first
        try:
            parsed = json.loads(stripped)
            if isinstance(parsed, dict):
                return parsed
        except Exception:
            pass

        # Fallback: scan all balanced {...} segments and parse from tail.
        segments: List[str] = []
        start = None
        depth = 0
        in_string = False
        escaped = False
        for idx, ch in enumerate(stripped):
            if in_string:
                if escaped:
                    escaped = False
                elif ch == "\\":
                    escaped = True
                elif ch == '"':
                    in_string = False
                continue

            if ch == '"':
                in_string = True
                continue
            if ch == "{":
                if depth == 0:
                    start = idx
                depth += 1
                continue
            if ch == "}" and depth > 0:
                depth -= 1
                if depth == 0 and start is not None:
                    segments.append(stripped[start:idx + 1])
                    start = None

        for seg in reversed(segments):
            try:
                parsed = json.loads(seg)
                if isinstance(parsed, dict):
                    return parsed
            except Exception:
                continue
        return {}

    @staticmethod
    def _model_validate(output_model: Any, payload: Dict[str, Any]) -> Any:
        if output_model is None:
            return None
        if not payload:
            return None
        try:
            if hasattr(output_model, "model_validate"):
                return output_model.model_validate(payload)
            return output_model(**payload)
        except Exception as exc:
            logger.debug(
                "Output-model validation failed for %s: %s",
                getattr(output_model, "__name__", output_model),
                exc,
            )
            return None

    @staticmethod
    def _is_length_finish_error(exc: Exception) -> bool:
        text = str(exc or "")
        exc_name = exc.__class__.__name__
        return (
            "LengthFinishReasonError" in exc_name
            or "LengthFinishReasonError" in text
            or "length limit was reached" in text.lower()
        )

    @staticmethod
    def _is_response_format_parse_error(exc: Exception) -> bool:
        text = str(exc or "")
        lowered = text.lower()
        return (
            "_AnnotatedAlias.__init__() missing" in text
            or "Cannot instantiate typing_extensions.Required" in text
            or "typing_extensions.required" in lowered
            or "json_invalid" in text
        )

    def run(
        self,
        prompt: str,
        output_model: Optional[Any] = None,
        allow_tool_calls: bool = True,
    ) -> TaskOutput:
        # Root-cause fix:
        # OpenAI-compatible structured parse (`response_format`) is unstable across SDK/model variants
        # in our runtime stack (e.g., typing_extensions.Required / AnnotatedAlias transform failures).
        # Keep a single generic pathway: plain text generation + local JSON extraction + local pydantic validation.
        # Can be re-enabled for experiments by setting CATDT_ENABLE_SERVER_STRUCTURED_PARSE=1.
        enable_server_structured_parse = os.getenv("CATDT_ENABLE_SERVER_STRUCTURED_PARSE", "0").strip() in {
            "1", "true", "TRUE", "yes", "YES"
        }
        response_format = (
            output_model
            if (output_model is not None and not allow_tool_calls and enable_server_structured_parse)
            else None
        )
        before_sig = getattr(self.chat_agent, "_last_tool_call_signature", None)

        def _step_with_optional_tools(rf: Optional[Any]):
            if allow_tool_calls:
                return self.chat_agent.step(prompt, response_format=rf)
            tool_dict = getattr(self.chat_agent, "tool_dict", {}) or {}
            removed_tools = list(tool_dict.values())
            removed_names = [tool.get_function_name() for tool in removed_tools]
            if removed_names:
                self.chat_agent.remove_tools(removed_names)
            try:
                return self.chat_agent.step(prompt, response_format=rf)
            finally:
                if removed_tools:
                    self.chat_agent.add_tools(removed_tools)

        try:
            response = _step_with_optional_tools(response_format)
        except Exception as exc:
            # Structured parsing can fail when model output is truncated at max tokens.
            # Retry once without server-side response_format, then parse JSON locally.
            if (
                (not allow_tool_calls)
                and (output_model is not None)
                and (response_format is not None)
                and (
                    self._is_length_finish_error(exc)
                    or self._is_response_format_parse_error(exc)
                )
            ):
                response = _step_with_optional_tools(None)
            else:
                raise

        after_sig = getattr(self.chat_agent, "_last_tool_call_signature", None)
        tool_call_record = None
        if after_sig and after_sig != before_sig:
            tool_call_record = getattr(self.chat_agent, "_last_tool_call_record", None)
        content = ""
        if hasattr(response, "msgs") and response.msgs:
            content = getattr(response.msgs[-1], "content", "") or ""

        # If tool calls were made but the final response has no text content,
        # do a follow-up call WITHOUT tools to force the LLM to output JSON.
        if allow_tool_calls and not content.strip():
            logger.info("Tool-calling loop ended with empty content; requesting final JSON output.")
            follow_up = (
                "Tool calls are complete. Now output ONLY the final JSON object. "
                "No tool calls, no explanations — just the JSON."
            )
            tool_dict = getattr(self.chat_agent, "tool_dict", {}) or {}
            removed_tools = list(tool_dict.values())
            removed_names = [tool.get_function_name() for tool in removed_tools]
            if removed_names:
                self.chat_agent.remove_tools(removed_names)
            try:
                follow_up_response = self.chat_agent.step(follow_up, response_format=None)
                if hasattr(follow_up_response, "msgs") and follow_up_response.msgs:
                    content = getattr(follow_up_response.msgs[-1], "content", "") or ""
            except Exception as exc:
                logger.warning("Follow-up JSON request failed: %s", exc)
            finally:
                if removed_tools:
                    self.chat_agent.add_tools(removed_tools)

        payload = self._extract_json(content)
        parsed = self._model_validate(output_model, payload)
        return TaskOutput(
            raw=content,
            json_dict=payload,
            pydantic=parsed,
            tool_call_record=tool_call_record,
        )
