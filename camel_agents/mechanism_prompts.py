"""Backward-compatibility re-exports.

All mechanism agent prompts now live in ``camel_agents.prompts`` alongside
Agent1-7 prompts.  This module re-exports for any code that imported from here.
"""

from camel_agents.prompts import (  # noqa: F401
    MECHANISM_CONTEXT_ROUTING_AGENT,
    MECHANISM_SEARCH_TRIAGE_AGENT,
    TaskPrompts as MechanismTaskPrompts,
    get_agent_role as get_mechanism_agent_role,
    get_all_agent_roles as get_mechanism_agent_roles,
)

__all__ = [
    "MECHANISM_CONTEXT_ROUTING_AGENT",
    "MECHANISM_SEARCH_TRIAGE_AGENT",
    "MechanismTaskPrompts",
    "get_mechanism_agent_roles",
    "get_mechanism_agent_role",
]
