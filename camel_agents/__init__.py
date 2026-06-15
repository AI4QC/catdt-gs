"""CAMEL-based multi-agent package for CatDT.

Keep package import lightweight: delay heavy CAMEL/workflow imports until actually used.
"""

__all__ = [
    "AgentRole",
    "TaskPrompts",
    "get_all_agent_roles",
    "get_agent_role",
    "CatDTCamelWorkflow",
    "ReactionContext",
    "AtomAddSpec",
    "AtomModifySpec",
    "PathwayStepSpec",
    "PathwayDesign",
    "ValidationReport",
    "WorkflowState",
    "CatDTConfig",
    "EvolvablePolicy",
    "Task",
    "TaskOutput",
    "CamelWorkflowAgent",
    # Mechanism search schemas
    "MechanismContext",
    "MechanismSearchConfig",
    "CatalyticStateRecord",
    "ElementaryStepCandidate",
    "CandidatePathway",
    "MultiPathwaySearchResult",
    "MechanismTaskPrompts",
]


def __getattr__(name):
    if name in {"AgentRole", "TaskPrompts", "get_all_agent_roles", "get_agent_role"}:
        from .prompts import AgentRole, TaskPrompts, get_all_agent_roles, get_agent_role

        mapping = {
            "AgentRole": AgentRole,
            "TaskPrompts": TaskPrompts,
            "get_all_agent_roles": get_all_agent_roles,
            "get_agent_role": get_agent_role,
        }
        return mapping[name]

    if name == "CatDTCamelWorkflow":
        from .workflow import CatDTCamelWorkflow

        return CatDTCamelWorkflow

    if name in {
        "ReactionContext",
        "AtomAddSpec",
        "AtomModifySpec",
        "PathwayStepSpec",
        "PathwayDesign",
        "ValidationReport",
        "WorkflowState",
        "CatDTConfig",
    }:
        from .schemas import (
            AtomAddSpec,
            AtomModifySpec,
            CatDTConfig,
            PathwayDesign,
            PathwayStepSpec,
            ReactionContext,
            ValidationReport,
            WorkflowState,
        )

        mapping = {
            "ReactionContext": ReactionContext,
            "AtomAddSpec": AtomAddSpec,
            "AtomModifySpec": AtomModifySpec,
            "PathwayStepSpec": PathwayStepSpec,
            "PathwayDesign": PathwayDesign,
            "ValidationReport": ValidationReport,
            "WorkflowState": WorkflowState,
            "CatDTConfig": CatDTConfig,
        }
        return mapping[name]

    if name == "EvolvablePolicy":
        from .policy import EvolvablePolicy

        return EvolvablePolicy

    if name in {"Task", "TaskOutput", "CamelWorkflowAgent"}:
        from .runtime import CamelWorkflowAgent, Task, TaskOutput

        mapping = {
            "Task": Task,
            "TaskOutput": TaskOutput,
            "CamelWorkflowAgent": CamelWorkflowAgent,
        }
        return mapping[name]

    if name in {
        "MechanismContext",
        "MechanismSearchConfig",
        "CatalyticStateRecord",
        "ElementaryStepCandidate",
        "CandidatePathway",
        "MultiPathwaySearchResult",
    }:
        from .mechanism_schemas import (
            MechanismContext,
            MechanismSearchConfig,
            CatalyticStateRecord,
            ElementaryStepCandidate,
            CandidatePathway,
            MultiPathwaySearchResult,
        )

        mapping = {
            "MechanismContext": MechanismContext,
            "MechanismSearchConfig": MechanismSearchConfig,
            "CatalyticStateRecord": CatalyticStateRecord,
            "ElementaryStepCandidate": ElementaryStepCandidate,
            "CandidatePathway": CandidatePathway,
            "MultiPathwaySearchResult": MultiPathwaySearchResult,
        }
        return mapping[name]

    if name == "MechanismTaskPrompts":
        from .mechanism_prompts import MechanismTaskPrompts

        return MechanismTaskPrompts

    raise AttributeError(f"module 'camel_agents' has no attribute '{name}'")
