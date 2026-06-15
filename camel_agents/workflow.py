"""CatDT CAMEL workflow: orchestration + runtime, helper logic delegated to tools."""

from __future__ import annotations

import copy
import json
import logging
import os
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.io import read, write
from camel.memories import ChatHistoryMemory, ScoreBasedContextCreator
from camel.toolkits import FunctionTool

_current_file_dir = Path(__file__).parent.resolve()
_project_root = _current_file_dir.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

from camel_agents.tools import CatDTTools
from camel_agents.prompts import (
    AGENT1_STRUCTURE_INITIALIZER,
    AGENT2_ADSORPTION_PREDICTOR,
    AGENT3_RECONSTRUCTION_SIMULATOR,
    AGENT4_PATHWAY_DESIGNER,
    AGENT5_PATHWAY_VALIDATOR,
    AGENT6_NEB_RUNNER,
    AGENT7_REPORT_GENERATOR,
    MECHANISM_CONTEXT_ROUTING_AGENT,
    TaskPrompts,
)
from camel_agents.camel_model_backend import get_camel_model_backend
from core.pathway.pathway_predictor import CompletePathwayResult, PathwayStep, PathwayPredictor
from camel_agents.schemas import (
    AtomAddSpec,
    CatDTConfig,
    FacetResult,
    PathwayDesign,
    PathwayStepSpec,
    ReactionContext,
    ValidationReport,
    WorkflowState,
)
from camel_agents.runtime import CamelWorkflowAgent, Task, TaskOutput
from camel_agents.policy import EvolvablePolicy

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)



def _identify_top_two_layer_indices(
    atoms: Atoms,
    exclude_indices: Optional[List[int]] = None,
    layer_tol: float = 1.5,
) -> List[int]:
    """Return indices of atoms in the top two atomic layers of a slab.

    Layers are detected by z-coordinate clustering: atoms whose z is within
    ``layer_tol`` Å of the same value belong to the same layer. The two
    highest-z clusters are returned.

    The default ``layer_tol`` of 1.5 Å is chosen to handle stepped surfaces
    like fcc(211), where the surface has multiple staggered z-sublayers
    within a single physical "layer". A typical fcc interlayer spacing is
    ~2 Å, so 1.5 Å correctly merges staggered sub-rows while keeping true
    layers separate.

    Args:
        atoms: slab (may include adsorbate atoms).
        exclude_indices: atoms to exclude from layer detection (e.g. adsorbate
            atoms, which sit above the surface and would pollute the top
            layer). Returned indices never include these.
        layer_tol: Å tolerance for layer grouping.
    """
    exclude = set(exclude_indices or [])
    zs = [(i, float(atoms.positions[i, 2])) for i in range(len(atoms)) if i not in exclude]
    if not zs:
        return []
    zs.sort(key=lambda x: -x[1])
    top_z = zs[0][1]
    top_layer = [i for i, z in zs if top_z - z < layer_tol]
    remaining = [(i, z) for i, z in zs if i not in set(top_layer)]
    if remaining:
        second_z = remaining[0][1]
        second_layer = [i for i, z in remaining if second_z - z < layer_tol]
    else:
        second_layer = []
    return sorted(top_layer + second_layer)


class CatDTCamelWorkflow:
    """CatDT CAMEL workflow orchestrator composed from dedicated modules."""

    def __init__(
        self,
        tools: Optional[CatDTTools] = None,
        max_pathway_iters: int = 10,
        llm: Optional[Any] = None,
        initialize_agents: bool = True,
    ):
        self.tools = tools or CatDTTools(output_base_dir="output/catdt_workflow")
        self._max_pathway_iters = min(max(1, int(max_pathway_iters)), 10)
        self._last_iteration_count = 0

        self.agent1 = None
        self.agent2 = None
        self.agent3 = None
        self.agent4 = None
        self.agent5 = None
        self.agent6 = None
        self.agent7 = None
        self.mechanism_agent = None

        self.llm = llm
        self.camel_tools: Dict[str, FunctionTool] = {}
        self.current_strategy = "balanced"
        self.evolvable_policy: Optional[EvolvablePolicy] = None
        self._pathway_energy_helper: Optional[PathwayPredictor] = None

        if initialize_agents:
            if self.llm is None:
                self._init_llm()
            self._init_agents()

    def _make_memory(self, token_limit: int = 12000, window_size: int = 50) -> ChatHistoryMemory:
        token_counter = getattr(self.llm, "token_counter", None)
        if token_counter is None:
            raise RuntimeError("Workflow memory requires an initialized model backend with a token_counter.")
        context_creator = ScoreBasedContextCreator(
            token_counter=token_counter,
            token_limit=token_limit,
        )
        return ChatHistoryMemory(context_creator=context_creator, window_size=window_size)

    def _init_llm(self) -> None:
        try:
            from dotenv import load_dotenv
            load_dotenv()
        except Exception:
            pass

        llm_temperature = float(os.getenv("OPENAI_TEMPERATURE", "0.2"))
        llm_max_tokens = int(os.getenv("OPENAI_MAX_TOKENS", "4096"))
        llm_top_p = float(os.getenv("OPENAI_TOP_P", "0.95"))
        llm_presence_penalty = os.getenv("OPENAI_PRESENCE_PENALTY")
        llm_frequency_penalty = os.getenv("OPENAI_FREQUENCY_PENALTY")
        llm_seed = os.getenv("OPENAI_SEED")

        self.llm = get_camel_model_backend(
            model=os.getenv("OPENAI_MODEL", "gpt-5.4"),
            base_url=os.getenv("OPENAI_BASE_URL"),
            api_key=os.getenv("OPENAI_API_KEY"),
            temperature=llm_temperature,
            max_tokens=llm_max_tokens,
            top_p=llm_top_p,
            presence_penalty=(float(llm_presence_penalty) if llm_presence_penalty is not None else None),
            frequency_penalty=(float(llm_frequency_penalty) if llm_frequency_penalty is not None else None),
            seed=(int(llm_seed) if llm_seed is not None and str(llm_seed).strip() else None),
        )

    def _init_camel_tools(self) -> None:
        if self.tools is None:
            raise RuntimeError("CatDTTools must be initialized before CAMEL tool registry setup.")
        if not hasattr(self.tools, "to_camel_tools"):
            raise RuntimeError("CatDTTools must implement to_camel_tools() for CAMEL workflow.")

        registry = self.tools.to_camel_tools()
        if not isinstance(registry, dict) or not registry:
            raise RuntimeError("CatDTTools.to_camel_tools() returned an empty/invalid registry.")

        required = {
            "generate_surfaces",
            "predict_adsorption_sites",
            "simulate_surface_reconstruction",
            "run_neb_for_steps",
            "compute_adsorption_energies",
            "run_kmc_simulation",
            "generate_final_report",
            "retrieve_agent45_memento_cases",
            "record_agent45_memento_case",
        }
        missing = sorted(k for k in required if k not in registry)
        if missing:
            raise RuntimeError(f"CAMEL tool registry missing required entries: {missing}")

        self.camel_tools = registry

    def _init_agents(self) -> None:
        if self.llm is None:
            self._init_llm()

        self._init_camel_tools()
        default_step_timeout = float(os.getenv("CATDT_AGENT_STEP_TIMEOUT_SEC", "600"))
        agent45_step_timeout = float(os.getenv("CATDT_AGENT45_STEP_TIMEOUT_SEC", "1800"))
        default_message_window = int(os.getenv("CATDT_AGENT_MESSAGE_WINDOW", "20"))
        agent45_message_window = int(os.getenv("CATDT_AGENT45_MESSAGE_WINDOW", "8"))
        default_token_limit = int(os.getenv("CATDT_AGENT_TOKEN_LIMIT", "200000"))
        agent45_token_limit = int(os.getenv("CATDT_AGENT45_TOKEN_LIMIT", "200000"))
        agent45_memory_token_limit = int(os.getenv("CATDT_AGENT45_MEMORY_TOKEN_LIMIT", "6000"))
        agent45_memory_window = int(os.getenv("CATDT_AGENT45_MEMORY_WINDOW", "12"))
        agent45_summarize_threshold = int(
            os.getenv("CATDT_AGENT45_SUMMARIZE_THRESHOLD", str(max(12, agent45_message_window * 2)))
        )
        agent45_memory_agent4 = self._make_memory(
            token_limit=agent45_memory_token_limit,
            window_size=agent45_memory_window,
        )
        agent45_memory_agent5 = self._make_memory(
            token_limit=agent45_memory_token_limit,
            window_size=agent45_memory_window,
        )

        self.agent1 = CamelWorkflowAgent(
            role_name=AGENT1_STRUCTURE_INITIALIZER.role,
            goal=AGENT1_STRUCTURE_INITIALIZER.goal,
            backstory=AGENT1_STRUCTURE_INITIALIZER.backstory,
            model_backend=self.llm,
            memory=self._make_memory(),
            tools=[self.camel_tools["generate_surfaces"]],
            step_timeout=default_step_timeout,
            message_window_size=default_message_window,
            token_limit=default_token_limit,
        )
        self.agent2 = CamelWorkflowAgent(
            role_name=AGENT2_ADSORPTION_PREDICTOR.role,
            goal=AGENT2_ADSORPTION_PREDICTOR.goal,
            backstory=AGENT2_ADSORPTION_PREDICTOR.backstory,
            model_backend=self.llm,
            memory=self._make_memory(),
            tools=[self.camel_tools["predict_adsorption_sites"]],
            step_timeout=default_step_timeout,
            message_window_size=default_message_window,
            token_limit=default_token_limit,
        )
        self.agent3 = CamelWorkflowAgent(
            role_name=AGENT3_RECONSTRUCTION_SIMULATOR.role,
            goal=AGENT3_RECONSTRUCTION_SIMULATOR.goal,
            backstory=AGENT3_RECONSTRUCTION_SIMULATOR.backstory,
            model_backend=self.llm,
            memory=self._make_memory(),
            tools=[self.camel_tools["simulate_surface_reconstruction"]],
            step_timeout=default_step_timeout,
            message_window_size=default_message_window,
            token_limit=default_token_limit,
        )
        agent4_tools: List[FunctionTool] = [
            self.camel_tools["suggest_staged_positions"],
            self.camel_tools["validate_proposed_positions"],
        ]
        agent5_tools: List[FunctionTool] = []  # Agent5 uses pre-check results only, no tool calls
        self.agent4 = CamelWorkflowAgent(
            role_name=AGENT4_PATHWAY_DESIGNER.role,
            goal=AGENT4_PATHWAY_DESIGNER.goal,
            backstory=AGENT4_PATHWAY_DESIGNER.backstory,
            model_backend=self.llm,
            memory=agent45_memory_agent4,
            tools=agent4_tools,
            step_timeout=agent45_step_timeout,
            message_window_size=agent45_message_window,
            summarize_threshold=agent45_summarize_threshold,
            token_limit=agent45_token_limit,
        )
        self.agent5 = CamelWorkflowAgent(
            role_name=AGENT5_PATHWAY_VALIDATOR.role,
            goal=AGENT5_PATHWAY_VALIDATOR.goal,
            backstory=AGENT5_PATHWAY_VALIDATOR.backstory,
            model_backend=self.llm,
            memory=agent45_memory_agent5,
            tools=agent5_tools,
            step_timeout=agent45_step_timeout,
            message_window_size=agent45_message_window,
            summarize_threshold=agent45_summarize_threshold,
            token_limit=agent45_token_limit,
        )
        self.agent6 = CamelWorkflowAgent(
            role_name=AGENT6_NEB_RUNNER.role,
            goal=AGENT6_NEB_RUNNER.goal,
            backstory=AGENT6_NEB_RUNNER.backstory,
            model_backend=self.llm,
            memory=self._make_memory(),
            tools=[self.camel_tools["compute_adsorption_energies"], self.camel_tools["run_kmc_simulation"]],
            step_timeout=default_step_timeout,
            message_window_size=default_message_window,
            token_limit=default_token_limit,
        )
        self.agent7 = CamelWorkflowAgent(
            role_name=AGENT7_REPORT_GENERATOR.role,
            goal=AGENT7_REPORT_GENERATOR.goal,
            backstory=AGENT7_REPORT_GENERATOR.backstory,
            model_backend=self.llm,
            memory=self._make_memory(),
            tools=[self.camel_tools["generate_final_report"]],
            step_timeout=default_step_timeout,
            message_window_size=default_message_window,
            token_limit=default_token_limit,
        )
        # Mechanism search agent (between Agent3 and Agent4/5)
        mechanism_tools: List[FunctionTool] = [
            self.camel_tools["initialize_mechanism_context"],
            self.camel_tools["generate_candidate_steps"],
            self.camel_tools["run_mechanism_search"],
            self.camel_tools["extract_pathway_shortlist"],
        ]
        self.mechanism_agent = CamelWorkflowAgent(
            role_name=MECHANISM_CONTEXT_ROUTING_AGENT.role,
            goal=MECHANISM_CONTEXT_ROUTING_AGENT.goal,
            backstory=MECHANISM_CONTEXT_ROUTING_AGENT.backstory,
            model_backend=self.llm,
            memory=self._make_memory(),
            tools=mechanism_tools,
            step_timeout=float(os.getenv("CATDT_MECHANISM_AGENT_TIMEOUT_SEC", "600")),
            message_window_size=default_message_window,
            token_limit=default_token_limit,
        )

    def _run_task(self, task: Task, use_memory: bool = False) -> TaskOutput:
        _ = use_memory  # memory is managed by ChatAgent instances directly
        allow_tool_calls = True if task.allow_tool_calls is None else bool(task.allow_tool_calls)
        output: TaskOutput
        retry_done = False
        api_retries = 0
        max_api_retries = 5
        while True:
            try:
                output = task.agent.run(
                    task.description,
                    output_model=task.output_pydantic,
                    allow_tool_calls=bool(allow_tool_calls),
                )
                break
            except Exception as exc:
                err_text = str(exc)

                is_api_error = any(
                    marker in err_text
                    for marker in ["503", "502", "429", "temporarily unavailable", "rate limit", "overloaded", "算力紧张"]
                )
                if is_api_error and api_retries < max_api_retries:
                    api_retries += 1
                    wait_sec = min(30 * api_retries, 120)
                    logger.warning(
                        "API error (attempt %d/%d), retrying in %ds: %s",
                        api_retries, max_api_retries, wait_sec, err_text[:200],
                    )
                    import time
                    time.sleep(wait_sec)
                    continue

                recoverable_tool_transcript_error = (
                    allow_tool_calls
                    and ("tool_use" in err_text)
                    and ("tool_result" in err_text)
                    and (not retry_done)
                )
                if recoverable_tool_transcript_error:
                    retry_done = True
                    logger.warning(
                        "Recoverable tool transcript error for task; resetting agent dialogue and retrying once: %s",
                        err_text,
                    )
                    try:
                        chat_agent = getattr(task.agent, "chat_agent", None)
                        if chat_agent is not None and hasattr(chat_agent, "reset"):
                            chat_agent.reset()
                    except Exception as reset_exc:
                        logger.warning("Failed to reset chat agent after recoverable error: %s", reset_exc)
                    continue
                raise
        if task.result_handler:
            task.result_handler(output)
        return output

    def _persist_prompt_snapshot(
        self,
        state: WorkflowState,
        iteration: int,
        agent_name: str,
        prompt_text: str,
    ) -> None:
        try:
            surface = self._get_surface()
            if surface is None:
                return
            pathway_dir = surface.agent_dir("agent45")
            prompts_dir = surface.agent_dir("agent45", subdir="prompts")

            iter_file = prompts_dir / f"iter_{iteration:02d}_{agent_name}_prompt.txt"
            iter_file.write_text(str(prompt_text or ""), encoding="utf-8")

            merged_file = pathway_dir / f"{agent_name}_prompts_full.txt"
            with open(merged_file, "a", encoding="utf-8") as fh:
                fh.write(f"\n\n{'=' * 90}\n")
                fh.write(f"iteration={iteration}; agent={agent_name}\n")
                fh.write(f"{'=' * 90}\n")
                fh.write(str(prompt_text or ""))
                fh.write("\n")
        except Exception as exc:
            logger.debug("Failed to persist prompt snapshot for %s iter %d: %s", agent_name, iteration, exc)

    def _export_pathway_iteration_visuals(self, state: WorkflowState, iteration: int) -> None:
        """Dump every Stage B step's reactant/product as VASP + PNG.

        Resilient against per-step failures: VASP writes are attempted first
        (cheap, reliable); PNG rendering (tachyon headless) is attempted
        after and a failure in one endpoint does not stop subsequent writes.
        """
        if not state.step_structures:
            return

        surface = self._get_surface()
        if surface is None:
            return
        viz_dir = surface.agent_dir(
            "agent45", subdir=f"iter_{iteration:02d}_step_visualizations"
        )

        try:
            from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
            visualizer = CatalystSurfaceVisualizer(quality="low", auto_expand=False)
        except Exception as exc:
            logger.warning("Iter %d viz: visualizer init failed (%s). VASPs will still be dumped.", iteration, exc)
            visualizer = None

        vasp_ok = 0
        png_ok = 0
        png_fail = 0
        for step_idx, step in enumerate(state.step_structures, start=1):
            step_name = str(step.get("name", f"step_{step_idx}"))
            safe_name = re.sub(r"[^\w\-]+", "_", step_name)

            for endpoint in ("reactant", "product"):
                atoms = step.get(endpoint)
                if not isinstance(atoms, Atoms):
                    continue

                vasp_path = viz_dir / f"step{step_idx:02d}_{safe_name}_{endpoint}.vasp"
                try:
                    write(vasp_path, atoms)
                    vasp_ok += 1
                except Exception as exc:
                    logger.warning("Iter %d step %d %s VASP write failed: %s", iteration, step_idx, endpoint, exc)
                    continue

                if visualizer is None:
                    continue
                png_path = viz_dir / f"step{step_idx:02d}_{safe_name}_{endpoint}.png"
                try:
                    # Render from the written VASP to avoid stale per-atom arrays
                    # on in-memory Atoms objects after staged atom insert/remove.
                    viz_atoms = read(vasp_path)
                    visualizer.visualize_structure(
                        viz_atoms,
                        output_file=str(png_path),
                        title=f"Iter {iteration} Step {step_idx} {endpoint.capitalize()}",
                        show=False,
                    )
                    png_ok += 1
                except Exception as exc:
                    png_fail += 1
                    logger.warning(
                        "Iter %d step %d %s PNG render failed: %s",
                        iteration, step_idx, endpoint, exc,
                    )

        logger.info(
            "Iter %d viz: wrote %d VASP files, %d PNGs (%d render failures) → %s",
            iteration, vasp_ok, png_ok, png_fail, viz_dir,
        )

    def _summarize_iteration_history(self, state: WorkflowState, max_items: int = 5) -> str:
        if not state.iteration_history:
            return ""

        lines: List[str] = []
        for item in state.iteration_history[-max_items:]:
            idx = item.get("iteration", "?")
            status = item.get("status", "UNKNOWN")
            issue_count = item.get("issue_count", 0)
            issue_preview = item.get("issue_preview", "")
            lines.append(
                f"iter={idx}; status={status}; issue_count={issue_count}; issue_preview={issue_preview}"
            )

        return "\n".join(lines)

    @staticmethod
    def _truncate_prompt_block(text: Any, max_chars: int, block_name: str) -> str:
        content = str(text or "").strip()
        if max_chars <= 0 or len(content) <= max_chars:
            return content

        keep = max(256, int(max_chars) - 128)
        dropped = len(content) - keep
        return f"{content[:keep]}\n...(truncated {dropped} chars from {block_name})"

    def _reset_agent45_dialogue(self) -> None:
        for agent_name, agent in (("agent4", self.agent4), ("agent5", self.agent5)):
            if agent is None:
                continue
            try:
                chat_agent = getattr(agent, "chat_agent", None)
                if chat_agent is not None and hasattr(chat_agent, "reset"):
                    chat_agent.reset()
            except Exception as exc:
                logger.warning("Failed to reset %s dialogue: %s", agent_name, exc)

    # === Mechanism search stage helpers ===

    @staticmethod
    def _should_enable_mechanism_stage(config) -> bool:
        """Check if the mechanism search stage should run."""
        return bool(getattr(config, "enable_mechanism_search", False))

    def _run_mechanism_search_stage(
        self,
        state: WorkflowState,
        config,
    ) -> Optional[Dict[str, Any]]:
        """Run mechanism search via Agent M1 (a real CAMEL agent with tools).

        The agent analyses the user's reaction description, decides the
        exploration mode, calls mechanism tools, and returns evaluated pathways.
        """
        if self.mechanism_agent is None:
            logger.warning("Mechanism agent not initialized, skipping")
            return {"error": "mechanism_agent not initialized"}

        surface = self._get_surface()
        if surface is None:
            logger.warning("Mechanism agent: surface not set, skipping")
            return {"error": "surface not set"}

        # Build context for the agent prompt
        initial_state = ""
        target_state = ""
        known_intermediates: List[str] = []

        if state.reaction_context:
            intermediates = list(state.reaction_context.intermediates or [])
            if len(intermediates) >= 2:
                initial_state = intermediates[0]
                target_state = intermediates[-1]
                known_intermediates = intermediates[1:-1]

        if not initial_state or not target_state:
            return {"error": "Cannot determine initial/target from reaction context"}

        # Determine fixed site from current adsorbate position
        fixed_site_str = "null"
        if surface.adsorbate_indices and surface.reconstructed_surface_path:
            try:
                ads_atoms = read(surface.reconstructed_surface_path)
                ads_pos = ads_atoms.get_positions()[surface.adsorbate_indices]
                fixed_site = ads_pos[np.argmin(ads_pos[:, 2])].tolist()
                fixed_site_str = str([round(x, 2) for x in fixed_site])
            except Exception:
                pass

        surface_facet = ""
        if surface.surface_path:
            import re as _re
            m = _re.search(r"(\d{3,4})", str(surface.surface_path))
            if m:
                surface_facet = m.group(1)

        exploration_mode = getattr(config, "mechanism_exploration_mode", "agent_guided")
        known_str = ", ".join(known_intermediates) if known_intermediates else "(none — recommend pathways)"

        prompt = f"""Analyse this catalytic reaction and run mechanism search.

Reaction: {state.reaction_description}
Initial state: {initial_state}
Target state: {target_state}
Surface: {surface.surface_path or "unknown"}  Facet: {surface_facet}
Known intermediates: {known_str}
Fixed adsorption site: {fixed_site_str}
User-requested exploration mode: {exploration_mode}

Call tools in this order:
1. initialize_mechanism_context(initial_state="{initial_state}", target_state="{target_state}", \
surface_id="surface", surface_facet="{surface_facet}", \
surface_structure_path="{surface.reconstructed_surface_path or surface.surface_path or ""}", \
clean_slab_path="{surface.clean_slab_path or ""}", \
known_intermediates={known_intermediates}, \
fixed_adsorption_site={fixed_site_str})
2. run_mechanism_search(exploration_mode="{exploration_mode}", max_depth={config.mechanism_max_depth}, \
beam_width={config.mechanism_beam_width}, delta_keep={config.mechanism_delta_keep}, \
delta_prune={config.mechanism_delta_prune}, max_state_evaluations={config.mechanism_max_state_evaluations})
3. extract_pathway_shortlist(search_result=<result from step 2>, \
output_base_dir="{state.output_base_dir}", run_id="{state.run_id}")

After all tools complete, output a JSON summary of the retained pathways."""

        try:
            task = Task(
                description=prompt,
                expected_output="Mechanism search completed",
                agent=self.mechanism_agent,
                allow_tool_calls=True,
            )
            output = self._run_task(task)

            # Extract results from tool call records
            # The tools have already stored results in self.tools._mechanism_context etc.
            search_result = getattr(self.tools, "_mechanism_last_search_result", None)
            shortlist = getattr(self.tools, "_mechanism_last_shortlist", None)

            if search_result and shortlist:
                return {
                    "status": "OK",
                    "context": self.tools._mechanism_context,
                    "search_result": search_result,
                    "shortlist": shortlist,
                }

            # If agent didn't call all tools, fall back to direct tool calls
            logger.info("Mechanism agent did not complete all tool calls, running direct pipeline")
            return self._run_mechanism_search_direct(state, config)

        except Exception as exc:
            logger.warning("Mechanism agent failed: %s; falling back to direct pipeline", exc)
            return self._run_mechanism_search_direct(state, config)

    def _run_mechanism_search_direct(
        self,
        state: WorkflowState,
        config,
    ) -> Optional[Dict[str, Any]]:
        """Direct (non-agent) mechanism search pipeline as fallback."""
        try:
            from camel_agents.mechanism_search import MechanismStageRunner

            runner = MechanismStageRunner(tools=self.tools, use_agents=False)
            return runner.run_after_agent3(
                state=state, config=config, surface=self._get_surface(),
            )
        except Exception as exc:
            logger.warning("Direct mechanism search failed: %s", exc, exc_info=True)
            return {"error": str(exc)}

    def _load_shortlist_into_state(
        self,
        state: WorkflowState,
        mechanism_output: Dict[str, Any],
    ) -> None:
        """Store mechanism search results into WorkflowState.

        Critically: updates ``state.reaction_context.intermediates`` with the
        best pathway so that Agent4/5 can consume it directly.
        """
        state.mechanism_context = mechanism_output.get("context")
        search_result = mechanism_output.get("search_result", {})
        shortlist = mechanism_output.get("shortlist", {})

        state.mechanism_search_result = search_result

        all_pathways = search_result.get("pathways", [])
        state.candidate_pathways = all_pathways
        state.retained_pathways = [p for p in all_pathways if p.get("is_retained", True)]
        state.pruned_pathways = [p for p in all_pathways if not p.get("is_retained", True)]
        state.pathway_manifests = shortlist.get("manifest", {})

        # Update reaction_context.intermediates with the best pathway
        # so Agent4/5 can consume it directly via existing _run_pathway_iteration
        if state.retained_pathways and state.reaction_context:
            best_pathway = state.retained_pathways[0]
            best_intermediates = [
                s.get("species_label", "")
                for s in best_pathway.get("states", [])
                if s.get("species_label")
            ]
            if len(best_intermediates) >= 2:
                state.reaction_context.intermediates = best_intermediates
                logger.info(
                    "Updated intermediates from mechanism search: %s",
                    best_intermediates,
                )

    def _build_agent45_transition_signature(self, state: WorkflowState) -> List[str]:
        intermediates = self._get_working_intermediate_sequence(state)
        if len(intermediates) < 2:
            return []
        return [
            f"{str(intermediates[i]).strip()}->{str(intermediates[i + 1]).strip()}".lower()
            for i in range(len(intermediates) - 1)
        ]

    def _build_agent45_memento_query(self, state: WorkflowState, iteration: int) -> str:
        surface = self._get_surface()
        reaction_type = str(state.reaction_context.reaction_type if state.reaction_context else "").strip()
        constraints = str(state.reaction_context.constraints if state.reaction_context else "").strip()
        transitions = self._build_agent45_transition_signature(state)
        feedback = str(state.last_feedback or "").strip()
        if len(feedback) > 240:
            feedback = feedback[:240] + "..."
        query = (
            f"iteration={int(iteration)}; "
            f"reaction_type={reaction_type}; "
            f"transitions={transitions}; "
            f"constraints={constraints}; "
            f"feedback={feedback}"
        )
        # Add surface context so parametric retriever can match similar surfaces
        if surface is not None and surface.surface_indices and state.step_structures:
            first_atoms = (state.step_structures[0] or {}).get("reactant")
            if first_atoms is not None:
                import numpy as _np
                from collections import Counter
                all_surf = list(surface.surface_indices)
                positions = first_atoms.get_positions()
                syms = first_atoms.get_chemical_symbols()
                surf_z = positions[all_surf, 2]
                top_z = float(_np.max(surf_z))
                top2_mask = surf_z >= (top_z - 4.0)
                top2_indices = [all_surf[i] for i in range(len(all_surf)) if top2_mask[i]]
                elem_counts = dict(Counter(syms[i] for i in top2_indices))
                cell = first_atoms.cell
                query += (
                    f"; surface: cell=[{cell[0][0]:.1f},{cell[1][1]:.1f},{cell[2][2]:.1f}] "
                    f"top_z={top_z:.1f} composition={elem_counts}"
                )
        return query

    def _retrieve_agent45_memento_context(self, state: WorkflowState, iteration: int) -> str:
        if self.tools is None:
            state.memento_retrieval = {}
            state.memento_query_text = ""
            return "(memento unavailable: tools not initialized)"

        query_text = self._build_agent45_memento_query(state=state, iteration=iteration)
        state.memento_query_text = query_text
        top_k = max(1, int(os.getenv("CATDT_AGENT45_MEMENTO_TOP_K", "6")))
        min_score = float(os.getenv("CATDT_AGENT45_MEMENTO_MIN_SCORE", "0.05"))
        include_negative = str(os.getenv("CATDT_AGENT45_MEMENTO_INCLUDE_NEGATIVE", "1")).strip().lower() not in {
            "0",
            "false",
            "no",
        }
        try:
            retrieval = self.tools.retrieve_agent45_memento_cases(
                query_text=query_text,
                top_k=top_k,
                include_negative=include_negative,
                min_score=min_score,
            )
            state.memento_retrieval = dict(retrieval or {})
            return str(state.memento_retrieval.get("prompt_block", "") or "(no memento cases)")
        except Exception as exc:
            logger.warning("Agent45 memento retrieval failed in iteration %d: %s", iteration, exc)
            state.memento_retrieval = {
                "query": query_text,
                "cases": [],
                "positive_cases": [],
                "negative_cases": [],
                "prompt_block": f"(memento retrieval failed: {exc})",
            }
            return str(state.memento_retrieval.get("prompt_block", ""))

    def _retrieve_agent45_knowledge_context(self, state: WorkflowState, iteration: int) -> str:
        if self.tools is None:
            state.knowledge_retrieval = {}
            return "(knowledge bank unavailable: tools not initialized)"

        query_text = state.memento_query_text or self._build_agent45_memento_query(state=state, iteration=iteration)
        top_k = max(1, int(os.getenv("CATDT_AGENT45_KNOWLEDGE_TOP_K", "4")))
        min_score = float(os.getenv("CATDT_AGENT45_KNOWLEDGE_MIN_SCORE", "0.03"))
        try:
            retrieval = self.tools.retrieve_agent45_knowledge_items(
                query_text=query_text,
                top_k=top_k,
                min_score=min_score,
            )
            state.knowledge_retrieval = dict(retrieval or {})
            return str(state.knowledge_retrieval.get("prompt_block", "") or "(no knowledge items)")
        except Exception as exc:
            logger.warning("Agent45 knowledge retrieval failed in iteration %d: %s", iteration, exc)
            state.knowledge_retrieval = {
                "query": query_text,
                "items": [],
                "prompt_block": f"(knowledge retrieval failed: {exc})",
            }
            return str(state.knowledge_retrieval.get("prompt_block", ""))

    def _retrieve_agent45_skill_context(self, state: WorkflowState, iteration: int) -> str:
        if self.tools is None:
            state.skill_retrieval = {}
            return "(skill bank unavailable: tools not initialized)"

        query_text = state.memento_query_text or self._build_agent45_memento_query(state=state, iteration=iteration)
        top_k = max(1, int(os.getenv("CATDT_AGENT45_SKILL_TOP_K", "4")))
        min_score = float(os.getenv("CATDT_AGENT45_SKILL_MIN_SCORE", "0.03"))
        try:
            retrieval = self.tools.retrieve_agent45_skill_items(
                query_text=query_text,
                top_k=top_k,
                min_score=min_score,
            )
            state.skill_retrieval = dict(retrieval or {})
            return str(state.skill_retrieval.get("prompt_block", "") or "(no skill items)")
        except Exception as exc:
            logger.warning("Agent45 skill retrieval failed in iteration %d: %s", iteration, exc)
            state.skill_retrieval = {
                "query": query_text,
                "items": [],
                "prompt_block": f"(skill retrieval failed: {exc})",
            }
            return str(state.skill_retrieval.get("prompt_block", ""))

    def _persist_agent45_memento_case(
        self,
        state: WorkflowState,
        iteration: int,
        reward: Optional[float] = None,
    ) -> None:
        if self.tools is None or state.reaction_context is None or state.validation_report is None:
            return

        surface = self._get_surface()
        design_outline: List[Dict[str, Any]] = []
        if state.pathway_design and state.pathway_design.steps:
            for spec in state.pathway_design.steps:
                step_outline: Dict[str, Any] = {
                    "step_name": spec.step_name,
                    "reactant_formula": spec.reactant_formula,
                    "product_formula": spec.product_formula,
                    "atoms_to_add_count": len(spec.atoms_to_add),
                    "atoms_to_remove_count": len(spec.atoms_to_remove),
                }
                # Record actual positions Agent4 chose for added atoms
                if spec.atoms_to_add:
                    step_outline["atoms_to_add"] = [
                        {
                            "element": a.element or a.species or a.symbol,
                            "position": a.position,
                            "reactant_position": a.reactant_position,
                            "product_position": a.product_position,
                        }
                        for a in spec.atoms_to_add
                    ]
                design_outline.append(step_outline)

        # Build surface top-2-layer snapshot (same for all steps, stored once)
        import numpy as _np
        surface_snapshot: Optional[Dict[str, Any]] = None
        if surface is not None and surface.surface_indices and state.step_structures:
            # Use the first step's reactant to get surface coordinates
            first_atoms = (state.step_structures[0] or {}).get("reactant")
            if first_atoms is not None:
                all_surf = list(surface.surface_indices)
                positions = first_atoms.get_positions()
                syms = first_atoms.get_chemical_symbols()
                surf_z = positions[all_surf, 2]
                surface_top_z = float(_np.max(surf_z))
                # Top 2 layers: atoms within 4.0 Å of surface_top_z
                top2_mask = surf_z >= (surface_top_z - 4.0)
                top2_indices = [all_surf[i] for i in range(len(all_surf)) if top2_mask[i]]
                cell = first_atoms.cell
                surface_snapshot = {
                    "cell": [round(float(cell[i][i]), 3) for i in range(3)],
                    "surface_top_z": round(surface_top_z, 3),
                    "top2_layers": [
                        {"elem": syms[i], "xyz": [round(float(c), 3) for c in positions[i]]}
                        for i in top2_indices
                    ],
                }

        # Build per-step adsorbate geometry snapshot from step_structures
        # Stores core (original) and staged (Agent4-added) atoms separately
        step_geometry_snapshot: List[Dict[str, Any]] = []
        for step in (state.step_structures or []):
            snap: Dict[str, Any] = {"name": step.get("name", "")}
            for endpoint in ("reactant", "product"):
                atoms_obj = step.get(endpoint)
                if atoms_obj is None:
                    continue
                ads_indices = list(step.get(f"{endpoint}_adsorbate_indices", []))
                staged_indices = list(step.get(f"{endpoint}_staged_indices", []))
                syms = atoms_obj.get_chemical_symbols()
                pos = atoms_obj.get_positions()
                staged_set = set(staged_indices)
                core_ads = [
                    {"elem": syms[i], "xyz": [round(float(c), 3) for c in pos[i]]}
                    for i in ads_indices if i not in staged_set and 0 <= i < len(atoms_obj)
                ]
                staged_atoms = []
                core_positions = _np.array([a["xyz"] for a in core_ads]) if core_ads else _np.zeros((0, 3))
                for si in staged_indices:
                    if si < 0 or si >= len(atoms_obj):
                        continue
                    sp = pos[si]
                    info: Dict[str, Any] = {
                        "elem": syms[si],
                        "xyz": [round(float(c), 3) for c in sp],
                    }
                    if len(core_positions) > 0:
                        dists = _np.linalg.norm(core_positions - sp, axis=1)
                        info["dist_to_nearest_ads"] = round(float(_np.min(dists)), 2)
                    staged_atoms.append(info)
                snap[endpoint] = {
                    "core": core_ads,
                    "staged": staged_atoms,
                }
            step_geometry_snapshot.append(snap)

        neb_summary: Dict[str, Any] = {}
        for key, value in (state.neb_results or {}).items():
            if isinstance(value, dict):
                neb_summary[str(key)] = {
                    "error": value.get("error"),
                    "converged": value.get("converged"),
                }
            else:
                neb_summary[str(key)] = {
                    "converged": getattr(value, "converged", None),
                    "activation_energy_forward": getattr(value, "activation_energy_forward", None),
                    "reaction_energy": getattr(value, "reaction_energy", None),
                }

        try:
            self.tools.record_agent45_memento_case(
                run_id=state.run_id,
                reaction_description=state.reaction_description,
                reaction_type=state.reaction_context.reaction_type,
                intermediates=list(state.reaction_context.intermediates or []),
                transition_signature=self._build_agent45_transition_signature(state),
                iteration=int(iteration),
                validation_status=state.validation_report.status,
                issues=list(state.validation_report.issues or []),
                feedback=str(state.validation_report.feedback or ""),
                design_outline=design_outline,
                energy_gate_report=dict(state.energy_gate_report or {}),
                neb_summary=neb_summary,
                reward=reward,
                extra_query=str(state.memento_query_text or ""),
                step_geometry_snapshot=step_geometry_snapshot,
                surface_snapshot=surface_snapshot,
            )
        except Exception as exc:
            logger.warning("Failed to persist Agent45 memento case at iteration %d: %s", iteration, exc)

    def _strategy_hint(self) -> str:
        hints = {
            "balanced": "平衡探索与稳健性，优先最小必要改动与可计算性。",
            "conservative": "优先保持几何连续与低风险映射，减少激进结构修改。",
            "exploratory": "允许在满足可计算前提下尝试替代步骤设计。",
        }
        return hints.get(self.current_strategy, hints["balanced"])

    def _create_reaction_context_task(self, state: WorkflowState) -> Task:
        description = TaskPrompts.reaction_context_parsing(state.reaction_description)

        def _apply_reaction_context(output: TaskOutput):
            payload = output.json_dict or {}
            context = output.pydantic
            if context is None:
                try:
                    context = ReactionContext(**payload)
                except Exception:
                    raise RuntimeError(
                        "Agent4 reaction-context parsing failed: LLM did not return valid ReactionContext JSON."
                    )

            # ── Universal consistency check ──
            # The intermediates chain MUST start with initial_adsorbate.
            # If the LLM omitted it (e.g. parsed intermediates as downstream
            # products only), prepend initial_adsorbate so subsequent stages
            # have a valid starting point. This is a structural invariant of
            # the Stage A / Stage B pipeline, independent of chemistry.
            initial = (context.initial_adsorbate or "").strip()
            chain = list(context.intermediates or [])
            if initial and (not chain or chain[0] != initial):
                if initial in chain:
                    # Initial appears mid-chain — truncate to start at it
                    idx = chain.index(initial)
                    chain = chain[idx:]
                    logger.warning(
                        "reaction_context: intermediates did not start at initial "
                        "'%s'; truncated chain to start there: %s",
                        initial, chain,
                    )
                else:
                    chain = [initial] + chain
                    logger.warning(
                        "reaction_context: intermediates did not include initial "
                        "'%s'; prepended: %s",
                        initial, chain,
                    )
                context.intermediates = chain

            state.reaction_context = context
            state.intermediates = list(context.intermediates)

        return Task(
            description=description,
            expected_output="JSON with reaction_type, initial_adsorbate, intermediates, constraints",
            agent=self.agent4,
            output_pydantic=ReactionContext,
            result_handler=_apply_reaction_context,
        )

    @staticmethod
    def _parse_pathway_step_formulas(step_name: str) -> Tuple[str, str]:
        raw = str(step_name or "").strip()
        if not raw:
            return "", ""
        if "_to_" in raw:
            left, right = raw.split("_to_", 1)
            return left.strip(), right.strip()
        if "->" in raw:
            left, right = raw.split("->", 1)
            return left.strip(), right.strip()
        return "", ""

    @staticmethod
    def _coerce_step_atoms_to_add(raw_step: Dict[str, Any]) -> List[Dict[str, Any]]:
        normalized: List[Dict[str, Any]] = []

        raw_adds = raw_step.get("atoms_to_add", [])
        if isinstance(raw_adds, list):
            for item in raw_adds:
                if not isinstance(item, dict):
                    continue

                species = str(
                    item.get("species") or item.get("element") or item.get("symbol") or ""
                ).strip()
                if not species:
                    continue

                entry: Dict[str, Any] = {
                    "species": species,
                    "reason": str(item.get("reason") or "").strip(),
                }

                for key in ("position", "reactant_position", "product_position"):
                    value = item.get(key)
                    if isinstance(value, list) and len(value) == 3:
                        try:
                            entry[key] = [float(value[0]), float(value[1]), float(value[2])]
                        except Exception:
                            continue

                if "position" not in entry:
                    if "reactant_position" in entry:
                        entry["position"] = list(entry["reactant_position"])
                    elif "product_position" in entry:
                        entry["position"] = list(entry["product_position"])

                normalized.append(entry)

        if normalized:
            return normalized

        staged_atoms = raw_step.get("staged_atoms", [])
        if not isinstance(staged_atoms, list):
            return normalized

        for item in staged_atoms:
            if not isinstance(item, dict):
                continue

            species = str(
                item.get("species") or item.get("element") or item.get("symbol") or ""
            ).strip()
            position = item.get("position")
            if not species or not (isinstance(position, list) and len(position) == 3):
                continue

            try:
                pos = [float(position[0]), float(position[1]), float(position[2])]
            except Exception:
                continue

            side = str(item.get("side") or "").strip().lower()
            entry = {
                "species": species,
                "position": pos,
                "reason": f"recovered_from_staged_atoms:{side or 'unspecified'}",
            }
            if side == "reactant":
                entry["reactant_position"] = list(pos)
            elif side == "product":
                entry["product_position"] = list(pos)

            normalized.append(entry)

        return normalized

    def _coerce_pathway_design_payload(
        self,
        state: WorkflowState,
        payload: Dict[str, Any],
    ) -> Dict[str, Any]:
        payload = dict(payload or {})
        context_seq = (
            list(state.reaction_context.intermediates)
            if state.reaction_context is not None
            else list(state.intermediates or [])
        )

        pathway_name = str(payload.get("pathway_name") or "").strip() or "repaired_pathway"
        description_text = str(payload.get("description") or "").strip()
        if not description_text:
            summary = payload.get("summary")
            summary_bits: List[str] = []
            if isinstance(summary, dict):
                validated = summary.get("validated_steps")
                failed = summary.get("failed_steps")
                if validated is not None:
                    summary_bits.append(f"validated_steps={validated}")
                if failed is not None:
                    summary_bits.append(f"failed_steps={failed}")
            status = str(payload.get("status") or "").strip()
            description_parts = [part for part in [pathway_name, status] if part]
            if summary_bits:
                description_parts.append(", ".join(summary_bits))
            description_text = " | ".join(description_parts) or "Recovered Agent4 pathway design"

        overall_reaction = str(payload.get("overall_reaction") or "").strip()
        if not overall_reaction:
            if len(context_seq) >= 2:
                overall_reaction = " -> ".join(context_seq)
            else:
                overall_reaction = str(state.reaction_description or pathway_name)

        normalized_steps: List[Dict[str, Any]] = []
        raw_steps = payload.get("steps", [])
        if isinstance(raw_steps, list):
            for idx, item in enumerate(raw_steps):
                if not isinstance(item, dict):
                    continue

                step_name = str(item.get("step_name") or item.get("name") or "").strip()
                reactant_formula = str(item.get("reactant_formula") or "").strip()
                product_formula = str(item.get("product_formula") or "").strip()

                parsed_reactant, parsed_product = self._parse_pathway_step_formulas(step_name)
                if not reactant_formula:
                    reactant_formula = parsed_reactant
                if not product_formula:
                    product_formula = parsed_product

                if idx < len(context_seq) - 1:
                    reactant_formula = reactant_formula or str(context_seq[idx])
                    product_formula = product_formula or str(context_seq[idx + 1])

                if not step_name and reactant_formula and product_formula:
                    step_name = f"{reactant_formula}_to_{product_formula}"

                confidence = item.get("confidence", 0.0)
                try:
                    confidence_value = float(confidence or 0.0)
                except Exception:
                    confidence_value = 0.0
                if confidence_value <= 0.0:
                    validation = item.get("validation")
                    if isinstance(validation, dict):
                        interp = str(validation.get("interpolation") or "").strip().lower()
                        ok = bool(validation.get("ok"))
                        if not ok or interp == "collision":
                            confidence_value = 0.2
                        elif interp == "warning":
                            confidence_value = 0.6
                        else:
                            confidence_value = 0.8

                raw_modify = item.get("atoms_to_modify", [])
                atoms_to_modify = [
                    dict(mod)
                    for mod in raw_modify
                    if isinstance(raw_modify, list) and isinstance(mod, dict) and "index" in mod
                ]
                raw_remove = item.get("atoms_to_remove", [])

                normalized_steps.append(
                    {
                        "step_name": step_name or f"step_{idx + 1}",
                        "reactant_formula": reactant_formula or "*UNKNOWN",
                        "product_formula": product_formula or "*UNKNOWN",
                        "reaction_type": str(item.get("reaction_type") or "").strip() or "elementary_transition",
                        "atoms_to_add": self._coerce_step_atoms_to_add(item),
                        "atoms_to_modify": atoms_to_modify,
                        "atoms_to_remove": list(raw_remove) if isinstance(raw_remove, list) else [],
                        "confidence": confidence_value,
                    }
                )

        confidence_value = payload.get("confidence", 0.0)
        try:
            top_confidence = float(confidence_value or 0.0)
        except Exception:
            top_confidence = 0.0
        if top_confidence <= 0.0 and normalized_steps:
            top_confidence = float(
                sum(float(step.get("confidence", 0.0) or 0.0) for step in normalized_steps)
                / len(normalized_steps)
            )

        return {
            "pathway_name": pathway_name,
            "description": description_text,
            "overall_reaction": overall_reaction,
            "steps": normalized_steps,
            "confidence": top_confidence,
        }

    def _create_pathway_design_task(
        self,
        state: WorkflowState,
        structure: Atoms,
        preplaced_products: Optional[Dict[str, Dict[str, Any]]] = None,
        extra_feedback: str = "",
    ) -> Task:
        surface = self._get_surface()
        reaction_context = {
            "reaction_type": state.reaction_context.reaction_type if state.reaction_context else "",
            "intermediates": state.reaction_context.intermediates if state.reaction_context else [],
            "constraints": state.reaction_context.constraints if state.reaction_context else "",
        }

        # Surface info: include compact descriptors + full POSCAR text (LLM cannot read local files).
        cell = structure.cell
        ads_for_top_z = list(surface.adsorbate_indices) if surface is not None else []
        surface_top_z = self._estimate_surface_top_z(structure, ads_for_top_z)
        surface_element = structure[0].symbol if len(structure) > 0 else "unknown"
        surface_poscar = self._atoms_to_poscar_text(structure)
        n_surf = len(surface.surface_indices) if surface is not None else 0
        surface_info = (
            f"cell_a={cell[0][0]:.3f}, cell_b={cell[1][1]:.3f}, cell_c={cell[2][2]:.3f} Å\n"
            f"surface_top_z={surface_top_z:.3f} Å\n"
            f"n_surface_atoms={n_surf}\n"
            f"surface_element={surface_element}\n"
            f"## SURFACE_POSCAR\n{surface_poscar}"
        )

        # Add top-layer surface atom positions so LLM can find adsorption sites
        if surface is not None and surface.surface_indices:
            import numpy as np
            all_surf = list(surface.surface_indices)
            positions = structure.get_positions()
            surf_z = positions[all_surf, 2]
            # Top layer: atoms within 2.0 Å of surface_top_z
            # (needs 2.0 Å because VSSR-MC reconstruction can add atoms
            #  ~1.5 Å above the original slab top layer)
            top_layer_mask = surf_z >= (surface_top_z - 2.0)
            top_layer_indices = [all_surf[i] for i in range(len(all_surf)) if top_layer_mask[i]]
            if top_layer_indices:
                surface_info += f"\n\ntop_layer_surface_atoms (use these to find hollow/bridge adsorption sites):"
                for idx in top_layer_indices:
                    x, y, z = positions[idx]
                    surface_info += f"\n  {idx}: {structure[idx].symbol}  {x:.4f}  {y:.4f}  {z:.4f}"

        # Per-step adsorbate coordinates (compact format from build_steps_payload)
        baseline_steps = "(none)"
        if state.tool_baseline_steps:
            baseline_steps = self.tools.build_steps_payload(self, state.tool_baseline_steps)

        # Build feedback: Agent5 text + previous round's adsorbate coords
        feedback = state.last_feedback or ""
        if str(extra_feedback or "").strip():
            feedback = (
                f"{feedback}\n\n【系统校验反馈】\n{str(extra_feedback).strip()}".strip()
            )
        if feedback and state.step_structures:
            prev_adsorbate_summary = self.tools.build_steps_payload(self, state.step_structures)
            feedback = f"{feedback}\n\n【前一轮各步吸附分子坐标】\n{prev_adsorbate_summary}"

        feedback_limit = int(os.getenv("CATDT_AGENT45_FEEDBACK_MAX_CHARS", "6000"))
        history_limit = int(os.getenv("CATDT_AGENT45_HISTORY_MAX_CHARS", "1200"))
        memento_limit = int(os.getenv("CATDT_AGENT45_MEMENTO_MAX_CHARS", "2000"))
        knowledge_limit = int(os.getenv("CATDT_AGENT45_KNOWLEDGE_MAX_CHARS", "2000"))
        skill_limit = int(os.getenv("CATDT_AGENT45_SKILL_MAX_CHARS", "2000"))

        feedback = self._truncate_prompt_block(feedback, feedback_limit, "agent4_feedback")
        history_text = self._truncate_prompt_block(
            self._summarize_iteration_history(state),
            history_limit,
            "agent4_history",
        )
        memento_context = self._truncate_prompt_block(
            (state.memento_retrieval or {}).get("prompt_block", ""),
            memento_limit,
            "agent4_memento",
        )
        knowledge_context = self._truncate_prompt_block(
            (state.knowledge_retrieval or {}).get("prompt_block", ""),
            knowledge_limit,
            "agent4_knowledge",
        )
        skill_context = self._truncate_prompt_block(
            (state.skill_retrieval or {}).get("prompt_block", ""),
            skill_limit,
            "agent4_skill",
        )

        description = TaskPrompts.pathway_design(
            reaction_context=reaction_context,
            surface_info=surface_info,
            baseline_steps=baseline_steps,
            memento_context=memento_context,
            knowledge_context=knowledge_context,
            skill_context=skill_context,
            feedback=feedback,
            history=history_text,
        )
        description += f"\n\n【当前策略提示】{self._strategy_hint()}"

        def _apply_pathway_design(output: TaskOutput):
            payload = output.json_dict or {}
            design = output.pydantic
            if design is None:
                if (not payload) and getattr(output, "raw", None):
                    raw_text = str(output.raw or "").strip()
                    candidate = raw_text
                    fenced = re.search(r"```json\s*(\{.*?\})\s*```", raw_text, flags=re.DOTALL | re.IGNORECASE)
                    if fenced:
                        candidate = fenced.group(1)
                    else:
                        loose = re.search(r"\{.*\}", raw_text, flags=re.DOTALL)
                        if loose:
                            candidate = loose.group(0)
                    try:
                        payload = json.loads(candidate)
                    except Exception:
                        payload = payload or {}

                try:
                    design = PathwayDesign(**self._coerce_pathway_design_payload(state, payload))
                except Exception as exc:
                    raw_text = str(output.raw or "").strip()
                    repair_prompt = (
                        "Convert the draft below into ONE valid PathwayDesign JSON object.\n"
                        "Requirements: JSON only, no markdown, no explanations.\n"
                        "Fields required: pathway_name, description, overall_reaction, steps, confidence.\n"
                        "Each step requires: step_name, reactant_formula, product_formula, reaction_type, "
                        "atoms_to_add, atoms_to_modify, atoms_to_remove, confidence.\n"
                        "Each atoms_to_add item requires: species, position [x,y,z], reason.\n\n"
                        "DRAFT_START\n"
                        f"{raw_text[:4000]}\n"
                        "DRAFT_END"
                    )
                    try:
                        repair_output = self.agent4.run(
                            repair_prompt,
                            output_model=PathwayDesign,
                            allow_tool_calls=False,
                        )
                        repaired_design = repair_output.pydantic
                        if repaired_design is None:
                            repair_payload = repair_output.json_dict or {}
                            if (not repair_payload) and getattr(repair_output, "raw", None):
                                repair_raw = str(repair_output.raw or "").strip()
                                repair_candidate = repair_raw
                                repair_fenced = re.search(
                                    r"```json\s*(\{.*?\})\s*```",
                                    repair_raw,
                                    flags=re.DOTALL | re.IGNORECASE,
                                )
                                if repair_fenced:
                                    repair_candidate = repair_fenced.group(1)
                                else:
                                    repair_loose = re.search(r"\{.*\}", repair_raw, flags=re.DOTALL)
                                    if repair_loose:
                                        repair_candidate = repair_loose.group(0)
                                try:
                                    repair_payload = json.loads(repair_candidate)
                                except Exception:
                                    repair_payload = repair_payload or {}
                            repaired_design = PathwayDesign(**self._coerce_pathway_design_payload(state, repair_payload))
                        design = repaired_design
                    except Exception as repair_exc:
                        raise RuntimeError(
                            "Agent4 pathway-design parsing failed: LLM did not return valid PathwayDesign JSON. "
                            f"parse_error={exc}; repair_error={repair_exc}; raw_output={raw_text}"
                        )

            design = self._align_pathway_design_to_context(state, design)
            state.pathway_design = design

        return Task(
            description=description,
            expected_output="JSON pathway design with step-level atom operations",
            agent=self.agent4,
            output_pydantic=PathwayDesign,
            allow_tool_calls=True,
            result_handler=_apply_pathway_design,
        )

    def _create_validation_task(self, state: WorkflowState) -> Task:
        programmatic_check = self._agent5_programmatic_validation(state.step_structures)
        energy_gate = dict(state.energy_gate_report or {})
        design_steps = state.pathway_design.steps if state.pathway_design else []
        structures_payload = self._build_agent5_structures_payload(state)
        fatal_issues = list(programmatic_check.get("fatal_issues", []))
        warning_issues = list(programmatic_check.get("warning_issues", []))
        energy_fatal = list(energy_gate.get("fatal_issues", []))
        energy_warning = list(energy_gate.get("warning_issues", []))
        for item in energy_fatal:
            if item not in fatal_issues:
                fatal_issues.append(item)
        for item in energy_warning:
            if item not in warning_issues:
                warning_issues.append(item)

        precheck_summary = str(programmatic_check.get("summary", "(none)") or "(none)")
        energy_summary = str(energy_gate.get("summary", "") or "").strip()
        if energy_summary:
            precheck_summary = f"{precheck_summary}\n[UMA Energy Gate]\n{energy_summary}"

        step_briefs: List[str] = []
        for idx, struct_step in enumerate(state.step_structures, start=1):
            spec = design_steps[idx - 1] if (idx - 1) < len(design_steps) else None
            adds = len(spec.atoms_to_add) if spec else 0
            rems = len(spec.atoms_to_remove) if spec else 0
            step_briefs.append(
                f"step={idx}:{struct_step.get('name', '')}; "
                f"reactant={struct_step.get('reactant_formula', '')}; "
                f"product={struct_step.get('product_formula', '')}; "
                f"adds={adds}; removes={rems}"
            )

        memento_limit = int(os.getenv("CATDT_AGENT45_MEMENTO_MAX_CHARS", "2000"))
        knowledge_limit = int(os.getenv("CATDT_AGENT45_KNOWLEDGE_MAX_CHARS", "2000"))
        skill_limit = int(os.getenv("CATDT_AGENT45_SKILL_MAX_CHARS", "2000"))
        memento_context = self._truncate_prompt_block(
            (state.memento_retrieval or {}).get("prompt_block", ""),
            memento_limit,
            "agent5_memento",
        )
        knowledge_context = self._truncate_prompt_block(
            (state.knowledge_retrieval or {}).get("prompt_block", ""),
            knowledge_limit,
            "agent5_knowledge",
        )
        skill_context = self._truncate_prompt_block(
            (state.skill_retrieval or {}).get("prompt_block", ""),
            skill_limit,
            "agent5_skill",
        )

        description = TaskPrompts.pathway_validation(
            steps_details="\n".join(step_briefs) if step_briefs else "(none)",
            structures_payload=structures_payload,
            precheck_summary=precheck_summary,
            precheck_issues=(
                "\n".join(
                    [f"[FATAL] {item}" for item in fatal_issues]
                    + [f"[WARNING] {item}" for item in warning_issues]
                )
                if (fatal_issues or warning_issues)
                else "(none)"
            ),
            memento_context=memento_context,
            knowledge_context=knowledge_context,
            skill_context=skill_context,
        )

        def _apply_validation(output: TaskOutput):
            payload = output.json_dict or {}
            # LLM may return issues as List[dict] instead of List[str]; normalize
            if "issues" in payload and isinstance(payload["issues"], list):
                normalized = []
                for item in payload["issues"]:
                    if isinstance(item, dict):
                        # Convert structured issue dict to readable string
                        parts = []
                        if "severity" in item:
                            parts.append(f"[{item['severity']}]")
                        if "step" in item:
                            parts.append(f"Step {item['step']}:")
                        if "issue" in item:
                            parts.append(str(item["issue"]))
                        elif "description" in item:
                            parts.append(str(item["description"]))
                        elif "message" in item:
                            parts.append(str(item["message"]))
                        normalized.append(" ".join(parts) if parts else str(item))
                    else:
                        normalized.append(str(item))
                payload["issues"] = normalized
            if "feedback" in payload and not isinstance(payload.get("feedback"), str):
                try:
                    payload["feedback"] = json.dumps(payload["feedback"], ensure_ascii=False)
                except Exception:
                    payload["feedback"] = str(payload.get("feedback"))
            report = output.pydantic
            if report is None:
                try:
                    report = ValidationReport(**payload)
                except Exception:
                    raw_text = str(output.raw or "").strip()
                    raise RuntimeError(
                        "Agent5 validation parsing failed: LLM did not return valid ValidationReport JSON. "
                        f"raw_output={raw_text}"
                    )
            # Extra safety: normalize report.issues if pydantic auto-parsed but issues contain non-str
            if report.issues:
                report.issues = [str(i) if not isinstance(i, str) else i for i in report.issues]

            combined_fatal_issues = list(fatal_issues)
            combined_warning_issues = list(warning_issues)

            issues = list(report.issues or [])
            for item in combined_fatal_issues + combined_warning_issues:
                if item not in issues:
                    issues.append(item)

            llm_status = str(report.status or "FAIL").upper()
            llm_fail = llm_status != "PASS"
            status = "FAIL" if combined_fatal_issues else "PASS"
            if combined_fatal_issues:
                gate_note = str(programmatic_check.get("feedback", "")).strip() or (
                    "Programmatic validation found endpoint issues that block NEB mapping; "
                    "please revise atom additions/removals while preserving immutable core coordinates."
                )
                if energy_summary:
                    gate_note = f"{gate_note}\nUMA energy gate detected endpoint instability; revise staged coordinates."
                feedback_text = str(report.feedback or "").strip()
                report.feedback = f"{feedback_text}\n{gate_note}".strip()
            elif llm_fail:
                advisory = (
                    "LLM returned non-PASS but no programmatic fatal issue was found; "
                    "treating these as advisory feedback."
                )
                feedback_text = str(report.feedback or "").strip()
                report.feedback = f"{feedback_text}\n{advisory}".strip()
                issues = [f"[ADVISORY] {item}" for item in issues]

            report.status = status
            report.issues = issues
            state.validation_report = report
            # Prefix feedback with format reminder so Agent4 doesn't drift
            if status == "FAIL":
                format_reminder = (
                    "IMPORTANT: Re-generate ALL atoms_to_add with corrected [x,y,z] positions. "
                    "atoms_to_modify and atoms_to_remove MUST remain empty arrays []. "
                    "Do NOT attempt to modify existing atoms — only add new ones.\n\n"
                )
                report.feedback = format_reminder + str(report.feedback or "")
            state.last_feedback = report.feedback

        return Task(
            description=description,
            expected_output="JSON validation report",
            agent=self.agent5,
            output_pydantic=ValidationReport,
            allow_tool_calls=False,
            result_handler=_apply_validation,
        )

    def _validate_pathway_with_agent5(self, state: WorkflowState) -> Dict[str, Any]:
        report = state.validation_report
        if not report:
            return {"status": "FAIL", "issues": ["No validation report"], "feedback": ""}
        return {
            "status": report.status,
            "issues": report.issues,
            "feedback": report.feedback,
        }

    def _find_neb_result_for_transition(
        self,
        state: WorkflowState,
        reactant: str,
        product: str,
        step_index: int,
        step_name: Optional[str] = None,
    ) -> Any:
        if not state.neb_results:
            return None

        candidates: List[str] = []

        def _add(name: str) -> None:
            if name and name not in candidates:
                candidates.append(name)

        if step_name:
            _add(step_name)

        reactant_raw = str(reactant or "").strip()
        product_raw = str(product or "").strip()
        reactant_clean = reactant_raw.lstrip("*")
        product_clean = product_raw.lstrip("*")

        _add(f"{reactant_raw}_to_{product_raw}")
        _add(f"{reactant_raw}->{product_raw}")
        _add(f"{reactant_clean}_to_{product_clean}")
        _add(f"{reactant_clean}->{product_clean}")

        reactant_key = self._canonical_species_label(reactant_raw)
        product_key = self._canonical_species_label(product_raw)
        if reactant_key and product_key:
            _add(f"{reactant_key}_to_{product_key}")
            _add(f"{reactant_key}->{product_key}")

        for key in candidates:
            if key in state.neb_results:
                return state.neb_results[key]

        if step_index < len(state.step_structures):
            struct_step_name = state.step_structures[step_index].get("name")
            if struct_step_name in state.neb_results:
                return state.neb_results[struct_step_name]

        neb_values = list(state.neb_results.values())
        if len(neb_values) == len(state.step_structures) and step_index < len(neb_values):
            return neb_values[step_index]

        return None

    @staticmethod
    def _extract_neb_activation_energy(neb_value: Any) -> Optional[float]:
        value = None
        if isinstance(neb_value, dict):
            value = neb_value.get("activation_energy_forward", None)
        else:
            value = getattr(neb_value, "activation_energy_forward", None)
        try:
            if value is None:
                return None
            f = float(value)
            if not np.isfinite(f):
                return None
            return f
        except Exception:
            return None

    def _collect_neb_quality_metrics(self, neb_results: Dict[str, Any]) -> Dict[str, float]:
        total = len(neb_results or {})
        if total <= 0:
            return {
                "total": 0.0,
                "converged_ratio": 0.0,
                "barrier_count": 0.0,
                "reasonable_barrier_ratio": 0.0,
            }

        converged = 0
        barriers: List[float] = []

        for value in (neb_results or {}).values():
            if isinstance(value, dict):
                converged_flag = value.get("converged", None)
                if converged_flag is None:
                    if "error" not in value:
                        converged += 1
                elif bool(converged_flag):
                    converged += 1
            else:
                converged_flag = getattr(value, "converged", None)
                if converged_flag is None or bool(converged_flag):
                    converged += 1

            ea = self._extract_neb_activation_energy(value)
            if ea is not None:
                barriers.append(float(ea))

        min_barrier = float(os.getenv("CATDT_AGENT45_BARRIER_MIN_EV", "0.05"))
        max_barrier = float(os.getenv("CATDT_AGENT45_BARRIER_MAX_EV", "3.0"))
        reasonable_count = sum(1 for ea in barriers if min_barrier <= ea <= max_barrier)
        barrier_count = len(barriers)

        return {
            "total": float(total),
            "converged_ratio": float(converged) / float(total),
            "barrier_count": float(barrier_count),
            "reasonable_barrier_ratio": (float(reasonable_count) / float(barrier_count)) if barrier_count > 0 else 0.0,
        }

    def _compute_iteration_reward(self, validation_passed: bool, neb_results: Dict[str, Any]) -> float:
        if not validation_passed:
            return 0.0
        if not neb_results:
            return 1.0

        metrics = self._collect_neb_quality_metrics(neb_results)
        converged_ratio = float(metrics.get("converged_ratio", 0.0))
        barrier_count = int(metrics.get("barrier_count", 0.0))
        reasonable_ratio = float(metrics.get("reasonable_barrier_ratio", 0.0)) if barrier_count > 0 else 0.0

        if converged_ratio <= 0.0 or reasonable_ratio <= 0.0:
            return 0.0

        reward = 0.4 * converged_ratio + 0.6 * reasonable_ratio
        return round(max(0.0, min(1.0, float(reward))), 6)

    def _run_pathway_iteration(self, state: WorkflowState, base_structure: Optional[Atoms], config: Optional[CatDTConfig]) -> None:
        if base_structure is None or config is None:
            return

        self._reset_agent45_dialogue()
        if self.tools is not None:
            self.tools._agent45_site_retry_state = {}
        state.neb_results = {}
        state.energy_gate_report = {}
        preplaced_products_cache: Optional[Dict[str, Dict[str, Any]]] = None

        for i in range(self._max_pathway_iters):
            self._last_iteration_count = i + 1
            iteration_num = i + 1
            logger.info("Pathway iteration %d/%d", i + 1, self._max_pathway_iters)

            if preplaced_products_cache is None:
                preplaced_products_cache = self._llm_driven_preplace_products(
                    state=state,
                    base_structure=base_structure,
                    iteration=i + 1,
                    config=config,
                )
            preplaced_products = preplaced_products_cache

            if self.agent4 is None:
                raise RuntimeError("Agent4 is required but not initialized.")
            if self.agent5 is None:
                raise RuntimeError("Agent5 is required but not initialized.")

            _ = self._retrieve_agent45_memento_context(state=state, iteration=iteration_num)
            _ = self._retrieve_agent45_knowledge_context(state=state, iteration=iteration_num)
            _ = self._retrieve_agent45_skill_context(state=state, iteration=iteration_num)

            # ── Sequential step-by-step: build → Agent4/5 → relax → chain ──
            intermediates = self._get_working_intermediate_sequence(state)
            n_steps = max(len(intermediates) - 1, 0)
            surface = self._get_surface()

            # Initialize step-0 reactant from the FIRST Stage A clean intermediate
            # (not from `surface.endpoint`, which has been progressively updated
            # by Stage A's relaxations all the way to the last intermediate —
            # its indices and atom count no longer match base_structure).
            first_key = self._canonical_species_label(intermediates[0]) if intermediates else ""
            first_entry = preplaced_products_cache.get(first_key) if preplaced_products_cache else None
            if first_entry and isinstance(first_entry.get("structure"), Atoms):
                current_structure = first_entry["structure"].copy()
                _s_idx = list(first_entry.get("surface_indices", []))
                _a_idx = list(first_entry.get("adsorbate_indices", []))
            else:
                # Fallback (shouldn't happen after Stage A): use base_structure
                # and whatever the surface currently holds.
                current_structure = base_structure.copy()
                _s_idx = list(surface.surface_indices) if surface else []
                _a_idx = list(surface.adsorbate_indices) if surface else []
            self._set_surface_indices(current_structure, _s_idx, _a_idx)
            # Also re-sync surface.endpoint back to the first intermediate so
            # downstream tools that read surface.adsorbate_indices get the
            # correct step-0 state, not the residual last-Stage-A endpoint.
            if surface is not None and first_entry is not None:
                try:
                    surface.update_endpoint(current_structure, _s_idx, _a_idx)
                except Exception as exc:
                    logger.warning("Failed to re-sync surface endpoint to step-0 intermediate: %s", exc)

            all_step_structures: List[Dict[str, Any]] = []
            for step_idx in range(n_steps):
                step_label = f"{intermediates[step_idx]}_to_{intermediates[step_idx + 1]}"
                logger.info(
                    "[iter %d] Building step %d/%d: %s",
                    iteration_num, step_idx + 1, n_steps, step_label,
                )

                # 1. Build single-step baseline (reactant = current_structure)
                baseline_entry = self.tools.build_single_step_baseline(
                    workflow=self,
                    state=state,
                    current_structure=current_structure,
                    preplaced_products=preplaced_products,
                    step_idx=step_idx,
                )
                state.tool_baseline_steps = [baseline_entry]
                self.tools._agent45_baseline_steps = [baseline_entry]

                # Temporarily scope intermediates to this step so _build_step_structures
                # alignment logic matches Agent4's single-step output.
                _saved_intermediates = list(state.intermediates or [])
                state.intermediates = [intermediates[step_idx], intermediates[step_idx + 1]]
                if state.reaction_context is not None:
                    state.reaction_context.intermediates = state.intermediates

                # 2. Agent4 designs staging for this single step
                agent4_retry_feedback = ""
                max_agent4_attempts = 2
                step_structure = None

                for agent4_attempt in range(max_agent4_attempts):
                    design_task = self._create_pathway_design_task(
                        state,
                        current_structure,
                        preplaced_products=preplaced_products,
                        extra_feedback=agent4_retry_feedback,
                    )
                    self._persist_prompt_snapshot(
                        state=state,
                        iteration=iteration_num,
                        agent_name=f"agent4_step{step_idx + 1:02d}",
                        prompt_text=design_task.description,
                    )
                    try:
                        self._run_task(design_task, use_memory=True)
                    except RuntimeError as exc:
                        err_text = str(exc)
                        design_issue = any(
                            marker in err_text
                            for marker in [
                                "pathway design has",
                                "formula matching failed",
                                "returned empty pathway design",
                            ]
                        )
                        if design_issue and (agent4_attempt + 1 < max_agent4_attempts):
                            agent4_retry_feedback = (
                                f"SYSTEM CHECK FAILED: {err_text}\n"
                                "RULES: "
                                "1) Use the exact reaction pathway given — step count must match exactly. "
                                "2) reactant_formula/product_formula must match the given pathway order. "
                                "3) ONLY use atoms_to_add with [x,y,z] positions. atoms_to_modify and atoms_to_remove MUST be empty arrays []."
                            )
                            logger.warning(
                                "Agent4 design failed in iter %d step %d attempt %d: %s",
                                iteration_num, step_idx + 1, agent4_attempt + 1, err_text,
                            )
                            continue
                        raise

                    if not state.pathway_design or not state.pathway_design.steps:
                        if agent4_attempt + 1 < max_agent4_attempts:
                            agent4_retry_feedback = (
                                "SYSTEM CHECK FAILED: pathway_design is empty. "
                                "Output a complete PathwayDesign JSON with atoms_to_add."
                            )
                            continue
                        raise RuntimeError("Agent4 returned empty pathway design.")

                    try:
                        built = self._build_step_structures(
                            state,
                            current_structure,
                            state.pathway_design.steps,
                            preplaced_products=preplaced_products,
                        )
                        if built:
                            step_structure = built[0]
                        break
                    except RuntimeError as exc:
                        err_text = str(exc)
                        design_issue = any(
                            marker in err_text
                            for marker in [
                                "requires atoms_to_add",
                                "endpoint element multisets differ",
                                "changed immutable",
                                "attempts to remove immutable core atoms",
                            ]
                        )
                        if design_issue and (agent4_attempt + 1 < max_agent4_attempts):
                            agent4_retry_feedback = (
                                f"SYSTEM CHECK FAILED: {err_text}\n"
                                "RULES: "
                                "1) Existing atom coordinates are READ-ONLY. "
                                "2) ONLY use atoms_to_add to balance elements. "
                                "3) Each step's reactant/product must have identical element multisets."
                            )
                            continue
                        raise
                else:
                    raise RuntimeError(f"Agent4 design could not pass checks for step {step_idx + 1}.")

                if step_structure is None:
                    raise RuntimeError(f"No step structure built for step {step_idx + 1}")

                # 3. Agent5 validates this single step
                state.step_structures = [step_structure]
                state.step_structures = self._agent4_postprocess_step_structures(
                    state=state, base_structure=current_structure,
                    step_structures=state.step_structures,
                )

                energy_gate_dir = surface.agent_dir(
                    "agent45",
                    subdir=f"iter_{iteration_num:02d}_energy_gate_step{step_idx + 1:02d}",
                )
                try:
                    step_energy_gate = self.tools.run_agent45_energy_gate(
                        workflow=self,
                        step_structures=state.step_structures,
                        output_dir=str(energy_gate_dir),
                    )
                except Exception as exc:
                    logger.warning("Energy gate failed for step %d: %s", step_idx + 1, exc)
                    step_energy_gate = {"status": "FAIL", "fatal_issues": [str(exc)], "warning_issues": []}

                validation_task = self._create_validation_task(state)
                self._persist_prompt_snapshot(
                    state=state, iteration=iteration_num,
                    agent_name=f"agent5_step{step_idx + 1:02d}",
                    prompt_text=validation_task.description,
                )
                self._run_task(validation_task, use_memory=True)

                # 4. Relax endpoints for this step
                relax_dir = surface.agent_dir(
                    "agent6", subdir=f"pre_neb_relax/step{step_idx + 1:02d}",
                )
                try:
                    state.step_structures = self._relax_neb_endpoints(
                        step_structures=state.step_structures,
                        output_dir=str(relax_dir),
                        fmax=0.05,
                        max_steps=200,
                    )
                except Exception as exc:
                    logger.warning("Endpoint relax failed for step %d: %s", step_idx + 1, exc)

                relaxed_step = state.step_structures[0]
                all_step_structures.append(relaxed_step)

                # Checkpoint after each pathway step
                self._auto_checkpoint(f"pathway_step{step_idx + 1:02d}")

                # 5. Chain: next step's reactant = CLEAN product from Stage A.
                # Per Phase-2 spec: staging atoms (added on this step's endpoints
                # for NEB element-count balance) must NOT propagate to the next
                # step. The next reactant is the clean intermediate produced by
                # _llm_driven_preplace_products, not the staged relaxed product.
                if step_idx + 1 < n_steps:
                    next_reactant_label = intermediates[step_idx + 1]
                    next_key = self._canonical_species_label(next_reactant_label)
                    next_entry = preplaced_products_cache.get(next_key) if preplaced_products_cache else None
                    if not (next_entry and isinstance(next_entry.get("structure"), Atoms)):
                        raise RuntimeError(
                            f"[iter {iteration_num}] Step {step_idx + 1}: "
                            f"no Stage A clean intermediate for '{next_reactant_label}'. "
                            "Stage A (_llm_driven_preplace_products) must produce "
                            "every intermediate before Stage B can consume them."
                        )
                    current_structure = next_entry["structure"].copy()
                    clean_surf = list(next_entry.get("surface_indices", []))
                    clean_ads = list(next_entry.get("adsorbate_indices", []))
                    self._update_indices(state, current_structure, clean_surf, clean_ads)
                    logger.info(
                        "[iter %d] Step %d done. Next reactant <- Stage A clean '%s' (%d atoms, ads=%d)",
                        iteration_num, step_idx + 1, next_reactant_label,
                        len(current_structure), len(clean_ads),
                    )

                # Restore full intermediates for next iteration
                state.intermediates = _saved_intermediates
                if state.reaction_context is not None:
                    state.reaction_context.intermediates = _saved_intermediates

            # ── All steps built. Collect results. ──
            state.step_structures = all_step_structures
            state.energy_gate_report = {}  # individual gates already ran

            self._export_pathway_iteration_visuals(state, i + 1)

            validation_passed = bool(
                state.validation_report and state.validation_report.status.upper() == "PASS"
            )

            state.neb_results = {}
            if validation_passed and config.calculate_barriers and state.step_structures and self.tools is not None:
                neb_dir = surface.agent_dir("agent6")
                try:
                    state.neb_results = self.tools.run_neb_for_steps(
                        steps=state.step_structures,
                        output_dir=str(neb_dir),
                        n_frames=max(config.neb_n_frames, 7),
                        fmax=0.1,
                        max_steps=max(config.neb_max_steps, 200),
                    )
                except Exception as exc:
                    logger.warning("NEB run failed in iteration %d: %s", i + 1, exc)
                    state.neb_results = {}
            elif validation_passed and config.calculate_barriers and state.step_structures:
                logger.warning("Skipping NEB run in iteration %d because tools are not initialized", i + 1)

            iter_reward = self._compute_iteration_reward(
                validation_passed=validation_passed,
                neb_results=state.neb_results or {},
            )

            self._persist_agent45_memento_case(
                state=state,
                iteration=iteration_num,
                reward=iter_reward,
            )

            if state.validation_report:
                preview = ""
                if state.validation_report.issues:
                    preview = "; ".join(state.validation_report.issues[:2])
                state.iteration_history.append(
                    {
                        "iteration": i + 1,
                        "status": state.validation_report.status,
                        "issue_count": len(state.validation_report.issues or []),
                        "issue_preview": preview,
                    }
                )

            if validation_passed:
                logger.info("Pathway validated at iteration %d", i + 1)
                break

    def _relax_neb_endpoints(
        self,
        step_structures: List[Dict[str, Any]],
        output_dir: str,
        fmax: float = 0.03,
        max_steps: int = 100,
        enforce_connectivity: bool = True,
    ) -> List[Dict[str, Any]]:
        """Relax adsorbate atoms (surface fixed) on NEB endpoints before NEB.

        This ensures both endpoints are true local energy minima, which is critical
        for obtaining meaningful NEB barriers — especially when co-adsorbates have
        been placed at surface sites by Agent4.
        """
        import ase.optimize as ase_opt
        from ase.constraints import FixAtoms
        from ase.io import write as ase_write

        predictor = self.tools.get_shared_fairchem_predictor(
            cache_key="uma_shared",
            model_name="uma-s-1p1",
            use_gpu=True,
            device="cuda",
            work_subdir="_shared_fairchem_global",
            keep_files=False,
            verbose=False,
        )
        predictor._load_model()
        calc = predictor._calculator

        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # Top-2 layer indices are frozen at end of agent3 in full pipeline
        # runs. In focused single-step/replay tests they may be unavailable;
        # fall back to adsorbate-only relaxation instead of aborting.
        surface = self._get_surface()
        frozen_top_layer = set(surface.top_layer_indices) if surface is not None and surface.top_layer_indices else set()
        if not frozen_top_layer:
            logger.warning(
                "_relax_neb_endpoints: surface.top_layer_indices unavailable; "
                "falling back to adsorbate/staged-only relaxation."
            )

        for step_idx, step in enumerate(step_structures):
            name = step.get("name", f"step{step_idx:02d}")

            for side in ("reactant", "product"):
                atoms = step.get(side)
                if not isinstance(atoms, Atoms):
                    continue

                ads_indices = {i for i in set(step.get(f"{side}_adsorbate_indices", [])) if 0 <= i < len(atoms)}
                staged_indices = {i for i in set(step.get(f"{side}_staged_indices", [])) if 0 <= i < len(atoms)}

                # Only keep indices that are still in range for this snapshot
                # (Stage B may have appended staging atoms beyond the surface range).
                top2_surface = {i for i in frozen_top_layer if 0 <= i < len(atoms)}

                if top2_surface:
                    # Full pipeline: staged atoms + top surface layers move.
                    # The adsorbate core stays at its Stage A local minimum.
                    free_indices = staged_indices | top2_surface
                else:
                    # Replay/unit-test mode without a registered surface.
                    free_indices = ads_indices | staged_indices

                if not free_indices:
                    continue

                original = atoms.copy()
                relax_free_for_connectivity = ads_indices | staged_indices
                if enforce_connectivity:
                    connectivity_before = self._free_atom_connectivity_signature(
                        original,
                        relax_free_for_connectivity,
                        staged_indices=staged_indices,
                    )
                    cross_connectivity_before = self._free_atom_cross_connectivity_signature(
                        original,
                        relax_free_for_connectivity,
                        staged_indices=staged_indices,
                    )
                else:
                    connectivity_before = None
                    cross_connectivity_before = None
                a = atoms.copy()
                a.pbc = [True, True, True]
                a.calc = calc
                fixed_indices = [i for i in range(len(a)) if i not in free_indices]
                logger.info(
                    "Pre-NEB relax %s/%s: %d free atoms (ads=%d, staged=%d, top2_surface=%d)",
                    name, side, len(free_indices), len(ads_indices), len(staged_indices), len(top2_surface),
                )
                if fixed_indices:
                    a.set_constraint(FixAtoms(indices=fixed_indices))

                log_file = str(out_path / f"{name}_{side}_relax.log")
                try:
                    e_before = a.get_potential_energy()
                    # Use a conservative optimizer step so staged atoms do not
                    # immediately collapse into new bonds before we can capture
                    # a lower-energy connectivity-preserving snapshot.
                    opt_cls = getattr(ase_opt, "FIRE")
                    try:
                        opt = opt_cls(a, logfile=log_file, maxstep=0.05)
                    except TypeError:
                        opt = opt_cls(a, logfile=log_file)
                    last_safe = original.copy()
                    try:
                        last_safe.set_constraint()
                    except Exception:
                        pass
                    safe_snapshot_steps = 0
                    connectivity_changed = False

                    if enforce_connectivity and connectivity_before is not None:
                        def _capture_last_safe_snapshot() -> None:
                            nonlocal last_safe, safe_snapshot_steps, connectivity_changed
                            connectivity_now = self._free_atom_connectivity_signature(
                                a,
                                relax_free_for_connectivity,
                                staged_indices=staged_indices,
                            )
                            cross_connectivity_now = self._free_atom_cross_connectivity_signature(
                                a,
                                relax_free_for_connectivity,
                                staged_indices=staged_indices,
                            )
                            new_cross_bonds = set(cross_connectivity_now) - set(cross_connectivity_before or ())
                            if connectivity_now != connectivity_before or new_cross_bonds:
                                connectivity_changed = True
                                raise StopIteration("free-atom connectivity changed during pre-NEB relax")
                            last_safe = a.copy()
                            try:
                                last_safe.set_constraint()
                            except Exception:
                                pass
                            safe_snapshot_steps = int(getattr(opt, "nsteps", 0))

                        if hasattr(opt, "attach"):
                            opt.attach(_capture_last_safe_snapshot, interval=1)

                    try:
                        opt.run(fmax=fmax, steps=max_steps)
                    except Exception as exc:
                        msg = str(exc)
                        if isinstance(exc, StopIteration) or "StopIteration" in msg:
                            logger.info(
                                "Pre-NEB relax stopped early for %s/%s: %s",
                                name,
                                side,
                                msg,
                            )
                        else:
                            raise
                    finally:
                        if enforce_connectivity and connectivity_before is not None and hasattr(opt, "detach"):
                            try:
                                opt.detach(_capture_last_safe_snapshot)
                            except Exception:
                                pass

                    n_steps = int(getattr(opt, "nsteps", 0))

                    selected = a
                    if enforce_connectivity and connectivity_before is not None:
                        selected_connectivity = self._free_atom_connectivity_signature(
                            selected,
                            relax_free_for_connectivity,
                            staged_indices=staged_indices,
                        )
                        selected_cross_connectivity = self._free_atom_cross_connectivity_signature(
                            selected,
                            relax_free_for_connectivity,
                            staged_indices=staged_indices,
                        )
                        new_cross_bonds = set(selected_cross_connectivity) - set(cross_connectivity_before or ())
                        if selected_connectivity != connectivity_before or new_cross_bonds or connectivity_changed:
                            selected = last_safe.copy()
                            chosen_label = (
                                f"last connectivity-safe snapshot (step {safe_snapshot_steps})"
                                if safe_snapshot_steps > 0
                                else "original endpoint"
                            )
                            logger.warning(
                                "Pre-NEB relax changed connectivity for %s/%s; keeping %s. before=%s after=%s",
                                name,
                                side,
                                chosen_label,
                                connectivity_before,
                                self._free_atom_connectivity_signature(
                                    a,
                                    relax_free_for_connectivity,
                                    staged_indices=staged_indices,
                                ),
                            )
                    e_after = e_before
                    try:
                        selected_for_energy = selected.copy()
                        selected_for_energy.calc = calc
                        e_after = selected_for_energy.get_potential_energy()
                    except Exception:
                        pass

                    logger.info(
                        "Pre-NEB relax %s/%s: E %.3f → %.3f eV (dE=%+.3f) steps=%d",
                        name, side, e_before, e_after, e_after - e_before, n_steps,
                    )

                    ase_write(str(out_path / f"{name}_{side}_relaxed.vasp"), selected)

                    # Update step structure in-place (remove constraint for NEB)
                    try:
                        selected.set_constraint()
                    except Exception:
                        pass
                    step[side] = selected
                    # Ensure indices survive relaxation in atoms.info
                    selected.info["surface_indices"] = [
                        i for i in range(len(selected))
                        if i not in set(step.get(f"{side}_adsorbate_indices", []))
                        and i not in set(step.get(f"{side}_staged_indices", []))
                    ]
                    selected.info["adsorbate_indices"] = list(step.get(f"{side}_adsorbate_indices", []))
                    selected.info["staged_indices"] = list(step.get(f"{side}_staged_indices", []))

                except Exception as exc:
                    logger.warning("Pre-NEB relax failed for %s/%s: %s", name, side, exc)

        return step_structures

    def _free_atom_connectivity_signature(
        self,
        atoms: Atoms,
        free_indices: List[int] | set[int],
        staged_indices: Optional[List[int] | set[int]] = None,
    ) -> Tuple[Tuple[int, int], ...]:
        valid = sorted({int(i) for i in free_indices if 0 <= int(i) < len(atoms)})
        staged = {int(i) for i in (staged_indices or set()) if 0 <= int(i) < len(atoms)}
        symbols = atoms.get_chemical_symbols()
        signature: List[Tuple[int, int]] = []

        for pos, i in enumerate(valid):
            if symbols[i] == "Cu":
                continue
            for j in valid[pos + 1:]:
                if symbols[j] == "Cu":
                    continue
                cutoff = 1.15 * (
                    covalent_radii[atomic_numbers[symbols[i]]] + covalent_radii[atomic_numbers[symbols[j]]]
                )
                if atoms.get_distance(i, j, mic=True) <= cutoff:
                    if staged and ((i in staged) != (j in staged)):
                        continue
                    signature.append((i, j))

        return tuple(signature)

    def _free_atom_cross_connectivity_signature(
        self,
        atoms: Atoms,
        free_indices: List[int] | set[int],
        staged_indices: List[int] | set[int],
    ) -> Tuple[Tuple[int, int], ...]:
        valid = sorted({int(i) for i in free_indices if 0 <= int(i) < len(atoms)})
        staged = {int(i) for i in staged_indices if 0 <= int(i) < len(atoms)}
        if not staged:
            return tuple()
        symbols = atoms.get_chemical_symbols()
        signature: List[Tuple[int, int]] = []

        for pos, i in enumerate(valid):
            if symbols[i] == "Cu":
                continue
            for j in valid[pos + 1:]:
                if symbols[j] == "Cu" or ((i in staged) == (j in staged)):
                    continue
                cutoff = 1.15 * (
                    covalent_radii[atomic_numbers[symbols[i]]] + covalent_radii[atomic_numbers[symbols[j]]]
                )
                if atoms.get_distance(i, j, mic=True) <= cutoff:
                    signature.append((i, j))

        return tuple(signature)

    def _compute_episode_reward(self, state: WorkflowState) -> float:
        validation_score = 1.0 if (state.validation_report and state.validation_report.status.upper() == "PASS") else 0.0

        neb_metrics = self._collect_neb_quality_metrics(state.neb_results or {})
        neb_total = int(neb_metrics.get("total", 0.0))
        if neb_total > 0:
            neb_score = (
                0.6 * float(neb_metrics.get("converged_ratio", 0.0))
                + 0.4 * float(neb_metrics.get("reasonable_barrier_ratio", 0.0))
            )
        else:
            neb_score = 0.0

        target_steps = max(len((state.reaction_context.intermediates if state.reaction_context else [])) - 1, 1)
        step_score = min(len(state.step_structures), target_steps) / target_steps

        reward = 0.5 * validation_score + 0.3 * neb_score + 0.2 * step_score
        return round(float(reward), 6)

    # ==================================================================
    # LLM-driven preplacement
    # ==================================================================

    def _llm_driven_preplace_products(
        self,
        state: WorkflowState,
        base_structure: Atoms,
        iteration: int,
        config: Optional[CatDTConfig] = None,
    ) -> Dict[str, Any]:
        """Build each intermediate structure sequentially using Agent4 LLM.

        intermediates[0] uses base_structure directly (from Agent2/3).
        Each subsequent intermediate is built by calling Agent4 to decide
        atom additions/removals, then relaxing with UMA.
        """
        from core.pathway.adsorbate_resolver import resolve_adsorbate
        from collections import Counter

        intermediates = self._get_working_intermediate_sequence(state)
        if not intermediates:
            return {}

        placed: Dict[str, Dict[str, Any]] = {}

        surface = self._get_surface()
        if surface is None:
            raise RuntimeError(
                "_llm_driven_preplace_products requires a CatalystSurface. "
                "Call run_on_surface() or run() first."
            )

        surf_idx = list(surface.surface_indices)
        ads_idx = list(surface.adsorbate_indices)

        # intermediates[0]: use base_structure as-is (this is the agent3 output)
        base_copy = base_structure.copy()
        self._set_surface_indices(base_copy, surf_idx, ads_idx)
        first_key = self._canonical_species_label(intermediates[0])
        # First intermediate is already registered as an IntermediateAdsorbate at
        # the end of agent3 (see _freeze_surface_after_agent3). Keep legacy
        # placed-dict entry for downstream consumers that still read it.
        placed[first_key] = {
            "label": intermediates[0],
            "structure": base_copy,
            "source": "agent3_input",
            "path": surface.reconstructed_surface_path or "",
            "surface_indices": surf_idx,
            "adsorbate_indices": ads_idx,
            "fixed_atom_count": int(len(base_copy)),
        }

        if len(intermediates) <= 1:
            return placed

        iter_dir = surface.agent_dir(
            "agent45", subdir=f"iter_{iteration:02d}_tool_preplacement",
        )

        # Save clean surface for reference
        clean_surface = self._build_clean_surface_from_known_indices(base_structure, state)
        write(iter_dir / "clean_surface_known_indices.vasp", clean_surface)

        from camel_agents.surface_model import StructureEndpoint
        prev_ep = StructureEndpoint.from_parts(base_copy, surf_idx, ads_idx)

        for idx, label in enumerate(intermediates[1:], start=1):
            if self._contains_gas_phase_hint(label):
                continue
            key = self._canonical_species_label(label)
            prev_label = intermediates[idx - 1]

            # Get target molecule from DB
            try:
                target_atoms, target_binding = resolve_adsorbate(label)
            except Exception as exc:
                logger.warning("Cannot resolve adsorbate '%s': %s", label, exc)
                continue

            # Format current adsorbate coordinates for prompt
            prev_syms = prev_ep.atoms.get_chemical_symbols()
            prev_pos = prev_ep.atoms.get_positions()
            current_ads_text = f"Species: {prev_label}\n"
            for i in prev_ep.adsorbate_indices:
                if 0 <= i < len(prev_ep):
                    current_ads_text += f"  [{i}] {prev_syms[i]}  ({prev_pos[i][0]:.4f}, {prev_pos[i][1]:.4f}, {prev_pos[i][2]:.4f})\n"

            # Format target molecule coordinates
            t_syms = target_atoms.get_chemical_symbols()
            t_pos = target_atoms.get_positions()
            target_mol_text = f"Species: {label}  Formula: {target_atoms.get_chemical_formula('hill')}  Binding: {target_binding}\n"
            for j, (s, p) in enumerate(zip(t_syms, t_pos)):
                target_mol_text += f"  [{j}] {s}  ({p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f})\n"

            # Compute element delta
            prev_ads_syms = [prev_syms[i] for i in prev_ep.adsorbate_indices if 0 <= i < len(prev_ep.atoms)]
            prev_counts = Counter(prev_ads_syms)
            target_counts = Counter(t_syms)
            delta_parts = []
            for elem in sorted(set(list(prev_counts.keys()) + list(target_counts.keys()))):
                d = target_counts.get(elem, 0) - prev_counts.get(elem, 0)
                if d != 0:
                    delta_parts.append(f"{d:+d}{elem}")
            element_delta = ", ".join(delta_parts) if delta_parts else "none"

            # Surface info
            cell = prev_ep.atoms.cell
            surface_top_z = self._estimate_surface_top_z(prev_ep.atoms, prev_ep.adsorbate_indices)
            surface_info = f"surface_top_z={surface_top_z:.3f} Å, cell=({cell[0][0]:.2f}, {cell[1][1]:.2f}, {cell[2][2]:.2f})"

            # ── Agent4 design + relax + Agent5 validation loop ──
            max_stageA_attempts = 3
            stageA_feedback = ""

            for stageA_attempt in range(max_stageA_attempts):
                # Create preplacement task for Agent4
                extra_feedback = (
                    f"\n\n[PREVIOUS ATTEMPT FAILED]\n{stageA_feedback}\n"
                    "Fix the issues above. Place atoms closer to the surface "
                    "and ensure they bond to the adsorbate backbone."
                ) if stageA_feedback else ""

                description = TaskPrompts.preplacement_design(
                    current_label=prev_label,
                    target_label=label,
                    current_adsorbate_coords=current_ads_text,
                    target_molecule_coords=target_mol_text,
                    surface_info=surface_info,
                    element_delta=element_delta,
                ) + extra_feedback

                # Temporarily set baseline for suggest_staged_positions tool
                baseline_entry = {
                    "name": f"{prev_label}_to_{label}",
                    "reactant": prev_ep.atoms.copy(),
                    "product": prev_ep.atoms.copy(),
                    "reactant_adsorbate_indices": prev_ep.adsorbate_indices,
                    "product_adsorbate_indices": prev_ep.adsorbate_indices,
                    "reactant_formula": prev_label,
                    "product_formula": label,
                }
                state.tool_baseline_steps = [baseline_entry]
                if self.tools is not None:
                    self.tools._agent45_baseline_steps = [baseline_entry]

                # Scope intermediates to this pair
                _saved = list(state.intermediates or [])
                state.intermediates = [prev_label, label]
                if state.reaction_context is not None:
                    state.reaction_context.intermediates = state.intermediates

                # Reset Agent4 dialogue for fresh context each step
                try:
                    chat_agent = getattr(self.agent4, "chat_agent", None)
                    if chat_agent is not None and hasattr(chat_agent, "reset"):
                        chat_agent.reset()
                except Exception:
                    pass

                # Call Agent4 LLM
                logger.info(
                    "[stageA] Building %s from %s (delta=%s, attempt %d/%d)",
                    label, prev_label, element_delta,
                    stageA_attempt + 1, max_stageA_attempts,
                )
                self._persist_prompt_snapshot(
                    state, iteration,
                    f"preplacement_{self._species_safe_name(label)}_attempt{stageA_attempt + 1}",
                    description,
                )

                def _apply_preplacement_design(output: "TaskOutput"):
                    payload = output.json_dict or {}
                    design_obj = output.pydantic
                    if design_obj is None:
                        try:
                            design_obj = PathwayDesign(**self._coerce_pathway_design_payload(state, payload))
                        except Exception:
                            design_obj = None
                    state.pathway_design = design_obj

                preplacement_task = Task(
                    description=description,
                    expected_output="Atom placement JSON",
                    agent=self.agent4,
                    output_pydantic=PathwayDesign,
                    result_handler=_apply_preplacement_design,
                )

                try:
                    self._run_task(preplacement_task, use_memory=True)
                except Exception as exc:
                    logger.warning("[stageA] Agent4 failed for %s (attempt %d): %s", label, stageA_attempt + 1, exc)
                    stageA_feedback = f"Agent4 call failed: {exc}"
                    continue

                # Apply Agent4's atom operations to prev_ep.atoms
                design = state.pathway_design
                new_structure = prev_ep.atoms.copy()
                new_ads_indices = list(prev_ep.adsorbate_indices)

                if design and design.steps:
                    step_spec = design.steps[0]

                    # Remove atoms first (highest index first)
                    remove_indices = []
                    for item in (step_spec.atoms_to_remove or []):
                        idx = None
                        if isinstance(item, dict):
                            idx = item.get("index")
                        elif isinstance(item, int):
                            idx = item
                        if idx is not None and isinstance(idx, int) and 0 <= idx < len(new_structure):
                            if idx in set(new_ads_indices):
                                remove_indices.append(idx)
                    for idx in sorted(set(remove_indices), reverse=True):
                        del new_structure[idx]
                        new_ads_indices = [j if j < idx else j - 1
                                           for j in new_ads_indices if j != idx]

                    # Add atoms
                    for atom_spec in (step_spec.atoms_to_add or []):
                        elements = self._resolve_addition_elements(atom_spec)
                        pos = atom_spec.position
                        if pos and len(pos) == 3 and elements:
                            for elem in elements:
                                new_structure.append(Atom(elem, position=[float(pos[0]), float(pos[1]), float(pos[2])]))
                                new_ads_indices.append(len(new_structure) - 1)

                new_ads_indices, removed_for_formula = self._enforce_stageA_target_formula(
                    new_structure,
                    new_ads_indices,
                    label,
                )
                if removed_for_formula:
                    logger.info(
                        "[stageA] Formula guard removed excess atoms for %s: %s",
                        label,
                        removed_for_formula,
                    )

                surface_indices = [i for i in range(len(new_structure)) if i not in set(new_ads_indices)]
                self._set_surface_indices(new_structure, surface_indices, new_ads_indices)

                # Save and relax
                pre_relax_path = str(iter_dir / f"{self._species_safe_name(label)}_pre_relax.vasp")
                write(pre_relax_path, new_structure)

                relax_output = str(iter_dir / f"{self._species_safe_name(label)}_product.vasp")
                relaxed_energy: Optional[float] = None
                if self.tools is not None:
                    try:
                        relax_fmax = float(getattr(config, "stagea_relax_fmax", 0.05))
                        relax_steps = int(getattr(config, "stagea_relax_steps", 200))
                        result = self.tools.relax_adsorbate_on_surface(
                            structure_path=pre_relax_path,
                            adsorbate_indices=new_ads_indices,
                            output_path=relax_output,
                            fmax=relax_fmax,
                            max_steps=relax_steps,
                            free_surface_indices=list(surface.top_layer_indices),
                        )
                        new_structure = read(relax_output)
                        self._set_surface_indices(new_structure, surface_indices, new_ads_indices)
                        relaxed_energy = float(result.get("energy", 0.0))
                        logger.info(
                            "[stageA] Relaxed %s: E=%.3f fmax=%.3f (top-2 free)",
                            label, relaxed_energy, result.get("fmax", 0.0),
                        )
                    except Exception as exc:
                        logger.warning("[stageA] Relax failed for %s: %s", label, exc)
                        write(relax_output, new_structure)

                # ── Agent5-style validation: check structural quality ──
                valid, feedback = self._validate_stageA_intermediate(
                    new_structure, new_ads_indices, label,
                )
                if valid:
                    logger.info("[stageA] Validation PASSED for %s", label)
                    break
                else:
                    logger.warning(
                        "[stageA] Validation FAILED for %s (attempt %d/%d): %s",
                        label, stageA_attempt + 1, max_stageA_attempts, feedback,
                    )
                    stageA_feedback = feedback
                    # Restore intermediates before retry
                    state.intermediates = _saved
                    if state.reaction_context is not None:
                        state.reaction_context.intermediates = _saved
            else:
                # All attempts exhausted — accept with warning
                logger.warning(
                    "[stageA] All %d attempts failed validation for %s. "
                    "Accepting last attempt with known quality issues.",
                    max_stageA_attempts, label,
                )

            # Restore intermediates after the loop
            state.intermediates = _saved
            if state.reaction_context is not None:
                state.reaction_context.intermediates = _saved

            placed[key] = {
                "label": label,
                "structure": new_structure,
                "source": "llm_preplacement",
                "path": relax_output,
                "surface_indices": surface_indices,
                "adsorbate_indices": new_ads_indices,
                "fixed_atom_count": int(len(new_structure)),
            }

            # Register clean Stage A product into the persistent store.
            try:
                surface.add_intermediate(
                    formula=label,
                    atoms=new_structure,
                    adsorbate_indices=new_ads_indices,
                    source="agent45_stageA",
                    parent_formula=prev_label,
                    energy=relaxed_energy,
                    file_path=relax_output,
                )
            except Exception as exc:
                logger.warning("[stageA] Failed to register '%s' into surface.intermediate_adsorbates: %s", label, exc)

            # Chain: next intermediate starts from this relaxed structure
            new_surf = [i for i in range(len(new_structure)) if i not in set(new_ads_indices)]
            prev_ep = StructureEndpoint.from_parts(new_structure, new_surf, new_ads_indices)

            # Write preplaced structure back to surface so checkpoint captures it
            self._update_indices(state, new_structure, new_surf, new_ads_indices)
            self._auto_checkpoint(f"preplacement_{self._species_safe_name(label)}")

        # Dump the Stage A chain as JSON for inspection
        try:
            chain_dump_path = surface.write_json(
                "agent45",
                "stageA_intermediate_chain.json",
                surface.intermediates_summary(),
            )
            logger.info("[stageA] Intermediate chain dumped: %s", chain_dump_path)
        except Exception as exc:
            logger.warning("[stageA] Failed to dump intermediate chain: %s", exc)

        return placed

    def _enforce_stageA_target_formula(
        self,
        atoms: Atoms,
        ads_indices: List[int],
        label: str,
    ) -> Tuple[List[int], List[Dict[str, Any]]]:
        """Trim Stage A adsorbate atoms to the target adsorbate formula.

        Stage A builds the next covalently bound intermediate, not the
        equal-element NEB endpoint. Leaving excess H on radicals turns
        *C3H7/*C3H6 back into saturated propane and corrupts staging.
        """
        target_tokens = self._extract_species_tokens(label)
        if not target_tokens:
            return sorted({i for i in ads_indices if 0 <= i < len(atoms)}), []

        target_counts = Counter(target_tokens)
        valid_ads = sorted({i for i in ads_indices if 0 <= i < len(atoms)})
        current_counts = Counter(atoms[i].symbol for i in valid_ads)
        remove_plan: List[int] = []
        positions = atoms.get_positions()

        heavy_ads = [i for i in valid_ads if atoms[i].symbol != "H"]
        anchor = (
            np.mean([positions[i] for i in heavy_ads], axis=0)
            if heavy_ads
            else np.mean([positions[i] for i in valid_ads], axis=0)
        )

        for element, count in current_counts.items():
            excess = int(count - target_counts.get(element, 0))
            if excess <= 0:
                continue
            candidates = [i for i in valid_ads if atoms[i].symbol == element]
            candidates.sort(
                key=lambda i: (
                    float(np.linalg.norm(positions[i] - anchor)),
                    float(positions[i, 2]),
                ),
                reverse=True,
            )
            remove_plan.extend(candidates[:excess])

        if not remove_plan:
            return valid_ads, []

        removed: List[Dict[str, Any]] = []
        for idx in sorted(set(remove_plan)):
            if 0 <= idx < len(atoms):
                p = atoms.positions[idx]
                removed.append(
                    {
                        "index": int(idx),
                        "symbol": str(atoms[idx].symbol),
                        "position": [float(p[0]), float(p[1]), float(p[2])],
                    }
                )

        new_ads, _ = self._delete_atoms_with_index_tracking(
            structure=atoms,
            adsorbate_indices=valid_ads,
            staged_indices=[],
            delete_indices=remove_plan,
        )
        return new_ads, removed

    @staticmethod
    def _validate_stageA_intermediate(
        atoms: Atoms,
        ads_indices: List[int],
        label: str,
    ) -> Tuple[bool, str]:
        """Validate a Stage A intermediate after UMA relax.

        Fully general — works for any adsorbate molecule on any surface.
        Uses a bond-graph approach: every adsorbate atom must be reachable
        (via covalent bonds) from at least one surface-bound atom. Atoms
        that are neither bonded to other adsorbate atoms nor to the surface
        are "dangling" and indicate a structural problem.

        Checks:
          1. At least one adsorbate atom is chemically bound to the surface
          2. All adsorbate atoms form a single connected component (including
             bonds to the surface as anchors)
          3. No same-element homonuclear diatomic molecules formed in gas
             phase (e.g. H-H < 0.85 Å, O-O < 1.30 Å)

        Returns (ok, feedback_text).
        """
        from ase.data import covalent_radii, atomic_numbers
        import numpy as np

        if not ads_indices:
            return True, ""

        syms = atoms.get_chemical_symbols()
        pos = atoms.get_positions()
        n = len(atoms)
        valid_ads = [i for i in ads_indices if 0 <= i < n]
        surf_set = set(range(n)) - set(valid_ads)

        if not valid_ads or not surf_set:
            return True, ""

        def _cov_cut(a, b):
            return 1.3 * (covalent_radii[atomic_numbers[syms[a]]] +
                          covalent_radii[atomic_numbers[syms[b]]])

        issues = []

        # Build bond graph among adsorbate atoms
        ads_local = {gi: li for li, gi in enumerate(valid_ads)}
        adj = {li: set() for li in range(len(valid_ads))}
        for li in range(len(valid_ads)):
            for lj in range(li + 1, len(valid_ads)):
                gi, gj = valid_ads[li], valid_ads[lj]
                d = float(np.linalg.norm(pos[gi] - pos[gj]))
                if d < _cov_cut(gi, gj):
                    adj[li].add(lj)
                    adj[lj].add(li)

        # 1. Surface binding: which ads atoms are bonded to a surface atom?
        surface_anchors: set = set()  # local indices of ads atoms bound to surface
        for li, gi in enumerate(valid_ads):
            for si in surf_set:
                d = float(np.linalg.norm(pos[gi] - pos[si]))
                if d < _cov_cut(gi, si):
                    surface_anchors.add(li)
                    break

        if not surface_anchors:
            issues.append(
                "No adsorbate atom is chemically bound to the surface. "
                "The entire adsorbate may have desorbed."
            )

        # 2. Connectivity: BFS from surface-anchored atoms through ads bond graph.
        #    Any atom NOT reachable = dangling (floating in vacuum).
        reachable: set = set()
        if surface_anchors:
            stack = list(surface_anchors)
            while stack:
                x = stack.pop()
                if x in reachable:
                    continue
                reachable.add(x)
                stack.extend(adj[x] - reachable)

        dangling = [valid_ads[li] for li in range(len(valid_ads)) if li not in reachable]
        if dangling:
            dangling_desc = ", ".join(
                f"{syms[i]}[{i}]" for i in dangling
            )
            issues.append(
                f"Dangling atoms (not connected to any surface-bound atom via "
                f"covalent bonds): {dangling_desc}. These atoms are floating "
                f"in vacuum and need to be repositioned."
            )

        # 3. Gas-phase homonuclear diatomic detection (H₂, O₂, N₂, etc.)
        #    Check if any pair of same-element ads atoms are at their
        #    gas-phase diatomic bond length (much shorter than typical
        #    surface bonds).
        gas_diatomic_cutoffs = {
            "H": 0.85,   # H₂ = 0.74 Å
            "O": 1.30,   # O₂ = 1.21 Å
            "N": 1.15,   # N₂ = 1.10 Å
            "F": 1.50,   # F₂ = 1.42 Å
            "Cl": 2.10,  # Cl₂ = 1.99 Å
        }
        for elem, cutoff in gas_diatomic_cutoffs.items():
            elem_ads = [i for i in valid_ads if syms[i] == elem]
            for ii in range(len(elem_ads)):
                for jj in range(ii + 1, len(elem_ads)):
                    d = float(np.linalg.norm(pos[elem_ads[ii]] - pos[elem_ads[jj]]))
                    if d < cutoff:
                        issues.append(
                            f"Gas-phase {elem}₂ detected: {elem}[{elem_ads[ii]}]-"
                            f"{elem}[{elem_ads[jj]}] d={d:.3f}Å (cutoff {cutoff}Å)"
                        )

        if issues:
            feedback = (
                f"[stageA validation FAIL for '{label}'] "
                + " | ".join(issues)
                + " — Reposition atoms so that every adsorbate atom is connected "
                "to the surface-bound backbone via covalent bonds."
            )
            return False, feedback
        return True, ""

    def _species_safe_name(self, species: str) -> str:
        import re
        base = self._canonical_species_label(species) or str(species or "").strip()
        return re.sub(r"[^A-Za-z0-9._-]+", "_", base) or "species"

    # ==================================================================
    # Multi-facet Wulff-weighted pipeline
    # ==================================================================

    # ── CatalystSurface population helpers ──────────────────────────────

    def _get_surface(self) -> Optional[Any]:
        return getattr(self, "_catalyst_surface", None)

    def _update_indices(self, state: WorkflowState, atoms: Atoms,
                        surface_indices: List[int], adsorbate_indices: List[int]) -> None:
        """Write surface / adsorbate indices onto the current CatalystSurface.

        CatalystSurface is the single source of truth for structural state.
        ``state`` is kept as a parameter for signature compatibility but is no
        longer written to (physical fields have been removed from WorkflowState).
        """
        self._set_surface_indices(atoms, surface_indices, adsorbate_indices)
        surface = self._get_surface()
        if surface is None:
            raise RuntimeError(
                "_update_indices: CatalystSurface not attached — call run_on_surface()."
            )
        surface.update_endpoint(atoms, surface_indices, adsorbate_indices)

    def _auto_checkpoint(self, stage_name: str) -> None:
        surface = self._get_surface()
        if surface is None or not surface.output_dir:
            return
        try:
            ckpt_path = surface.agent_dir("checkpoints") / f"checkpoint_after_{stage_name}.pkl"
            surface.save_checkpoint(ckpt_path)
            logger.info("Auto-checkpoint saved: %s", ckpt_path)
        except Exception as exc:
            logger.warning("Auto-checkpoint failed after %s: %s", stage_name, exc)

    def _write_agent2_to_surface(self, state: WorkflowState) -> None:
        surface = self._get_surface()
        if surface is None:
            return
        from camel_agents.surface_model import Agent2Result
        try:
            ads_atoms = read(surface.surface_with_adsorbate_path) if surface.surface_with_adsorbate_path else None
        except Exception:
            ads_atoms = None
        surface.agent2 = Agent2Result(
            best_structure=ads_atoms or Atoms(),
            adsorption_energy=0.0,
            surface_indices=list(surface.surface_indices or []),
            adsorbate_indices=list(surface.adsorbate_indices or []),
            best_config_path=surface.surface_with_adsorbate_path,
        )
        if ads_atoms is not None:
            surface.update_endpoint(ads_atoms, list(surface.surface_indices or []), list(surface.adsorbate_indices or []))
        logger.info("CatalystSurface: Agent2 result attached")
        self._auto_checkpoint("agent2")

    def _write_agent3a_to_surface(self, state: WorkflowState) -> None:
        surface = self._get_surface()
        if surface is None:
            return
        from camel_agents.surface_model import Agent3aResult, Agent3Result
        try:
            slab = read(surface.clean_slab_path) if surface.clean_slab_path else None
        except Exception:
            slab = None
        a3a = Agent3aResult(
            relaxed_slab=slab or Atoms(),
            relaxed_slab_path=surface.clean_slab_path,
            surface_indices=list(surface.surface_indices or []),
        )
        if surface.agent3 is None:
            surface.agent3 = Agent3Result()
        surface.agent3.agent3a = a3a
        logger.info("CatalystSurface: Agent3a result attached")
        self._auto_checkpoint("agent3a")

    def _write_agent3b_to_surface(self, state: WorkflowState, mc_result_pkl: str) -> None:
        surface = self._get_surface()
        if surface is None:
            return
        from camel_agents.surface_model import Agent3bResult, Agent3Result
        try:
            recon = read(surface.reconstructed_surface_path) if surface.reconstructed_surface_path else None
        except Exception:
            recon = None

        pre_mc_surf = list(surface.agent2.surface_indices) if surface.agent2 else []
        pre_mc_ads = list(surface.agent2.adsorbate_indices) if surface.agent2 else []
        pre_mc_count = len(pre_mc_surf) + len(pre_mc_ads)
        mc_added = [i for i in (surface.surface_indices or []) if i >= pre_mc_count]

        a3b = Agent3bResult(
            reconstructed_surface=recon or Atoms(),
            post_mc_surface_indices=list(surface.surface_indices or []),
            post_mc_adsorbate_indices=list(surface.adsorbate_indices or []),
            pre_mc_surface_indices=pre_mc_surf,
            pre_mc_adsorbate_indices=pre_mc_ads,
            mc_added_adatom_indices=mc_added,
            result_pkl_path=mc_result_pkl,
            reconstructed_path=surface.reconstructed_surface_path,
        )
        if surface.agent3 is None:
            surface.agent3 = Agent3Result()
        surface.agent3.agent3b = a3b
        if recon is not None:
            surface.update_endpoint(recon, list(surface.surface_indices or []), list(surface.adsorbate_indices or []))

        # ── Freeze the surface and register the first intermediate ──
        self._freeze_surface_after_agent3(state, recon)

        logger.info("CatalystSurface: Agent3b result attached (mc_added=%d adatoms)", len(mc_added))
        self._auto_checkpoint("agent3b")

    def _freeze_surface_after_agent3(
        self,
        state: WorkflowState,
        reconstructed_atoms: Optional[Atoms],
    ) -> None:
        """Called at the end of agent3: freeze surface_indices + first_adsorbate_indices
        + top_layer_indices, and register the initial adsorbate as the first
        entry in surface.intermediate_adsorbates.

        After this call, the canonical catalyst atoms are locked; subsequent
        code must use ``surface.first_adsorbate_indices`` and
        ``surface.top_layer_indices`` as the source of truth.
        """
        surface = self._get_surface()
        if surface is None or reconstructed_atoms is None:
            return

        if surface.is_frozen:
            return  # already frozen (e.g. resumed run)

        ads_idx = list(surface.adsorbate_indices or surface.endpoint.adsorbate_indices or [])

        top_layer = _identify_top_two_layer_indices(
            reconstructed_atoms,
            exclude_indices=ads_idx,
        )

        surface.freeze(
            top_layer_indices=top_layer,
            first_adsorbate_indices=ads_idx,
        )

        # Register the initial adsorbate as the first intermediate in the chain
        initial_label = ""
        if state.reaction_context is not None:
            initial_label = str(
                getattr(state.reaction_context, "initial_adsorbate", "") or ""
            ).strip()
        if not initial_label:
            syms = reconstructed_atoms.get_chemical_symbols()
            initial_label = "*" + "".join(syms[i] for i in ads_idx if 0 <= i < len(syms))

        if initial_label and not surface.has_intermediate(initial_label):
            try:
                surface.register_intermediate_file(
                    formula=initial_label,
                    atoms=reconstructed_atoms,
                    adsorbate_indices=ads_idx,
                    source="agent3b",
                    agent_key="agent3",
                )
            except Exception as exc:
                logger.warning("Failed to register initial intermediate '%s': %s", initial_label, exc)

        logger.info(
            "Surface frozen: n_surface=%d, n_first_adsorbate=%d, n_top_layer=%d, initial='%s'",
            len(surface.surface_indices),
            len(surface.first_adsorbate_indices),
            len(surface.top_layer_indices),
            initial_label,
        )

    def _write_agent6_to_surface(self, state: WorkflowState, neb_summary: Dict, pathway_result: Any) -> None:
        surface = self._get_surface()
        if surface is None:
            return
        from camel_agents.surface_model import Agent6Result, CompletePathwayInfo

        barriers_fwd: Dict[str, float] = {}
        barriers_rev: Dict[str, float] = {}
        converged_map: Dict[str, bool] = {}
        for step_name, r in neb_summary.items():
            if isinstance(r, dict) and "error" not in r:
                ea_f = r.get("activation_energy_forward")
                ea_r = r.get("activation_energy_reverse")
                if ea_f is not None:
                    barriers_fwd[step_name] = float(ea_f)
                if ea_r is not None:
                    barriers_rev[step_name] = float(ea_r)
                converged_map[step_name] = bool(r.get("converged", False))

        max_b = max(barriers_fwd.values()) if barriers_fwd else None
        rds = max(barriers_fwd, key=barriers_fwd.get) if barriers_fwd else None

        surface.agent6 = Agent6Result(
            neb_results=neb_summary,
            barriers_forward=barriers_fwd,
            barriers_reverse=barriers_rev,
            converged=converged_map,
            rate_determining_step=rds,
            max_barrier=max_b,
        )
        surface.pathway_info = CompletePathwayInfo(
            surface_formula=getattr(pathway_result, "surface_formula", ""),
            intermediates=list(state.intermediates or []),
            adsorbate_energies=dict(state.adsorbate_energies or {}),
            adsorption_energies=dict(state.adsorption_energies or {}),
            rate_determining_step=rds,
            max_barrier=max_b,
            overall_reaction_energy=getattr(pathway_result, "overall_reaction_energy", None),
        )
        surface.intermediates = list(state.intermediates or [])
        if state.reaction_context:
            surface.reaction_context = state.reaction_context
        logger.info("CatalystSurface: Agent6 result attached (max_barrier=%s)", max_b)
        self._auto_checkpoint("agent6")

    def _write_agent7_to_surface(self, state: WorkflowState, kmc_pkl: Optional[str], tof: float, production_rates: Dict) -> None:
        surface = self._get_surface()
        if surface is None:
            return
        from camel_agents.surface_model import Agent7Result
        surface.agent7 = Agent7Result(
            tof=tof if tof else None,
            production_rates=dict(production_rates or {}),
            catmap_result_path=kmc_pkl,
        )
        logger.info("CatalystSurface: Agent7 result attached (tof=%s)", tof)
        self._auto_checkpoint("agent7")

    def _run_single_facet_pipeline(
        self,
        state: WorkflowState,
        config: CatDTConfig,
        surface_path: str,
        facet_id: str,
        agent3_adsorbate_elements_for_mc: List[str],
        agent3_mc_temperature: float,
    ) -> Dict[str, Any]:
        """Run Agent3a→Agent2→Agent3b→Agent4/5→NEB→CatMAP for one facet.

        Mutates *state* in place. Returns a result dict with pathway_result_pkl,
        neb_summary, kmc_result_pkl, tof, and production_rates.
        """
        is_multi = config.multi_facet
        step_prefix = f"facet_{facet_id}/" if is_multi else ""

        # Resume from surface if agents already completed
        surface = self._get_surface()
        _done_3a = surface is not None and surface.agent3 is not None and surface.agent3.agent3a is not None
        _done_2 = surface is not None and surface.agent2 is not None
        _done_3b = surface is not None and surface.agent3 is not None and surface.agent3.agent3b is not None

        # Restore state from surface
        if surface is not None and surface.state is not None:
            for attr in vars(surface.state):
                if not attr.startswith("_"):
                    val = getattr(surface.state, attr, None)
                    if val is not None:
                        setattr(state, attr, val)

        if not _done_3a:
            slab_atoms = read(surface_path)
            surface.surface_path = surface_path
            surface.clean_slab_path = surface_path
            self._update_indices(state, slab_atoms, list(range(len(slab_atoms))), [])
            write(surface_path, slab_atoms)
        state.current_facet_id = facet_id

        # ---- Agent2 output handler (closure over state) ----
        def _apply_agent2_output(_output: TaskOutput):
            record = getattr(_output, "tool_call_record", None)
            if record is None:
                raise RuntimeError("Agent2 did not call `predict_adsorption_sites`.")
            tool_name = getattr(record, "tool_name", None)
            args = getattr(record, "args", {}) or {}
            result = getattr(record, "result", None)
            if tool_name not in {"predict_adsorption_sites", "predict_adsorption_sites_tool"}:
                raise RuntimeError(f"Agent2 called unexpected tool: {tool_name}")

            expected_surface = str(Path(surface.surface_path).resolve())
            actual_surface_raw = str(args.get("surface_path", "")).strip()
            if not actual_surface_raw:
                raise RuntimeError("Agent2 tool args missing `surface_path`.")
            actual_surface = str(Path(actual_surface_raw).resolve())
            if actual_surface != expected_surface:
                raise RuntimeError(
                    f"Agent2 used wrong surface_path: {actual_surface_raw} (expected {surface.surface_path})"
                )
            if str(args.get("adsorbate_smi", "")) != str(state.reaction_context.initial_adsorbate):
                raise RuntimeError("Agent2 used wrong `adsorbate_smi`.")
            if int(args.get("num_sites", -1)) != int(config.num_adsorption_sites):
                raise RuntimeError("Agent2 used wrong `num_sites`.")
            if str(args.get("run_id", "")) != str(config.run_id):
                raise RuntimeError("Agent2 used wrong `run_id`.")

            ads_result_pkl = str(result or "").strip()
            if not ads_result_pkl:
                raise RuntimeError("Agent2 tool returned empty result path.")
            # Reject obviously non-path results (e.g. tool error messages that
            # got returned as strings) before passing to Path.exists(), which
            # raises OSError on long strings.
            if len(ads_result_pkl) > 4096 or ads_result_pkl.startswith("Tool execution failed"):
                raise RuntimeError(f"Agent2 tool returned error instead of result path: {ads_result_pkl[:300]}")
            if not Path(ads_result_pkl).exists():
                raise RuntimeError(f"Agent2 tool result does not exist: {ads_result_pkl}")

            ads_result = self.tools._load_result_from_pickle(ads_result_pkl)
            best_adsorption_config_path = Path(ads_result_pkl).parent / "best_adsorption_config.vasp"
            surface.surface_with_adsorbate_path = str(best_adsorption_config_path)

            ads_atoms = None
            try:
                best_result = getattr(ads_result, "best_result", None)
                final_structure = getattr(best_result, "final_structure", None)
                if isinstance(final_structure, Atoms):
                    ads_atoms = final_structure.copy()
            except Exception:
                ads_atoms = None
            if ads_atoms is None:
                ads_atoms = read(surface.surface_with_adsorbate_path)

            surface_indices, adsorbate_indices = self._infer_surface_adsorbate_indices(
                ads_atoms, fallback_surface_count=surface.n_surface_atoms,
            )
            if not adsorbate_indices and surface.n_surface_atoms is not None:
                adsorbate_indices = list(range(min(surface.n_surface_atoms, len(ads_atoms)), len(ads_atoms)))
                surface_indices = [i for i in range(len(ads_atoms)) if i not in set(adsorbate_indices)]

            self._update_indices(state, ads_atoms, surface_indices, adsorbate_indices)
            write(surface.surface_with_adsorbate_path, ads_atoms)
            self._ensure_clean_slab_consistency(state, ads_atoms, surface_indices, Path(ads_result_pkl).parent)

        # ---- Step 3a: Clean-slab relaxation (NOT MC) ----
        # The VSSR-MC framework (Du et al. 2023 Nat. Comp. Sci.) samples
        # **adsorbate coverage/identity on virtual sites above a pristine
        # surface**. For a CLEAN slab with no chemical adsorbates, there
        # is nothing for VSSR-MC to sample — the framework's "adsorbate"
        # always refers to adatoms sitting above the pristine surface
        # (e.g. GaN Ga adlayer, Si DAS adatoms, CO on Cu), never the
        # substrate top-layer atoms. Running VSSR-MC on a clean slab
        # either adds spurious adatoms (semi-grand) or stalls in
        # prepare_canonical (canonical without physical adsorbates).
        #
        # Instead, we simply UMA-relax the clean slab here. The real MC
        # reconstruction happens at step 3b, with the actual chemical
        # adsorbate from Agent2 driving the virtual-site sampling.
        if config.mc_two_step_reconstruction and config.mc_total_sweeps > 0 and not _done_3a:
            logger.info(
                "[%s] Step 3a: UMA-relaxing clean slab (no MC — clean "
                "surfaces have no adsorbates to sample)",
                facet_id,
            )
            from core.pathway.fairchem_predictor import FairchemPredictor
            fairchem_root_path = str(Path(state.output_base_dir).parent.parent / "deps/fairchem")
            if not Path(fairchem_root_path).exists():
                fairchem_root_path = str(_project_root / "deps/fairchem")
            uma_model_path = str(_project_root / "deps/fairchem_models/uma-s-1p1.pt")

            relax_dir = surface.agent_dir(
                "agent3", subdir=f"{step_prefix.rstrip('/') or 'default'}_03a_clean_relax" if step_prefix else "03a_clean_relax",
            )
            fairchem_relax = FairchemPredictor(
                fairchem_root=fairchem_root_path,
                model_name="uma-s-1p1",
                model_path=uma_model_path,
                use_gpu=True,
                work_dir=str(relax_dir / "uma_relax"),
            )
            clean_slab_atoms = read(surface.surface_path)
            t_relax = time.time()
            relax_result = fairchem_relax.predict_energy(
                clean_slab_atoms, relax=True, fmax=0.05, max_steps=200,
            )
            relaxed_path = str(relax_dir / "relaxed_clean_slab.vasp")
            write(relaxed_path, clean_slab_atoms)
            surface.surface_path = relaxed_path
            surface.clean_slab_path = relaxed_path
            self._update_indices(state, clean_slab_atoms, list(range(len(clean_slab_atoms))), [])
            logger.info(
                "[%s] Step 3a complete: UMA-relaxed clean slab %s (%d atoms, "
                "E=%.3f eV, converged=%s, %.1fs)",
                facet_id, clean_slab_atoms.get_chemical_formula(),
                len(clean_slab_atoms), float(relax_result.energy),
                getattr(relax_result, "converged", True),
                time.time() - t_relax,
            )
            self._write_agent3a_to_surface(state)

        # ---- Step 2: Adsorption ----
        if not _done_2:
            agent2_task = Task(
                description=TaskPrompts.agent2_place_adsorbate(
                    surface_path=surface.surface_path or "",
                    initial_adsorbate=state.reaction_context.initial_adsorbate,
                    num_sites=config.num_adsorption_sites,
                    run_id=config.run_id,
                    workflow_step=f"{step_prefix}02_adsorption",
                    llm_review_context=config.llm_adsorption_review_context if config.enable_llm_adsorption_review else None,
                ),
                expected_output="Adsorption completed",
                agent=self.agent2,
                result_handler=_apply_agent2_output,
            )
            self._run_task(agent2_task)
            self._write_agent2_to_surface(state)
        else:
            logger.info("[%s] Skipping Step 2 (surface.agent2 exists)", facet_id)

        # ---- Step 3b: Reconstruction with adsorbate ----
        if not _done_3b:
            agent3b_ads_atoms = read(surface.surface_with_adsorbate_path)
            agent3b_surface_indices = _identify_top_two_layer_indices(
                agent3b_ads_atoms,
                exclude_indices=surface.adsorbate_indices or [],
            )
            logger.info(
                "[%s] Step 3b: Reconstructing with adsorbate "
                "(top-2 layers = %d, adsorbate = %d / %d total)",
                facet_id, len(agent3b_surface_indices),
                len(surface.adsorbate_indices or []), len(agent3b_ads_atoms),
            )
            agent3b_step_name = "03b_reconstruction_adsorbate" if config.mc_two_step_reconstruction else "03_reconstruction"
            mc_result_pkl = self.tools.simulate_surface_reconstruction(
                surface_with_adsorbate_path=surface.surface_with_adsorbate_path,
                adsorbates_elements_for_mc=agent3_adsorbate_elements_for_mc,
                temperature_k=agent3_mc_temperature,
                total_sweeps=config.mc_total_sweeps,
                run_id=config.run_id,
                workflow_step=f"{step_prefix}{agent3b_step_name}",
                clean_slab_path=surface.clean_slab_path,
                surface_indices=agent3b_surface_indices,
                adsorbate_indices=surface.adsorbate_indices,
                num_adsorbates_override=getattr(config, "mc_num_adsorbates_override", None),
                chem_pots=getattr(config, "mc_chem_pots", None),
                use_seed_for_virtual_sites=getattr(config, "mc_use_seed_for_virtual_sites", False),
                adsorbate_exclusion_radius_A=getattr(
                    config, "mc_adsorbate_exclusion_radius_A", None
                ),
                existing_atom_exclusion_radius_A=getattr(
                    config, "mc_existing_atom_exclusion_radius_A", None
                ),
            )
            if not Path(mc_result_pkl).exists():
                raise RuntimeError(f"Agent3b tool result does not exist: {mc_result_pkl}")
            reconstructed_path = Path(mc_result_pkl).parent / "lowest_energy_reconstructed_structure.vasp"
            surface.reconstructed_surface_path = str(reconstructed_path)
            recon_atoms = read(surface.reconstructed_surface_path)
            reference_atoms = None
            try:
                reference_atoms = read(surface.surface_with_adsorbate_path) if surface.surface_with_adsorbate_path else None
            except Exception:
                reference_atoms = None
            surface_indices, adsorbate_indices = self._infer_adsorbate_indices_after_reconstruction(
                reference_atoms=reference_atoms,
                reference_adsorbate_indices=surface.adsorbate_indices,
                reconstructed_atoms=recon_atoms,
                fallback_surface_count=surface.n_surface_atoms,
            )
            self._update_indices(state, recon_atoms, surface_indices, adsorbate_indices)
            write(surface.reconstructed_surface_path, recon_atoms)
            base_structure = read(surface.reconstructed_surface_path)
            self._set_surface_indices(base_structure, surface.surface_indices or [], surface.adsorbate_indices or [])
            self._write_agent3b_to_surface(state, mc_result_pkl)
        else:
            logger.info("[%s] Skipping Step 3b (surface.agent3.agent3b exists)", facet_id)
            base_structure = surface.atoms.copy()
            self._set_surface_indices(base_structure, list(surface.surface_indices), list(surface.adsorbate_indices))
            # Resumed run: ensure freeze-state is present (idempotent)
            if not surface.is_frozen:
                self._freeze_surface_after_agent3(state, base_structure)

        # All agent3 completion paths converge here: surface must be frozen.
        if surface is not None and not surface.is_frozen:
            raise RuntimeError(
                "CatalystSurface was not frozen at end of agent3. "
                "This is a programming error — _freeze_surface_after_agent3 should have been called."
            )

        if getattr(config, "stop_after_agent3b", False):
            logger.info("[%s] stop_after_agent3b=True; returning after Agent3b", facet_id)
            return {
                "pathway_result_pkl": None,
                "kmc_result_pkl": None,
                "neb_summary": {},
                "tof": 0.0,
                "production_rates": {},
            }

        # ---- Mechanism search (optional) ----
        if self._should_enable_mechanism_stage(config):
            mechanism_output = self._run_mechanism_search_stage(state, config)
            if mechanism_output and mechanism_output.get("status") == "OK":
                self._load_shortlist_into_state(state, mechanism_output)
                logger.info("[%s] Mechanism search: %d retained pathways", facet_id, len(state.retained_pathways))
            else:
                logger.warning("[%s] Mechanism search failed, using default pathway", facet_id)

        # ---- Agent4/5 pathway iteration ----
        self._run_pathway_iteration(state, base_structure, config)

        # ---- NEB retry if needed ----
        if config.calculate_barriers and state.step_structures and not state.neb_results:
            logger.warning("[%s] No NEB results, retrying step-by-step", facet_id)
            neb_dir = surface.agent_dir(
                "agent6", subdir=f"{step_prefix.rstrip('/')}" if step_prefix else None,
            )
            retry_results: Dict[str, Any] = {}
            for step in state.step_structures:
                try:
                    step_result = self.tools.run_neb_for_steps(
                        steps=[step],
                        output_dir=str(neb_dir / step["name"]),
                        n_frames=7, fmax=0.1, max_steps=200,
                    )
                    retry_results.update(step_result)
                except Exception as exc:
                    logger.warning("Per-step NEB retry failed for %s: %s", step.get("name"), exc)
                    retry_results[step.get("name", "unknown_step")] = {"error": str(exc)}
            state.neb_results = retry_results

        # ---- Build NEB summary ----
        neb_summary: Dict[str, Any] = {}
        for step_name, neb_value in state.neb_results.items():
            if isinstance(neb_value, dict):
                neb_summary[step_name] = {"error": neb_value.get("error", str(neb_value))}
            else:
                neb_summary[step_name] = {
                    "activation_energy_forward": getattr(neb_value, "activation_energy_forward", None),
                    "activation_energy_reverse": getattr(neb_value, "activation_energy_reverse", None),
                    "reaction_energy": getattr(neb_value, "reaction_energy", None),
                    "converged": getattr(neb_value, "converged", None),
                }
        if config.calculate_barriers:
            neb_dir = surface.agent_dir(
                "agent6", subdir=f"{step_prefix.rstrip('/')}" if step_prefix else None,
            )
            with open(neb_dir / "neb_results_summary.json", "w", encoding="utf-8") as f:
                json.dump(neb_summary, f, ensure_ascii=False, indent=2)

        # ---- Adsorption energies + pathway result ----
        working_intermediates = self._get_working_intermediate_sequence(state)
        adsorption_dir = surface.agent_dir(
            "agent45", subdir=f"{step_prefix.rstrip('/')}_adsorption_energies" if step_prefix else "adsorption_energies",
        )

        # Pick a "fixed site" for adsorption-energy computation: the lowest-z
        # position among the FIRST adsorbate's atoms. Use first_adsorbate_indices
        # (frozen at agent3b), NOT the current surface.adsorbate_indices — the
        # latter has been updated by Stage A/B and may contain indices that are
        # out of range for base_structure.
        fixed_site = None
        first_ads = list(surface.first_adsorbate_indices or surface.adsorbate_indices or [])
        valid_first_ads = [i for i in first_ads if 0 <= i < len(base_structure)]
        if valid_first_ads:
            ads_pos = base_structure.get_positions()[valid_first_ads]
            fixed_site = ads_pos[np.argmin(ads_pos[:, 2])].tolist()

        adsorption_intermediates: List[str] = []
        for species in working_intermediates:
            if self._is_adsorbate_species(species) and species not in adsorption_intermediates:
                adsorption_intermediates.append(species)
        if not adsorption_intermediates and working_intermediates:
            adsorption_intermediates = list(dict.fromkeys(working_intermediates))

        adsorption_result = self.tools.compute_adsorption_energies(
            surface_path=surface.reconstructed_surface_path,
            intermediates=adsorption_intermediates,
            output_dir=str(adsorption_dir),
            num_sites=config.num_pathway_sites,
            max_steps=config.adsorption_relax_max_steps,
            fixed_site=fixed_site,
        )
        state.adsorbate_energies = adsorption_result["adsorbate_energies"]
        state.adsorption_energies = adsorption_result["adsorption_energies"]
        state.gas_frequencies = adsorption_result.get("gas_frequencies", {}) or {}

        # Build transition entries from the FULL intermediates chain — not from
        # state.pathway_design.steps which only holds the LAST Stage B step
        # (the sequential loop overwrites it each iteration).
        transition_entries: List[Dict[str, str]] = []
        if len(working_intermediates) >= 2:
            for i in range(len(working_intermediates) - 1):
                transition_entries.append({
                    "step_name": f"{working_intermediates[i]}_to_{working_intermediates[i + 1]}",
                    "reactant": working_intermediates[i],
                    "product": working_intermediates[i + 1],
                })

        steps: List[PathwayStep] = []
        for i, entry in enumerate(transition_entries):
            reactant = str(entry.get("reactant", "")).strip()
            product = str(entry.get("product", "")).strip()
            step_name = str(entry.get("step_name", f"{reactant}_to_{product}"))
            reactant_energy = state.adsorbate_energies.get(reactant)
            product_energy = state.adsorbate_energies.get(product)
            reaction_energy = self._compute_reaction_energy_with_correction(
                reactant, product, reactant_energy, product_energy,
            )
            neb = self._find_neb_result_for_transition(state, reactant, product, i, step_name=step_name)
            activation_energy = None
            ts_energy = None
            neb_rxn_energy = None
            if neb is not None and not isinstance(neb, dict):
                activation_energy = getattr(neb, "activation_energy_forward", None)
                if getattr(neb, "transition_state_index", None) is not None:
                    ts_energy = neb.energies[neb.transition_state_index]
                # Raw NEB reaction energy (reactant and product endpoints share
                # the same atom count — staged atoms included — so this is the
                # correct ΔE for CatMAP's A_s -> B_s + H_s splits).
                neb_rxn_energy = getattr(neb, "reaction_energy", None)
            elif isinstance(neb, dict):
                neb_rxn_energy = neb.get("reaction_energy")
            steps.append(PathwayStep(
                step_index=i, name=step_name,
                reactant_adsorbate=reactant, product_adsorbate=product,
                reactant_energy=reactant_energy or 0.0, product_energy=product_energy or 0.0,
                reaction_energy=reaction_energy or 0.0,
                activation_energy=activation_energy, transition_state_energy=ts_energy,
                neb_reaction_energy=float(neb_rxn_energy) if neb_rxn_energy is not None else None,
            ))

        rate_determining_step = None
        max_barrier = None
        valid_barriers = [s for s in steps if s.activation_energy is not None]
        if valid_barriers:
            rate_determining_step = max(valid_barriers, key=lambda s: s.activation_energy)
            rate_determining_step.is_rate_determining = True
            max_barrier = rate_determining_step.activation_energy

        overall_reaction_energy = 0.0
        if working_intermediates:
            first_e = state.adsorbate_energies.get(working_intermediates[0])
            last_e = state.adsorbate_energies.get(working_intermediates[-1])
            if first_e is not None and last_e is not None:
                overall_reaction_energy = last_e - first_e
            elif steps:
                overall_reaction_energy = float(sum(s.reaction_energy for s in steps))

        pathway_dir = surface.agent_dir(
            "agent45", subdir=f"{step_prefix.rstrip('/')}" if step_prefix else None,
        )
        pathway_result = CompletePathwayResult(
            surface_formula=base_structure.get_chemical_formula(),
            adsorbates=working_intermediates or state.reaction_context.intermediates,
            adsorbate_energies=state.adsorbate_energies,
            adsorption_energies=state.adsorption_energies,
            steps=steps, rate_determining_step=rate_determining_step,
            max_barrier=max_barrier, overall_reaction_energy=overall_reaction_energy,
            output_dir=str(pathway_dir),
            gas_frequencies=dict(getattr(state, "gas_frequencies", {}) or {}),
        )
        pathway_pkl = self.tools._save_result_to_pickle(pathway_result, pathway_dir / "complete_pathway_result.pkl")

        # ---- KMC ----
        kmc_result_pkl = None
        tof = 0.0
        production_rates: Dict[str, float] = {}
        if config.kmc_temperature_k and config.gas_pressures:
            try:
                kmc_result_pkl = self.tools.run_kmc_simulation(
                    pathway_result_path=pathway_pkl,
                    temperature_k=config.kmc_temperature_k,
                    pressures=config.gas_pressures,
                    run_id=config.run_id,
                    workflow_step=f"{step_prefix}05_kmc_simulation",
                )
                kmc_result = self.tools._load_result_from_pickle(kmc_result_pkl)
                # SinglePointResult exposes TOF via get_turnover_frequency();
                # older KMC backends expose it as `.tof`. Prefer the callable
                # when available so CatMAP production rates are surfaced
                # correctly.
                tof_attr = getattr(kmc_result, "tof", None)
                if callable(getattr(kmc_result, "get_turnover_frequency", None)):
                    try:
                        tof = float(kmc_result.get_turnover_frequency())
                    except Exception:
                        tof = float(tof_attr or 0.0)
                else:
                    tof = float(tof_attr or 0.0)
                production_rates = dict(getattr(kmc_result, "production_rates", {}) or {})
            except Exception as exc:
                logger.warning("[%s] KMC simulation failed: %s", facet_id, exc)

        # ---- Populate CatalystSurface: Agent6 + Agent7 ----
        self._write_agent6_to_surface(state, neb_summary, pathway_result)
        self._write_agent7_to_surface(state, kmc_result_pkl, tof, production_rates)

        return {
            "pathway_result_pkl": pathway_pkl,
            "kmc_result_pkl": kmc_result_pkl,
            "neb_summary": neb_summary,
            "tof": tof,
            "production_rates": production_rates,
        }

    def _run_multi_facet_pipeline(
        self,
        state: WorkflowState,
        config: CatDTConfig,
        agent3_adsorbate_elements_for_mc: List[str],
        agent3_mc_temperature: float,
    ) -> None:
        """Run the per-facet pipeline for each top-N Wulff surface and aggregate."""
        for facet_id, surface_path in state.all_surface_paths.items():
            miller_index = state.all_miller_indices.get(facet_id)
            area_fraction = state.all_area_fractions.get(facet_id, 0.0)
            logger.info(
                "===== Multi-facet: processing facet %s (area=%.1f%%) =====",
                facet_id, area_fraction * 100,
            )

            # Deep copy state so per-facet mutations are isolated
            facet_state = copy.deepcopy(state)
            # Clear per-facet mutable fields to avoid stale data
            facet_state.step_structures = []
            facet_state.tool_baseline_steps = []
            facet_state.neb_results = {}
            facet_state.pathway_design = None
            facet_state.validation_report = None
            facet_state.energy_gate_report = {}
            facet_state.last_feedback = ""
            facet_state.iteration_history = []

            try:
                result = self._run_single_facet_pipeline(
                    state=facet_state,
                    config=config,
                    surface_path=surface_path,
                    facet_id=facet_id,
                    agent3_adsorbate_elements_for_mc=agent3_adsorbate_elements_for_mc,
                    agent3_mc_temperature=agent3_mc_temperature,
                )
                state.facet_results.append(FacetResult(
                    facet_id=facet_id,
                    miller_index=miller_index,
                    area_fraction=area_fraction,
                    surface_path=surface_path,
                    pathway_result_pkl=result.get("pathway_result_pkl"),
                    kmc_result_pkl=result.get("kmc_result_pkl"),
                    tof=result.get("tof", 0.0),
                    production_rates=result.get("production_rates", {}),
                    neb_summary=result.get("neb_summary", {}),
                ))
            except Exception as exc:
                logger.error("Facet %s failed: %s", facet_id, exc, exc_info=True)
                state.facet_results.append(FacetResult(
                    facet_id=facet_id,
                    miller_index=miller_index,
                    area_fraction=area_fraction,
                    surface_path=surface_path,
                    error=str(exc),
                ))

        self._aggregate_wulff_tof(state)

    @staticmethod
    def _aggregate_wulff_tof(state: WorkflowState) -> None:
        """Compute Wulff-weighted TOF: TOF_total = Σ (a_i / Σa_ok) × TOF_i."""
        successful = [f for f in state.facet_results if f.error is None and f.tof > 0]
        if not successful:
            state.aggregated_tof = 0.0
            return

        total_area = sum(f.area_fraction for f in successful)
        if total_area <= 0:
            state.aggregated_tof = 0.0
            return

        state.aggregated_tof = sum(
            (f.area_fraction / total_area) * f.tof for f in successful
        )

        all_species: set = set()
        for f in successful:
            all_species.update(f.production_rates.keys())
        state.aggregated_production_rates = {
            sp: sum(
                (f.area_fraction / total_area) * f.production_rates.get(sp, 0.0)
                for f in successful
            )
            for sp in all_species
        }
        logger.info(
            "Wulff-weighted TOF = %.4e (from %d/%d facets)",
            state.aggregated_tof, len(successful), len(state.facet_results),
        )

    def run_on_surface(self, surface: "CatalystSurface", user_input_config: dict) -> "CatalystSurface":
        """Run the full pipeline with CatalystSurface as the central data carrier.

        Each agent stage writes its result to surface.agentX in real time.
        """
        from camel_agents.surface_model import Agent1Result

        surface.output_dir = str(user_input_config.get("output_base_dir", ""))
        surface.reaction_description = str(user_input_config.get("reaction_description", ""))

        if surface.agent1 is None and surface.clean_slab is not None:
            surface.agent1 = Agent1Result(
                clean_slab=surface.clean_slab,
                slab_path=user_input_config.get("initial_surface_path"),
                surface_indices=list(surface.surface_indices),
            )

        self._catalyst_surface = surface
        if self.tools is not None:
            self.tools._current_surface = surface
        try:
            self.run(user_input_config)
        finally:
            self._catalyst_surface = None
            if self.tools is not None:
                self.tools._current_surface = None

        ckpt_path = Path(surface.output_dir) / "surface_checkpoint.pkl"
        try:
            surface.save_checkpoint(ckpt_path)
            logger.info("Saved CatalystSurface checkpoint to %s", ckpt_path)
        except Exception as exc:
            logger.warning("Failed to save checkpoint: %s", exc)

        return surface

    def run(self, user_input_config: dict) -> dict:
        config = CatDTConfig.create_from_args(**user_input_config)
        if not config.reaction_description:
            raise ValueError("reaction_description is required")
        if not config.bulk_structure_path and not config.initial_surface_path:
            raise ValueError("bulk_structure_path or initial_surface_path is required")

        desired_output_dir = Path(config.output_base_dir).resolve()
        current_output_dir = None
        try:
            current_output_dir = Path(getattr(self.tools, "output_base_dir", "")).resolve()
        except Exception:
            current_output_dir = None

        # Sync electrochemical / model params into existing tools instance
        if self.tools is not None:
            self.tools.mc_energy_model = getattr(config, "mc_energy_model", "CHGNet")
            self.tools.potential_she = getattr(config, "potential_she", None)
            self.tools.ph = getattr(config, "ph", None)
            self.tools.adsorption_backend = getattr(config, "adsorption_backend", "adsorbml")
            self.tools.adsorbml_num_sites = getattr(config, "adsorbml_num_sites", 20)
            self.tools.adsorbml_placement_mode = getattr(config, "adsorbml_placement_mode", "random_site_heuristic_placement")
            self.tools.adsorbml_interstitial_gap = getattr(config, "adsorbml_interstitial_gap", 0.1)
            self.tools.adsorbml_relax_steps = getattr(config, "adsorbml_relax_steps", 200)
            self.tools.adsorbml_relax_fmax = getattr(config, "adsorbml_relax_fmax", 0.05)
            if getattr(config, "mp_api_key", None):
                self.tools.mp_api_key = config.mp_api_key
                os.environ["MP_API_KEY"] = config.mp_api_key
            # Force re-creation of DT instance to pick up new params
            self.tools.dt_instance = None

        if self.tools is None or current_output_dir != desired_output_dir:
            self.tools = CatDTTools(
                output_base_dir=config.output_base_dir,
                mc_energy_model=getattr(config, "mc_energy_model", "CHGNet"),
                potential_she=getattr(config, "potential_she", None),
                ph=getattr(config, "ph", None),
                mp_api_key=getattr(config, "mp_api_key", None),
                adsorption_backend=getattr(config, "adsorption_backend", "adsorbml"),
                adsorbml_num_sites=getattr(config, "adsorbml_num_sites", 20),
                adsorbml_placement_mode=getattr(config, "adsorbml_placement_mode", "random_site_heuristic_placement"),
                adsorbml_interstitial_gap=getattr(config, "adsorbml_interstitial_gap", 0.1),
                adsorbml_relax_steps=getattr(config, "adsorbml_relax_steps", 200),
                adsorbml_relax_fmax=getattr(config, "adsorbml_relax_fmax", 0.05),
            )
            self._init_camel_tools()
            if all(agent is not None for agent in [self.agent1, self.agent2, self.agent3, self.agent4, self.agent5, self.agent6, self.agent7]):
                self._init_agents()
            # Re-link surface reference after tools recreated
            self.tools._current_surface = self._get_surface()

        policy_path = Path(config.output_base_dir) / "evolvability_policy.json"
        self.evolvable_policy = EvolvablePolicy(policy_path=policy_path)
        self.current_strategy = self.evolvable_policy.choose_strategy()
        logger.info("Evolvable strategy selected: %s", self.current_strategy)

        surface = self._get_surface()
        # Auto-create a default CatalystSurface if the caller invoked run() directly
        # (i.e. without run_on_surface). Physical state lives on surface now, so it
        # must exist for any pipeline method that reads it.
        if surface is None:
            from camel_agents.surface_model import CatalystSurface, StructureEndpoint
            from ase import Atoms
            default_atoms = Atoms()
            default_ep = StructureEndpoint.from_parts(default_atoms, [], [])
            surface = CatalystSurface(
                surface_id=str(config.run_id or "run"),
                endpoint=default_ep,
                output_dir=str(Path(config.output_base_dir) / (config.run_id or "")),
            )
            self._catalyst_surface = surface
            if self.tools is not None:
                self.tools._current_surface = surface

        if surface.state is not None and surface.surface_path:
            state = surface.state
            state.run_id = config.run_id
            state.output_base_dir = config.output_base_dir
            state.reaction_description = config.reaction_description
            logger.info("Restored WorkflowState from surface checkpoint (surface_path=%s)", surface.surface_path)
        else:
            state = WorkflowState(
                run_id=config.run_id,
                output_base_dir=config.output_base_dir,
                reaction_description=config.reaction_description,
            )
            surface.state = state

        def _apply_agent1_output(_output: TaskOutput):
            if config.initial_surface_path:
                surface.surface_path = config.initial_surface_path
                slab_atoms = read(config.initial_surface_path)
            else:
                try:
                    record = getattr(_output, "tool_call_record", None)
                    if record is None:
                        raise RuntimeError("Agent1 did not call `generate_surfaces`.")
                    tool_name = getattr(record, "tool_name", None)
                    args = getattr(record, "args", {}) or {}
                    result = getattr(record, "result", None)
                    if tool_name not in {"generate_surfaces", "generate_surfaces_tool"}:
                        raise RuntimeError(f"Agent1 called unexpected tool: {tool_name}")

                    expected_bulk = str(Path(config.bulk_structure_path).resolve())
                    actual_bulk_raw = str(args.get("bulk_structure_path", "")).strip()
                    if not actual_bulk_raw:
                        raise RuntimeError("Agent1 tool args missing `bulk_structure_path`.")
                    actual_bulk = str(Path(actual_bulk_raw).resolve())
                    if actual_bulk != expected_bulk:
                        raise RuntimeError(
                            f"Agent1 used wrong bulk_structure_path: {actual_bulk_raw} (expected {config.bulk_structure_path})"
                        )
                    if int(args.get("top_n_surfaces", -1)) != int(config.top_n_surfaces):
                        raise RuntimeError("Agent1 used wrong `top_n_surfaces`.")
                    if str(args.get("run_id", "")) != str(config.run_id):
                        raise RuntimeError("Agent1 used wrong `run_id`.")
                    if str(args.get("workflow_step", "")) != "01_surfaces":
                        raise RuntimeError("Agent1 used wrong `workflow_step`.")

                    surff_result_pkl = str(result or "").strip()
                    if not surff_result_pkl:
                        raise RuntimeError("Agent1 tool returned empty result path.")
                    if not Path(surff_result_pkl).exists():
                        raise RuntimeError(f"Agent1 tool result does not exist: {surff_result_pkl}")

                    surff_result = self.tools._load_result_from_pickle(surff_result_pkl)
                    top_surface = surff_result.get_top_n(1)[0]
                    miller_str = "".join(map(str, top_surface.miller_index))
                    surface.surface_path = str(Path(surff_result_pkl).parent / f"slab_{miller_str}.vasp")
                    slab_atoms = read(surface.surface_path)

                    # Store all top-N surfaces for multi-facet mode
                    n_facets = config.multi_facet_top_n if config.multi_facet else 1
                    for surf in surff_result.get_top_n(n_facets):
                        ms = "".join(map(str, surf.miller_index))
                        slab_path = str(Path(surff_result_pkl).parent / f"slab_{ms}.vasp")
                        if Path(slab_path).exists():
                            state.all_surface_paths[ms] = slab_path
                            state.all_area_fractions[ms] = surf.area_fraction
                            state.all_miller_indices[ms] = surf.miller_index

                except Exception as exc:
                    raise RuntimeError(f"Agent1 surface generation failed: {exc}") from exc

            surface.clean_slab_path = surface.surface_path
            self._update_indices(state, slab_atoms, list(range(len(slab_atoms))), [])
            write(surface.surface_path, slab_atoms)

        _skip_agent1 = surface is not None and surface.agent1 is not None and surface.surface_path
        if _skip_agent1:
            logger.info("Skipping Agent1 (restored from surface checkpoint, surface_path=%s)", surface.surface_path)
        else:
            agent1_task = Task(
                description=TaskPrompts.agent1_initialize_surface(
                    reaction_description=config.reaction_description,
                    bulk_structure_path=config.bulk_structure_path,
                    run_id=config.run_id,
                    top_n_surfaces=config.top_n_surfaces,
                    workflow_step="01_surfaces",
                    initial_surface_path=config.initial_surface_path,
                ),
                expected_output="Surface initialized",
                agent=self.agent1,
                result_handler=_apply_agent1_output,
            )
            self._run_task(agent1_task)

        if getattr(state, 'reaction_context', None) and hasattr(state.reaction_context, 'intermediates') and len(state.reaction_context.intermediates) >= 2:
            logger.info("Skipping reaction context parsing (restored from checkpoint, %d intermediates)", len(state.reaction_context.intermediates))
        else:
            reaction_task = self._create_reaction_context_task(state)
            self._run_task(reaction_task)
            if not state.reaction_context:
                raise RuntimeError("Agent4 reaction-context task did not produce valid context.")
            self._normalize_reaction_context(state)

        if len(state.reaction_context.intermediates) < 2:
            raise RuntimeError("Agent4 reaction-context must include at least two intermediates.")

        # --- Determine MC parameters (shared by all facets) ---
        agent3_adsorbate_elements_for_mc = config.adsorbate_elements_for_mc
        if not agent3_adsorbate_elements_for_mc:
            slab = read(surface.clean_slab_path)
            symbols = slab.get_chemical_symbols()
            if not symbols:
                raise ValueError("Could not infer adsorbate_elements_for_mc from clean slab")
            main_symbol = max(set(symbols), key=symbols.count)
            agent3_adsorbate_elements_for_mc = [main_symbol]
        agent3_mc_temperature = config.mc_temperature_k if config.mc_temperature_k is not None else 500.0

        # ================================================================
        # Branch: multi-facet vs single-facet pipeline
        # ================================================================
        if config.multi_facet and len(state.all_surface_paths) > 1:
            logger.info(
                "Multi-facet mode: processing %d facets %s",
                len(state.all_surface_paths), list(state.all_surface_paths.keys()),
            )
            self._run_multi_facet_pipeline(
                state=state,
                config=config,
                agent3_adsorbate_elements_for_mc=agent3_adsorbate_elements_for_mc,
                agent3_mc_temperature=agent3_mc_temperature,
            )

            # Use the best facet's pathway for the final report
            best_facet = max(
                (f for f in state.facet_results if f.error is None),
                key=lambda f: f.area_fraction,
                default=None,
            )
            pathway_pkl = best_facet.pathway_result_pkl if best_facet else None
            kmc_result_pkl = best_facet.kmc_result_pkl if best_facet else None

            results = {
                "multi_facet": True,
                "facet_results": [
                    {
                        "facet_id": f.facet_id,
                        "miller_index": f.miller_index,
                        "area_fraction": f.area_fraction,
                        "tof": f.tof,
                        "production_rates": f.production_rates,
                        "neb_summary": f.neb_summary,
                        "error": f.error,
                    }
                    for f in state.facet_results
                ],
                "aggregated_tof": state.aggregated_tof,
                "aggregated_production_rates": state.aggregated_production_rates,
                "pathway_result_pkl": pathway_pkl,
                "selected_strategy": self.current_strategy,
                "pathway_iterations": self._last_iteration_count,
            }

        else:
            # Single-facet mode (original behavior)
            facet_result = self._run_single_facet_pipeline(
                state=state,
                config=config,
                surface_path=surface.surface_path,
                facet_id="primary",
                agent3_adsorbate_elements_for_mc=agent3_adsorbate_elements_for_mc,
                agent3_mc_temperature=agent3_mc_temperature,
            )

            if getattr(config, "stop_after_agent3b", False):
                return {
                    "stopped_after": "agent3b",
                    "selected_strategy": self.current_strategy,
                    "surface_with_adsorbate_path": surface.surface_with_adsorbate_path,
                    "reconstructed_surface_path": surface.reconstructed_surface_path,
                    "agent3b_checkpoint": str(
                        Path(surface.output_dir)
                        / "checkpoints"
                        / "checkpoint_after_agent3b.pkl"
                    ),
                }

            pathway_pkl = facet_result["pathway_result_pkl"]
            kmc_result_pkl = facet_result.get("kmc_result_pkl")

            results = {
                "pathway_result_pkl": pathway_pkl,
                "neb_results_summary": facet_result["neb_summary"],
                "selected_strategy": self.current_strategy,
                "pathway_iterations": self._last_iteration_count,
            }
            if kmc_result_pkl:
                results["kmc_result_pkl"] = kmc_result_pkl

        # ================================================================
        # Final report + policy update (shared by both modes)
        # ================================================================
        final_report = self.tools.generate_final_report(
            run_id=config.run_id,
            surff_result_path=None,
            adsorbdiff_result_path=None,
            vssr_mc_result_path=None,
            pathway_result_path=pathway_pkl,
            kmc_result_path=kmc_result_pkl,
        )
        results["final_report_path"] = final_report

        reward = self._compute_episode_reward(state)
        if self.evolvable_policy is not None:
            self.evolvable_policy.update(strategy=self.current_strategy, reward=reward)

        results["episode_reward"] = reward
        results["validation_status"] = state.validation_report.status if state.validation_report else "UNKNOWN"

        return results


    def __getattr__(self, name: str) -> Any:
        tools = self.__dict__.get("tools")
        if tools is not None:
            try:
                return getattr(tools, name)
            except AttributeError:
                pass
        raise AttributeError(f"{self.__class__.__name__} has no attribute {name}")


def main(user_input_config: dict):
    workflow = CatDTCamelWorkflow()
    return workflow.run(user_input_config)


__all__ = [
    "CatDTCamelWorkflow",
    "EvolvablePolicy",
    "get_camel_model_backend",
    "main",
    "CatalystSurface",
]
