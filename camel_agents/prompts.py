"""
CatDT Agent Prompts (CAMEL-ready)

Design principles:
1. Prompts are fully generic — no reaction-specific rules or examples.
2. Agent4/Agent5 collaborate via multi-round feedback (up to N rounds).
3. No hardcoded element add/remove rules — LLM decides based on input.
4. Full structure coordinates given once; subsequent rounds only show adsorbate coords.

Backup of previous Chinese version: prompts_zh_backup.py
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class AgentRole:
    name: str
    role: str
    goal: str
    backstory: str

    def to_dict(self) -> Dict[str, str]:
        return {
            "role": self.role,
            "goal": self.goal,
            "backstory": self.backstory,
        }


class TaskPrompts:
    """任务提示词模板"""

    @staticmethod
    def reaction_context_parsing(reaction_description: str) -> str:
        return f"""[TASK] Parse a catalytic reaction description into a deterministic pathway.

[INPUT REACTION DESCRIPTION]
{reaction_description}

[HOW TO EXTRACT initial_adsorbate]

initial_adsorbate is the adsorbed species that carries the main molecular
skeleton of the reaction onto the surface. Pick it with this procedure:

1. Scan all elementary reactions top-to-bottom. Classify each reaction as
   either:
   (a) a "bonding-partner activation" step: its left-hand side contains only
       bare gas-phase molecules and empty sites, and its right-hand side is
       a set of atomic single-element adsorbates acting as bonding partners
       for the main chain; or
   (b) a "main-skeleton" step: at least one left-hand-side species is a
       polyatomic adsorbate that will be transformed into the next main
       intermediate.

2. initial_adsorbate is the main-skeleton adsorbate on the LHS of the
   FIRST main-skeleton step (step 1b). If that species appears in gas phase
   form in the reactant, normalize it to its surface-adsorbed label.

3. Normalize the label: use the leading-asterisk form (e.g. "XY*" → "*XY").

[HOW TO BUILD intermediates]

intermediates is the ordered list of main-skeleton adsorbed species that the
reactant passes through on its way to the terminal product. Rules:

1. intermediates[0] MUST equal initial_adsorbate.

2. intermediates[-1] is the terminal main-skeleton product (adsorbed or
   desorbing into gas phase) of the main linear pathway.

3. For each step i: intermediates[i+1] is the main-skeleton species that
   results from one elementary reaction consuming intermediates[i] (add or
   remove a bonding partner). Follow the chain deterministically — pick one
   branch, do not include parallel side branches.

4. Do NOT list bare atomic bonding partners as standalone intermediates;
   they are consumed in each main-skeleton step but are not pathway endpoints.

5. Do NOT list side-branch products that split off from the main skeleton,
   unless the user explicitly requested them as the terminal product.

6. If the user's reaction description already enumerates the main-skeleton
   chain explicitly (e.g. "A -> B -> C -> D"), use that chain verbatim and
   only normalize the labels.

[OTHER FIELDS]

* reaction_type = one short sentence describing the chemistry class.
* constraints = any explicit user constraints (operating conditions, surface
  type, excluded branches, etc.). Keep it short.

[OUTPUT FORMAT — JSON ONLY, NO PROSE]
{{
  "reaction_type": "...",
  "initial_adsorbate": "*X",
  "intermediates": ["*X", "*Y", "*Z", "..."],
  "constraints": "..."
}}
"""

    @staticmethod
    def preplacement_design(
        current_label: str,
        target_label: str,
        current_adsorbate_coords: str,
        target_molecule_coords: str,
        surface_info: str,
        element_delta: str,
    ) -> str:
        """Agent4 preplacement prompt: build next intermediate from current structure."""
        return f"""[TASK] Build the next reaction intermediate on the surface.

You are transforming the current surface adsorbate into the next intermediate
in a reaction pathway. The current adsorbate is already relaxed on the surface.
You must add or remove atoms to convert it into the target intermediate.

═══════════════════════════════════════════════════════════
SECTION 1 — WHAT TO DO
═══════════════════════════════════════════════════════════
Current adsorbate: {current_label}
Target adsorbate:  {target_label}
Element changes needed: {element_delta}

For each atom to ADD:
  • Look at the target molecule's 3D structure to understand WHERE that atom
    should go relative to the existing adsorbate atoms.
  • Use ``suggest_staged_positions_tool`` to get candidate surface sites near
    the adsorbate.
  • Place the new atom at a chemically reasonable position: at bond distance
    from the atom it should bond to, oriented like in the target molecule.

For each atom to REMOVE:
  • Compare current adsorbate elements with target molecule elements.
  • Identify which atom(s) in the current adsorbate should be removed.
  • Remove the atom furthest from the adsorbate center of the matching element.
  • Specify the atom INDEX (from the current adsorbate coordinates above).

═══════════════════════════════════════════════════════════
SECTION 2 — REFERENCE STRUCTURES
═══════════════════════════════════════════════════════════
[Current adsorbate on surface]
{current_adsorbate_coords}

[Target molecule (from database, gas-phase geometry)]
{target_molecule_coords}

[Surface information]
{surface_info}

═══════════════════════════════════════════════════════════
SECTION 3 — WORKFLOW (MUST follow all steps)
═══════════════════════════════════════════════════════════
1. Examine the element delta. Decide which atoms to add/remove.
2. For atoms to add: call ``suggest_staged_positions_tool`` to get candidate
   positions near the adsorbate, then choose positions consistent with the
   target molecule's bonding geometry.
3. For atoms to remove: specify indices of atoms to remove.
4. Call ``validate_proposed_positions_tool`` to check your proposal.
5. Output the final JSON.

IMPORTANT: After outputting the JSON, the system will automatically apply your
changes and relax the structure using UMA. You do NOT need to call a relax tool.

═══════════════════════════════════════════════════════════
SECTION 4 — OUTPUT FORMAT
═══════════════════════════════════════════════════════════
{{
  "pathway_name": "{current_label}_to_{target_label}",
  "description": "Transform {current_label} to {target_label}",
  "overall_reaction": "{current_label} -> {target_label}",
  "steps": [
    {{
      "step_name": "{current_label}_to_{target_label}",
      "reactant_formula": "{current_label}",
      "product_formula": "{target_label}",
      "reaction_type": "preplacement",
      "atoms_to_add": [
        {{"species": "H", "position": [x, y, z], "reason": "add H at bond distance from C"}}
      ],
      "atoms_to_modify": [],
      "atoms_to_remove": [
        {{"index": 110, "reason": "remove O atom to match target formula"}}
      ],
      "confidence": 0.0
    }}
  ],
  "confidence": 0.0
}}"""

    @staticmethod
    def pathway_design(
        reaction_context: Dict[str, Any],
        surface_info: str,
        baseline_steps: str,
        memento_context: str = "",
        knowledge_context: str = "",
        skill_context: str = "",
        feedback: str = "",
        history: str = "",
    ) -> str:
        """Agent4 路径端点修正提示词。

        Args:
            reaction_context: 反应上下文 dict (reaction_type, intermediates, constraints)
            surface_info: 表面概要信息（表面顶部 z 坐标、晶胞参数等）
            baseline_steps: 每步吸附分子坐标与元素差异（compact 格式）
            feedback: 上一轮 Agent5 反馈
            history: 迭代历史摘要
        """
        feedback_section = f"\n【上一轮反馈（含前一轮吸附分子坐标）】\n{feedback}\n" if feedback else ""
        history_section = f"\n【迭代历史】\n{history}\n" if history else ""
        memento_section = memento_context if memento_context else "(no retrieved memento cases)"
        knowledge_section = knowledge_context if knowledge_context else "(no retrieved knowledge items)"
        skill_section = skill_context if skill_context else "(no retrieved skill items)"

        return f"""[TASK] Agent4 — Place Co-adsorbate Atoms for ONE Reaction Step

You are processing ONE elementary step in a sequential NEB pipeline.
The reactant structure is the relaxed product from the previous step (or the
initial adsorbed surface for the first step). The product structure comes from
tool-based adsorbate preplacement. Both share the SAME surface geometry.

Your job: equalise element counts between reactant and product by adding
co-adsorbate atoms at physically correct surface adsorption sites so that NEB
can compute a meaningful activation barrier.

═══════════════════════════════════════════════════════════
SECTION 1 — CO-ADSORBATE STRATEGY
═══════════════════════════════════════════════════════════
NEB requires IDENTICAL element multisets on both endpoints. The element_delta
tells you exactly what to add:
• element_delta = +X → add X to the **reactant** (atom joins the adsorbate in the product)
• element_delta = -X → add X to the **product** (atom departs the adsorbate)
• element_delta = none → no atoms needed, endpoints already balanced

Co-adsorbates MUST sit at the nearest-neighbor surface adsorption site
(hollow/bridge) in the bonding direction from the anchor atom. Typical distance:
1.5–4.5 Å from the bonding anchor. Both endpoints must be true local energy minima.

Example: *HCOO* + H* → *HCOOH* + *
  product has +1H vs reactant → add 1H to reactant as co-adsorbate at a surface
  hollow site near the HCOO* adsorbate.

═══════════════════════════════════════════════════════════
SECTION 2 — PLACEMENT RULES
═══════════════════════════════════════════════════════════
1. USE [hint] POSITIONS — each step with element_delta ≠ none shows [hint] lines
   with up to 3 candidate co-adsorption sites sorted nearest-first.
   Site 1 is the default. Copy the [x,y,z] directly.

2. INCOMING atoms (element_delta = +X, added to reactant):
   • Place at site 1 on first attempt; switch to site 2/3 only after failure feedback.
   • Site 1 is the nearest-neighbor surface site in the bonding direction.
   • Height above surface: use the [hint] z-coordinate (already correct).
   • Distance from bonding anchor atom: 1.5–4.5 Å.
   • Must NOT overlap with existing atoms (> 0.8 Å from all atoms).

3. OUTGOING atoms (element_delta = -X, added to product):
   • Same rules. For multi-atom fragments (e.g., OH): place anchor at [hint] site,
     H at ~0.98 Å above.

4. HARD LIMITS:
   • Position within 1.5 Å of a [hint] site
   • Distance from bonding anchor: 1.5–5.0 Å
   • All atom–atom distances > 0.8 Å
   • ABOVE top surface z (never inside the slab)
   • Co-adsorbate position displaced ≥ 2.0 Å from its bonded position on the other endpoint

═══════════════════════════════════════════════════════════
SECTION 3 — INPUT DATA (this step only)
═══════════════════════════════════════════════════════════
[Reaction Context]
reaction_type: {reaction_context.get("reaction_type", "")}
intermediates: {reaction_context.get("intermediates", [])}
constraints: {reaction_context.get("constraints", "")}

[Surface Information]
{surface_info}

[Current Step — Baseline Structure, Element Delta & Co-adsorption Hints]
{baseline_steps}

[Memento CaseBank]
{memento_section}

[KnowledgeBank]
{knowledge_section}

[SkillBank]
{skill_section}
{feedback_section}{history_section}
═══════════════════════════════════════════════════════════
SECTION 4 — WORKFLOW
═══════════════════════════════════════════════════════════
1. Read the element_delta for this step. If none → output empty atoms_to_add.
2. Propose [x,y,z] positions based on the [hint] lines and your chemical knowledge.
3. Call ``validate_proposed_positions_tool(steps=[...])`` to check your proposal.
4. If validation fails, call ``suggest_staged_positions_tool(requests=[...])``
   for alternatives, then re-validate.
5. Output the PathwayDesign JSON.

You MUST call ``validate_proposed_positions_tool`` before outputting the final JSON.

═══════════════════════════════════════════════════════════
SECTION 5 — HARD CONSTRAINTS
═══════════════════════════════════════════════════════════
1. Existing atom coordinates are READ-ONLY.
2. Output ONLY atoms_to_add. atoms_to_modify and atoms_to_remove MUST be [].
3. Every step with element_delta ≠ none MUST have atoms_to_add with [x,y,z].
4. Steps with element_delta = none → atoms_to_add = [].

═══════════════════════════════════════════════════════════
SECTION 6 — OUTPUT FORMAT
═══════════════════════════════════════════════════════════
{{
  "pathway_name": "...",
  "description": "...",
  "overall_reaction": "...",
  "steps": [
    {{
      "step_name": "...",
      "reactant_formula": "...",
      "product_formula": "...",
      "reaction_type": "...",
      "atoms_to_add": [
        {{
          "species": "H",
          "position": [x, y, z],
          "reason": "H co-adsorbate at hollow site 2.8 Å from HCOO centroid"
        }}
      ],
      "atoms_to_modify": [],
      "atoms_to_remove": [],
      "confidence": 0.0
    }}
  ],
  "confidence": 0.0
}}"""

    @staticmethod
    def pathway_validation(
        steps_details: str,
        structures_payload: str,
        precheck_summary: str = "",
        precheck_issues: str = "",
        memento_context: str = "",
        knowledge_context: str = "",
        skill_context: str = "",
    ) -> str:
        """Agent5 端点验证提示词。

        Args:
            steps_details: 每步摘要（step/formula/adds/removes）
            structures_payload: 每步吸附分子坐标（compact 格式）
            precheck_summary: 程序预检概览
            precheck_issues: 程序预检问题明细（fatal/warning）
        """
        precheck_text = precheck_summary if precheck_summary else "(none)"
        precheck_issue_text = precheck_issues if precheck_issues else "(none)"
        memento_section = memento_context if memento_context else "(no retrieved memento cases)"
        knowledge_section = knowledge_context if knowledge_context else "(no retrieved knowledge items)"
        skill_section = skill_context if skill_context else "(no retrieved skill items)"

        return f"""[TASK] Agent5 — Validate ONE Reaction Step's NEB Endpoints

You validate Agent4's co-adsorbate placement for ONE elementary step in a
sequential NEB pipeline. The reactant is the relaxed product from the previous
step, ensuring surface continuity. Focus on whether the added co-adsorbate
atoms are at physically reasonable positions.

═══════════════════════════════════════════════════════════
SECTION 1 — PASS / FAIL RULES
═══════════════════════════════════════════════════════════
• Programmatic pre-check "fatal" issues → MUST FAIL
• Programmatic pre-check "warning" issues → do NOT affect PASS/FAIL
• No fatal issues → MUST PASS
• Interpolation warnings remain advisory, but interpolation findings marked fatal
  by the programmatic pre-check MUST FAIL.

═══════════════════════════════════════════════════════════
SECTION 2 — INPUT DATA (this step only)
═══════════════════════════════════════════════════════════
[Step Summary]
{steps_details}

[Structure Coordinates]
{structures_payload}

[Programmatic Pre-check Summary]
{precheck_text}

[Pre-check Issues (fatal vs warning)]
{precheck_issue_text}

[Memento CaseBank]
{memento_section}

[KnowledgeBank]
{knowledge_section}

[SkillBank]
{skill_section}

═══════════════════════════════════════════════════════════
SECTION 3 — REVIEW CHECKLIST
═══════════════════════════════════════════════════════════
1. ELEMENT CONSISTENCY: reactant and product must have identical element multisets → mismatch = FAIL
2. CO-ADSORBATE QUALITY:
   • 2.5–5.0 Å from main adsorbate centroid (> 5.0 Å = FAIL)
   • At or near a surface hollow/bridge site, 0.5–2.5 Å above top surface layer
   • Not inside the slab (FAIL if below surface)
3. ATOM DISTANCES: all pair distances > 0.8 Å (overlap = FAIL)
4. ADSORBATE FORCES: high forces (> 10 eV/Å) suggest bad geometry → FAIL with specific coordinates

Do NOT override programmatic pre-check results. 0 fatal issues → MUST PASS.

═══════════════════════════════════════════════════════════
SECTION 4 — OUTPUT
═══════════════════════════════════════════════════════════
When FAIL: feedback must be ACTIONABLE — which atom, what's wrong, suggested fix.
When PASS with warnings: list warnings in issues, status stays PASS.

{{
  "status": "PASS" or "FAIL",
  "issues": ["..."],
  "feedback": "..."
}}"""

    # ----- Mechanism search prompts -----

    @staticmethod
    def mechanism_recommend_pathways(
        initial_state: str,
        target_state: str,
        surface_desc: str,
        constraints: str = "",
        n_pathways: int = 3,
    ) -> str:
        """Prompt LLM to recommend multiple competing reaction pathways."""
        constraint_line = f"\nConstraints: {constraints}" if constraints else ""
        return f"""You are a heterogeneous catalysis expert. Propose {n_pathways} distinct, chemically plausible reaction pathways for:

Reaction: {initial_state} → {target_state} on {surface_desc}
{constraint_line}

Requirements:
1. Each pathway is an ordered list of surface intermediate labels (e.g. *CO2, *COOH, *CO, *CHO, …).
2. Use standard adsorbate notation: *X for adsorbed species, X(g) for gas-phase species.
3. Each consecutive pair must differ by exactly one elementary step (single bond break/form, +H, -H, +OH, −OH, H₂O release, etc.).
4. The {n_pathways} pathways must represent genuinely different mechanisms (not minor variants of the same route).
5. Include both {initial_state} as the first entry and {target_state} as the last entry in every pathway.
6. Prioritise mechanisms that are well-documented in literature for this surface type.

Output ONLY a JSON array of arrays, no explanation:
[
  ["{initial_state}", "…", "…", "{target_state}"],
  ["{initial_state}", "…", "…", "{target_state}"]
]"""

    @staticmethod
    def mechanism_context_routing(
        reaction_description: str,
        surface_info: str,
        initial_state: str,
        target_state: str,
        known_intermediates: str,
        environment: str,
        care_domain_report: str = "",
    ) -> str:
        """Agent M1 task prompt."""
        return f"""[TASK] Agent M1 — Mechanism Context Routing & Candidate Expansion

[Reaction] {reaction_description}
[Surface] {surface_info}
[Initial] {initial_state}  [Target] {target_state}
[Known Intermediates] {known_intermediates}
[Environment] {environment}
[CARE Domain] {care_domain_report or "(not checked)"}

Workflow:
1. Call ``initialize_mechanism_context`` with the above.
2. Call ``generate_candidate_steps``.
3. Summarise the candidates and your chemical reasoning."""

    @staticmethod
    def mechanism_search_triage(
        candidate_summary: str,
        search_config: str,
        n_candidates: int,
    ) -> str:
        """Agent M2 task prompt."""
        return f"""[TASK] Agent M2 — Free-Energy Pruning & Pathway Export

[Candidates] {n_candidates} initial candidates
{candidate_summary}

[Config] {search_config}

Workflow:
1. Call ``run_mechanism_search``.
2. Review retained vs pruned pathways.
3. Call ``extract_pathway_shortlist``.
4. Output structured summary for Agent4/5."""

    # ----- Agent 1-3, 6-7 简短提示词 -----

    @staticmethod
    def agent1_initialize_surface(
        reaction_description: str,
        bulk_structure_path: str,
        run_id: str,
        top_n_surfaces: int,
        workflow_step: str,
        initial_surface_path: str = "",
    ) -> str:
        return f"""你是 Agent1（表面初始化）。

任务：
1) 若 `initial_surface_path` 非空，仅确认该路径可用于后续流程，不调用工具。
2) 若 `initial_surface_path` 为空，必须只调用一次工具 `generate_surfaces`，且参数必须与下述值完全一致，不得改写：
   - bulk_structure_path: {bulk_structure_path}
   - top_n_surfaces: {top_n_surfaces}
   - run_id: {run_id}
   - workflow_step: {workflow_step}
3) 禁止臆造任何路径或默认值。

上下文：
- reaction_description: {reaction_description}
- initial_surface_path: {initial_surface_path}

最后仅输出简短 JSON：{{"status":"OK"}}"""

    @staticmethod
    def agent2_place_adsorbate(
        surface_path: str,
        initial_adsorbate: str,
        num_sites: int,
        run_id: str,
        workflow_step: str,
        llm_review_context: Optional[str] = None,
    ) -> str:
        return f"""你是 Agent2（吸附位点预测）。

必须只调用一次工具 `predict_adsorption_sites`，且参数必须与下述值完全一致，不得改写：
- surface_path: {surface_path}
- adsorbate_smi: {initial_adsorbate}
- num_sites: {num_sites}
- run_id: {run_id}
- workflow_step: {workflow_step}
- llm_review_context: {llm_review_context}

禁止使用其他路径、其他吸附物或其他数值。
最后仅输出简短 JSON：{{"status":"OK"}}"""

    @staticmethod
    def agent3_reconstruct_surface(
        surface_with_adsorbate_path: str,
        adsorbates_elements_for_mc: List[str],
        temperature_k: float,
        total_sweeps: int,
        run_id: str,
        workflow_step: str,
        clean_slab_path: Optional[str] = None,
        surface_indices: Optional[List[int]] = None,
        adsorbate_indices: Optional[List[int]] = None,
    ) -> str:
        return f"""你是 Agent3（表面重构）。

必须只调用一次工具 `simulate_surface_reconstruction`，且参数必须与下述值完全一致，不得改写：
- surface_with_adsorbate_path: {surface_with_adsorbate_path}
- adsorbates_elements_for_mc: {adsorbates_elements_for_mc}
- temperature_k: {temperature_k}
- total_sweeps: {total_sweeps}
- run_id: {run_id}
- workflow_step: {workflow_step}
- clean_slab_path: {clean_slab_path}
- surface_indices: {surface_indices}
- adsorbate_indices: {adsorbate_indices}

禁止使用其他参数组合。
最后仅输出简短 JSON：{{"status":"OK"}}"""

    @staticmethod
    def agent6_run_neb(step_structures: List[Dict[str, Any]]) -> str:
        step_names = [s.get("name", "unknown") for s in step_structures]
        return f"对以下步骤执行 NEB 并返回能垒摘要：{step_names}"

    @staticmethod
    def agent7_generate_report(
        pathway_name: str,
        reaction_type: str,
        adsorption_energies: Dict[str, float],
        neb_results: Dict[str, Any],
    ) -> str:
        return (
            "整合工作流结果并生成报告。"
            f" pathway_name={pathway_name}, reaction_type={reaction_type}, "
            f"adsorption_keys={list(adsorption_energies.keys())}, neb_steps={list(neb_results.keys())}"
        )


# ----- Agent Role Definitions -----

AGENT1_STRUCTURE_INITIALIZER = AgentRole(
    name="Agent1 Structure Initializer",
    role="Catalyst Surface Initializer",
    goal="Initialize a simulation-ready catalyst surface",
    backstory="Specialized in bulk-to-surface conversion and structural normalization.",
)

AGENT2_ADSORPTION_PREDICTOR = AgentRole(
    name="Agent2 Adsorption Predictor",
    role="Adsorption Site Planner",
    goal="Produce stable initial adsorption configurations",
    backstory="Specialized in adsorption placement, ranking and metadata continuity.",
)

AGENT3_RECONSTRUCTION_SIMULATOR = AgentRole(
    name="Agent3 Reconstruction Simulator",
    role="Surface Reconstruction Sampler",
    goal="Sample reconstructed low-energy surfaces",
    backstory="Specialized in thermodynamic sampling and robust structure handoff.",
)

AGENT4_PATHWAY_DESIGNER = AgentRole(
    name="Agent4 Pathway Designer",
    role="Reaction Pathway Designer",
    goal="Generate universal, executable NEB-ready pathway steps",
    backstory="General-purpose pathway reasoner for arbitrary catalytic reactions.",
)

AGENT5_PATHWAY_VALIDATOR = AgentRole(
    name="Agent5 Pathway Validator",
    role="Pathway Validation Auditor",
    goal="Validate pathway quality and return actionable feedback",
    backstory="Focuses on mapping consistency and geometric computability without reaction-specific hardcoding.",
)

AGENT6_NEB_RUNNER = AgentRole(
    name="Agent6 NEB Runner",
    role="Barrier Computation Specialist",
    goal="Execute barrier calculations and summarize reliability",
    backstory="Specialized in NEB execution, convergence review and barrier extraction.",
)

AGENT7_REPORT_GENERATOR = AgentRole(
    name="Agent7 Report Generator",
    role="Workflow Reporter",
    goal="Aggregate all outputs into final report artifacts",
    backstory="Specialized in traceable multi-stage synthesis and reporting.",
)

# --- Mechanism search agents (between Agent3 and Agent4/5) ---

MECHANISM_CONTEXT_ROUTING_AGENT = AgentRole(
    name="MechanismContextRouter",
    role="Mechanism Search Coordinator",
    goal=(
        "Analyse the user's reaction specification, determine the exploration strategy, "
        "call mechanism search tools, and produce evaluated candidate pathways."
    ),
    backstory=(
        "You are the mechanism search coordinator in the CatDT workflow. "
        "You have access to these tools:\n"
        "  • initialize_mechanism_context — set up the search context\n"
        "  • generate_candidate_steps — enumerate candidate reactions\n"
        "  • run_mechanism_search — execute beam search with UMA energy pruning\n"
        "  • extract_pathway_shortlist — export retained pathways for Agent4/5\n\n"
        "Your workflow:\n"
        "1. Call initialize_mechanism_context with the reaction info.\n"
        "2. Decide exploration_mode based on user request:\n"
        "   - If user specified intermediates → 'agent_guided' (fast, evaluate the given path)\n"
        "   - If user asked for systematic/thorough exploration → 'systematic'\n"
        "   - If user wants both → 'both_sequential'\n"
        "   - Default when only initial+target given → 'agent_guided' (LLM recommends paths)\n"
        "3. Call run_mechanism_search with the chosen mode.\n"
        "4. Call extract_pathway_shortlist to prepare Agent4/5 handoff.\n"
        "5. Output a brief JSON summary of the results."
    ),
)

MECHANISM_SEARCH_TRIAGE_AGENT = AgentRole(
    name="MechanismSearchTriage",
    role="Free-Energy Pruning & Pathway Export Coordinator",
    goal=(
        "Evaluate candidate pathways via free-energy estimation, prune unpromising "
        "branches, select the shortlist of high-value pathways, and export them "
        "to structure files for Agent4/5 consumption."
    ),
    backstory=(
        "You are the mechanism search triage coordinator. You receive candidate "
        "states and steps from Agent M1 and decide which branches to keep. "
        "Your workflow: "
        "(1) call the search engine tool to run beam search with pruning, "
        "(2) review the retained pathways and explain why each was kept, "
        "(3) invoke the materialization & export tool to write structure files, "
        "(4) prepare the handoff payload for Agent4/5. "
        "You NEVER estimate energies directly — you call tools. "
        "You NEVER construct atomic structures — tools do that."
    ),
)


def get_agent_role(agent_name: str) -> Optional[AgentRole]:
    roles = {
        "agent1": AGENT1_STRUCTURE_INITIALIZER,
        "agent2": AGENT2_ADSORPTION_PREDICTOR,
        "agent3": AGENT3_RECONSTRUCTION_SIMULATOR,
        "agent4": AGENT4_PATHWAY_DESIGNER,
        "agent5": AGENT5_PATHWAY_VALIDATOR,
        "agent6": AGENT6_NEB_RUNNER,
        "agent7": AGENT7_REPORT_GENERATOR,
        "mechanism_context_router": MECHANISM_CONTEXT_ROUTING_AGENT,
        "mechanism_search_triage": MECHANISM_SEARCH_TRIAGE_AGENT,
    }
    return roles.get(agent_name.lower())


def get_all_agent_roles() -> Dict[str, AgentRole]:
    return {
        "agent1": AGENT1_STRUCTURE_INITIALIZER,
        "agent2": AGENT2_ADSORPTION_PREDICTOR,
        "agent3": AGENT3_RECONSTRUCTION_SIMULATOR,
        "agent4": AGENT4_PATHWAY_DESIGNER,
        "agent5": AGENT5_PATHWAY_VALIDATOR,
        "agent6": AGENT6_NEB_RUNNER,
        "agent7": AGENT7_REPORT_GENERATOR,
        "mechanism_context_router": MECHANISM_CONTEXT_ROUTING_AGENT,
        "mechanism_search_triage": MECHANISM_SEARCH_TRIAGE_AGENT,
    }
