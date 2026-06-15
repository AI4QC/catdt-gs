"""
CatDT Agent Prompts (CAMEL-ready)

设计原则：
1. 提示词完全通用，不包含针对某一反应的规则或示例。
2. Agent4/Agent5 通过多轮反馈协作（最多 N 轮）逐步修正。
3. 不预设固定元素增删规则，所有增删改由 LLM 基于输入自行判断。
4. 完整结构坐标仅在首次输入一次，后续只给吸附分子坐标。
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
        return f"""【任务】解析催化反应上下文

【输入反应描述】
{reaction_description}

【要求】
1. 提取 reaction_type: 反应类型概述（通用文本）。
2. 提取 initial_adsorbate: 初始主吸附态标签。
3. 提取 intermediates: 有序中间体序列（必须含起点和终点，且顺序可执行）。
4. 提取 constraints: 用户显式约束与计算约束摘要。
5. 若用户未给完整中间体序列，请基于一般催化机理补全"最小可执行主路径"。
6. 禁止引入与用户目标无关的支路。
7. 仅输出 JSON，不要额外解释。

【输出格式】
{{
  "reaction_type": "...",
  "initial_adsorbate": "...",
  "intermediates": ["...", "..."],
  "constraints": "..."
}}
"""

    @staticmethod
    def pathway_design(
        reaction_context: Dict[str, Any],
        surface_info: str,
        baseline_steps: str,
        staging_plan_text: str = "",
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
            staging_plan_text: 程序预计算的暂驻方案（元素差异 + 建议坐标）
            feedback: 上一轮 Agent5 反馈
            history: 迭代历史摘要
        """
        feedback_section = f"\n【上一轮反馈（含前一轮吸附分子坐标）】\n{feedback}\n" if feedback else ""
        history_section = f"\n【迭代历史】\n{history}\n" if history else ""
        staging_section = staging_plan_text if staging_plan_text else "(无需暂驻)"
        memento_section = memento_context if memento_context else "(no retrieved memento cases)"
        knowledge_section = knowledge_context if knowledge_context else "(no retrieved knowledge items)"
        skill_section = skill_context if skill_context else "(no retrieved skill items)"

        return f"""【任务】Agent4 生成 NEB 端点修正方案

你的任务是为每一步反应补齐缺失的原子，使 reactant 和 product 两端的元素组成完全一致（NEB 要求）。

【核心原理】
工具已为每步生成了 reactant 和 product 的吸附质结构。由于相邻步骤的吸附质组成不同（如 CO→CHO 多了一个 H），你需要在缺少原子的一端补上对应原子。
- 若 product 比 reactant 多 H：在 reactant 端的吸附质附近放一个 H
- 若 reactant 比 product 多 OH：在 product 端的吸附质附近放一个解离的 O 和 H
- 新增原子的位置必须参考已有吸附质坐标，放在合理的化学距离处

【原子放置指南——关键，必须严格遵守】
1. 新增原子必须放在该端已有吸附质原子附近（距最近吸附质原子 1.0-2.0 Å）
2. 确定位置的方法：
   a) 在 baseline 坐标中找到与新增原子化学上最相关的吸附质原子（锚点）
   b) 将新增原子放在锚点上方或旁边 1.0-2.0 Å 处
   c) 例如：若锚点原子在 [x, y, z]，新增原子可放在 [x, y, z+1.5]
3. 绝对禁止将新增原子放在距所有吸附质原子 >3.0 Å 的位置——这会导致 FAIL
4. NEB 要求同一原子在 reactant 与 product 之间的位移尽可能小（<2 Å 理想，<3 Å 可接受）
5. 多个新增原子之间距离应 >0.9 Å
6. 所有坐标使用笛卡尔坐标（Å），与 baseline 一致
7. 不要把新增原子放到表面层或远离吸附位点的空位上

【反应上下文】
reaction_type: {reaction_context.get("reaction_type", "")}
intermediates: {reaction_context.get("intermediates", [])}
constraints: {reaction_context.get("constraints", "")}

【催化表面信息】
{surface_info}

【每步基线结构与元素差异】
{baseline_steps}

【元素差异汇总】
{staging_section}

【三层记忆说明】
1) CaseBank：保存历史成功/失败案例，提供具体反例与可复用片段。
2) KnowledgeBank：从多案例抽象出的通用规则与约束，减少重复错误。
3) SkillBank：在规则基础上沉淀的可执行策略模板，提升跨反应可迁移性。
目的：你必须综合三层信息，在保持通用性的前提下生成当前步骤修正方案。

【历史案例记忆（Memento 检索）】
{memento_section}

【KnowledgeBank 检索结果（Top-K）】
{knowledge_section}

【SkillBank 检索结果（Top-K）】
{skill_section}
{feedback_section}{history_section}
【硬约束】
1. 工具阶段已有原子坐标绝对只读，不得改写或重排。
2. 仅输出新增原子：使用 atoms_to_add；atoms_to_modify 和 atoms_to_remove 必须为空。
3. 每个有元素差异的步骤，必须提供 atoms_to_add 及其精确的 [x, y, z] 坐标。
4. 坐标必须参考 baseline 中已有吸附质的坐标来确定（在其附近）。
5. 若某步无元素差异（element_delta=none），atoms_to_add 为空数组。
6. 不要输出推理过程；仅输出一个 JSON 对象。

【输出格式】
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
          "reason": "placed near C at [cx, cy, cz], distance ~1.1 Å"
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
        staging_plan_text: str = "",
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
            staging_plan_text: 预计算暂驻方案参考
        """
        precheck_text = precheck_summary if precheck_summary else "(none)"
        precheck_issue_text = precheck_issues if precheck_issues else "(none)"
        memento_section = memento_context if memento_context else "(no retrieved memento cases)"
        knowledge_section = knowledge_context if knowledge_context else "(no retrieved knowledge items)"
        skill_section = skill_context if skill_context else "(no retrieved skill items)"

        return f"""【任务】Agent5 审核 Agent4 端点

你的任务是判断 Agent4 生成的端点结构是否可用于 NEB 计算。

【PASS/FAIL 判定规则——严格遵守】
- 程序预检中标记为 "fatal" 的问题 → 必须 FAIL
- 程序预检中标记为 "warning" 的问题 → 不影响 PASS/FAIL 判定
- 如果没有任何 fatal 问题 → 必须 PASS
- 跨步骤核心连续性不匹配（pymatgen sort 导致）→ 这是已知行为，属于 warning，不应判 FAIL
- 插值路径如果是 warning 级别 → 不应判 FAIL；如果程序预检标为 fatal → 必须 FAIL

【步骤摘要】
{steps_details}

【每步结构文本】
{structures_payload}

【程序预检摘要】
{precheck_text}

【程序预检问题明细（区分 fatal 和 warning）】
{precheck_issue_text}

【三层记忆说明】
1) CaseBank：历史具体案例，用于对照相似失败/成功模式。
2) KnowledgeBank：跨案例抽象规则，用于判断是否违反普适约束。
3) SkillBank：可执行审查策略模板，用于提高审核一致性与可操作反馈质量。
目的：你需要同时参考三层记忆做审核，不得只依赖单一来源。

【历史案例记忆（Memento 检索）】
{memento_section}

【KnowledgeBank 检索结果（Top-K）】
{knowledge_section}

【SkillBank 检索结果（Top-K）】
{skill_section}

【重点审查项——新增原子位置】
- 新增（staged）原子必须在距最近吸附质原子 3.0 Å 以内
- 若程序预检报告 staged 原子距离过远（>3.0 Å）→ 这是 fatal 问题，必须 FAIL
- FAIL 时必须在 feedback 中给出该原子应移到的具体坐标（参考预检建议位置）

【反馈要求】
- 如果 FAIL：在 feedback 中给出具体修改建议（哪个步骤、哪个原子、移到哪里，附带具体坐标）
- 如果 PASS 但有 warning：在 issues 中列出 warning，但 status 仍为 PASS
- feedback 应该可操作：告诉 Agent4 具体的坐标调整方向和目标坐标

【输出格式】
{{
  "status": "PASS" 或 "FAIL",
  "issues": ["..."],
  "feedback": "..."
}}"""

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


def get_agent_role(agent_name: str) -> Optional[AgentRole]:
    roles = {
        "agent1": AGENT1_STRUCTURE_INITIALIZER,
        "agent2": AGENT2_ADSORPTION_PREDICTOR,
        "agent3": AGENT3_RECONSTRUCTION_SIMULATOR,
        "agent4": AGENT4_PATHWAY_DESIGNER,
        "agent5": AGENT5_PATHWAY_VALIDATOR,
        "agent6": AGENT6_NEB_RUNNER,
        "agent7": AGENT7_REPORT_GENERATOR,
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
    }
