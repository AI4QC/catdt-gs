# CatDT Agent 提示词系统重构说明

## 概述

本次重构将分散在各个文件中的 LLM 提示词集中到 `camel_agents/prompts.py`，实现了：
1. **提示词专业化** - 每个 Agent 的角色定义和任务提示词都经过优化
2. **集中管理** - 所有提示词统一存放，便于维护和版本控制
3. **代码解耦** - 业务逻辑与提示词分离，提高代码可读性

## 文件变更

### 1. 新增文件：`camel_agents/prompts.py`

集中存放所有 Agent 提示词，包含：

#### Agent 角色定义 (`AgentRole` 数据类)

| Agent | 名称 | 角色 | 目标 |
|-------|------|------|------|
| Agent1 | Structure Initializer | 催化剂表面结构初始化专家 | 从体相结构生成或加载稳定的催化剂表面 |
| Agent2 | Adsorption Site Predictor | 吸附位点预测与吸附质放置专家 | 在给定表面上预测最优吸附位点 |
| Agent3 | Surface Reconstruction Simulator | 表面重构与蒙特卡洛模拟专家 | 运行 VSSR-MC 表面重构模拟 |
| Agent4 | Reaction Pathway Designer | 催化反应路径设计专家 | 解析反应描述，设计 NEB 就绪的反应路径 |
| Agent5 | Pathway Validator | 反应路径验证与质量评估专家 | 验证反应路径的化学合理性和 NEB 结果 |
| Agent6 | NEB Runner | NEB 计算执行专家 | 运行 NEB 计算获取反应能垒和过渡态 |
| Agent7 | Report Generator | 催化计算报告生成与可视化专家 | 整合所有计算结果，生成专业报告 |

#### 任务提示词模板 (`TaskPrompts` 类)

- `reaction_context_parsing()` - 解析反应描述
- `pathway_design()` - 设计反应路径
- `pathway_validation()` - 验证反应路径
- `energy_profile_estimation()` - 估算能线图
- `neb_results_validation()` - 验证 NEB 结果
- `agent1_surface_initialization()` - Agent1 表面初始化任务
- `agent2_adsorption_prediction()` - Agent2 吸附预测任务
- `agent3_surface_reconstruction()` - Agent3 表面重构任务
- `agent6_neb_execution()` - Agent6 NEB 执行任务
- `agent7_report_generation()` - Agent7 报告生成任务

#### 默认反应描述

- `DEFAULT_HER_DESCRIPTION` - 氢气析出反应 (HER)
- `DEFAULT_NRR_DESCRIPTION` - 氮气还原反应 (NRR)
- `DEFAULT_CO2RR_DESCRIPTION` - 二氧化碳还原反应 (CO2RR)

#### 系统提示词

- `SYSTEM_PROMPTS` - 不同领域专家的系统提示词

### 2. 修改文件：`camel_agents/__init__.py`

导出 prompts 模块，便于其他代码使用：

```python
from .prompts import (
    AgentRole,
    TaskPrompts,
    get_all_agent_roles,
    get_agent_role,
    DEFAULT_HER_DESCRIPTION,
    DEFAULT_NRR_DESCRIPTION,
)
```

### 3. 修改文件：`camel_agents/catdt_gas_solid_workflow.py`

**变更内容：**
- 导入 prompts 模块中的角色定义和任务提示词
- 修改 `_init_agents()` 方法，使用 `AgentRole.to_dict()` 初始化 Agent
- 修改 `_create_reaction_context_task()` 方法，使用 `TaskPrompts.reaction_context_parsing()`
- 修改 `_create_pathway_design_task()` 方法，使用 `TaskPrompts.pathway_design()`
- 修改 `_create_validation_task()` 方法，使用 `TaskPrompts.pathway_validation()` 或 `TaskPrompts.neb_results_validation()`

**代码示例：**
```python
# 修改前
self.agent4 = Agent(
    role="Reaction Pathway Designer",
    goal="Parse reaction description once and design NEB-ready pathway structures",
    backstory="You must not re-identify surfaces/adsorbates; use provided indices.",
    llm=self.llm,
    verbose=True,
)

# 修改后
self.agent4 = Agent(
    **AGENT4_PATHWAY_DESIGNER.to_dict(),
    llm=self.llm,
    verbose=True,
)
```

### 4. 修改文件：`camel_agents/pathway_agents_only.py`

**变更内容：**
- 导入 prompts 模块
- 修改 `_init_agents()` 方法使用专业角色定义
- 修改 `parse_reaction_context()` 使用 `TaskPrompts.reaction_context_parsing()`
- 修改 `design_pathway()` 使用 `TaskPrompts.pathway_design()`
- 修改 `validate_pathway()` 使用 `TaskPrompts.pathway_validation()`
- 修改 `estimate_energy_profile()` 使用 `TaskPrompts.energy_profile_estimation()`
- 修改 `get_default_her_description()` 和 `get_default_nrr_description()` 使用 prompts 中的默认描述

### 5. 修改文件：`camel_agents/resumable_workflow.py`

**变更内容：**
- 导入 `TaskPrompts`（为后续扩展预留）

## 提示词改进亮点

### 1. Agent 4 - 反应路径设计专家

**改进前：**
```python
role="Reaction Pathway Designer"
goal="Parse reaction description once and design NEB-ready pathway structures"
backstory="You must not re-identify surfaces/adsorbates; use provided indices."
```

**改进后：**
```python
role="催化反应路径设计专家"
goal="解析反应描述，设计 NEB 就绪的反应路径，包括中间体和过渡态结构"
backstory="""你是一位催化机理和反应路径设计的专家，擅长：
- 解析自然语言描述的反应，提取反应类型、中间体序列和约束条件
- 设计符合化学直觉的反应路径（加氢、脱氢、C-O/C-N 键断裂/形成等）
- 为每一步反应指定原子级别的修改（添加/移除原子）
- 保持所有中间体在同一反应位点（same-site constraint）
- 使用提供的表面/吸附质索引，绝不重新识别原子
你设计的路径将直接用于 NEB 计算，必须保证每一步的化学合理性。"""
```

### 2. 任务提示词结构化

**改进前：**
```python
description = f"""Parse the reaction description once and output JSON with:
- reaction_type
- initial_adsorbate
...
"""
```

**改进后：**
```python
description = f"""【任务】解析催化反应描述

【输入描述】
{reaction_description}

【解析要求】
请仔细解析上述反应描述，提取以下信息并以 JSON 格式返回：

1. **reaction_type** (string): 
   - 反应类型，如 "CO 加氢制甲烷"、"氮气还原"...
   
2. **initial_adsorbate** (string):
   - 初始吸附质，格式为 "*化学式"...

3. **intermediates** (array of strings):
   - 完整的中间体序列...

4. **constraints** (string):
   - 用户指定的任何约束条件...

【输出格式】
只返回以下 JSON 格式，不要有其他文字：
{{...}}

【重要提示】
- 中间体序列必须保持化学逻辑的一致性
- 如为加氢反应，中间体应逐步增加 H 原子
...
"""
```

## 使用方式

### 在其他模块中使用提示词

```python
from camel_agents.prompts import (
    AGENT4_PATHWAY_DESIGNER,
    TaskPrompts,
    DEFAULT_HER_DESCRIPTION,
)

# 使用 Agent 角色定义
agent = Agent(
    **AGENT4_PATHWAY_DESIGNER.to_dict(),
    llm=llm,
    verbose=True,
)

# 使用任务提示词
description = TaskPrompts.pathway_design(
    reaction_context={...},
    adsorbate_info=...,
    feedback=...,
)

# 使用默认描述
her_desc = DEFAULT_HER_DESCRIPTION
```

### 修改提示词

只需编辑 `camel_agents/prompts.py` 文件，所有使用到该提示词的地方会自动更新。

## 后续建议

1. **版本控制** - 建议对 prompts.py 进行版本控制，记录每次修改
2. **A/B 测试** - 可以通过修改提示词对比不同版本的效果
3. **扩展性** - 可以方便地添加新的 Agent 角色和任务提示词
4. **多语言支持** - 可以添加英文版的提示词，通过参数切换

## 文件清单

- ✅ `camel_agents/prompts.py` - 新增
- ✅ `camel_agents/__init__.py` - 修改
- ✅ `camel_agents/catdt_gas_solid_workflow.py` - 修改
- ✅ `camel_agents/pathway_agents_only.py` - 修改
- ✅ `camel_agents/resumable_workflow.py` - 修改
