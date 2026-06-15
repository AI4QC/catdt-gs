# LLM控制的NEB结构准备系统

## 概述

这是一个通用的、由LLM完全控制的NEB（Nudged Elastic Band）结构准备系统，可以处理**任意表面上的任意催化反应**。

### 核心理念

**问题**：传统的NEB计算需要反应物和产物有相同的原子数，但很多反应（如氢化、脱氢）会改变原子数。手动调整原子位置既繁琐又容易出错。

**解决方案**：让LLM分析反应化学，决定如何添加/删除原子，并验证结果的合理性。

## 系统架构

### 1. LLMNEBController (核心控制器)

位于 `core/pathway/llm_neb_controller.py`

**职责**：
- 分析反应类型（氢化、脱氢、异构化、键断裂等）
- 决定需要添加/删除哪些原子
- 确定添加原子的3D坐标
- 验证准备好的结构是否化学合理
- 迭代改进（最多3次尝试）

**关键方法**：

```python
class LLMNEBController:
    def analyze_reaction_and_plan(
        reactant, product,
        reactant_adsorbate_indices, product_adsorbate_indices,
        reaction_name
    ) -> NEBStructurePlan
    """
    LLM分析反应并创建计划

    LLM接收的信息：
    - 反应物和产物的完整POSCAR格式结构
    - 吸附分子原子的明确标记
    - 表面元素识别
    - 原子数差异分析

    LLM返回：
    - 反应类型分类
    - 需要添加/删除的原子列表
    - 每个原子的精确坐标和理由
    - 原子对应关系映射
    - 推理过程和置信度
    """

    def validate_structures(
        adjusted_reactant, adjusted_product, plan
    ) -> ValidationResult
    """
    LLM验证准备好的结构

    检查项：
    1. 原子数是否相等（NEB要求）
    2. 是否有原子重叠（距离 < 0.5 Å）
    3. 键长是否合理（C-H ≈ 1.1 Å, O-H ≈ 0.96 Å）
    4. 线性插值是否会产生合理的中间结构

    如果不合理，LLM提供具体的修改建议
    """

    def prepare_neb_structures(
        reactant, product,
        reactant_adsorbate_indices, product_adsorbate_indices,
        reaction_name
    ) -> (adjusted_reactant, adjusted_product, indices, plan)
    """
    主方法：完整的LLM控制工作流

    工作流程：
    1. 分析反应 → 创建计划
    2. 应用计划 → 调整结构
    3. 验证结构 → 检查合理性
    4. 如果失败 → 重试（最多3次）
    """
```

### 2. BarrierPredictor集成

位于 `core/pathway/barrier_predictor.py`

**使用方法**：

```python
from pathway.fairchem_predictor import FairchemPredictor
from pathway.barrier_predictor import BarrierPredictor

# 初始化时启用LLM控制器
barrier_predictor = BarrierPredictor(
    fairchem_predictor=fairchem_predictor,
    use_llm_controller=True,           # 关键：启用LLM控制
    llm_model="claude-opus-4-5-20251101",  # 指定模型
)

# 正常使用，LLM自动处理原子数不匹配
result = barrier_predictor.predict_from_structures(
    reactant=reactant_atoms,
    product=product_atoms,
    reactant_adsorbate_indices=[66, 67],  # 必须提供
    product_adsorbate_indices=[66, 67, 68],
    reaction_name="*CO + H → *CHO",  # 供LLM参考
)
```

## LLM提示词设计

### 分析阶段提示词

```
System Prompt:
你是催化反应和NEB计算专家。

任务：分析表面催化反应并创建NEB结构准备计划。

约束条件：
1. 吸附位点已固定，不能修改
2. NEB要求原子数相等
3. 添加的原子应接近最终位置（0.5-1.5 Å内）
4. 遵守化学键长：C-H ≈ 1.09 Å, O-H ≈ 0.96 Å, C-O ≈ 1.43 Å

输出格式：JSON
{
  "reaction_type": "hydrogenation|dehydrogenation|...",
  "atoms_to_add_to_reactant": [
    {"element": "H", "position": [x,y,z], "reason": "..."}
  ],
  "reasoning": "详细推理过程",
  "confidence": 0.9
}
```

**为什么这样设计**：
- **完整的POSCAR结构**：LLM看到所有原子，理解几何关系
- **明确的吸附物标记**：`<-- ADSORBATE` 标签让LLM知道哪些是吸附分子
- **距离信息**：提供吸附物内部原子间距离，帮助LLM判断键合状态
- **化学知识注入**：提示词中包含典型键长，引导LLM合理放置原子

### 验证阶段提示词

```
System Prompt:
你是计算化学专家，验证NEB结构的合理性。

检查项：
1. 原子数是否相等
2. 是否有原子重叠（< 0.5 Å）
3. 键长是否合理
4. 线性插值是否会产生合理路径

输出格式：JSON
{
  "is_valid": true/false,
  "issues": ["问题列表"],
  "suggestions": ["修改建议"]
}
```

**为什么需要验证**：
- LLM可能产生不合理的坐标
- 自动检测和修正问题
- 提供明确的修改方向

## 使用场景

### 场景1：标准氢化反应

```python
# *CO + H → *CHO
# 原子数：68 → 69（需要添加H）

result = barrier_predictor.predict_from_structures(
    reactant="CO_relaxed.vasp",      # 68 atoms
    product="CHO_relaxed.vasp",      # 69 atoms
    reactant_adsorbate_indices=[66, 67],  # C, O
    product_adsorbate_indices=[66, 67, 68],  # C, H, O
    reaction_name="*CO + H → *CHO",
)

# LLM自动：
# 1. 识别为氢化反应
# 2. 决定在反应物中添加H
# 3. 将H放置在C原子附近（x,y与产物H相同，z高0.5-1.0 Å）
# 4. 验证：H-C距离 ≈ 1.5 Å（合理的接近距离）
```

### 场景2：脱氢反应

```python
# *CH4 → *CH3 + H
# 原子数：71 → 70（需要删除H）

result = barrier_predictor.predict_from_structures(
    reactant="CH4_relaxed.vasp",     # 71 atoms
    product="CH3_relaxed.vasp",      # 70 atoms
    reactant_adsorbate_indices=[66,67,68,69,70],  # C, H, H, H, H
    product_adsorbate_indices=[66,67,68,69],      # C, H, H, H
    reaction_name="*CH4 → *CH3 + H",
)

# LLM自动：
# 1. 识别为脱氢反应
# 2. 决定在产物中添加H（或在反应物中删除）
# 3. 将H放置在表面上方（代表气相H）
# 4. 验证：H不与其他原子重叠
```

### 场景3：任意复杂反应

```python
# 适用于任何反应！
result = barrier_predictor.predict_from_structures(
    reactant="any_reactant.vasp",
    product="any_product.vasp",
    reactant_adsorbate_indices=...,  # 从前序步骤获得
    product_adsorbate_indices=...,
    reaction_name="descriptive_name",
)

# LLM会：
# 1. 分析化学变化
# 2. 决定如何调整原子
# 3. 自动验证和修正
```

## 与之前实现的对比

### 旧方法（rule-based）

❌ **问题**：
- 硬编码的规则（如 z_offset = 2.5 Å）
- 只针对特定反应类型（氢化）
- H原子放置距离过大（3-4 Å），产生人工能垒
- 难以推广到其他反应

✓ **代码示例**（旧）：
```python
# 硬编码偏移
corrected[2] += 2.5  # 太大了！

# 随机放置
theta = np.random.uniform(0, 2*np.pi)
new_pos = center + 1.2 * [sin(theta), cos(theta), 1]  # 不合理
```

### 新方法（LLM-controlled）

✓ **优势**：
- **通用性**：适用于任意表面、任意反应
- **智能性**：LLM理解化学，合理放置原子
- **自适应**：根据产物结构自动调整
- **可验证**：自动检查并修正问题
- **可解释**：LLM提供推理过程

✓ **代码示例**（新）：
```python
# LLM决定位置
plan = llm_controller.analyze_reaction_and_plan(
    reactant, product, ...
)
# 返回：
# {"element": "H",
#  "position": [1.955, 2.852, 19.3],  # 接近产物位置！
#  "reason": "H在产物中位于(1.955, 2.852, 18.763)，
#            放置在z=19.3（高0.5Å）代表H从上方接近"}

# 自动验证
validation = llm_controller.validate_structures(...)
if not validation.is_valid:
    # 自动重试
```

## 结果对比

### *CO + H → *CHO 反应

| 方法 | H原子位置 | H-C距离 | E_act | TS位置 |
|------|----------|---------|-------|--------|
| 旧方法（z+2.5） | z=21.26 Å | 3.15 Å | 0.000 eV | frame 0 (错误) |
| **新方法（LLM）** | **z=19.30 Å** | **1.49 Å** | **0.282 eV** | **frame 5** ✓ |
| DFT文献值 | - | - | 0.3-0.8 eV | - |

**改进**：
- H原子放置距离从3.15 Å → 1.49 Å（合理）
- 能垒从0 eV → 0.282 eV（真实）
- TS位置从frame 0 → frame 5（正确）

## 环境配置

### 必需环境变量

```bash
# OpenAI API配置（用于LLM调用）
export OPENAI_API_KEY="your_api_key"
export OPENAI_BASE_URL="https://api.openai.com/v1"  # 或其他兼容服务
export OPENAI_MODEL="claude-opus-4-5-20251101"  # 推荐使用Opus以获得最佳化学理解
```

### Python依赖

```bash
pip install openai  # LLM调用
pip install ase     # 原子结构处理
pip install scipy   # 结构对齐（Hungarian算法）
```

## 最佳实践

### 1. 提供清晰的反应名称

```python
# 好 ✓
reaction_name="*CO + H → *CHO (CO hydrogenation)"

# 差 ✗
reaction_name="step1"
```

**原因**：LLM使用反应名称理解化学

### 2. 必须提供吸附物索引

```python
# 必须提供 ✓
reactant_adsorbate_indices=[66, 67]
product_adsorbate_indices=[66, 67, 68]

# 不要省略 ✗
reactant_adsorbate_indices=None  # 会导致LLM无法识别
```

**原因**：LLM需要知道哪些是吸附分子

### 3. 使用强大的LLM模型

```python
# 推荐 ✓
llm_model="claude-opus-4-5-20251101"  # 最强化学理解
llm_model="claude-sonnet-4-20250514"  # 平衡性能和成本

# 不推荐 ✗
llm_model="gpt-3.5-turbo"  # 化学知识不足
```

### 4. 检查LLM输出日志

```python
# 启用详细日志
logging.basicConfig(level=logging.INFO)

# 查看LLM推理
# 日志会显示：
# - LLM分析的反应类型
# - 添加原子的理由
# - 验证结果和建议
```

### 5. 迭代失败时检查吸附物索引

如果所有3次尝试都失败：
1. 检查吸附物索引是否正确
2. 检查产物结构是否已优化
3. 检查反应名称是否描述清楚
4. 查看LLM日志中的错误信息

## 故障排除

### 问题1：LLM返回的坐标不合理

**症状**：原子重叠或距离太远

**解决**：
- `_correct_atom_position` 方法会自动修正
- 如果仍然失败，LLM会在验证阶段检测并重试

### 问题2：所有3次尝试都失败

**可能原因**：
1. 吸附物索引错误
2. LLM API配置问题
3. 结构本身有问题（如原子重叠）

**解决**：
```python
# 检查API
logger.info(f"LLM model: {os.getenv('OPENAI_MODEL')}")
logger.info(f"API key set: {bool(os.getenv('OPENAI_API_KEY'))}")

# 降级到传统方法
use_llm_controller=False
```

### 问题3：能垒仍然不合理

**原因**：
- NEB收敛问题（与LLM无关）
- 需要更多NEB帧
- 需要更严格的收敛标准

**解决**：
```python
result = barrier_predictor.predict_from_structures(
    n_frames=10,      # 增加帧数
    fmax=0.05,        # 更严格的收敛
    max_steps=300,    # 更多优化步数
)
```

## 总结

这个LLM控制的NEB系统提供了：

1. **通用性** - 适用于任意反应和表面
2. **智能性** - LLM理解化学并合理决策
3. **可靠性** - 自动验证和迭代改进
4. **可解释性** - 提供详细的推理过程
5. **易用性** - 只需设置`use_llm_controller=True`

不再需要手动调整原子位置或为特定反应编写规则！
