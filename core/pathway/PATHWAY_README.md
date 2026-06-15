# Pathway Module - 反应路径能量与能垒预测

## 概述

pathway 模块提供了完整的催化反应路径分析功能，包括：

1. **自由能预测** (`fairchem_predictor.py`)：使用 Fairchem 模型预测吸附能
2. **能垒预测** (`barrier_predictor.py`)：使用 NEB 方法计算反应能垒
3. **完整路径分析** (`pathway_predictor.py`)：整合自由能和能垒，生成完整的能量剖面图

这些模块基于 [Fairchem](https://github.com/facebookresearch/fairchem) 库，支持使用 UMA、eSCN 等先进的 ML 势能模型。

## 安装依赖

### 1. 克隆 Fairchem 库

```bash
cd /home/zhilong/workspace/catdt/deps
git clone https://github.com/facebookresearch/fairchem.git
```

### 2. 安装 Fairchem

```bash
cd fairchem
pip install -e packages/fairchem-core
```

### 3. 验证安装

```python
from fairchem.core import pretrained_mlip
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
print("✓ Fairchem installed successfully")
```

注意：首次使用 UMA 模型需要登录 Hugging Face：
```bash
huggingface-cli login
```

## 模块 1: FairchemPredictor - 自由能预测

### 功能

- 预测单个结构的能量
- 预测吸附能（E_ads = E(ads+surf) - E(surf) - E(ads_ref)）
- 批量预测反应路径上多个吸附质的能量

### 使用示例

```python
import sys
sys.path.insert(0, 'core')
from pathway.fairchem_predictor import FairchemPredictor

# 初始化预测器
predictor = FairchemPredictor(
    fairchem_root="deps/fairchem",
    model_name="uma-s-1p1",  # 或 "uma-m-1p1" (更准确但更慢)
    use_gpu=True,
)

# 1. 预测单个结构的能量
result = predictor.predict_energy(
    structure="surface.vasp",
    relax=True,
    fmax=0.05,  # 力收敛标准 (eV/Å)
)
print(f"Energy: {result.energy:.3f} eV")

# 2. 预测吸附能
ads_result = predictor.predict_adsorption_energy(
    surface="clean_surface.vasp",
    adsorbate="*CO",  # SMILES 格式或 Atoms 对象
    num_sites=10,     # 尝试的吸附位点数
)
print(f"Adsorption energy: {ads_result.adsorption_energy:.3f} eV")

# 3. 预测反应路径能量
pathway_result = predictor.predict_pathway_energies(
    surface="surface.vasp",
    adsorbates=["*O", "*OH", "*OOH"],
    num_sites=5,
)
print(pathway_result.summary())
```

### 可用模型

| 模型 | 参数量 | 速度 | 精度 | 推荐用途 |
|------|--------|------|------|----------|
| `uma-s-1p1` | 6.6M / 150M | 快 | 高 | 快速筛选、大规模计算 |
| `uma-m-1p1` | 50M / 1.4B | 慢 | 更高 | 精确计算、最终验证 |

### 原子参考能量

模块内置了常见元素的参考能量（用于计算吸附能）：

```python
ATOMIC_REFERENCE_ENERGIES = {
    "H": -3.477 eV,
    "C": -7.282 eV,
    "N": -8.083 eV,
    "O": -7.204 eV,
    "F": -4.891 eV,
    "S": -4.659 eV,
}
```

## 模块 2: BarrierPredictor - 能垒预测

### 功能

- 使用 NEB (Nudged Elastic Band) 方法计算反应能垒
- 支持从结构文件或反应式自动生成 NEB 路径
- 使用 CatTSunami 的 AutoFrame 自动生成初始/最终帧

### 使用示例

```python
import sys
sys.path.insert(0, 'core')
from pathway.fairchem_predictor import FairchemPredictor
from pathway.barrier_predictor import BarrierPredictor

# 初始化（共享 Fairchem 模型）
fairchem = FairchemPredictor(
    fairchem_root="deps/fairchem",
    model_name="uma-s-1p1",
)

barrier_predictor = BarrierPredictor(
    fairchem_predictor=fairchem,
)

# 1. 从反应物和产物结构预测能垒
result = barrier_predictor.predict_from_structures(
    reactant="reactant.vasp",
    product="product.vasp",
    n_frames=10,       # NEB 帧数
    fmax=0.05,         # 力收敛标准
)

print(f"Activation energy (forward): {result.activation_energy_forward:.3f} eV")
print(f"Activation energy (reverse): {result.activation_energy_reverse:.3f} eV")
print(f"Reaction energy: {result.reaction_energy:.3f} eV")

# 2. 从反应式预测能垒（需要 CatTSunami 数据库）
result = barrier_predictor.predict_from_reaction(
    reaction_str="*CH -> *C + *H",
    surface="Pt_111.vasp",
    n_frames=10,
)
```

### NEB 参数说明

- `n_frames`: NEB 路径的帧数（包括反应物和产物）。通常 10-15 帧足够
- `fmax`: 力收敛标准。0.05 eV/Å 是常用值
- `max_steps`: 最大优化步数。NEB 通常需要 200-500 步
- `k`: 弹簧常数 (eV/Å²)。默认 1.0
- `climb`: 是否使用爬升像方法（推荐）

### NEB 优化策略

BarrierPredictor 使用两阶段优化：

1. **阶段 1**：不使用爬升像，快速收敛到大致路径
2. **阶段 2**：启用爬升像，精确定位过渡态

## 模块 3: PathwayPredictor - 完整路径分析

### 功能

- 整合自由能和能垒预测
- 自动识别速率决定步骤（RDS）
- 生成完整的能量剖面图

### 使用示例

```python
import sys
sys.path.insert(0, 'core')
from pathway.pathway_predictor import PathwayPredictor

# 初始化预测器
predictor = PathwayPredictor(
    fairchem_root="deps/fairchem",
    model_name="uma-s-1p1",
    use_gpu=True,
)

# 预测完整反应路径
result = predictor.predict_pathway(
    surface="Pt_111.vasp",
    adsorbates=["*O", "*OH", "*OOH", "*"],  # 按反应顺序
    calculate_barriers=True,  # 是否计算能垒
    num_sites=5,              # 每个吸附质尝试的位点数
    n_frames=10,              # NEB 帧数
)

# 打印结果
print(result.summary())

# 绘制能量剖面图
predictor.plot_energy_profile(
    result,
    output="energy_profile.png",
    show_barriers=True,
)
```

### 输出示例

```
Complete Reaction Pathway Analysis
================================================================================
Surface: Pt36
Number of intermediates: 4
Number of steps: 3
Overall reaction energy: -0.523 eV

================================================================================
Adsorption Energies:
--------------------------------------------------------------------------------
  *O         : E_ads =  -1.234 eV, E_total =   -345.123 eV
  *OH        : E_ads =  -0.876 eV, E_total =   -344.965 eV
  *OOH       : E_ads =  -0.543 eV, E_total =   -344.732 eV
  *          : E_ads =   0.000 eV, E_total =   -343.600 eV

================================================================================
Reaction Steps:
--------------------------------------------------------------------------------

Step 1: *O_to_*OH [RDS]
  *O -> *OH
  Reaction energy (ΔE): +0.158 eV
  Activation energy (E_act): 0.876 eV
  Transition state energy: -344.247 eV

Step 2: *OH_to_*OOH
  *OH -> *OOH
  Reaction energy (ΔE): +0.233 eV
  Activation energy (E_act): 0.654 eV

Step 3: *OOH_to_*
  *OOH -> *
  Reaction energy (ΔE): +1.132 eV
  Activation energy (E_act): 0.432 eV

================================================================================
Rate-Determining Step:
  Step 1: *O_to_*OH
  Barrier: 0.876 eV

Maximum barrier: 0.876 eV
```

## 测试

运行测试脚本：

```bash
# 快速测试（不包括能垒计算）
python test/test_pathway_predictor.py

# 完整测试（包括能垒计算，需要较长时间）
RUN_SLOW_TESTS=1 python test/test_pathway_predictor.py
```

测试将创建：
- 测试用的 Pt(111) 表面
- 多个吸附质的能量预测
- 能量剖面图

## 性能优化

### 1. 使用 Turbo 模式（适用于 MD 和长优化）

```python
predictor = FairchemPredictor(
    fairchem_root="deps/fairchem",
    model_name="uma-s-1p1",
    inference_mode="turbo",  # 1.5-2x 加速
)
```

### 2. 多 GPU 并行

```python
# 需要安装 Ray: pip install ray
predictor = FairchemPredictor(
    fairchem_root="deps/fairchem",
    model_name="uma-s-1p1",
    inference_mode="turbo",
    workers=4,  # 使用 4 个 GPU
)
```

### 3. 减少计算量

- 使用 `uma-s-1p1` 而不是 `uma-m-1p1`
- 减少 `num_sites`（位点数）
- 使用较大的 `fmax`（如 0.1 而不是 0.05）
- 减少 NEB `n_frames`（如 5-7 帧用于快速估计）

## 常见问题

### Q1: 如何处理自定义吸附质？

可以直接传入 ASE Atoms 对象：

```python
from ase import Atoms

custom_adsorbate = Atoms('COOH',
                         positions=[...],
                         cell=[...])

result = predictor.predict_adsorption_energy(
    surface="surface.vasp",
    adsorbate=custom_adsorbate,
)
```

### Q2: 如何处理异常检测（解吸、解离等）？

FairchemPredictor 自动检测异常：

```python
ads_result = predictor.predict_adsorption_energy(...)

if ads_result.is_anomalous:
    print("Warning: Anomaly detected!")
    # 重新计算或使用其他位点
```

### Q3: 如何保存和恢复计算？

所有计算结果都保存在 `output_dir` 中：

```python
result = predictor.predict_pathway(
    surface="surface.vasp",
    adsorbates=["*O", "*OH"],
    output_dir="my_calculation",
)

# 结构保存在 my_calculation/*.vasp
# 轨迹保存在 my_calculation/*.traj
```

### Q4: NEB 计算失败怎么办？

尝试：
1. 增加 `max_steps`
2. 增大 `fmax`（降低收敛标准）
3. 减少 `n_frames`
4. 检查反应物和产物是否已充分优化

## 参考文献

1. **Fairchem**: Chanussot et al. (2021). Open Catalyst 2020 (OC20) Dataset. *ACS Catalysis*, 11(10), 6059-6072.

2. **UMA Models**: Barroso et al. (2025). Universal Models for Atoms. *Nature Communications*.

3. **CatTSunami**: Esterhuizen et al. (2022). Comprehensive Catalyst Screening via NEB Calculations. *Nature Catalysis*.

4. **NEB Method**: Henkelman & Jónsson (2000). Improved Tangent Estimate in the NEB Method. *Journal of Chemical Physics*, 113(22), 9978-9985.

## 维护信息

- **创建日期**: 2026-01-19
- **作者**: Claude
- **依赖**: Fairchem >= 2.0, ASE >= 3.22
- **测试状态**: 已测试
