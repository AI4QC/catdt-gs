# Pathway 模块安装和使用指南

## 安装步骤

### 1. 安装 Fairchem 库

```bash
# 进入 deps 目录
cd /home/zhilong/workspace/catdt/deps

# 如果还没有克隆 fairchem（已经克隆过可跳过）
git clone https://github.com/facebookresearch/fairchem.git

# 安装 fairchem
cd fairchem
pip install -e packages/fairchem-core
```

### 2. 配置 Hugging Face 访问（使用 UMA 模型需要）

UMA 模型需要 Hugging Face 账户和访问权限：

```bash
# 登录 Hugging Face
huggingface-cli login

# 访问 https://huggingface.co/facebook/UMA 请求访问权限
```

如果无法获取 UMA 访问权限，可以使用其他模型（见下文）。

### 3. 验证安装

```bash
cd /home/zhilong/workspace/catdt
python -c "
import sys
sys.path.insert(0, 'core')
from pathway.fairchem_predictor import FairchemPredictor
print('✓ Pathway module imported successfully')
"
```

## 替代模型选项

如果无法使用 UMA 模型，可以使用以下替代方案：

### 选项 1: 使用 CP-MACE（推荐）

CP-MACE 已经在项目的其他模块中使用过，不需要 Hugging Face 认证：

```python
# 使用 core/reconstruction/cp_mace_predictor.py 中的 CP-MACE
# 需要先训练或加载预训练模型
```

### 选项 2: 使用 CHGNet

CHGNet 是一个开源的通用势能模型：

```python
from chgnet.model import CHGNet

# 修改 fairchem_predictor.py 使用 CHGNet
# 或直接使用 CHGNet 作为 ASE calculator
```

### 选项 3: 使用 Fairchem 的其他开源模型

Fairchem 还提供了其他不需要特殊权限的模型（需要查看 Fairchem v1 文档）。

## 快速开始示例

### 示例 1: 使用 CP-MACE 进行能量预测

这个示例展示如何结合现有的 CP-MACE predictor 使用：

```python
import sys
sys.path.insert(0, 'core')

from reconstruction.cp_mace_predictor import CPMACEPredictor
from pathway.barrier_predictor import BarrierPredictor
from ase.io import read

# 1. 使用 CP-MACE 预测器
cp_mace = CPMACEPredictor(
    cp_mace_root="deps/CP-MACE",
    use_gpu=True,
)

# 2. 加载预训练模型
model_paths = ["path/to/model1.model", "path/to/model2.model"]

# 3. 模拟和优化结构
sim_result = cp_mace.simulate(
    structure="surface_with_adsorbate.xyz",
    model_paths=model_paths,
    target_potential=-3.36,
    temperature=300.0,
    steps=1000,
)

print(f"Final energy: {sim_result.energy_history[-1]:.3f} eV")
```

### 示例 2: 不依赖 Fairchem 的简化版本

如果只需要基本的能量和能垒预测，可以使用任何 ASE calculator：

```python
from ase.calculators.emt import EMT  # 示例：使用 EMT
from ase.optimize import BFGS
from ase.mep import NEB
from ase.io import read

# 使用任何 ASE calculator
calc = EMT()

# 优化结构
atoms = read("structure.vasp")
atoms.calc = calc
opt = BFGS(atoms)
opt.run(fmax=0.05)

# NEB 计算
reactant = read("reactant.vasp")
product = read("product.vasp")
images = [reactant]
for i in range(8):
    images.append(reactant.copy())
images.append(product)

neb = NEB(images)
neb.interpolate()

for image in images:
    image.calc = calc

opt = BFGS(neb)
opt.run(fmax=0.05)
```

## 测试说明

由于 UMA 模型需要认证，测试脚本可能无法直接运行。您可以：

### 1. 在获取 UMA 访问权限后测试

```bash
# 完成 Hugging Face 登录后
python test/test_pathway_predictor.py
```

### 2. 使用其他模型测试代码结构

修改 `test_pathway_predictor.py` 中的模型名称：

```python
# 在测试脚本中将
model_name="uma-s-1p1"

# 改为其他可用的模型
model_name="your_model_name"
```

### 3. 单元测试（不需要下载模型）

创建 mock 测试来验证代码结构：

```python
# test_pathway_structure.py
import sys
sys.path.insert(0, 'core')

from pathway.fairchem_predictor import FairchemPredictor
from pathway.barrier_predictor import BarrierPredictor
from pathway.pathway_predictor import PathwayPredictor

# 测试类是否可以正确初始化（不加载模型）
print("Testing module imports...")
print("✓ FairchemPredictor imported")
print("✓ BarrierPredictor imported")
print("✓ PathwayPredictor imported")
print("\nAll modules can be imported successfully!")
```

## 实际使用流程

### 完整工作流程示例

```python
import sys
sys.path.insert(0, 'core')
from pathway.pathway_predictor import PathwayPredictor

# 假设您已经：
# 1. 完成了 Hugging Face 登录
# 2. 获得了 UMA 模型访问权限
# 3. 准备好了催化剂表面结构

predictor = PathwayPredictor(
    fairchem_root="deps/fairchem",
    model_name="uma-s-1p1",
    use_gpu=True,
    work_dir="output/my_pathway_analysis",
)

# 分析 OER 反应路径
result = predictor.predict_pathway(
    surface="my_catalyst_surface.vasp",
    adsorbates=["*", "*O", "*OH", "*OOH"],  # OER 中间体
    calculate_barriers=True,
    num_sites=10,
    n_frames=10,
)

# 查看结果
print(result.summary())

# 生成能量剖面图
predictor.plot_energy_profile(
    result,
    output="oer_energy_profile.png"
)
```

## 常见问题解决

### Q1: 无法下载 UMA 模型

**解决方案**:
1. 访问 https://huggingface.co/facebook/UMA
2. 点击 "Request Access"
3. 等待审批（通常几小时到几天）
4. 使用 `huggingface-cli login` 登录
5. 重新运行代码

### Q2: 计算太慢

**解决方案**:
- 使用较小的模型 (`uma-s-1p1` 而不是 `uma-m-1p1`)
- 启用 turbo 模式：`inference_mode="turbo"`
- 减少位点数：`num_sites=3`
- 使用较宽松的收敛标准：`fmax=0.1`
- 减少 NEB 帧数：`n_frames=5`

### Q3: GPU 内存不足

**解决方案**:
- 使用 CPU：`use_gpu=False`
- 减小批次大小
- 使用较小的模型
- 逐个计算而不是批量计算

### Q4: 想使用其他模型

**解决方案**:
修改 `fairchem_predictor.py` 的 `_load_model` 方法以支持其他模型或 calculator。

## 下一步

1. **获取 UMA 访问权限** - 这是使用最先进模型的推荐方式
2. **准备您的催化剂系统** - 使用前面的模块（SurFF, AdsorbDiff, VSSR-MC）
3. **定义反应路径** - 确定要研究的反应中间体
4. **运行分析** - 使用 pathway 模块预测能量和能垒
5. **解释结果** - 识别速率决定步骤，优化催化剂

## 技术支持

如果遇到问题：
1. 检查 README 和文档
2. 查看 Fairchem 官方文档
3. 检查错误日志
4. 尝试简化的测试用例

## 参考资源

- Fairchem 文档: https://github.com/facebookresearch/fairchem
- UMA 论文: https://arxiv.org/abs/...
- CatTSunami: https://github.com/facebookresearch/fairchem/tree/main/src/fairchem/applications/cattsunami
