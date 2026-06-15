# Energy Diagram Plotter - User Guide

## 简介

`energy_diagram_plotter.py` 是一个专业的自由能图可视化工具，结合了 pMuTT 和 catplot 的优点，提供出版级质量的能量图。

### 核心特性

✅ **灵活的输入格式**
- 简单列表：`[0.0, 1.2, -0.5]`
- 字典格式：`{'label': 'TS1', 'energy': 1.2, 'is_ts': True}`
- EnergyStep 对象

✅ **智能能垒拟合**
- 使用 catplot 的 spline 插值（比 pMuTT 更平滑）
- 自动检测过渡态
- 可调节峰宽和平台长度

✅ **专业的视觉效果**
- 4种预设配色方案（Material、Nature、Science、Elegant）
- 阴影效果增强立体感
- LaTeX 格式化学式（自动上下标）
- 高分辨率导出（默认 300 DPI）

✅ **出版级质量**
- 遵循 Nature/Science 期刊图表规范
- 可自定义所有样式参数
- 支持 PNG/PDF/SVG 等格式

## 快速开始

### 示例 1：最简单的用法

```python
from energy_diagram_plotter import quick_plot

# 一行代码绘制能量图
quick_plot(
    energies=[0.0, 1.2, -0.5, 0.8, -2.5],
    labels=['H2 + O2', 'TS1', 'H2O2', 'TS2', '2 H2O'],
    output='my_diagram.png'
)
```

### 示例 2：自定义样式

```python
from energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

# 创建自定义样式
style = DiagramStyle(
    color_scheme='nature',          # Nature 期刊配色
    figsize=(12, 7),                # 图片尺寸
    line_width=3.0,                 # 线宽
    use_shadow=True,                # 启用阴影
    show_energies=True,             # 显示能量值
    ylabel_units='kJ/mol'           # 能量单位
)

# 创建绘图器
plotter = EnergyDiagramPlotter(style=style)

# 绘制
plotter.plot(
    energies=[0.0, 1.2, -0.5, 0.8, -2.5],
    labels=['H₂ + O₂', 'TS1', 'H₂O₂', 'TS2', '2 H₂O']
)

# 保存
plotter.save('custom_diagram.png', dpi=300)
plotter.export_data('data.csv')  # 同时导出数据
```

### 示例 3：使用字典输入（最灵活）

```python
from energy_diagram_plotter import EnergyDiagramPlotter

# 定义能量步骤（字典格式）
steps = [
    {'label': 'CO* + O*', 'energy': 0.0, 'is_ts': False},
    {'label': 'TS-COOH', 'energy': 0.95, 'is_ts': True},
    {'label': 'CO-O*', 'energy': 0.15, 'is_ts': False, 'color': '#FF6B6B'},
    {'label': 'TS-CO2', 'energy': 0.85, 'is_ts': True},
    {'label': 'CO₂ + 2*', 'energy': -2.35, 'is_ts': False},
]

plotter = EnergyDiagramPlotter()
plotter.plot(steps)
plotter.save('co_oxidation.png')
```

## 输入格式详解

### 格式 1：能量列表（最简单）

```python
energies = [0.0, 1.2, -0.5, 0.8, -2.5]
labels = ['A', 'TS1', 'B', 'TS2', 'C']

# 脚本会自动检测 TS（相邻能量的峰值）
plotter.plot(energies, labels)

# 或者手动指定 TS
is_ts = [False, True, False, True, False]
plotter.plot(energies, labels, is_ts)
```

### 格式 2：字典列表（推荐）

```python
steps = [
    {
        'label': 'Reactants',       # 标签
        'energy': 0.0,              # 能量值
        'is_ts': False,             # 是否是过渡态
        'color': '#1976D2',         # 自定义颜色（可选）
        'line_style': '-'           # 线型（可选）
    },
    {
        'label': 'TS1',
        'energy': 1.2,
        'is_ts': True
    },
    # ...
]

plotter.plot(steps)
```

### 格式 3：EnergyStep 对象（最灵活）

```python
from energy_diagram_plotter import EnergyStep

steps = [
    EnergyStep(label='H₂ + O₂', energy=0.0, is_ts=False),
    EnergyStep(label='TS1', energy=1.2, is_ts=True, color='#FF0000'),
    EnergyStep(label='H₂O', energy=-2.5, is_ts=False),
]

plotter.plot(steps)
```

## 化学式格式化

脚本会自动将化学式转换为 LaTeX 格式：

| 输入 | 输出 | 效果 |
|------|------|------|
| `H2O` | `H$_2$O` | H₂O |
| `CO2` | `CO$_2$` | CO₂ |
| `H+` | `H$^+$` | H⁺ |
| `O2-` | `O$_2^-$` | O₂⁻ |
| `1/2 O2` | `$\frac{1}{2}$ O$_2$` | ½ O₂ |

**支持的特殊字符**：
- 下标：`H2O` → H₂O
- 上标（电荷）：`Fe3+` → Fe³⁺
- 分数：`1/2 O2` → ½ O₂
- 点号（吸附位）：`CO*` → CO*

## 样式配置

### 预设配色方案

#### 1. Material Design（默认）
```python
style = DiagramStyle(color_scheme='material')
```
- 主色：蓝色 `#1976D2`
- TS色：红色 `#D32F2F`
- 适合：演示、报告

#### 2. Nature 期刊
```python
style = DiagramStyle(color_scheme='nature')
```
- 主色：深蓝 `#0C4B8E`
- TS色：暗红 `#C3423F`
- 适合：Nature 系列期刊投稿

#### 3. Science 期刊
```python
style = DiagramStyle(color_scheme='science')
```
- 主色：深蓝 `#003C71`
- TS色：红色 `#B31B1B`
- 适合：Science 系列期刊投稿

#### 4. Elegant（优雅）
```python
style = DiagramStyle(color_scheme='elegant')
```
- 主色：深灰蓝 `#34495E`
- TS色：柔和红 `#E74C3C`
- 适合：论文、学术海报

### 自定义样式参数

```python
style = DiagramStyle(
    # === 颜色设置 ===
    color_scheme='material',        # 基础配色
    primary_color='#1976D2',        # 主要颜色（非TS步骤）
    ts_color='#D32F2F',             # 过渡态颜色
    background_color='#FFFFFF',     # 背景颜色

    # === 线条属性 ===
    line_width=2.5,                 # 主线宽度
    ts_line_width=2.0,              # TS线宽度
    platform_length=1.0,            # 平台长度
    ts_peak_width=0.8,              # TS峰宽度

    # === 阴影效果 ===
    use_shadow=True,                # 是否使用阴影
    shadow_color='#CCCCCC',         # 阴影颜色
    shadow_offset=0.05,             # 阴影偏移（y轴范围的百分比）

    # === 插值设置 ===
    interp_method='spline',         # 'spline' 或 'quadratic'
    interp_points=100,              # 每段的插值点数

    # === 标签设置 ===
    label_fontsize=12,              # 状态标签字体大小
    energy_fontsize=10,             # 能量值字体大小
    show_energies=True,             # 是否显示能量值
    energy_format='.2f',            # 能量值格式（2位小数）

    # === 坐标轴设置 ===
    xlabel='Reaction Coordinate',   # x轴标签
    ylabel='Free Energy',           # y轴标签
    ylabel_units='kcal/mol',        # y轴单位
    grid=False,                     # 是否显示网格

    # === 图片设置 ===
    figsize=(10, 6),                # 图片尺寸（英寸）
    dpi=300                         # 分辨率
)
```

## 高级功能

### 1. 多条反应路径对比

```python
import matplotlib.pyplot as plt
from energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# 路径 1
plotter1 = EnergyDiagramPlotter(style=DiagramStyle(color_scheme='material'))
plotter1.ax = axes[0]
plotter1.plot(energies_path1, labels_path1)
axes[0].set_title('Path 1: Direct Mechanism')

# 路径 2
plotter2 = EnergyDiagramPlotter(style=DiagramStyle(color_scheme='nature'))
plotter2.ax = axes[1]
plotter2.plot(energies_path2, labels_path2)
axes[1].set_title('Path 2: Sequential Mechanism')

plt.tight_layout()
plt.savefig('comparison.png', dpi=300)
```

### 2. 设置参考点

```python
# 将第 3 个状态设为 0 能量
plotter.plot(energies, labels, reference_index=2)
```

### 3. 精细控制插值

```python
style = DiagramStyle(
    interp_method='spline',         # 更平滑（推荐）
    # interp_method='quadratic',    # 更尖锐
    interp_points=200,              # 更多点 = 更平滑
    platform_length=1.2,            # 更长的平台
    ts_peak_width=0.6               # 更窄的峰（更尖锐）
)
```

### 4. 导出多种格式

```python
plotter.save('diagram.png', dpi=300)    # 高分辨率 PNG
plotter.save('diagram.pdf')             # 矢量 PDF
plotter.save('diagram.svg')             # 矢量 SVG
plotter.export_data('data.csv')         # 数据 CSV
```

## 实际案例

### 案例 1：CO 氧化机理

```python
from energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

# Pt(111) 表面 CO 氧化
steps = [
    {'label': 'CO* + O*', 'energy': 0.00, 'is_ts': False},
    {'label': 'TS-CO-O', 'energy': 0.95, 'is_ts': True},
    {'label': 'CO-O*', 'energy': 0.15, 'is_ts': False},
    {'label': 'TS-CO₂', 'energy': 0.85, 'is_ts': True},
    {'label': 'CO₂ + 2*', 'energy': -2.35, 'is_ts': False},
]

style = DiagramStyle(
    color_scheme='nature',
    ylabel='Free Energy',
    ylabel_units='eV',
    figsize=(10, 6),
    use_shadow=True,
    show_energies=True
)

plotter = EnergyDiagramPlotter(style=style)
plotter.plot(steps)
plotter.save('co_oxidation_pt111.png', dpi=300)
```

### 案例 2：氢析出反应（HER）

```python
# Volmer-Heyrovsky 机理
steps = [
    {'label': 'H⁺ + e⁻', 'energy': 0.0, 'is_ts': False},
    {'label': 'TS-Volmer', 'energy': 0.3, 'is_ts': True},
    {'label': 'H*', 'energy': -0.1, 'is_ts': False},
    {'label': 'TS-Heyrovsky', 'energy': 0.5, 'is_ts': True},
    {'label': 'H₂', 'energy': -0.8, 'is_ts': False},
]

style = DiagramStyle(
    color_scheme='science',
    ylabel='Free Energy',
    ylabel_units='eV vs. RHE',
    platform_length=1.5,
    ts_peak_width=0.6
)

plotter = EnergyDiagramPlotter(style=style)
plotter.plot(steps)
plotter.save('her_mechanism.png', dpi=300)
```

### 案例 3：CO₂ 还原

```python
# CO₂ → CO 机理
co2_reduction = [
    {'label': 'CO₂ + *', 'energy': 0.0, 'is_ts': False},
    {'label': 'TS-COOH', 'energy': 1.2, 'is_ts': True},
    {'label': 'COOH*', 'energy': 0.5, 'is_ts': False},
    {'label': 'TS-CO', 'energy': 0.8, 'is_ts': True},
    {'label': 'CO* + H₂O', 'energy': -0.3, 'is_ts': False},
]

style = DiagramStyle(
    color_scheme='elegant',
    ylabel='Free Energy',
    ylabel_units='eV',
    line_width=3.0,
    use_shadow=True,
    shadow_offset=0.08
)

plotter = EnergyDiagramPlotter(style=style)
plotter.plot(co2_reduction)
plotter.save('co2_reduction.png', dpi=300)
```

## 与 pMuTT 集成

可以直接从 pMuTT 的 Reaction 对象提取能量：

```python
from pmutt.reaction import Reaction
from energy_diagram_plotter import EnergyDiagramPlotter

# 假设你有 pMuTT Reaction 对象
reactions = [rxn1, rxn2, rxn3]

# 提取能量和标签
energies = []
labels = []
is_ts = []

for rxn in reactions:
    # 反应物
    E_reactants = rxn.get_H(T=500, units='kcal/mol', state='reactants')
    energies.append(E_reactants)
    labels.append(rxn.get_species_str(state='reactants'))
    is_ts.append(False)

    # 过渡态（如果有）
    if rxn.transition_state:
        E_ts = rxn.get_H(T=500, units='kcal/mol', state='transition state')
        energies.append(E_ts)
        labels.append('TS')
        is_ts.append(True)

    # 产物
    E_products = rxn.get_H(T=500, units='kcal/mol', state='products')
    energies.append(E_products)
    labels.append(rxn.get_species_str(state='products'))
    is_ts.append(False)

# 绘图
plotter = EnergyDiagramPlotter()
plotter.plot(energies, labels, is_ts)
plotter.save('from_pmutt.png')
```

## 常见问题

### Q1: 如何调整 TS 峰的尖锐程度？

```python
style = DiagramStyle(
    ts_peak_width=0.5,     # 较小的值 → 更尖锐
    interp_method='quadratic'  # 使用二次插值更尖锐
)
```

### Q2: 如何去掉能量数值标签？

```python
style = DiagramStyle(show_energies=False)
```

### Q3: 如何更改能量单位？

```python
style = DiagramStyle(
    ylabel='Gibbs Free Energy',
    ylabel_units='kJ/mol'  # 或 'eV', 'kcal/mol', 'Ha' 等
)
```

### Q4: 如何自定义颜色？

```python
# 方法 1：修改预设配色
style = DiagramStyle(
    color_scheme='material',
    primary_color='#FF5722',  # 覆盖主色
    ts_color='#9C27B0'        # 覆盖TS色
)

# 方法 2：为每个步骤指定颜色
steps = [
    {'label': 'A', 'energy': 0.0, 'color': '#FF0000'},
    {'label': 'TS', 'energy': 1.2, 'is_ts': True, 'color': '#00FF00'},
    # ...
]
```

### Q5: 为什么能垒曲线看起来不平滑？

增加插值点数：
```python
style = DiagramStyle(interp_points=200)  # 默认100
```

### Q6: 如何设置图片尺寸适合期刊要求？

```python
# Nature: 单栏 89 mm, 双栏 183 mm
# 转换：1 inch = 25.4 mm
style = DiagramStyle(
    figsize=(7.2, 4.5),  # 183 mm 宽
    dpi=300              # Nature 要求 ≥300 DPI
)
```

## API 参考

### 主要类

#### `EnergyDiagramPlotter`
主绘图类。

**方法**：
- `__init__(style=None)`: 初始化
- `plot(energies, labels, is_ts_list, reference_index, show)`: 绘制图表
- `save(filename, dpi, **kwargs)`: 保存图片
- `export_data(filename)`: 导出数据到 CSV

#### `DiagramStyle`
样式配置类，包含所有可自定义参数。

#### `EnergyStep`
单个能量步骤的数据类。

**属性**：
- `label`: 标签
- `energy`: 能量值
- `is_ts`: 是否为过渡态
- `color`: 自定义颜色
- `line_style`: 线型

### 便捷函数

#### `quick_plot()`
一行代码绘图。

```python
quick_plot(energies, labels, is_ts_list, output, color_scheme, **style_kwargs)
```

#### `format_chemical_formula()`
格式化化学式为 LaTeX。

```python
formatted = format_chemical_formula('H2O')  # → 'H$_2$O'
```

## 性能优化

对于大量数据点：
```python
style = DiagramStyle(
    interp_points=50,   # 减少插值点（默认100）
    dpi=150             # 降低分辨率（默认300）
)
```

## 已知限制

1. **LaTeX 上下标**：需要 matplotlib 支持，某些系统可能显示异常
2. **多条路径**：当前版本需要手动创建子图
3. **动画**：不支持动态能量图

## 贡献和反馈

如有问题或建议，请在 CatDT 仓库提交 Issue。

## License

MIT License

---

**更新日期**: 2026-01-18
**版本**: 1.0.0
