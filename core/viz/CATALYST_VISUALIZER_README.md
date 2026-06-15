# Catalyst Surface Visualizer - 使用说明

高质量催化剂表面可视化工具，支持多种渲染后端（Tachyon、OSPRay、OpenGL、Anari、matplotlib）。

## 主要特性

- 🎨 **多种渲染器**：支持 5 种渲染后端，从快速预览到出版物级质量
- 🔬 **自动扩胞**：小结构（<15原子）自动扩展为超胞，便于观察
- 📹 **轨迹动画**：支持生成高质量 GIF 动画
- 🎯 **灵活配置**：可自定义视角、原子大小、质量等级
- 🧊 **晶格显示**：可选显示晶胞边界
- 🌈 **标准配色**：使用 Jmol 配色方案

---

## 安装依赖

### 基础依赖（必需）
```bash
pip install ase matplotlib pillow numpy
```

### OVITO 依赖（高质量渲染）
```bash
pip install ovito
```

**注意**：OVITO 是使用 Tachyon、OSPRay、OpenGL、Anari 渲染器所必需的。如果只使用 matplotlib，则不需要安装。

---

## 快速开始

### 1. 基础使用（默认设置）
```bash
# 使用默认 Tachyon 渲染器，自动扩胞
python catalyst_surface_visualizer.py structure.xyz
```

**输出**：`structure.png`（2048x2048，高质量）

### 2. 指定渲染器
```bash
# 使用 OSPRay 渲染器
python catalyst_surface_visualizer.py structure.xyz --renderer ospray

# 使用 matplotlib（无需 OVITO）
python catalyst_surface_visualizer.py structure.xyz --renderer matplotlib
```

### 3. 调整质量
```bash
# 超高质量（4096x4096）
python catalyst_surface_visualizer.py structure.xyz --quality ultra

# 快速预览（800x800）
python catalyst_surface_visualizer.py structure.xyz --quality low
```

### 4. 轨迹动画
```bash
# 生成 GIF 动画
python catalyst_surface_visualizer.py trajectory.traj --fps 2
```

---

## 渲染器对比

| 渲染器 | 速度 | 质量 | 阴影 | 环境光遮蔽 | 推荐场景 |
|--------|------|------|------|------------|----------|
| **tachyon** ⭐ | 中等 | 高 | ✅ | ✅ | 默认选择，平衡质量和速度 |
| **ospray** | 慢 | 极高 | ✅ | ✅ | 出版物图片，最佳质量 |
| **opengl** | 快 | 中 | ❌ | ❌ | 快速预览，实时渲染 |
| **anari** | 中等 | 高 | ✅ | ✅ | 硬件加速（需 NVIDIA GPU） |
| **matplotlib** | 快 | 中低 | ⚠️ | ❌ | 不安装 OVITO 时使用 |

---

## 完整参数说明

### 必需参数
```bash
input                 # 输入文件（xyz, traj, cif, POSCAR 等）
```

### 输出选项
```bash
-o, --output FILE     # 输出文件名（默认：自动生成）
```

### 渲染选项
```bash
--renderer RENDERER   # 渲染器选择
                      # 可选：tachyon, ospray, opengl, anari, matplotlib
                      # 默认：tachyon

--quality QUALITY     # 渲染质量
                      # 可选：low, medium, high, ultra
                      # 默认：high

--elevation DEGREE    # 视角仰角（度）
                      # 默认：30（30度斜上方往下看）

--azimuth DEGREE      # 视角方位角（度）
                      # 默认：45

--scale FLOAT         # 原子半径缩放因子
                      # 默认：1.0
                      # 建议范围：0.5-2.0

--background COLOR    # 背景颜色
                      # 默认：white
```

### 结构操作
```bash
--show-cell           # 显示晶胞边界
                      # 默认：不显示

--no-expand           # 禁用自动扩胞
                      # 默认：启用自动扩胞

--expand-threshold N  # 自动扩胞阈值（原子数）
                      # 默认：15

--supercell NX,NY,NZ  # 超胞尺寸
                      # 默认：2,2,1
                      # 格式：逗号分隔，无空格
```

### 动画选项
```bash
--fps N               # GIF 帧率
                      # 默认：2
```

### 其他选项
```bash
--color-scheme SCHEME # 配色方案
                      # 可选：jmol, cpk
                      # 默认：jmol

--show                # 显示图像（仅 matplotlib）
```

---

## 使用示例

### 示例 1：默认渲染
```bash
python catalyst_surface_visualizer.py CO_on_Cu111.xyz
```
**效果**：
- 使用 Tachyon 渲染器
- 如果原子数 < 15，自动扩展为 2x2x1 超胞
- 30° 俯视角度
- 高质量输出（2048x2048）

### 示例 2：出版物级图片
```bash
python catalyst_surface_visualizer.py surface.cif \
    --renderer ospray \
    --quality ultra \
    --scale 1.2 \
    -o publication_figure.png
```
**效果**：
- OSPRay 光线追踪
- 超高分辨率（4096x4096）
- 原子稍大（1.2倍）

### 示例 3：显示晶胞
```bash
python catalyst_surface_visualizer.py POSCAR \
    --show-cell \
    --no-expand
```
**效果**：
- 显示晶胞边界
- 不自动扩胞，保持原始结构

### 示例 4：自定义视角
```bash
python catalyst_surface_visualizer.py structure.xyz \
    --elevation 60 \
    --azimuth 30 \
    --scale 1.5
```
**效果**：
- 60° 仰角（更高角度）
- 30° 方位角
- 原子放大 1.5 倍

### 示例 5：快速预览
```bash
python catalyst_surface_visualizer.py structure.xyz \
    --renderer opengl \
    --quality low
```
**效果**：
- OpenGL 快速渲染
- 低分辨率（800x800）
- 适合快速检查结构

### 示例 6：小分子可视化
```bash
python catalyst_surface_visualizer.py CO2.xyz \
    --no-expand \
    --scale 2.0 \
    --quality high
```
**效果**：
- 禁用自动扩胞（小分子不需要）
- 原子放大 2 倍以便观察
- 高质量渲染

### 示例 7：轨迹动画
```bash
python catalyst_surface_visualizer.py MD_trajectory.traj \
    --fps 5 \
    --quality medium \
    -o animation.gif
```
**效果**：
- 5 帧/秒的流畅动画
- 中等质量（平衡速度和质量）
- 自动扩胞应用于所有帧

### 示例 8：大超胞
```bash
python catalyst_surface_visualizer.py small_unit.xyz \
    --supercell 3,3,2 \
    --expand-threshold 20
```
**效果**：
- 扩展为 3x3x2 超胞
- 当原子数 < 20 时触发扩胞

### 示例 9：不同渲染器对比
```bash
# 生成多个渲染器对比图
for renderer in tachyon ospray opengl anari; do
    python catalyst_surface_visualizer.py structure.xyz \
        --renderer $renderer \
        --no-expand \
        --quality high \
        -o comparison_${renderer}.png
done
```

### 示例 10：批量处理
```bash
# 批量渲染所有 xyz 文件
for file in *.xyz; do
    python catalyst_surface_visualizer.py "$file" \
        --quality high \
        -o "${file%.xyz}.png"
done
```

---

## 质量等级详解

### Low（快速预览）
- **分辨率**：800x800
- **采样**：4 samples/pixel
- **用途**：快速检查结构，测试参数
- **渲染时间**：~2-5秒

### Medium（日常使用）
- **分辨率**：1200x1200
- **采样**：8 samples/pixel
- **用途**：报告、展示
- **渲染时间**：~5-10秒

### High（推荐）⭐
- **分辨率**：2048x2048
- **采样**：16 samples/pixel
- **用途**：论文、海报
- **渲染时间**：~10-20秒

### Ultra（极致质量）
- **分辨率**：4096x4096
- **采样**：32 samples/pixel
- **用途**：期刊封面、高分辨率出版物
- **渲染时间**：~30-60秒

---

## 自动扩胞功能

### 工作原理
当结构原子数少于阈值（默认 15）时，自动扩展为超胞以便于观察表面结构。

### 配置选项
```bash
# 修改扩胞阈值
--expand-threshold 20

# 自定义超胞尺寸
--supercell 3,3,1

# 完全禁用扩胞
--no-expand
```

### 典型场景
- ✅ **小单胞**（2-6 原子）：自动扩展为 2x2x1
- ✅ **表面原子团**（<15 原子）：扩展后更直观
- ❌ **大超胞**（>15 原子）：不扩展，直接渲染
- ❌ **小分子**（使用 `--no-expand`）：保持原样

---

## 支持的文件格式

通过 ASE 支持多种格式：

| 格式 | 扩展名 | 说明 |
|------|--------|------|
| XYZ | `.xyz` | 通用原子坐标格式 |
| CIF | `.cif` | 晶体学信息文件 |
| POSCAR | `POSCAR`, `CONTCAR` | VASP 结构文件 |
| Trajectory | `.traj` | ASE 轨迹文件 |
| PDB | `.pdb` | 蛋白质数据库格式 |
| LAMMPS | `.lmp`, `.lammps` | LAMMPS 数据文件 |
| Gaussian | `.log`, `.out` | Gaussian 输出 |

**轨迹格式**：自动识别多帧文件，生成 GIF 动画

---

## 常见问题 (FAQ)

### Q1: 如何选择渲染器？
**A**:
- **日常使用**：`tachyon`（默认，质量好速度快）
- **最佳质量**：`ospray`（光线追踪，适合出版）
- **快速预览**：`opengl`（实时渲染）
- **无 OVITO**：`matplotlib`（无需额外依赖）

### Q2: 为什么图片里有晶格线条？
**A**: 默认情况下不显示晶格线条。如需显示，添加 `--show-cell` 参数。

### Q3: 如何控制原子大小？
**A**: 使用 `--scale` 参数：
```bash
# 原子变小
--scale 0.8

# 原子变大
--scale 1.5
```

### Q4: 扩胞后结构太大怎么办？
**A**:
```bash
# 使用 1x1x1（不扩展）
--supercell 1,1,1

# 或完全禁用
--no-expand
```

### Q5: 渲染很慢怎么办？
**A**:
- 降低质量：`--quality low` 或 `--quality medium`
- 使用更快的渲染器：`--renderer opengl`
- 减小结构：使用 `--no-expand`

### Q6: OpenGL 报错 "Qt platform plugin" 怎么办？
**A**:
```bash
# 设置环境变量
export QT_QPA_PLATFORM=offscreen
python catalyst_surface_visualizer.py structure.xyz --renderer opengl
```

### Q7: 如何调整视角？
**A**:
```bash
# 俯视（elevation 0-30）
--elevation 20 --azimuth 45

# 侧视（elevation 30-60）
--elevation 45 --azimuth 0

# 仰视（elevation 60-90）
--elevation 75 --azimuth 30
```

### Q8: 轨迹动画帧率如何选择？
**A**:
- **慢动画**：`--fps 1-2`（适合详细观察）
- **正常**：`--fps 5-10`（流畅展示）
- **快速**：`--fps 15-30`（快速浏览）

### Q9: 内存不足怎么办？
**A**:
- 使用 `--quality low` 或 `--quality medium`
- 对于轨迹，减少帧数或分段处理
- 禁用扩胞：`--no-expand`

### Q10: 如何生成透明背景？
**A**: 目前不支持透明背景。可以使用白色背景后期处理：
```bash
--background white
```

---

## 性能优化建议

### 1. 快速迭代
```bash
# 用 low 质量快速测试参数
python catalyst_surface_visualizer.py structure.xyz \
    --quality low --renderer opengl

# 确认效果后再用高质量渲染
python catalyst_surface_visualizer.py structure.xyz \
    --quality ultra --renderer ospray
```

### 2. 批量处理优化
```bash
# 使用并行处理（需要 GNU parallel）
parallel -j 4 python catalyst_surface_visualizer.py {} ::: *.xyz
```

### 3. 大轨迹处理
```bash
# 抽帧处理（每 5 帧取 1 帧）
python -c "
from ase.io import read, write
traj = read('large_traj.traj', index='::5')
write('sampled_traj.traj', traj)
"

# 然后渲染采样后的轨迹
python catalyst_surface_visualizer.py sampled_traj.traj
```

---

## 输出文件说明

### 单帧输出
```
structure.png
├── 格式：PNG
├── 分辨率：根据 quality 设置
├── 色深：24-bit RGB
└── 大小：~50KB - 500KB
```

### 轨迹输出
```
trajectory.gif
├── 格式：GIF
├── 帧数：等于轨迹帧数
├── 帧率：--fps 设置
└── 大小：~100KB - 2MB
```

---

## 技术细节

### 默认配置
- **渲染器**：Tachyon（高质量光线追踪）
- **视角**：elevation=30°, azimuth=45°（30度斜上方俯视）
- **质量**：high（2048x2048, 16 samples/pixel）
- **扩胞**：自动（<15原子时扩展为2x2x1）
- **晶格**：隐藏
- **配色**：Jmol 标准配色

### 坐标系统
- **X轴**：水平向右
- **Y轴**：垂直屏幕向外
- **Z轴**：垂直向上
- **原点**：结构质心

### 相机设置
- **投影**：透视投影（FOV = 35°）
- **距离**：自动调整（结构尺寸 × 3）
- **对焦**：始终对准结构中心

---

## 故障排除

### 安装问题
```bash
# ASE 安装失败
pip install --upgrade pip
pip install ase --no-cache-dir

# OVITO 安装失败（需要 Python 3.8-3.11）
pip install ovito --pre

# matplotlib 显示问题
pip install --upgrade matplotlib pillow
```

### 渲染问题
```bash
# OVITO 模块导入失败
python -c "import ovito; print(ovito.__version__)"

# 检查可用渲染器
python -c "
from ovito.vis import TachyonRenderer, OSPRayRenderer
print('Tachyon: OK')
print('OSPRay: OK')
"
```

### 文件读取问题
```bash
# 检查文件格式
python -c "from ase.io import read; atoms = read('structure.xyz'); print(len(atoms))"

# 转换文件格式
python -c "from ase.io import read, write; write('output.xyz', read('input.cif'))"
```

---

## 版本信息

**当前版本**：1.0.0
**更新日期**：2026-01-18
**兼容性**：
- Python 3.8+
- ASE 3.22+
- OVITO 3.9+
- matplotlib 3.5+

---

## 作者与致谢

**作者**：CatDT Team
**联系**：[项目仓库]

**依赖项目**：
- [ASE](https://wiki.fysik.dtu.dk/ase/) - 原子模拟环境
- [OVITO](https://www.ovito.org/) - 开放可视化工具
- [matplotlib](https://matplotlib.org/) - Python 绘图库

---

## 许可证

本工具基于依赖库的许可证使用。请查看各依赖项目的许可证要求。

---

## 更新日志

### v1.0.0 (2026-01-18)
- ✨ 初始版本发布
- ✅ 支持 5 种渲染器（Tachyon, OSPRay, OpenGL, Anari, matplotlib）
- ✅ 自动扩胞功能
- ✅ 轨迹动画支持
- ✅ 可配置晶格显示
- ✅ 多种质量等级

---

## 快速参考卡片

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Catalyst Surface Visualizer - 快速参考
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

基本用法：
  python catalyst_surface_visualizer.py structure.xyz

常用选项：
  --renderer tachyon|ospray|opengl|anari|matplotlib
  --quality low|medium|high|ultra
  --show-cell               显示晶格
  --no-expand              禁用扩胞
  --elevation 30           视角仰角
  --azimuth 45             视角方位角
  --scale 1.5              原子大小
  -o output.png            输出文件

轨迹动画：
  python catalyst_surface_visualizer.py traj.traj --fps 5

质量对比：
  low    : 800x800,  4 spp  (~3秒)
  medium : 1200x1200, 8 spp  (~7秒)
  high   : 2048x2048, 16 spp (~15秒) ⭐
  ultra  : 4096x4096, 32 spp (~45秒)

渲染器对比：
  tachyon    : 默认，高质量  ⭐
  ospray     : 最佳质量，慢
  opengl     : 快速预览
  anari      : 硬件加速
  matplotlib : 无需 OVITO

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```
