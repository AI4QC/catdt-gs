---
name: catdt-router
description: Intelligent routing skill for CatDT (Catalysis Digital Twin). Maps any catalysis-related user request to the correct module(s) — surface prediction, adsorption, reconstruction, NEB barriers, mechanism search (UniMech), microkinetics, active learning, visualization — with automatic parameter configuration and mandatory trajectory visualization.
---

# CatDT Intelligent Router

CatDT is a multi-scale heterogeneous catalysis simulation platform with 15+ standalone modules. This skill routes user requests to the correct module(s), configures parameters, generates runnable Python scripts, and enforces auto-visualization of all structural evolution.

## How to Use This Skill

1. Match user intent against the **Intent Catalog** (Section 1)
2. Read `config.json` for environment paths and default parameters
3. For API details, search `references/files.md` for the relevant module (use `## File:` markers)
4. For UniMech mechanism search, read `references/unimech.md`
5. For visualization APIs, read `references/visualization.md`
6. Adapt the matching template from `templates/` (Section 6)
7. **Always add auto-visualization** — see Section 5

## Environment Setup

All templates use `CATDT_ROOT` for portability:
```python
import os, sys
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)
```
Dependency paths (relative to CATDT_ROOT) and default parameters are in `config.json`.

---

## 1. Intent Catalog — What Users Might Ask

### 1.1 Surface & Structure Generation

| User Says (examples) | Action | Template |
|---|---|---|
| "预测Cu的暴露面" / "predict exposed facets of Pt" | SurFF Wulff prediction | `templates/surface.py` |
| "生成Cu(111) slab" / "build a Ni(100) surface" | Direct slab build via ASE, skip SurFF | `templates/surface.py` (set `ELEMENT`+`MILLER`) |
| "Cu和Ni合金的表面能" / "surface energy of PtRu alloy" | SurFF on alloy bulk | `templates/surface.py` |
| "TiO2(110)表面" / "SrTiO3(001) surface" | Oxide slab build | `templates/surface.py` (use `ase.build`) |
| "生成Wulff形貌" / "Wulff construction" | SurFF → Wulff shape | `templates/surface.py` |
| "有台阶的表面" / "stepped surface like Cu(211)" | Higher-index slab build | `templates/surface.py` |
| "掺杂表面" / "Ni-doped TiO2 surface" | Build slab → substitute atoms | Custom from `templates/surface.py` |

### 1.2 Adsorption

| User Says (examples) | Action | Template |
|---|---|---|
| "CO在Cu(111)上的吸附位点" / "adsorption sites for OH on Pt" | AdsorbDiff prediction | `templates/adsorption.py` |
| "吸附能是多少" / "binding energy of *O on Pd" | AdsorbDiff + energy extraction | `templates/adsorption.py` |
| "比较不同位点的吸附能" / "compare top, bridge, hollow" | AdsorbDiff with multiple sites | `templates/adsorption.py` (set `NUM_SITES` high) |
| "共吸附" / "co-adsorption of CO and O" | Sequential AdsorbDiff calls | Chain two `templates/adsorption.py` |
| "在重构后的表面上吸附" / "adsorb on reconstructed surface" | VSSR-MC → AdsorbDiff chain | `templates/reconstruction.py` → `templates/adsorption.py` |
| "不同覆盖度下的吸附" / "coverage-dependent adsorption" | Multiple AdsorbDiff runs | Loop over `templates/adsorption.py` |

**Critical**: Adsorbate names MUST use `*` prefix: `*CO`, `*OH`, `*H`, `*O`, `*OOH`. Never pass bare `CO`.

### 1.3 Surface Reconstruction

| User Says (examples) | Action | Template |
|---|---|---|
| "500K下Cu表面重构" / "surface reconstruction at 800K" | Thermal VSSR-MC | `templates/reconstruction.py` |
| "反应条件下的表面结构" / "surface under operating conditions" | VSSR-MC at reaction T | `templates/reconstruction.py` |
| "SMSI效应" / "strong metal-support interaction" | VSSR-MC on metal/oxide interface | `templates/reconstruction.py` |
| "退火模拟" / "annealing simulation" | VSSR-MC with annealing schedule | `templates/reconstruction.py` (set `perform_annealing=True`) |
| "正则系综MC" / "canonical ensemble MC" | Fixed adsorbate count | `templates/reconstruction.py` (set `CANONICAL=True`) |
| "电化学条件下的表面" / "surface at 1.23V vs SHE, pH=14" | Electrochemical VSSR-MC (Pourbaix) | `templates/reconstruction.py` (set `POTENTIAL_SHE`, `PH`) |
| "Pourbaix图" / "Pourbaix diagram" | ElectrochemicalSurfacePredictor | Search `references/files.md` for `electrochemical_surface_predictor` |
| "表面稳定性" / "surface stability at different T" | Multiple VSSR-MC at different T | Loop over `templates/reconstruction.py` |

### 1.4 Reaction Barriers & NEB

| User Says (examples) | Action | Template |
|---|---|---|
| "CO*→CHO*的能垒" / "activation energy for *O + *H → *OH" | NEB from species labels (need structures) | `templates/barrier_with_viz.py` |
| "我有初末态结构，算NEB" / "compute NEB from reactant/product" | Direct NEB from files | `templates/barrier_with_viz.py` |
| "过渡态结构" / "find transition state" | NEB → extract TS | `templates/barrier_with_viz.py` |
| "反应能垒太高了，检查一下" / "barrier seems wrong, debug" | Check endpoints, validate, re-run NEB | `templates/pathway_step.py` (or `EnhancedNEBValidator` directly) |
| "我只有反应物，自己构造产物再算能垒" / "design product from reactant + edit plan, validate, then NEB" | Agent4-style edit → Agent5 validation → NEB | `templates/pathway_step.py` (MODE="edit") |
| "校验一对NEB初末态是否合规" / "validate a NEB endpoint pair before running" | Programmatic + reaction-aware validation only | `templates/pathway_step.py` (MODE="endpoints", skip NEB section) |
| "弛豫这个结构" / "relax this structure" | FairchemPredictor relax | `templates/relaxation_with_viz.py` (MODE="relax") |
| "单点能计算" / "single point energy" | FairchemPredictor energy | `templates/relaxation_with_viz.py` (MODE="single_point") |
| "用更大的模型" / "use larger model for accuracy" | Switch to uma-m-1p1 | Set `MODEL="uma-m-1p1"` in template |

**Warning**: Staged-atom NEB (H placed on surface hollow near adsorbate) cannot reliably compute hydrogenation barriers. Only bond-breaking/dissociation steps give trustworthy Ea.

### 1.5 Mechanism & Pathway Search (UniMech)

| User Says (examples) | Action | Template |
|---|---|---|
| "CO加氢到甲醇的所有可能路径" / "find all pathways for CO2 → methanol" | UniMech systematic search | `templates/mechanism_search.py` |
| "比较不同反应机理" / "compare competing mechanisms" | UniMech both_sequential | `templates/mechanism_search.py` (mode=`both_sequential`) |
| "快速筛选反应路径" / "quick pathway screening" | UniMech agent_guided | `templates/mechanism_search.py` (mode=`agent_guided`) |
| "反应网络" / "build reaction network" | UniMech systematic | `templates/mechanism_search.py` (mode=`systematic`) |
| "找到速率决定步骤" / "identify rate-determining step" | UniMech + NEB on retained paths | Chain `templates/mechanism_search.py` → `templates/barrier_with_viz.py` |
| "PDH机理" / "propane dehydrogenation mechanism" | UniMech for C3H8 → C3H6 | `templates/mechanism_search.py` |
| "FT合成路径" / "Fischer-Tropsch pathway" | UniMech for CO+H2 → CnHm | `templates/mechanism_search.py` |

See `references/unimech.md` for full API — 3 exploration modes, beam search with sibling pruning, RDKit-based operations.

### 1.6 Microkinetics & KMC

| User Says (examples) | Action | Template |
|---|---|---|
| "计算TOF" / "turnover frequency" | CatMAP microkinetics | `templates/microkinetics.py` |
| "选择性分析" / "product selectivity" | CatMAP with multiple products | `templates/microkinetics.py` |
| "覆盖度" / "surface coverage under reaction conditions" | CatMAP steady-state | `templates/microkinetics.py` |
| "温度对活性的影响" / "activity vs temperature" | CatMAP sweep over T | Loop `templates/microkinetics.py` |
| "Arrhenius图" / "Arrhenius plot" | CatMAP at multiple T + plot | `templates/microkinetics.py` + matplotlib |
| "热力学校正" / "thermodynamic corrections, partition functions" | pMuTT | Search `references/files.md` for `pmutt_predictor` |
| "电化学电流密度" / "current density, Faradaic efficiency" | CatMAP electrochemical mode | `templates/microkinetics.py` (electrochemical) |

### 1.7 Molecular Dynamics

| User Says (examples) | Action | Template |
|---|---|---|
| "MD模拟" / "run molecular dynamics" | Langevin MD via FairChem | `templates/relaxation_with_viz.py` (MODE="md") |
| "看看这个表面在高温下稳定吗" / "thermal stability at 1000K" | MD at target T | `templates/relaxation_with_viz.py` (MODE="md") |
| "退火" / "simulated annealing" | MD with T schedule | Custom from `templates/relaxation_with_viz.py` |
| "扩散系数" / "diffusion coefficient" | MD + MSD analysis | `templates/relaxation_with_viz.py` + post-processing |

### 1.8 Active Learning & Model Training

| User Says (examples) | Action | Template |
|---|---|---|
| "训练MACE力场" / "train CP-MACE model" | CPMACEAdapter | Search `references/files.md` for `cp_mace_predictor` |
| "主动学习" / "active learning loop" | AL framework | Search `references/files.md` for `framework.py` (active_learning) |
| "constant-potential模拟" / "constant potential NEB" | CPMACEPredictor simulate | Search `references/files.md` for `cp_mace_predictor` |

### 1.9 Visualization & Analysis

| User Says (examples) | Action | Template |
|---|---|---|
| "画能量图" / "plot energy diagram" | EnergyDiagramPlotter | `templates/visualization.py` |
| "渲染这个结构" / "render structure as PNG" | CatalystSurfaceVisualizer | `templates/visualization.py` |
| "轨迹动画" / "animate trajectory as GIF" | visualize_trajectory | `templates/visualization.py` |
| "出一张发表级的图" / "publication-quality figure" | DiagramStyle with nature/science scheme | `templates/visualization.py` |
| "看看这个.traj文件" / "visualize this trajectory file" | Read .traj → GIF | `templates/visualization.py` |
| "NEB路径动画" / "animate NEB path" | visualize_neb_trajectory | See `references/visualization.md` |
| "MC重构过程动画" / "MC reconstruction animation" | visualize_mc_reconstruction | See `references/visualization.md` |
| "把output目录的所有结果整合成HTML报告" / "aggregate scattered .traj/.vasp/.json into one HTML report" | Walk dir → render GIFs+PNGs+energy diagrams → index HTML | `templates/report_aggregator.py` |
| "生成HTML报告" / "generate visualization summary" | CatalysisVisualizationManager | `templates/report_aggregator.py` (or see `references/visualization.md` for advanced) |

### 1.10 Full Pipeline & Compound Requests

| User Says (examples) | Action | Template |
|---|---|---|
| "从Cu bulk跑完整流程" / "end-to-end digital twin" | CatDTCamelWorkflow (needs LLM) | `templates/full_pipeline.py` |
| "CO氧化在Pt上的完整路径和能垒" / "full CO oxidation on Pt" | SurFF → AdsorbDiff → VSSR-MC → NEB chain | Chain multiple templates |
| "比较不同晶面的催化活性" / "compare activity across facets" | Multi-facet: SurFF top-N → NEB each | `templates/full_pipeline.py` (MULTI_FACET=True) |
| "PDH在Ni@TiOx上的研究" / "study PDH on supported catalyst" | Build interface → VSSR-MC → NEB | Chain templates |
| "帮我筛选催化剂" / "screen catalysts for HER" | Multiple surfaces × NEB → compare Ea | Script with loops |

### 1.11 Common Reaction Types — Quick Parameter Reference

| Reaction | Initial State | Target | Key Intermediates |
|---|---|---|---|
| CO oxidation | *CO + *O | CO2(g) | *CO, *O, *OCO |
| CO hydrogenation → methanol | *CO | CH3OH(g) | *CHO, *CH2O, *CH3O, *CH3OH |
| CO hydrogenation → methane | *CO | CH4(g) | *CHO, *CH2O, *CH3, *CH4 |
| CO2 reduction | *CO2 | *CO / CH4(g) / C2H4(g) | *COOH, *CO, *CHO, *C |
| PDH (propane dehydrogenation) | *C3H8 | C3H6(g) + H2(g) | *C3H7, *C3H6 |
| HER | *H | H2(g) | *H (Volmer-Heyrovsky/Tafel) |
| OER | *OH2 | O2(g) | *OH, *O, *OOH |
| ORR | *O2 | H2O(g) | *OOH, *O, *OH |
| NRR | *N2 | NH3(g) | *N2H, *N2H2, *NH, *NH2 |
| NH3 synthesis | N2(g) | NH3(g) | *N, *NH, *NH2, *NH3 |
| Water-gas shift | *CO + *H2O | CO2(g) + H2(g) | *COOH |
| Methane reforming | *CH4 + *H2O | *CO + H2(g) | *CH3, *CH2, *CH, *C, *OH |
| Ethanol dehydrogenation | *C2H5OH | CH3CHO(g) + H2(g) | *C2H5O, *CH3CHO |

---

## 2. Routing Decision Tree

```
User Request
│
├─ Already has structure files (reactant.vasp, product.vasp)?
│   ├─ Wants barrier/NEB → BarrierPredictor           [templates/barrier_with_viz.py]
│   ├─ Wants relaxation → FairchemPredictor            [templates/relaxation_with_viz.py]
│   ├─ Wants visualization → CatalystSurfaceVisualizer [templates/visualization.py]
│   └─ Wants MD → Langevin dynamics                    [templates/relaxation_with_viz.py MODE=md]
│
├─ Has surface + adsorbate name?
│   ├─ Wants adsorption sites → AdsorbDiff             [templates/adsorption.py]
│   ├─ Wants reconstruction → VSSR-MC                  [templates/reconstruction.py]
│   └─ Wants both → AdsorbDiff → VSSR-MC chain
│
├─ Has bulk structure / element name?
│   ├─ Wants surfaces → SurFF                          [templates/surface.py]
│   ├─ Wants full pipeline → CatDTCamelWorkflow        [templates/full_pipeline.py]
│   └─ Wants partial → SurFF → AdsorbDiff → ...chain
│
├─ Wants to explore mechanisms / reaction network?
│   └─ → UniMech mechanism search                      [templates/mechanism_search.py]
│       Mode: agent_guided (fast) | systematic (thorough) | both_sequential
│
├─ Has barriers, wants rates/kinetics?
│   └─ → CatMAP microkinetics                          [templates/microkinetics.py]
│
├─ Wants electrochemical analysis?
│   └─ → Set POTENTIAL_SHE + PH in reconstruction or full pipeline
│
├─ Wants energy diagram / visualization only?
│   └─ → Standalone visualization                      [templates/visualization.py]
│
└─ Unclear / complex multi-step request?
    └─ Break into steps, chain templates sequentially
```

---

## 3. Parameter Inference Rules

When user input is incomplete, infer parameters:

| User Provides | Inferred Action |
|---|---|
| Element name only ("Cu", "Pt", "铜") | `ase.build.bulk(element, crystal)` → SurFF |
| Element + Miller index ("Cu(111)") | Skip SurFF, `ase.build.fcc111(element, ...)` |
| Adsorbate without `*` ("CO", "OH") | Auto-prepend: "CO" → "*CO" |
| No temperature | Default 500 K (MC and KMC) |
| No NEB frames | Default 11 (optimal for CI-NEB) |
| No NEB fmax | Default 0.05 eV/A |
| No energy model for MC | Default "CHGNet" (CHGNetNFF) |
| No energy model for NEB/relax | Default "uma-s-1p1" |
| "electrochemical" / mentions potential or pH | Enable Pourbaix mode → require POTENTIAL_SHE + PH |
| No output directory | `output/{module}/{YYYYMMDD_HHMMSS}/` |
| "快速" / "quick" | Use agent_guided mode for mechanism; reduce MC sweeps |
| "精确" / "accurate" / "高精度" | Use uma-m-1p1; increase NEB frames to 15; increase MC sweeps |
| "发表级" / "publication quality" | Use ultra quality + nature/science color scheme for viz |

---

## 4. Module Import Quick Reference

```python
# ── Surface ──
from core.surface.surff_predictor import SurFFPredictor

# ── Adsorption ──
from core.reconstruction.adsorbdiff_predictor import AdsorbDiffPredictor

# ── Reconstruction ──
from core.reconstruction.vssr_mc_predictor import VSSRMCPredictor
from core.reconstruction.electrochemical_surface_predictor import ElectrochemicalSurfacePredictor
from core.reconstruction.cp_mace_predictor import CPMACEPredictor

# ── Energy & Relaxation ──
from core.pathway.fairchem_predictor import FairchemPredictor

# ── NEB Barriers ──
from core.pathway.barrier_predictor import BarrierPredictor

# ── Pathway Analysis ──
from core.pathway.pathway_predictor import PathwayPredictor

# ── Mechanism Search (UniMech) ──
from core.pathway.search_engine import MechanismSearchEngine
from core.pathway.candidate_generators import SystematicCRNExplorer, CandidateGeneratorRouter

# ── Microkinetics ──
from core.kmc.catmap_predictor import CatMAPPredictor
from core.kmc.pmutt_predictor import pMuTTPredictor

# ── Active Learning ──
from core.active_learning.framework import CPMACEAdapter, ALConfig

# ── Visualization ──
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter, quick_plot
from core.viz.visualization_manager import CatalysisVisualizationManager

# ── Tool Facade (all tools, no LLM) ──
from camel_agents.tools import CatDTTools

# ── Full Agent Pipeline (needs LLM API key) ──
from camel_agents.workflow import CatDTCamelWorkflow
from camel_agents.schemas import CatDTConfig
```

---

## 5. Auto-Visualization Policy (MANDATORY)

**Every computation producing structures or trajectories MUST generate visualization. Never return bare numerical results without visual output.**

| Computation | Must Produce | Implementation |
|---|---|---|
| Structure relaxation | Before/after PNG + trajectory GIF + energy convergence plot | `CatalystSurfaceVisualizer` + matplotlib |
| NEB barrier | Energy diagram PNG + NEB path GIF + TS structure PNG | `EnergyDiagramPlotter` + `CatalystSurfaceVisualizer` |
| MD simulation | Trajectory GIF + T(t) and E(t) evolution plot | `CatalystSurfaceVisualizer` + matplotlib |
| MC reconstruction | MC trajectory GIF + energy/acceptance history plot | `CatalystSurfaceVisualizer` + matplotlib |
| Adsorption sites | Top-N adsorption configuration PNGs | `CatalystSurfaceVisualizer` |
| Surface prediction | Top-N slab PNGs | `CatalystSurfaceVisualizer` |
| Mechanism search | Best pathway energy diagram PNG | `EnergyDiagramPlotter` |
| Microkinetics | Coverage bar chart + TOF summary | matplotlib |

See `references/visualization.md` for complete API and code recipes.

---

## 6. Template Reference

All templates are in the `templates/` directory. Each is a complete, runnable Python script.

| Template | Module(s) | When to Use |
|---|---|---|
| `templates/surface.py` | SurFF | Surface prediction, Wulff shape, slab generation |
| `templates/adsorption.py` | AdsorbDiff | Adsorption site prediction, binding energy |
| `templates/reconstruction.py` | VSSR-MC | Surface reconstruction (thermal + electrochemical) |
| `templates/barrier_with_viz.py` | BarrierPredictor | NEB activation barriers from existing reactant/product files |
| `templates/pathway_step.py` | EnhancedNEBValidator + BarrierPredictor | Agent4/5-style single-step: design product from edit plan, validate, then NEB |
| `templates/mechanism_search.py` | UniMech (CatDTTools) | Reaction pathway discovery, mechanism comparison |
| `templates/microkinetics.py` | CatMAP | TOF, selectivity, coverage, kinetics |
| `templates/relaxation_with_viz.py` | FairchemPredictor | Structure optimization, single-point energy, MD |
| `templates/full_pipeline.py` | CatDTCamelWorkflow | Full 8-agent digital twin pipeline |
| `templates/visualization.py` | Viz modules | Standalone rendering, energy diagrams, trajectory GIFs |
| `templates/report_aggregator.py` | CatalysisVisualizationManager | Walk an output dir, auto-generate GIFs/PNGs/energy figures, emit one HTML index |

**Usage pattern**: Copy template → modify the Configuration section at top → run.

---

## 7. Known Limitations

1. **Staged-atom NEB**: H on surface hollow near organic adsorbate is NOT a local minimum — NEB interior frames find 3-4 eV lower energies. Only dissociation/bond-breaking steps produce reliable barriers.
2. **AdsorbDiff `*` prefix**: Adsorbate names MUST start with `*`. Passing `CO` instead of `*CO` will fail silently or give wrong results.
3. **Full pipeline needs LLM**: `CatDTCamelWorkflow` requires `OPENAI_API_KEY` environment variable.
4. **GPU strongly recommended**: CPU mode works but 10-100x slower for all ML models (SurFF, AdsorbDiff, UMA, CHGNet).
5. **Electrochemical mode**: Requires `MP_API_KEY` for Materials Project Pourbaix data.
6. **OVITO for visualization**: `CatalystSurfaceVisualizer` requires OVITO Python package. Falls back to ASE if unavailable.
7. **RDKit for mechanism search**: `SystematicCRNExplorer` requires RDKit. Falls back to formula enumeration if unavailable.

---

## 8. Detailed Reference Files

| File | What It Contains |
|---|---|
| `config.json` | Configurable paths (deps, models), default parameters, visualization settings |
| `references/files.md` | Compressed source code of all 65 Python files — search with `## File: <path>` |
| `references/project-structure.md` | Directory tree with line counts |
| `references/summary.md` | Codebase statistics |
| `references/unimech.md` | UniMech mechanism search: 3 modes, beam search, data schemas, tool API |
| `references/visualization.md` | Visualization APIs, auto-viz recipes, output directory structure |
