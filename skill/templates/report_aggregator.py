"""CatDT Template: Aggregate Heterogeneous Outputs into a Single HTML Report

Scans an output directory for trajectory files (.traj), structure files
(.vasp / .cif / POSCAR*), and pickled CatDT results, then renders:

  - per-trajectory animation GIF (CatalystSurfaceVisualizer)
  - per-structure PNG (CatalystSurfaceVisualizer)
  - energy diagrams from any energy series found in JSON / pickle
  - a single index HTML page with publication-grade colour scheme

This template approximates Agent 7 (Orchestration & Visualization) for the
case where the user already has scattered results on disk and wants a tidy
report without re-running the full workflow.

Heuristics (intentionally simple — adapt for non-standard layouts):
  - Files matching `*neb*` or under directories containing "neb" → NEB path
  - Files matching `*mc*` or under directories containing "reconstruction"
    or "mc" → MC trajectory
  - All other .traj files → generic trajectory (treated as relaxation)
"""

import os, sys, json, pickle
from datetime import datetime
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────
CATDT_ROOT = os.environ.get("CATDT_ROOT", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CATDT_ROOT)

INPUT_DIR = "output"                                  # directory to scan recursively
OUTPUT_DIR = os.path.join(
    "output", "report", datetime.now().strftime("%Y%m%d_%H%M%S")
)
COLOR_SCHEME = "nature"                               # "material" | "nature" | "science" | "elegant"
QUALITY = "high"                                       # "low" | "medium" | "high" | "ultra"
RENDERER = "tachyon"                                  # OVITO renderer
MAX_TRAJECTORIES = 50                                 # safety cap
MAX_STRUCTURES = 50

# ── Execution ──────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

from ase.io import read

from core.viz.visualization_manager import CatalysisVisualizationManager

manager = CatalysisVisualizationManager(
    output_dir=OUTPUT_DIR,
    quality=QUALITY,
    renderer=RENDERER,
    color_scheme=COLOR_SCHEME,
)
viz_dir = manager.viz_dir

print(f"[Scan] Walking {INPUT_DIR} ...")
input_root = Path(INPUT_DIR).resolve()
if not input_root.is_dir():
    raise FileNotFoundError(f"INPUT_DIR not found: {input_root}")

trajectories: list[dict] = []
structures: list[dict] = []
energy_series: list[dict] = []

for path in sorted(input_root.rglob("*")):
    if not path.is_file():
        continue
    rel = path.relative_to(input_root)
    name_lc = path.name.lower()
    parent_lc = path.parent.name.lower()

    # Skip our own output dir to avoid infinite include.
    if Path(OUTPUT_DIR).resolve() in path.resolve().parents:
        continue

    if path.suffix == ".traj":
        if any(tag in name_lc or tag in parent_lc for tag in ("neb", "barrier", "pathway")):
            kind = "neb"
        elif any(tag in name_lc or tag in parent_lc for tag in ("mc", "reconstruction", "vssr")):
            kind = "mc"
        else:
            kind = "trajectory"
        trajectories.append({"path": path, "rel": rel, "kind": kind})
    elif path.suffix in (".vasp", ".cif") or path.name.startswith("POSCAR") or path.name.startswith("CONTCAR"):
        structures.append({"path": path, "rel": rel})
    elif path.suffix == ".json":
        # Try to find energy series inside.
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception:
            continue
        if isinstance(data, dict):
            for key in ("energies", "energy_history", "neb_energies", "mc_energy_history"):
                series = data.get(key)
                if isinstance(series, list) and len(series) >= 2 and all(isinstance(v, (int, float)) for v in series):
                    energy_series.append({"path": path, "rel": rel, "key": key, "values": series})

print(f"  trajectories: {len(trajectories)}  structures: {len(structures)}  energy series: {len(energy_series)}")

# Cap to keep the report sane.
trajectories = trajectories[:MAX_TRAJECTORIES]
structures = structures[:MAX_STRUCTURES]

viz_outputs: dict[str, str] = {}
sections: list[dict] = []

# ── 1. Trajectories → GIFs ─────────────────────────────────────
from core.viz.catalyst_surface_visualizer import CatalystSurfaceVisualizer

structure_viz = manager.structure_visualizer
for entry in trajectories:
    rel = entry["rel"]
    kind = entry["kind"]
    safe = str(rel).replace(os.sep, "_").replace(".", "_")
    gif_path = viz_dir / f"traj_{kind}_{safe}.gif"
    try:
        frames = read(str(entry["path"]), index=":")
        if not isinstance(frames, list):
            frames = [frames]
        if len(frames) < 2:
            png_path = viz_dir / f"struct_{safe}.png"
            structure_viz.visualize_structure(frames[0], output_file=str(png_path),
                                              title=str(rel), show=False)
            sections.append({"kind": "structure", "rel": rel, "img": png_path, "n_frames": 1})
            viz_outputs[f"struct_{safe}"] = str(png_path)
        else:
            structure_viz.visualize_trajectory(
                frames, output_file=str(gif_path), fps=3,
                titles=[f"{rel.name} F{i}" for i in range(len(frames))],
            )
            sections.append({"kind": kind, "rel": rel, "img": gif_path, "n_frames": len(frames)})
            viz_outputs[f"traj_{safe}"] = str(gif_path)
        print(f"  ✓ {rel} → {gif_path.name}")
    except Exception as exc:
        print(f"  ✗ {rel}: {exc}")

# ── 2. Standalone structures → PNGs ────────────────────────────
for entry in structures:
    rel = entry["rel"]
    safe = str(rel).replace(os.sep, "_").replace(".", "_")
    png_path = viz_dir / f"struct_{safe}.png"
    try:
        atoms = read(str(entry["path"]))
        structure_viz.visualize_structure(atoms, output_file=str(png_path),
                                          title=str(rel), show=False)
        sections.append({"kind": "structure", "rel": rel, "img": png_path, "n_frames": 1})
        viz_outputs[f"struct_{safe}"] = str(png_path)
        print(f"  ✓ {rel} → {png_path.name}")
    except Exception as exc:
        print(f"  ✗ {rel}: {exc}")

# ── 3. Energy series → diagrams ────────────────────────────────
from core.viz.energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

style = DiagramStyle(color_scheme=COLOR_SCHEME) if hasattr(DiagramStyle, "__init__") else None
for entry in energy_series:
    rel = entry["rel"]
    key = entry["key"]
    values = entry["values"]
    safe = str(rel).replace(os.sep, "_").replace(".", "_")
    png_path = viz_dir / f"energy_{key}_{safe}.png"
    try:
        plotter = EnergyDiagramPlotter(style=style) if style is not None else EnergyDiagramPlotter()
        plotter.plot(
            energies=values,
            labels=[f"F{i}" for i in range(len(values))],
        )
        plotter.save(str(png_path))
        sections.append({"kind": "energy", "rel": rel, "img": png_path,
                         "n_frames": len(values), "extra": key})
        viz_outputs[f"energy_{safe}"] = str(png_path)
        print(f"  ✓ energy[{key}] from {rel} → {png_path.name}")
    except Exception as exc:
        print(f"  ✗ energy[{key}] from {rel}: {exc}")

# ── 4. Render the HTML index ───────────────────────────────────
def _kind_label(k: str) -> str:
    return {"neb": "NEB Path", "mc": "MC Reconstruction",
            "trajectory": "Trajectory", "structure": "Structure",
            "energy": "Energy Diagram"}.get(k, k)


def _section_html(s: dict) -> str:
    img = Path(s["img"]).relative_to(viz_dir)
    nav_id = str(s["rel"]).replace(os.sep, "_").replace(".", "_")
    extra = f" — {s.get('extra')}" if s.get("extra") else ""
    n_frames = s.get("n_frames", 1)
    return (
        f'<div class="viz-section" id="sec_{nav_id}">'
        f'<h2>{_kind_label(s["kind"])} — {s["rel"]}{extra}</h2>'
        f'<p>Frames: {n_frames}</p>'
        f'<img src="{img}" alt="{s["rel"]}" loading="lazy">'
        f"</div>"
    )


nav_html = "".join(
    f'<li><a href="#sec_{str(s["rel"]).replace(os.sep, "_").replace(".", "_")}">'
    f'[{_kind_label(s["kind"])}] {s["rel"]}</a></li>'
    for s in sections
)

body = "\n".join(_section_html(s) for s in sections)

html_path = Path(OUTPUT_DIR) / "report.html"
html_path.write_text(
    f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>CatDT Aggregate Report — {COLOR_SCHEME}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
          max-width: 1200px; margin: 0 auto; padding: 20px; background: #fafafa; }}
  h1 {{ color: #1976D2; border-bottom: 3px solid #1976D2; padding-bottom: 10px; }}
  h2 {{ color: #388E3C; margin-top: 30px; }}
  nav {{ background: #fff; border: 1px solid #ddd; padding: 12px 20px; border-radius: 8px; }}
  nav ul {{ columns: 2; padding-left: 18px; }}
  nav a {{ color: #1976D2; text-decoration: none; }}
  nav a:hover {{ text-decoration: underline; }}
  .viz-section {{ background: #fff; padding: 20px; margin: 20px 0;
                  border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.08); }}
  img {{ max-width: 100%; border: 1px solid #eee; border-radius: 4px; }}
  .meta {{ color: #888; font-size: 0.9em; margin-bottom: 20px; }}
</style>
</head>
<body>
<h1>CatDT Aggregate Report</h1>
<p class="meta">
  Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} —
  scanned <code>{input_root}</code> —
  {len(trajectories)} trajectories, {len(structures)} structures,
  {len(energy_series)} energy series — colour scheme: <b>{COLOR_SCHEME}</b>.
</p>
<nav><h3>Contents</h3><ul>{nav_html}</ul></nav>
{body}
</body>
</html>
""",
    encoding="utf-8",
)

# Manifest for downstream consumption.
manifest_path = Path(OUTPUT_DIR) / "report_manifest.json"
manifest_path.write_text(json.dumps({
    "input_dir": str(input_root),
    "html": str(html_path),
    "viz_dir": str(viz_dir),
    "color_scheme": COLOR_SCHEME,
    "sections": [{"kind": s["kind"], "rel": str(s["rel"]),
                  "img": str(Path(s["img"]).relative_to(Path(OUTPUT_DIR)))}
                 for s in sections],
}, indent=2), encoding="utf-8")

print(f"\nReport: {html_path}")
print(f"Assets: {viz_dir}")
print(f"Manifest: {manifest_path}")
