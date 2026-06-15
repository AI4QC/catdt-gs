from __future__ import annotations

import json
import logging
import pickle
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from ase.io import read

from core.pathway.pathway_predictor import CompletePathwayResult, PathwayStep


logger = logging.getLogger(__name__)


_GAS_LABEL_RE = re.compile(r"^(?P<name>.+?)(?:\(g\)|_g)$")
_ITER_RE = re.compile(r"iter_(\d+)_energy_gate")


@dataclass
class WorkflowStepRecord:
    step_index: int
    step_name: str
    reactant_label: str
    product_label: str
    raw_reactant_energy: float
    raw_product_energy: float


def _is_gas_label(label: str) -> bool:
    return _GAS_LABEL_RE.match(str(label).strip()) is not None


def _strip_gas_suffix(label: str) -> str:
    match = _GAS_LABEL_RE.match(str(label).strip())
    return match.group("name") if match else str(label).strip()


def _normalize_surface_state_label(label: str) -> str:
    label = str(label).strip()
    if not label:
        return "*"
    if label == "*":
        return "*"
    if label.startswith("*"):
        return label
    return f"*{label}"


def _normalize_step_labels(reactant_raw: str, product_raw: str) -> Tuple[str, str]:
    reactant_raw = str(reactant_raw).strip()
    product_raw = str(product_raw).strip()

    reactant_is_gas = _is_gas_label(reactant_raw)
    product_is_gas = _is_gas_label(product_raw)

    if reactant_is_gas and not product_is_gas:
        return "*", _normalize_surface_state_label(product_raw)
    if product_is_gas and not reactant_is_gas:
        return _normalize_surface_state_label(reactant_raw), "*"

    return (
        _normalize_surface_state_label(_strip_gas_suffix(reactant_raw)),
        _normalize_surface_state_label(_strip_gas_suffix(product_raw)),
    )


def _iteration_sort_key(path: Path) -> Tuple[int, str]:
    match = _ITER_RE.search(path.parent.name)
    return (int(match.group(1)) if match else -1, str(path))


def _find_latest_energy_gate(run_dir: Path) -> Path:
    candidates = sorted(
        run_dir.glob("04_pathway/iter_*_energy_gate/agent45_energy_gate.json"),
        key=_iteration_sort_key,
    )
    if not candidates:
        raise FileNotFoundError(f"No agent45 energy gate JSON found under {run_dir}")
    return candidates[-1]


def _load_step_records(energy_gate_path: Path) -> List[WorkflowStepRecord]:
    payload = json.loads(energy_gate_path.read_text(encoding="utf-8"))
    steps = payload.get("steps", [])
    if not steps:
        raise ValueError(f"No steps found in {energy_gate_path}")

    records: List[WorkflowStepRecord] = []
    for ordinal, item in enumerate(steps):
        step_name = str(item.get("step_name", "")).strip()
        if "_to_" not in step_name:
            raise ValueError(f"Unsupported step name format: {step_name}")
        reactant_raw, product_raw = step_name.split("_to_", 1)
        reactant_label, product_label = _normalize_step_labels(reactant_raw, product_raw)
        endpoints = item.get("endpoints", {})
        reactant_energy = float(endpoints.get("reactant", {}).get("energy"))
        product_energy = float(endpoints.get("product", {}).get("energy"))

        records.append(
            WorkflowStepRecord(
                step_index=int(item.get("step_index", ordinal + 1)) - 1,
                step_name=step_name,
                reactant_label=reactant_label,
                product_label=product_label,
                raw_reactant_energy=reactant_energy,
                raw_product_energy=product_energy,
            )
        )
    return records


def _infer_surface_formula(run_dir: Path) -> Tuple[str, str]:
    candidates = [
        run_dir / "04_pathway" / "agent45_only_prepare" / "clean_surface.vasp",
        run_dir / "04_pathway" / "iter_01_tool_preplacement" / "clean_surface_known_indices.vasp",
    ]
    candidates.extend(sorted(run_dir.glob("04_pathway/**/*clean_surface*.vasp")))

    for candidate in candidates:
        if candidate.exists():
            atoms = read(str(candidate))
            return atoms.get_chemical_formula(), str(candidate)

    fallback_candidates = list(run_dir.glob("05_pre_neb_relax/*_reactant_relaxed.vasp"))
    fallback_candidates.extend(run_dir.glob("05_pre_neb_relax/*_product_relaxed.vasp"))
    if fallback_candidates:
        atoms = read(str(fallback_candidates[0]))
        return atoms.get_chemical_formula(), str(fallback_candidates[0])

    return run_dir.name, ""


def _extract_transition_state_from_traj(traj_path: Path) -> Dict[str, Any]:
    frames = read(str(traj_path), index=":")
    available: List[Tuple[int, float]] = []
    for frame_index, atoms in enumerate(frames):
        calc = getattr(atoms, "calc", None)
        if calc is None:
            continue
        results = getattr(calc, "results", {}) or {}
        energy = results.get("energy")
        if energy is None:
            continue
        available.append((frame_index, float(energy)))

    if not available:
        raise ValueError(f"No frame energies found in {traj_path}")

    interior = [(idx, energy) for idx, energy in available if 0 < idx < len(frames) - 1]
    source = interior if interior else available
    ts_index, ts_energy = max(source, key=lambda item: item[1])

    return {
        "traj_path": str(traj_path),
        "num_frames": len(frames),
        "num_energy_frames": len(available),
        "used_interior_only": bool(interior),
        "transition_state_index": ts_index,
        "transition_state_energy_raw": ts_energy,
    }


def _load_neb_endpoint_structures(traj_path: Path) -> Tuple[Any, Any]:
    frames = read(str(traj_path), index=":")
    if not frames:
        raise ValueError(f"No frames found in {traj_path}")
    return frames[0], frames[-1]


def _find_neb_traj(neb_dir: Path, step_name: str) -> Path:
    candidate = neb_dir / f"{step_name}_neb.traj"
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"NEB trajectory not found for step {step_name}: {candidate}")


def _build_adsorbate_energy_map(sequence_labels: Iterable[str], sequence_energies: Iterable[float]) -> Dict[str, float]:
    energies: Dict[str, float] = {}
    for label, energy in zip(sequence_labels, sequence_energies):
        energies[str(label)] = float(energy)
    return energies


def _build_adsorption_energy_map(adsorbate_energies: Dict[str, float]) -> Dict[str, float]:
    star_energy = adsorbate_energies.get("*")
    if star_energy is None:
        return {label: 0.0 for label in adsorbate_energies}
    return {label: float(energy - star_energy) for label, energy in adsorbate_energies.items()}


def build_complete_pathway_result_from_workflow(
    run_dir: str,
    *,
    output_dir: Optional[str] = None,
    endpoint_energy_getter: Optional[Callable[[str, str, Any], float]] = None,
    energy_mode: str = "absolute-continuity",
) -> Tuple[CompletePathwayResult, Dict[str, Any]]:
    run_path = Path(run_dir).resolve()
    energy_gate_path = _find_latest_energy_gate(run_path)
    neb_dir = run_path / "06_neb"
    records = _load_step_records(energy_gate_path)
    surface_formula, surface_formula_source = _infer_surface_formula(run_path)

    step_diagnostics: List[Dict[str, Any]] = []
    raw_step_entries: List[Dict[str, Any]] = []

    for index, record in enumerate(records):
        traj_path = _find_neb_traj(neb_dir, record.step_name)
        reactant_atoms, product_atoms = _load_neb_endpoint_structures(traj_path)

        if endpoint_energy_getter is not None:
            raw_reactant_energy = float(endpoint_energy_getter(record.step_name, "reactant", reactant_atoms))
            raw_product_energy = float(endpoint_energy_getter(record.step_name, "product", product_atoms))
            energy_source = "neb_endpoint_single_point"
        else:
            raw_reactant_energy = record.raw_reactant_energy
            raw_product_energy = record.raw_product_energy
            energy_source = "energy_gate"

        ts_info = _extract_transition_state_from_traj(traj_path)
        raw_ts_energy = float(ts_info["transition_state_energy_raw"])
        raw_barrier = max(0.0, raw_ts_energy - record.raw_reactant_energy)
        raw_step_entries.append(
            {
                "step_index": index,
                "step_name": record.step_name,
                "reactant_label": record.reactant_label,
                "product_label": record.product_label,
                "energy_source": energy_source,
                "energy_gate_reactant_energy": record.raw_reactant_energy,
                "energy_gate_product_energy": record.raw_product_energy,
                "raw_reactant_energy": raw_reactant_energy,
                "raw_product_energy": raw_product_energy,
                "step_free_energy_delta": raw_product_energy - raw_reactant_energy,
                "traj_path": ts_info["traj_path"],
                "transition_state_index": ts_info["transition_state_index"],
                "transition_state_energy_raw": raw_ts_energy,
                "activation_energy_forward": raw_barrier,
                "num_frames": ts_info["num_frames"],
                "num_energy_frames": ts_info["num_energy_frames"],
                "used_interior_only": ts_info["used_interior_only"],
            }
        )

    pathway_steps: List[PathwayStep] = []
    sequence_labels: List[str] = []
    sequence_energies: List[float] = []

    if energy_mode == "step-delta-relative":
        relative_product_energies: List[float] = [0.0] * len(raw_step_entries)
        relative_reactant_energies: List[float] = [0.0] * len(raw_step_entries)

        if raw_step_entries:
            current_product_energy = 0.0
            for rev_index in range(len(raw_step_entries) - 1, -1, -1):
                entry = raw_step_entries[rev_index]
                relative_product_energies[rev_index] = current_product_energy
                relative_reactant_energies[rev_index] = current_product_energy - float(entry["step_free_energy_delta"])
                current_product_energy = relative_reactant_energies[rev_index]

        for index, entry in enumerate(raw_step_entries):
            step_reactant_energy = relative_reactant_energies[index]
            step_product_energy = relative_product_energies[index]
            transition_state_energy = step_reactant_energy + float(entry["activation_energy_forward"])

            pathway_steps.append(
                PathwayStep(
                    step_index=index,
                    name=str(entry["step_name"]),
                    reactant_adsorbate=str(entry["reactant_label"]),
                    product_adsorbate=str(entry["product_label"]),
                    reactant_energy=step_reactant_energy,
                    product_energy=step_product_energy,
                    reaction_energy=float(entry["step_free_energy_delta"]),
                    activation_energy=float(entry["activation_energy_forward"]),
                    transition_state_energy=transition_state_energy,
                )
            )

            if not sequence_labels:
                sequence_labels.append(str(entry["reactant_label"]))
                sequence_energies.append(step_reactant_energy)
            sequence_labels.append(str(entry["product_label"]))
            sequence_energies.append(step_product_energy)

            entry["continuity_reactant_energy"] = step_reactant_energy
            entry["reactant_continuity_delta"] = 0.0

    else:
        for index, entry in enumerate(raw_step_entries):
            continuity_delta = 0.0
            step_reactant_energy = float(entry["raw_reactant_energy"])
            if index > 0 and raw_step_entries[index - 1]["product_label"] == entry["reactant_label"]:
                previous_product_energy = pathway_steps[index - 1].product_energy
                continuity_delta = abs(float(entry["raw_reactant_energy"]) - previous_product_energy)
                step_reactant_energy = previous_product_energy

            step_product_energy = float(entry["raw_product_energy"])
            transition_state_energy = step_reactant_energy + float(entry["activation_energy_forward"])

            pathway_steps.append(
                PathwayStep(
                    step_index=index,
                    name=str(entry["step_name"]),
                    reactant_adsorbate=str(entry["reactant_label"]),
                    product_adsorbate=str(entry["product_label"]),
                    reactant_energy=step_reactant_energy,
                    product_energy=step_product_energy,
                    reaction_energy=step_product_energy - step_reactant_energy,
                    activation_energy=float(entry["activation_energy_forward"]),
                    transition_state_energy=transition_state_energy,
                )
            )

            if not sequence_labels:
                sequence_labels.append(str(entry["reactant_label"]))
                sequence_energies.append(step_reactant_energy)
            sequence_labels.append(str(entry["product_label"]))
            sequence_energies.append(step_product_energy)

            entry["continuity_reactant_energy"] = step_reactant_energy
            entry["reactant_continuity_delta"] = continuity_delta

    step_diagnostics.extend(raw_step_entries)

    adsorbate_energies = _build_adsorbate_energy_map(sequence_labels, sequence_energies)
    adsorption_energies = _build_adsorption_energy_map(adsorbate_energies)

    rate_determining_step = max(pathway_steps, key=lambda step: step.activation_energy or float("-inf"), default=None)
    max_barrier = rate_determining_step.activation_energy if rate_determining_step is not None else None
    if rate_determining_step is not None:
        rate_determining_step.is_rate_determining = True

    result = CompletePathwayResult(
        surface_formula=surface_formula,
        adsorbates=sequence_labels,
        adsorbate_energies=adsorbate_energies,
        adsorption_energies=adsorption_energies,
        steps=pathway_steps,
        rate_determining_step=rate_determining_step,
        max_barrier=max_barrier,
        overall_reaction_energy=sequence_energies[-1] - sequence_energies[0],
        output_dir=output_dir or str(run_path / "07_posthoc_agent6_from_workflow"),
    )

    diagnostics = {
        "run_dir": str(run_path),
        "energy_gate_path": str(energy_gate_path),
        "neb_dir": str(neb_dir),
        "surface_formula": surface_formula,
        "surface_formula_source": surface_formula_source,
        "adsorbate_sequence": sequence_labels,
        "energy_mode": energy_mode,
        "step_diagnostics": step_diagnostics,
    }
    return result, diagnostics


def save_reconstructed_agent6_inputs(
    output_dir: str,
    pathway_result: CompletePathwayResult,
    diagnostics: Dict[str, Any],
    kmc_result: Optional[Any] = None,
) -> Dict[str, str]:
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    pathway_path = out_dir / "complete_pathway_result.pkl"
    with pathway_path.open("wb") as fh:
        pickle.dump(pathway_result, fh)

    summary = {
        **diagnostics,
        "surface_formula": pathway_result.surface_formula,
        "adsorbates": pathway_result.adsorbates,
        "adsorbate_energies": pathway_result.adsorbate_energies,
        "max_barrier": pathway_result.max_barrier,
        "overall_reaction_energy": pathway_result.overall_reaction_energy,
    }
    summary_path = out_dir / "workflow_agent6_input_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    saved = {
        "pathway_result": str(pathway_path),
        "summary": str(summary_path),
    }

    if kmc_result is not None:
        kmc_path = out_dir / "catmap_result.pkl"
        with kmc_path.open("wb") as fh:
            pickle.dump(kmc_result, fh)
        saved["kmc_result"] = str(kmc_path)

    return saved


def apply_barrier_overrides(
    pathway_result: CompletePathwayResult,
    diagnostics: Dict[str, Any],
    barrier_overrides: Dict[str, float],
) -> None:
    if not barrier_overrides:
        return

    for step in pathway_result.steps:
        if step.name in barrier_overrides:
            new_barrier = float(barrier_overrides[step.name])
            step.activation_energy = new_barrier
            step.transition_state_energy = float(step.reactant_energy + new_barrier)
            step.is_rate_determining = False

    step_diag_map = {
        str(item.get("step_name")): item
        for item in diagnostics.get("step_diagnostics", [])
    }
    for step_name, new_barrier in barrier_overrides.items():
        if step_name in step_diag_map:
            step_diag_map[step_name]["activation_energy_forward"] = float(new_barrier)
            continuity_reactant_energy = float(step_diag_map[step_name].get("continuity_reactant_energy", 0.0))
            step_diag_map[step_name]["transition_state_energy_overridden"] = continuity_reactant_energy + float(new_barrier)

    valid_barriers = [step for step in pathway_result.steps if step.activation_energy is not None]
    if valid_barriers:
        pathway_result.rate_determining_step = max(valid_barriers, key=lambda step: step.activation_energy)
        pathway_result.rate_determining_step.is_rate_determining = True
        pathway_result.max_barrier = float(pathway_result.rate_determining_step.activation_energy)
    else:
        pathway_result.rate_determining_step = None
        pathway_result.max_barrier = None

    diagnostics["barrier_overrides"] = {name: float(value) for name, value in barrier_overrides.items()}


def run_catmap_single_point_from_pathway(
    pathway_result: CompletePathwayResult,
    *,
    temperature: float,
    pressures: Dict[str, float],
    output_dir: Optional[str] = None,
    voltage: Optional[float] = None,
    pH: float = 0.0,
    predictor_factory: Optional[Callable[..., Any]] = None,
    mechanism_builder: Optional[Callable[[Any, CompletePathwayResult, Dict[str, float]], None]] = None,
) -> Any:
    work_dir = str(Path(output_dir).resolve()) if output_dir else tempfile.mkdtemp(prefix="catdt-agent6-catmap-")

    if predictor_factory is None:
        from core.kmc.catmap_predictor import CatMAPPredictor

        predictor_factory = CatMAPPredictor

    if mechanism_builder is None:
        from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin

        dt = GasSolidDigitalTwin.__new__(GasSolidDigitalTwin)
        dt.logger = logger
        mechanism_builder = dt._build_catmap_mechanism

    predictor = predictor_factory(work_dir=work_dir)
    mechanism_builder(predictor, pathway_result, dict(pressures))

    gas_names = list(getattr(predictor, "gases", {}).keys())
    if any(len(str(name)) == 1 for name in gas_names):
        predictor.set_thermo_modes(gas_mode="frozen_gas", adsorbate_mode="harmonic_adsorbate")

    predictor.set_conditions(
        temperature=temperature,
        pressures=dict(pressures),
        voltage=voltage,
        pH=float(pH),
    )
    return predictor.run_single_point(surface=pathway_result.surface_formula, output_dir=work_dir)
