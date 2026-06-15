from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, Iterable, Optional

from core.pathway.pathway_predictor import CompletePathwayResult


def _rebuild_adsorbate_maps(pathway_result: CompletePathwayResult) -> None:
    sequence_labels = []
    sequence_energies = []

    for step in pathway_result.steps:
        if not sequence_labels:
            sequence_labels.append(step.reactant_adsorbate)
            sequence_energies.append(float(step.reactant_energy))
        sequence_labels.append(step.product_adsorbate)
        sequence_energies.append(float(step.product_energy))

    adsorbate_energies = {
        str(label): float(energy)
        for label, energy in zip(sequence_labels, sequence_energies)
    }
    star_energy = adsorbate_energies.get("*", 0.0)
    adsorption_energies = {
        label: float(energy - star_energy)
        for label, energy in adsorbate_energies.items()
    }

    pathway_result.adsorbates = sequence_labels
    pathway_result.adsorbate_energies = adsorbate_energies
    pathway_result.adsorption_energies = adsorption_energies
    pathway_result.overall_reaction_energy = float(sequence_energies[-1] - sequence_energies[0])


def _update_rds(pathway_result: CompletePathwayResult) -> None:
    for step in pathway_result.steps:
        step.is_rate_determining = False

    valid = [step for step in pathway_result.steps if step.activation_energy is not None]
    if not valid:
        pathway_result.rate_determining_step = None
        pathway_result.max_barrier = None
        return

    pathway_result.rate_determining_step = max(valid, key=lambda step: float(step.activation_energy))
    pathway_result.rate_determining_step.is_rate_determining = True
    pathway_result.max_barrier = float(pathway_result.rate_determining_step.activation_energy)


def rebuild_pathway_with_step_parameter_overrides(
    base_pathway_result: CompletePathwayResult,
    *,
    reaction_energy_overrides: Optional[Dict[str, float]] = None,
    activation_energy_overrides: Optional[Dict[str, float]] = None,
) -> CompletePathwayResult:
    reaction_energy_overrides = reaction_energy_overrides or {}
    activation_energy_overrides = activation_energy_overrides or {}

    rebuilt = deepcopy(base_pathway_result)
    if not rebuilt.steps:
        return rebuilt

    current_product_energy = float(base_pathway_result.steps[-1].product_energy)

    for step in reversed(rebuilt.steps):
        reaction_energy = float(reaction_energy_overrides.get(step.name, step.reaction_energy))
        activation_energy = activation_energy_overrides.get(step.name, step.activation_energy)
        activation_energy = None if activation_energy is None else float(activation_energy)

        step.product_energy = current_product_energy
        step.reaction_energy = reaction_energy
        step.reactant_energy = current_product_energy - reaction_energy
        step.activation_energy = activation_energy
        step.transition_state_energy = (
            None if activation_energy is None else float(step.reactant_energy + activation_energy)
        )
        current_product_energy = float(step.reactant_energy)

    _rebuild_adsorbate_maps(rebuilt)
    _update_rds(rebuilt)
    return rebuilt


def find_smallest_shift_meeting_target(
    records: Iterable[Dict[str, Any]],
    *,
    target_tof: float,
    tof_key: str = "ch4_tof",
) -> Optional[Dict[str, Any]]:
    candidates = [
        record
        for record in records
        if record.get(tof_key) is not None and float(record[tof_key]) >= float(target_tof)
    ]
    if not candidates:
        return None

    return min(
        candidates,
        key=lambda record: (
            abs(float(record.get("delta_from_base", 0.0))),
            abs(float(record.get("candidate_value", 0.0))),
        ),
    )
