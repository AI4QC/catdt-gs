import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from ase import Atom, Atoms

from camel_agents.schemas import AtomAddSpec, PathwayStepSpec

logger = logging.getLogger(__name__)

DEFAULT_RENDERER = "tachyon"
DEFAULT_CLEARANCE = 1.1


def get_default_her_description() -> str:
    return (
        "Hydrogen evolution style pathway where adsorbed hydrogen species can appear sequentially "
        "and a second hydrogen addition may be required before molecular release."
    )


def get_default_nrr_description() -> str:
    return (
        "Nitrogen reduction style pathway where activated nitrogen receives multiple hydrogen equivalents "
        "and can ultimately form two NH3 products."
    )


def align_added_atoms_to_anchor(positions: List[List[float]], anchor_position: List[float]) -> List[List[float]]:
    if not positions:
        return []

    anchor = [float(anchor_position[0]), float(anchor_position[1]), float(anchor_position[2])]
    base = positions[0]
    aligned: List[List[float]] = []
    for pos in positions:
        dx = float(pos[0]) - float(base[0])
        dy = float(pos[1]) - float(base[1])
        dz = float(pos[2]) - float(base[2])
        aligned.append([anchor[0] + dx, anchor[1] + dy, float(base[2]) + dz])

    # First atom sits exactly above anchor in x/y, retains original z
    aligned[0][2] = float(positions[0][2])
    return aligned


def select_anchor_position(
    slab: Atoms,
    surface_indices: Optional[List[int]] = None,
    anchor_index: Optional[int] = None,
) -> List[float]:
    if len(slab) == 0:
        return [0.0, 0.0, 0.0]

    if surface_indices is None:
        surface_indices = list(range(len(slab)))

    if anchor_index is not None and 0 <= int(anchor_index) < len(slab):
        return slab.positions[int(anchor_index)].tolist()

    if not surface_indices:
        return slab.positions[0].tolist()

    top_idx = max(surface_indices, key=lambda idx: slab.positions[idx][2])
    return slab.positions[top_idx].tolist()


def enforce_minimum_height(
    positions: List[List[float]],
    surface_max_z: float,
    clearance: float = DEFAULT_CLEARANCE,
) -> List[List[float]]:
    min_allowed = float(surface_max_z) + float(clearance)
    adjusted: List[List[float]] = []
    for pos in positions:
        x, y, z = float(pos[0]), float(pos[1]), float(pos[2])
        if z < min_allowed:
            z = min_allowed
        adjusted.append([x, y, z])
    return adjusted


def _resolve_remove_indices(
    atoms_to_remove: List[Any],
    adsorbate_indices: List[int],
    use_adsorbate_local_indices: bool,
) -> List[int]:
    resolved: List[int] = []
    for raw in atoms_to_remove:
        idx_raw = raw.get("index") if isinstance(raw, dict) else raw
        try:
            idx = int(idx_raw)
        except Exception:
            continue

        if use_adsorbate_local_indices:
            if 0 <= idx < len(adsorbate_indices):
                resolved.append(adsorbate_indices[idx])
            elif idx in adsorbate_indices:
                resolved.append(idx)
        else:
            if idx >= 0:
                resolved.append(idx)

    return sorted(set(resolved))


def apply_step_modification(
    structure: Atoms,
    step: PathwayStepSpec,
    use_adsorbate_local_indices: bool = True,
) -> Tuple[Atoms, List[int]]:
    product = structure.copy()
    adsorbate_indices = list(product.info.get("adsorbate_indices", []))

    # Add atoms
    for add_spec in step.atoms_to_add:
        symbol = (add_spec.species or "").strip()[:2]
        if not symbol:
            continue
        symbol = symbol[0].upper() + (symbol[1].lower() if len(symbol) > 1 else "")
        if add_spec.position and len(add_spec.position) == 3:
            position = [float(add_spec.position[0]), float(add_spec.position[1]), float(add_spec.position[2])]
        else:
            position = [0.0, 0.0, float(product.positions[:, 2].max() + 1.5)]
        product.append(Atom(symbol, position=position))
        adsorbate_indices.append(len(product) - 1)

    # Remove atoms (protect surface-only structures in local-index mode)
    to_remove = _resolve_remove_indices(step.atoms_to_remove, adsorbate_indices, use_adsorbate_local_indices)
    if use_adsorbate_local_indices and not adsorbate_indices:
        to_remove = []

    if to_remove:
        remove_set = set(i for i in to_remove if 0 <= i < len(product))
        keep_mask = [i not in remove_set for i in range(len(product))]
        index_map: Dict[int, int] = {}
        next_idx = 0
        for old_idx, keep in enumerate(keep_mask):
            if keep:
                index_map[old_idx] = next_idx
                next_idx += 1
        product = product[keep_mask]
        adsorbate_indices = [index_map[i] for i in adsorbate_indices if i in index_map]

    ads_set = set(adsorbate_indices)
    product.info["adsorbate_indices"] = adsorbate_indices
    product.info["surface_indices"] = [i for i in range(len(product)) if i not in ads_set]

    return product, adsorbate_indices


def build_energy_diagram_inputs(
    intermediates: List[str],
    energies: List[float],
    barriers: List[float],
) -> Tuple[List[str], List[float], Dict[str, float]]:
    labels = list(intermediates)
    energy_list = [float(e) for e in energies]
    barrier_map: Dict[str, float] = {}

    for i in range(min(len(barriers), max(len(intermediates) - 1, 0))):
        key = f"{intermediates[i]} -> {intermediates[i + 1]}"
        barrier_map[key] = float(barriers[i])

    return labels, energy_list, barrier_map


@dataclass
class PathwayRunState:
    reaction_description: str
    anchor_position: List[float]
    surface_indices: List[int]
    adsorbate_indices: List[int]


class PathwayAgentsRunner:
    """Compatibility runner placeholder.

    Full pathway orchestration is now handled by `CatDTCamelWorkflow`.
    """

    def __init__(self, *args: Any, **kwargs: Any):
        _ = args, kwargs

    def run(self, state: PathwayRunState, base_structure: Atoms) -> PathwayRunState:
        _ = base_structure
        return state
