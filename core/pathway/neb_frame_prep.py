"""Shared endpoint preparation for interpolation-sensitive NEB workflows."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii


@dataclass
class PreparedEndpointPair:
    reactant: Atoms
    product: Atoms
    reactant_adsorbate_indices: List[int]
    product_adsorbate_indices: List[int]
    reordered: bool


def reorder_product_to_match_reactant(
    reactant: Atoms,
    product: Atoms,
    product_adsorbate_indices: Sequence[int] | None = None,
) -> tuple[Atoms, List[int], bool]:
    """Reorder product atoms to match reactant atom order when symbols match."""
    product_adsorbate_indices = list(product_adsorbate_indices or [])
    if len(reactant) != len(product):
        return product.copy(), product_adsorbate_indices, False

    reactant_symbols = reactant.get_chemical_symbols()
    product_symbols = product.get_chemical_symbols()
    if sorted(reactant_symbols) != sorted(product_symbols):
        return product.copy(), product_adsorbate_indices, False
    if reactant_symbols == product_symbols:
        copied = product.copy()
        copied.info.setdefault("adsorbate_indices", list(product_adsorbate_indices))
        return copied, list(product_adsorbate_indices), False

    reactant_positions = reactant.get_positions()
    product_positions = product.get_positions()

    try:
        from scipy.optimize import linear_sum_assignment  # type: ignore
    except Exception:
        linear_sum_assignment = None

    assignment: Dict[int, int] = {}
    grouped: Dict[str, Tuple[List[int], List[int]]] = {}
    for idx, sym in enumerate(reactant_symbols):
        grouped.setdefault(sym, ([], []))[0].append(idx)
    for idx, sym in enumerate(product_symbols):
        grouped.setdefault(sym, ([], []))[1].append(idx)

    for sym, (reactant_group, product_group) in grouped.items():
        if len(reactant_group) != len(product_group):
            return product.copy(), product_adsorbate_indices, False
        if not reactant_group:
            continue

        if linear_sum_assignment is not None:
            cost = np.zeros((len(reactant_group), len(product_group)), dtype=float)
            for i, reactant_idx in enumerate(reactant_group):
                for j, product_idx in enumerate(product_group):
                    cost[i, j] = float(
                        np.linalg.norm(reactant_positions[reactant_idx] - product_positions[product_idx])
                    )
            row_ind, col_ind = linear_sum_assignment(cost)
            for row, col in zip(row_ind.tolist(), col_ind.tolist()):
                assignment[reactant_group[row]] = product_group[col]
        else:
            remaining = list(product_group)
            for reactant_idx in reactant_group:
                chosen = min(
                    remaining,
                    key=lambda product_idx: float(
                        np.linalg.norm(product_positions[product_idx] - reactant_positions[reactant_idx])
                    ),
                )
                assignment[reactant_idx] = chosen
                remaining.remove(chosen)

    try:
        reorder_sequence = [assignment[idx] for idx in range(len(reactant_symbols))]
    except Exception:
        return product.copy(), product_adsorbate_indices, False

    if reorder_sequence == list(range(len(reactant_symbols))):
        copied = product.copy()
        copied.info.setdefault("adsorbate_indices", list(product_adsorbate_indices))
        return copied, list(product_adsorbate_indices), False

    reordered_product = product[reorder_sequence]
    old_to_new = {old: new for new, old in enumerate(reorder_sequence)}
    new_ads = sorted(old_to_new[idx] for idx in product_adsorbate_indices if idx in old_to_new)
    new_ads_set = set(new_ads)
    reordered_product.info["adsorbate_indices"] = list(new_ads)
    reordered_product.info["surface_indices"] = [
        idx for idx in range(len(reordered_product)) if idx not in new_ads_set
    ]
    return reordered_product, new_ads, True


def prepare_endpoint_pair_for_interpolation(
    reactant: Atoms,
    product: Atoms,
    reactant_adsorbate_indices: Sequence[int] | None = None,
    product_adsorbate_indices: Sequence[int] | None = None,
    reorder_product_atoms: bool = True,
) -> PreparedEndpointPair:
    """Prepare endpoints so interpolation follows the nearest periodic images per fragment."""
    prepared_reactant = reactant.copy()
    prepared_product = product.copy()

    reactant_ads = _valid_indices(
        prepared_reactant,
        reactant_adsorbate_indices or prepared_reactant.info.get("adsorbate_indices", []),
    )
    product_ads = _valid_indices(
        prepared_product,
        product_adsorbate_indices or prepared_product.info.get("adsorbate_indices", []),
    )

    reordered = False
    if reorder_product_atoms:
        prepared_product, product_ads, reordered = reorder_product_to_match_reactant(
            prepared_reactant,
            prepared_product,
            product_adsorbate_indices=product_ads,
        )

    reactant_graph = _build_adsorbate_graph(prepared_reactant, reactant_ads)
    product_graph = _build_adsorbate_graph(prepared_product, product_ads)

    reactant_components = _connected_components(reactant_ads, reactant_graph)
    product_components = _connected_components(product_ads, product_graph)

    _unwrap_components_inplace(prepared_reactant, reactant_components, reactant_graph)
    _unwrap_components_inplace(prepared_product, product_components, product_graph)

    product_symbols = prepared_product.get_chemical_symbols()
    for component in product_components:
        anchor_idx = _choose_component_anchor(component, product_symbols)
        nearest_anchor = _nearest_image_position(
            reference=prepared_reactant.positions[anchor_idx],
            moving=prepared_product.positions[anchor_idx],
            cell=np.array(prepared_product.cell, dtype=float),
            pbc=np.array(prepared_product.pbc, dtype=bool),
        )
        translation = nearest_anchor - prepared_product.positions[anchor_idx]
        if np.linalg.norm(translation) > 0:
            for atom_idx in component:
                prepared_product.positions[atom_idx] = prepared_product.positions[atom_idx] + translation

    _unwrap_components_inplace(prepared_product, product_components, product_graph)

    return PreparedEndpointPair(
        reactant=prepared_reactant,
        product=prepared_product,
        reactant_adsorbate_indices=reactant_ads,
        product_adsorbate_indices=product_ads,
        reordered=reordered,
    )


def interpolated_path_min_distances(
    reactant: Atoms,
    product: Atoms,
    adsorbate_indices: Sequence[int] | None = None,
    n_images: int = 7,
) -> tuple[float, int, float, int]:
    """Compute simple linear-interpolation distance diagnostics on prepared endpoints."""
    if len(reactant) != len(product):
        return -1.0, -1, -1.0, -1

    adsorbate_indices = _valid_indices(reactant, adsorbate_indices or [])
    frame_count = max(3, int(n_images))
    reactant_pos = np.array(reactant.get_positions(), dtype=float)
    product_pos = np.array(product.get_positions(), dtype=float)
    frame_atoms = reactant.copy()

    min_any_pair = float("inf")
    min_any_pair_frame = -1
    min_ads_surface = float("inf")
    min_ads_surface_frame = -1

    for frame_idx in range(frame_count):
        alpha = frame_idx / float(frame_count - 1)
        interp_pos = (1.0 - alpha) * reactant_pos + alpha * product_pos
        frame_atoms.set_positions(interp_pos)

        min_any = _min_distance_any_pair(frame_atoms)
        if 0 < min_any < min_any_pair:
            min_any_pair = min_any
            min_any_pair_frame = frame_idx

        min_ads_surf = _min_adsorbate_surface_distance(frame_atoms, adsorbate_indices)
        if 0 < min_ads_surf < min_ads_surface:
            min_ads_surface = min_ads_surf
            min_ads_surface_frame = frame_idx

    if min_any_pair == float("inf"):
        min_any_pair = -1.0
    if min_ads_surface == float("inf"):
        min_ads_surface = -1.0

    return float(min_any_pair), int(min_any_pair_frame), float(min_ads_surface), int(min_ads_surface_frame)


def _valid_indices(atoms: Atoms, indices: Iterable[int]) -> List[int]:
    return sorted({int(idx) for idx in indices if 0 <= int(idx) < len(atoms)})


def _build_adsorbate_graph(atoms: Atoms, adsorbate_indices: Sequence[int]) -> Dict[int, set[int]]:
    graph: Dict[int, set[int]] = {idx: set() for idx in adsorbate_indices}
    if len(adsorbate_indices) < 2:
        return graph

    symbols = atoms.get_chemical_symbols()
    for pos, idx_i in enumerate(adsorbate_indices):
        for idx_j in adsorbate_indices[pos + 1 :]:
            cutoff = _bond_cutoff(symbols[idx_i], symbols[idx_j])
            dist = float(atoms.get_distance(idx_i, idx_j, mic=True))
            if dist <= cutoff:
                graph[idx_i].add(idx_j)
                graph[idx_j].add(idx_i)
    return graph


def _bond_cutoff(symbol_a: str, symbol_b: str) -> float:
    try:
        radius_sum = covalent_radii[atomic_numbers[symbol_a]] + covalent_radii[atomic_numbers[symbol_b]]
    except Exception:
        radius_sum = 1.6
    return float(min(2.60, max(0.80, 1.30 * radius_sum)))


def _connected_components(indices: Sequence[int], graph: Dict[int, set[int]]) -> List[List[int]]:
    remaining = set(indices)
    components: List[List[int]] = []
    while remaining:
        start = min(remaining)
        queue = deque([start])
        seen = {start}
        remaining.remove(start)
        component = [start]
        while queue:
            current = queue.popleft()
            for neighbor in graph.get(current, set()):
                if neighbor in seen:
                    continue
                seen.add(neighbor)
                if neighbor in remaining:
                    remaining.remove(neighbor)
                queue.append(neighbor)
                component.append(neighbor)
        components.append(sorted(component))
    return components


def _unwrap_components_inplace(
    atoms: Atoms,
    components: Sequence[Sequence[int]],
    graph: Dict[int, set[int]],
) -> None:
    if not components or not np.any(atoms.pbc):
        return
    positions = np.array(atoms.get_positions(), dtype=float)
    cell = np.array(atoms.cell, dtype=float)
    if cell.shape != (3, 3):
        return
    try:
        inv_t = np.linalg.pinv(cell.T)
    except Exception:
        return

    symbols = atoms.get_chemical_symbols()
    for component in components:
        if len(component) < 2:
            continue
        anchor_idx = _choose_component_anchor(component, symbols)
        unwrapped: Dict[int, np.ndarray] = {anchor_idx: np.array(positions[anchor_idx], dtype=float)}
        queue = deque([anchor_idx])
        while queue:
            current = queue.popleft()
            current_frac = inv_t @ np.array(unwrapped[current], dtype=float)
            for neighbor in sorted(graph.get(current, set())):
                if neighbor in unwrapped:
                    continue
                neighbor_frac = inv_t @ np.array(positions[neighbor], dtype=float)
                delta = neighbor_frac - (inv_t @ np.array(positions[current], dtype=float))
                for dim in range(3):
                    if bool(atoms.pbc[dim]):
                        delta[dim] = float(delta[dim] - np.round(delta[dim]))
                unwrapped[neighbor] = cell.T @ (current_frac + delta)
                queue.append(neighbor)
        for atom_idx, cart_pos in unwrapped.items():
            positions[atom_idx] = cart_pos
    atoms.set_positions(positions)


def _choose_component_anchor(component: Sequence[int], symbols: Sequence[str]) -> int:
    non_h = [idx for idx in component if symbols[idx] != "H"]
    if non_h:
        return min(non_h)
    return min(component)


def _nearest_image_position(
    reference: np.ndarray,
    moving: np.ndarray,
    cell: np.ndarray,
    pbc: np.ndarray,
) -> np.ndarray:
    if cell.shape != (3, 3):
        return np.array(moving, dtype=float)
    try:
        inv_t = np.linalg.pinv(cell.T)
    except Exception:
        return np.array(moving, dtype=float)

    ref_frac = inv_t @ np.array(reference, dtype=float)
    moving_frac = inv_t @ np.array(moving, dtype=float)
    delta = ref_frac - moving_frac
    shift = np.zeros(3, dtype=float)
    for dim in range(3):
        if bool(pbc[dim]):
            shift[dim] = float(np.round(delta[dim]))
    return np.array(moving, dtype=float) + cell.T @ shift


def _min_distance_any_pair(structure: Atoms) -> float:
    if len(structure) < 2:
        return -1.0
    dists = structure.get_all_distances(mic=False)
    np.fill_diagonal(dists, np.inf)
    value = float(np.min(dists))
    return value if np.isfinite(value) else -1.0


def _min_adsorbate_surface_distance(structure: Atoms, adsorbate_indices: Sequence[int]) -> float:
    if len(structure) < 2 or not adsorbate_indices:
        return -1.0

    ads_set = set(adsorbate_indices)
    surface_indices = [idx for idx in range(len(structure)) if idx not in ads_set]
    if not surface_indices:
        return -1.0

    positions = structure.get_positions()
    min_dist = float("inf")
    for ads_idx in adsorbate_indices:
        for surf_idx in surface_indices:
            dist = float(np.linalg.norm(positions[ads_idx] - positions[surf_idx]))
            min_dist = min(min_dist, dist)
    return min_dist if np.isfinite(min_dist) else -1.0
