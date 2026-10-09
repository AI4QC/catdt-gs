"""Position calculations, distance metrics, atom manipulation, PBC handling, and staging."""

from __future__ import annotations

import re
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers, covalent_radii

from core.pathway.neb_frame_prep import (
    interpolated_path_min_distances,
    prepare_endpoint_pair_for_interpolation,
)

from camel_agents.schemas import (
    AtomAddSpec,
    PathwayStepSpec,
    WorkflowState,
)

from .common import logger

_current_file_dir = Path(__file__).parent.resolve()
_project_root = _current_file_dir.parent.parent

class WorkflowGeometryMixin:
    @staticmethod
    def _safe_position(candidate: Any, fallback: List[float], structure: Optional[Atoms] = None) -> List[float]:
        if isinstance(candidate, list) and len(candidate) == 3:
            try:
                position = [float(candidate[0]), float(candidate[1]), float(candidate[2])]
                if not all(np.isfinite(v) for v in position):
                    return list(fallback)

                if structure is not None and len(structure) > 0:
                    cell = np.array(structure.cell)
                    if cell.shape == (3, 3):
                        lengths = np.linalg.norm(cell, axis=1)
                        max_len = float(np.max(lengths)) if lengths.size else 0.0
                        if max_len > 0 and max(abs(v) for v in position) > (3.0 * max_len):
                            return list(fallback)

                        if bool(np.any(structure.pbc)):
                            try:
                                frac = np.linalg.pinv(cell.T) @ np.array(position, dtype=float)
                                for dim in range(3):
                                    if bool(structure.pbc[dim]):
                                        frac[dim] = float(frac[dim] % 1.0)
                                wrapped = cell.T @ frac
                                return [float(wrapped[0]), float(wrapped[1]), float(wrapped[2])]
                            except (ValueError, np.linalg.LinAlgError):
                                return list(fallback)
                return position
            except Exception:
                return list(fallback)
        return list(fallback)

    @staticmethod
    def _resolve_atom_index(idx_value: Any, adsorbate_indices: List[int], total_atoms: int) -> int:
        _ = total_atoms  # Pathway atom operations are restricted to adsorbate atoms only.
        try:
            idx = int(idx_value)
        except Exception:
            return -1

        if idx in adsorbate_indices:
            return idx
        if 0 <= idx < len(adsorbate_indices):
            return adsorbate_indices[idx]
        return -1

    @staticmethod
    def _minimum_distance(symbol_a: str, symbol_b: str, scale: float = 0.85) -> float:
        z_a = atomic_numbers.get(symbol_a)
        z_b = atomic_numbers.get(symbol_b)
        if z_a is None or z_b is None:
            return 1.05
        return max(1.05, scale * (covalent_radii[z_a] + covalent_radii[z_b]))

    def _resolve_addition_elements(self, atom_spec: AtomAddSpec) -> List[str]:
        candidates = [
            getattr(atom_spec, "species", ""),
            getattr(atom_spec, "element", ""),
            getattr(atom_spec, "symbol", ""),
        ]

        tokens: List[str] = []
        for raw in candidates:
            text = str(raw or "").strip()
            if not text:
                continue
            parsed = self._extract_species_tokens(text)
            if parsed:
                tokens.extend(parsed)
                continue
            canon = self._canonical_element(text)
            if canon in atomic_numbers:
                tokens.append(canon)

        return tokens

    def _step_formula_target_counts(self, step: PathwayStepSpec) -> Counter:
        reactant_tokens = Counter(self._extract_species_tokens(step.reactant_formula))
        product_tokens = Counter(self._extract_species_tokens(step.product_formula))
        if not reactant_tokens and not product_tokens:
            return Counter()

        target: Counter = Counter()
        for elem in sorted(set(reactant_tokens.keys()) | set(product_tokens.keys())):
            target[elem] = max(reactant_tokens.get(elem, 0), product_tokens.get(elem, 0))
        return target

    @staticmethod
    def _count_elements_on_indices(structure: Atoms, indices: List[int]) -> Counter:
        counts: Counter = Counter()
        for idx in indices:
            if 0 <= idx < len(structure):
                counts[structure[idx].symbol] += 1
        return counts

    def _delete_atoms_with_index_tracking(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        staged_indices: List[int],
        delete_indices: List[int],
    ) -> Tuple[List[int], List[int]]:
        to_delete = sorted({idx for idx in delete_indices if 0 <= idx < len(structure)}, reverse=True)
        if not to_delete:
            return list(adsorbate_indices), list(staged_indices)

        delete_set = set(to_delete)
        for idx in to_delete:
            del structure[idx]

        def _remap(old_idx: int) -> Optional[int]:
            if old_idx in delete_set:
                return None
            shift = sum(1 for d in to_delete if d < old_idx)
            new_idx = old_idx - shift
            if 0 <= new_idx < len(structure):
                return new_idx
            return None

        new_ads: List[int] = []
        for idx in adsorbate_indices:
            mapped = _remap(int(idx))
            if mapped is not None:
                new_ads.append(mapped)

        new_staged: List[int] = []
        for idx in staged_indices:
            mapped = _remap(int(idx))
            if mapped is not None:
                new_staged.append(mapped)

        return sorted(set(new_ads)), sorted(set(new_staged))

    def _enforce_formula_target_multiset(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        staged_indices: List[int],
        fixed_count: int,
        target_counts: Counter,
    ) -> Tuple[List[int], List[int]]:
        if not target_counts:
            return sorted(set(adsorbate_indices)), sorted(set(staged_indices))

        ads = sorted({idx for idx in adsorbate_indices if 0 <= idx < len(structure)})
        staged = sorted({idx for idx in staged_indices if 0 <= idx < len(structure)})
        if not ads:
            return ads, staged

        removable = [idx for idx in ads if idx >= max(0, fixed_count) or idx in set(staged)]
        removable_set = set(removable)
        core_refs = [idx for idx in ads if idx not in removable_set]

        positions = structure.get_positions()
        if core_refs:
            anchor = np.mean(positions[core_refs], axis=0)
        else:
            anchor = np.mean(positions[ads], axis=0)

        current_counts = self._count_elements_on_indices(structure, ads)
        remove_plan: List[int] = []

        for element, current in current_counts.items():
            target = int(target_counts.get(element, 0))
            excess = int(current - target)
            if excess <= 0:
                continue

            candidates = [idx for idx in removable if structure[idx].symbol == element]
            if not candidates:
                continue

            candidates.sort(
                key=lambda idx: (
                    float(np.linalg.norm(positions[idx] - anchor)),
                    float(positions[idx][2]),
                ),
                reverse=True,
            )
            remove_plan.extend(candidates[:excess])

        if remove_plan:
            ads, staged = self._delete_atoms_with_index_tracking(
                structure=structure,
                adsorbate_indices=ads,
                staged_indices=staged,
                delete_indices=remove_plan,
            )

        current_counts = self._count_elements_on_indices(structure, ads)
        add_offset = 0
        for element, target in sorted(target_counts.items()):
            missing = int(target - current_counts.get(element, 0))
            for _ in range(max(0, missing)):
                pos = self._get_unbonded_position(structure, ads)
                pos[0] += 0.20 * add_offset
                pos[1] -= 0.20 * add_offset
                structure.append(Atom(element, position=pos))
                new_idx = len(structure) - 1
                ads.append(new_idx)
                staged.append(new_idx)
                add_offset += 1

        return sorted(set(ads)), sorted(set(staged))

    def _relax_close_contacts(self, structure: Atoms, movable_indices: List[int], max_iter: int = 8) -> None:
        if not movable_indices:
            return

        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()
        movable = set(movable_indices)

        for _ in range(max_iter):
            changed = False
            for i in range(len(structure)):
                for j in range(i + 1, len(structure)):
                    min_dist = self._minimum_distance(symbols[i], symbols[j])
                    dist = float(np.linalg.norm(positions[j] - positions[i]))
                    if dist >= min_dist:
                        continue

                    if j in movable:
                        move_idx, anchor_idx = j, i
                    elif i in movable:
                        move_idx, anchor_idx = i, j
                    else:
                        continue

                    direction = positions[move_idx] - positions[anchor_idx]
                    norm = float(np.linalg.norm(direction))
                    if norm < 1e-8:
                        direction = np.array([0.0, 0.0, 1.0])
                        norm = 1.0
                    direction = direction / norm

                    positions[move_idx] = positions[anchor_idx] + direction * (min_dist + 0.05)
                    changed = True

            if not changed:
                break
        else:
            min_remaining = float("inf")
            for i in range(len(structure)):
                for j in range(i + 1, len(structure)):
                    min_remaining = min(
                        min_remaining,
                        float(np.linalg.norm(positions[j] - positions[i])),
                    )
            logger.warning(
                "_relax_close_contacts: max_iter=%d exhausted with min pair "
                "distance %.3f Å; close contacts may remain unresolved.",
                max_iter,
                min_remaining,
            )

        structure.set_positions(positions)

    def _get_pathway_energy_helper(self) -> Optional["PathwayPredictor"]:
        if self._pathway_energy_helper is not None:
            return self._pathway_energy_helper

        fairchem_root = _project_root / "deps" / "fairchem"
        if not fairchem_root.exists():
            return None

        try:
            from core.pathway.pathway_predictor import PathwayPredictor
            self._pathway_energy_helper = PathwayPredictor(
                fairchem_root=str(fairchem_root),
                use_gpu=True,
                verbose=False,
                keep_files=False,
            )
        except Exception as exc:
            logger.warning("Failed to initialize PathwayPredictor energy helper: %s", exc)
            self._pathway_energy_helper = None

        return self._pathway_energy_helper

    def _compute_reaction_energy_with_correction(
        self,
        reactant_label: str,
        product_label: str,
        reactant_energy: Optional[float],
        product_energy: Optional[float],
    ) -> float:
        if reactant_energy is None or product_energy is None:
            return 0.0

        corrected = float(product_energy - reactant_energy)
        helper = self._get_pathway_energy_helper()
        if helper is None:
            return corrected

        try:
            element_change = helper.calculate_element_change(reactant_label, product_label)
            if element_change:
                corrected += float(helper.calculate_gas_phase_correction(element_change))
        except Exception as exc:
            logger.debug("Reaction energy correction skipped for %s -> %s: %s", reactant_label, product_label, exc)

        return corrected

    @staticmethod
    def _all_pair_distances_mic(structure: Atoms) -> np.ndarray:
        """Full distance matrix; minimum-image convention when any PBC is set."""
        use_mic = bool(np.any(structure.pbc))
        return structure.get_all_distances(mic=use_mic)

    def _min_pair_distance(self, structure: Atoms, indices: List[int]) -> float:
        n = len(structure)
        out_of_range = [i for i in indices if not (0 <= i < n)]
        if out_of_range:
            raise IndexError(
                f"_min_pair_distance: indices {out_of_range} out of range for "
                f"structure with {n} atoms (valid 0..{n-1}). full indices={list(indices)}"
            )
        if len(indices) < 2:
            return -1.0
        dists = self._all_pair_distances_mic(structure)
        sub = dists[np.ix_(list(indices), list(indices))]
        iu = np.triu_indices(len(indices), k=1)
        vals = sub[iu]
        return float(np.min(vals)) if vals.size else -1.0

    @staticmethod
    def _min_distance_any_pair(structure: Atoms) -> float:
        if len(structure) < 2:
            return -1.0
        dists = WorkflowGeometryMixin._all_pair_distances_mic(structure)
        iu = np.triu_indices(len(structure), k=1)
        vals = dists[iu]
        return float(np.min(vals)) if vals.size else -1.0

    @staticmethod
    def _min_adsorbate_surface_distance(structure: Atoms, adsorbate_indices: List[int]) -> float:
        if len(structure) < 2 or not adsorbate_indices:
            return -1.0

        valid_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(structure)]
        if not valid_ads:
            return -1.0

        ads_set = set(valid_ads)
        surface_indices = [i for i in range(len(structure)) if i not in ads_set]
        if not surface_indices:
            return -1.0

        dists = WorkflowGeometryMixin._all_pair_distances_mic(structure)
        sub = dists[np.ix_(valid_ads, surface_indices)]
        return float(np.min(sub)) if sub.size else -1.0

    def _interpolated_path_min_distances(
        self,
        reactant: Atoms,
        product: Atoms,
        adsorbate_indices: List[int],
        n_images: int = 7,
    ) -> Tuple[float, int, float, int]:
        if len(reactant) != len(product):
            return -1.0, -1, -1.0, -1

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

            min_any = self._min_distance_any_pair(frame_atoms)
            if min_any > 0 and min_any < min_any_pair:
                min_any_pair = min_any
                min_any_pair_frame = frame_idx

            min_ads_surf = self._min_adsorbate_surface_distance(frame_atoms, adsorbate_indices)
            if min_ads_surf > 0 and min_ads_surf < min_ads_surface:
                min_ads_surface = min_ads_surf
                min_ads_surface_frame = frame_idx

        if min_any_pair == float("inf"):
            min_any_pair = -1.0
        if min_ads_surface == float("inf"):
            min_ads_surface = -1.0

        return float(min_any_pair), int(min_any_pair_frame), float(min_ads_surface), int(min_ads_surface_frame)

    def _unwrap_adsorbate_positions(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        reference_indices: Optional[List[int]] = None,
    ) -> np.ndarray:
        if not adsorbate_indices:
            return np.zeros((0, 3), dtype=float)

        positions = structure.get_positions()
        unwrapped = np.array([positions[idx] for idx in adsorbate_indices], dtype=float)

        if len(structure) == 0 or (not bool(np.any(structure.pbc))):
            return unwrapped

        cell = np.array(structure.cell, dtype=float)
        if cell.shape != (3, 3):
            return unwrapped

        try:
            inv_t = np.linalg.pinv(cell.T)
        except Exception:
            return unwrapped

        ref_pool = [idx for idx in (reference_indices or adsorbate_indices) if idx in set(adsorbate_indices)]
        if not ref_pool:
            ref_pool = list(adsorbate_indices)
        ref_idx = ref_pool[0]
        ref_frac = inv_t @ np.array(positions[ref_idx], dtype=float)

        for row, atom_idx in enumerate(adsorbate_indices):
            frac = inv_t @ np.array(positions[atom_idx], dtype=float)
            delta = frac - ref_frac
            for dim in range(3):
                if bool(structure.pbc[dim]):
                    delta[dim] = float(delta[dim] - np.round(delta[dim]))
            frac_aligned = ref_frac + delta
            unwrapped[row] = cell.T @ frac_aligned

        return unwrapped

    def _align_adsorbate_images_inplace(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        reference_indices: Optional[List[int]] = None,
    ) -> None:
        if not adsorbate_indices:
            return
        aligned = self._unwrap_adsorbate_positions(
            structure=structure,
            adsorbate_indices=adsorbate_indices,
            reference_indices=reference_indices,
        )
        for row, atom_idx in enumerate(adsorbate_indices):
            if 0 <= atom_idx < len(structure):
                structure.positions[atom_idx] = aligned[row]

    def _align_movable_adsorbates_to_core(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        core_reference_indices: List[int],
        movable_indices: List[int],
    ) -> None:
        """Align only movable adsorbate atoms to core reference periodic image."""
        if not adsorbate_indices or not movable_indices:
            return

        ordered_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(structure)]
        if not ordered_ads:
            return

        aligned = self._unwrap_adsorbate_positions(
            structure=structure,
            adsorbate_indices=ordered_ads,
            reference_indices=core_reference_indices or ordered_ads,
        )
        row_map = {idx: row for row, idx in enumerate(ordered_ads)}

        for idx in movable_indices:
            row = row_map.get(idx)
            if row is None:
                continue
            if 0 <= idx < len(structure):
                structure.positions[idx] = aligned[row]

    def _adsorbate_centroid_xy(self, structure: Atoms, adsorbate_indices: List[int]) -> Optional[np.ndarray]:
        if not adsorbate_indices:
            return None
        positions = self._unwrap_adsorbate_positions(structure, adsorbate_indices)
        if positions.size == 0:
            return None
        return np.mean(positions[:, :2], axis=0)

    def _estimate_surface_top_z(self, structure: Atoms, adsorbate_indices: List[int]) -> float:
        """Estimate exposed top-surface z from all non-adsorbate atoms.

        SMSI catalysts can have a thin TiOx overlayer on a much larger metal
        particle. Counting the dominant slab element would pick the buried
        metal top instead of the exposed oxide shell.
        """
        if len(structure) == 0:
            return 0.0

        ads_set = set(adsorbate_indices or [])
        surface_candidates = [i for i in range(len(structure)) if i not in ads_set]
        positions = structure.get_positions()

        if not surface_candidates:
            return float(np.max(positions[:, 2]))

        return float(np.max(positions[surface_candidates, 2]))

    @staticmethod
    def _contains_gas_phase_hint(*texts: Any) -> bool:
        for text in texts:
            t = str(text or "").lower()
            if "(g" in t or " gas" in t or "gas-phase" in t:
                return True
        return False

    def _match_adsorbate_atoms_by_element(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_ads: List[int],
        product_ads: List[int],
    ) -> List[Tuple[int, int]]:
        if not reactant_ads or not product_ads:
            return []

        r_pos = reactant.get_positions()
        p_pos = product.get_positions()
        r_sym = reactant.get_chemical_symbols()
        p_sym = product.get_chemical_symbols()

        pairs: List[Tuple[int, int]] = []
        used_product: set[int] = set()

        try:
            from scipy.optimize import linear_sum_assignment  # type: ignore
        except Exception:
            linear_sum_assignment = None

        by_elem: Dict[str, Tuple[List[int], List[int]]] = {}
        for idx in reactant_ads:
            by_elem.setdefault(r_sym[idx], ([], []))[0].append(idx)
        for idx in product_ads:
            by_elem.setdefault(p_sym[idx], ([], []))[1].append(idx)

        for elem, (r_list, p_list) in by_elem.items():
            if not r_list or not p_list:
                continue

            if linear_sum_assignment is not None:
                cost = np.zeros((len(r_list), len(p_list)), dtype=float)
                for i, r_idx in enumerate(r_list):
                    for j, p_idx in enumerate(p_list):
                        cost[i, j] = float(np.linalg.norm(r_pos[r_idx] - p_pos[p_idx]))
                row_ind, col_ind = linear_sum_assignment(cost)
                for ri, cj in zip(row_ind.tolist(), col_ind.tolist()):
                    r_idx = r_list[ri]
                    p_idx = p_list[cj]
                    if p_idx in used_product:
                        continue
                    pairs.append((r_idx, p_idx))
                    used_product.add(p_idx)
            else:
                for r_idx in r_list:
                    best_idx = None
                    best_dist = float("inf")
                    for p_idx in p_list:
                        if p_idx in used_product:
                            continue
                        dist = float(np.linalg.norm(r_pos[r_idx] - p_pos[p_idx]))
                        if dist < best_dist:
                            best_dist = dist
                            best_idx = p_idx
                    if best_idx is not None:
                        pairs.append((r_idx, best_idx))
                        used_product.add(best_idx)

        return pairs

    def _reorder_product_to_match_reactant(
        self,
        reactant: Atoms,
        product: Atoms,
        product_adsorbate_indices: List[int],
        product_staged_indices: Optional[List[int]] = None,
        atom_mapping_constraints: Optional[List[Dict[str, Any]]] = None,
    ) -> Tuple[Atoms, List[int], bool] | Tuple[Atoms, List[int], List[int], bool]:
        include_staged = product_staged_indices is not None
        staged_in = list(product_staged_indices or [])

        def _ret(atoms: Atoms, ads: List[int], staged: List[int], did_reorder: bool):
            if include_staged:
                return atoms, ads, staged, did_reorder
            return atoms, ads, did_reorder

        if len(reactant) != len(product):
            return _ret(product, list(product_adsorbate_indices), staged_in, False)

        reactant_symbols = reactant.get_chemical_symbols()
        product_symbols = product.get_chemical_symbols()
        product_positions = product.get_positions()
        reactant_positions = reactant.get_positions()
        if sorted(reactant_symbols) != sorted(product_symbols):
            return _ret(product, list(product_adsorbate_indices), staged_in, False)
        constraints = self._valid_atom_mapping_constraints(
            reactant=reactant,
            product=product,
            atom_mapping_constraints=atom_mapping_constraints,
        )

        if reactant_symbols == product_symbols and not constraints:
            diffs = product_positions - reactant_positions
            cell = reactant.get_cell()
            if bool(np.any(reactant.pbc)) and float(np.linalg.norm(cell[0])) > 0.1:
                frac = np.linalg.solve(cell.T, diffs.T).T
                frac -= np.round(frac)
                diffs = (cell.T @ frac.T).T
            displacements = np.linalg.norm(diffs, axis=1)
            if np.max(displacements) < 3.0:
                return _ret(product, list(product_adsorbate_indices), staged_in, False)

        try:
            from scipy.optimize import linear_sum_assignment  # type: ignore
        except Exception:
            linear_sum_assignment = None

        assignment: Dict[int, int] = {ridx: pidx for ridx, pidx in constraints}
        constrained_r = set(assignment)
        constrained_p = set(assignment.values())
        by_symbol: Dict[str, Tuple[List[int], List[int]]] = {}
        for ridx, sym in enumerate(reactant_symbols):
            by_symbol.setdefault(sym, ([], []))[0].append(ridx)
        for pidx, sym in enumerate(product_symbols):
            by_symbol.setdefault(sym, ([], []))[1].append(pidx)

        for sym, (r_group, p_group) in by_symbol.items():
            if len(r_group) != len(p_group):
                return _ret(product, list(product_adsorbate_indices), staged_in, False)
            if not r_group:
                continue
            r_group = [idx for idx in r_group if idx not in constrained_r]
            p_group = [idx for idx in p_group if idx not in constrained_p]
            if len(r_group) != len(p_group):
                return _ret(product, list(product_adsorbate_indices), staged_in, False)
            if not r_group:
                continue

            if linear_sum_assignment is not None:
                cost = np.zeros((len(r_group), len(p_group)), dtype=float)
                cell = reactant.get_cell()
                pbc = reactant.pbc
                use_mic = bool(np.any(pbc)) and float(np.linalg.norm(cell[0])) > 0.1
                for i, ridx in enumerate(r_group):
                    for j, pidx in enumerate(p_group):
                        diff = product_positions[pidx] - reactant_positions[ridx]
                        if use_mic:
                            frac = np.linalg.solve(cell.T, diff)
                            frac -= np.round(frac)
                            diff = cell.T @ frac
                        cost[i, j] = float(np.linalg.norm(diff))
                row_ind, col_ind = linear_sum_assignment(cost)
                for rr, cc in zip(row_ind.tolist(), col_ind.tolist()):
                    assignment[r_group[rr]] = p_group[cc]
            else:
                remaining = list(p_group)
                for ridx in r_group:
                    chosen = min(
                        remaining,
                        key=lambda pidx: float(np.linalg.norm(product_positions[pidx] - reactant_positions[ridx])),
                    )
                    assignment[ridx] = chosen
                    remaining.remove(chosen)

        try:
            reorder_sequence = [assignment[idx] for idx in range(len(reactant_symbols))]
        except Exception:
            return _ret(product, list(product_adsorbate_indices), staged_in, False)

        if reorder_sequence == list(range(len(reactant_symbols))):
            return _ret(product, list(product_adsorbate_indices), staged_in, False)

        reordered_product = product[reorder_sequence]
        old_to_new = {old: new for new, old in enumerate(reorder_sequence)}
        new_ads = [old_to_new[idx] for idx in product_adsorbate_indices if idx in old_to_new]
        new_ads = sorted(new_ads)
        new_staged = [old_to_new[idx] for idx in staged_in if idx in old_to_new]
        new_staged = sorted(new_staged)
        new_ads_set = set(new_ads)
        reordered_product.info["adsorbate_indices"] = list(new_ads)
        reordered_product.info["staged_indices"] = list(new_staged)
        reordered_product.info["surface_indices"] = [
            i for i in range(len(reordered_product))
            if i not in new_ads_set and i not in set(new_staged)
        ]
        reordered_product.info["_reorder_sequence"] = list(reorder_sequence)
        if atom_mapping_constraints:
            reordered_product.info["atom_mapping_constraints"] = [
                {
                    **dict(item),
                    "product_index_after_reorder": old_to_new.get(int(item.get("product_index", -1)), None),
                }
                for item in atom_mapping_constraints
                if isinstance(item, dict)
            ]
        return _ret(reordered_product, new_ads, new_staged, True)

    @staticmethod
    def _valid_atom_mapping_constraints(
        reactant: Atoms,
        product: Atoms,
        atom_mapping_constraints: Optional[List[Dict[str, Any]]],
    ) -> List[Tuple[int, int]]:
        constraints: List[Tuple[int, int]] = []
        used_r: set[int] = set()
        used_p: set[int] = set()
        r_symbols = reactant.get_chemical_symbols()
        p_symbols = product.get_chemical_symbols()

        for item in atom_mapping_constraints or []:
            if not isinstance(item, dict):
                continue
            try:
                ridx = int(item.get("reactant_index"))
                pidx = int(item.get("product_index"))
            except Exception:
                continue
            if not (0 <= ridx < len(reactant) and 0 <= pidx < len(product)):
                continue
            if ridx in used_r or pidx in used_p:
                continue
            if r_symbols[ridx] != p_symbols[pidx]:
                continue
            expected = str(item.get("element", "") or "").strip()
            if expected and expected != r_symbols[ridx]:
                continue
            constraints.append((ridx, pidx))
            used_r.add(ridx)
            used_p.add(pidx)
        return constraints

    def _get_unbonded_position(self, structure: Atoms, adsorbate_indices: List[int]) -> List[float]:
        positions = structure.get_positions()
        if adsorbate_indices:
            ref = np.mean(positions[adsorbate_indices], axis=0)
            surface_z = self._estimate_surface_top_z(structure, adsorbate_indices)
            z_target = float(np.clip(max(ref[2] + 0.6, surface_z + 1.0), surface_z + 0.9, surface_z + 3.2))
            return [float(ref[0] + 1.0), float(ref[1] + 0.6), z_target]

        ref = positions[np.argmax(positions[:, 2])]
        surface_z = float(np.max(positions[:, 2]))
        return [float(ref[0] + 1.0), float(ref[1] + 1.0), float(surface_z + 1.2)]

    def _get_bonded_position(self, structure: Atoms, adsorbate_indices: List[int], element: str) -> List[float]:
        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()

        if adsorbate_indices:
            anchor_idx = min(adsorbate_indices, key=lambda idx: positions[idx][2])
            ref = positions[anchor_idx]
            anchor_symbol = symbols[anchor_idx]
        else:
            anchor_idx = int(np.argmax(positions[:, 2]))
            ref = positions[anchor_idx]
            anchor_symbol = symbols[anchor_idx]

        bond_dist = max(1.0, self._minimum_distance(anchor_symbol, element, scale=0.9))
        return [float(ref[0]), float(ref[1]), float(ref[2] + bond_dist)]

    def _find_surface_adsorption_site(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        min_dist_from_ads: float = 3.5,
        max_dist_from_ads: float = 5.5,
        offset: int = 0,
    ) -> Optional[List[float]]:
        """Find a hollow/bridge adsorption site on the surface near the adsorbate.

        Uses adsorbate atom positions to find the nearest surface metal atoms
        (by 3D distance), then builds a hollow/bridge site from them.
        Works correctly on stepped/kinked surfaces where z-layer heuristics fail.
        """
        if len(structure) == 0:
            return None

        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()

        ads_set = set(adsorbate_indices or [])
        heavy_ads = [i for i in adsorbate_indices if symbols[i] not in ("H",) and 0 <= i < len(structure)]
        ref_ads = heavy_ads if heavy_ads else [i for i in adsorbate_indices if 0 <= i < len(structure)]
        if not ref_ads:
            return None
        ref_pos = np.mean(positions[ref_ads], axis=0)

        surface_indices = [i for i in range(len(structure)) if i not in ads_set]
        if not surface_indices:
            return None

        # Sort ALL surface metal atoms by 3D distance to adsorbate reference
        dist_list = []
        for idx in surface_indices:
            d = float(np.linalg.norm(positions[idx] - ref_pos))
            dist_list.append((idx, d))
        dist_list.sort(key=lambda x: x[1])

        candidates = [(idx, d) for idx, d in dist_list if min_dist_from_ads <= d <= max_dist_from_ads]
        if not candidates:
            candidates = [(idx, d) for idx, d in dist_list if 1.5 <= d <= 6.5]
        if not candidates:
            candidates = [(idx, d) for idx, d in dist_list if d > 1.0]
        if not candidates:
            return None

        group_size = min(3, len(candidates))
        start = (offset * max(group_size, 2)) % len(candidates)
        site_atoms = []
        for i in range(group_size):
            site_atoms.append(candidates[(start + i) % len(candidates)])

        site_center = np.mean([positions[idx, :2] for idx, _ in site_atoms], axis=0)
        site_z = float(np.mean([positions[idx, 2] for idx, _ in site_atoms]))

        z = float(site_z + 1.7 + 0.2 * offset)
        return [float(site_center[0]), float(site_center[1]), z]

    def enumerate_surface_sites_near_anchor(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        anchor_pos: np.ndarray,
        element: str,
        min_dist: float = 2.5,
        max_dist: float = 4.0,
    ) -> List[Dict[str, Any]]:
        """Enumerate ALL surface adsorption sites near a given anchor atom.

        Uses pymatgen ``AdsorbateSiteFinder`` (in an isolated subprocess) for
        robust site identification on any surface type — metals, oxides, alloys,
        stepped surfaces — then filters by distance from the anchor atom.

        The pymatgen subprocess is the ONLY site-enumeration path. If it fails,
        this method raises RuntimeError (retries are handled inside the
        subprocess call; no heuristic fallback).

        Parameters
        ----------
        structure : Atoms
            Slab + adsorbate structure.
        adsorbate_indices : list of int
            Indices of the main adsorbate atoms.
        anchor_pos : ndarray, shape (3,)
            3D position of the anchor atom (the atom the co-adsorbate will bond
            to / dissociate from in the other endpoint).
        element : str
            Element to place (reserved for future site-height hinting; currently
            not consumed by the pymatgen worker, which treats sites as generic).
        min_dist, max_dist : float
            Distance range from anchor_pos (Å).

        Returns
        -------
        list of dict
            Each entry: ``{"position": [x,y,z], "dist_from_anchor": float, "site_label": str}``
            Sorted by distance from anchor (closest first). Empty list if the
            subprocess succeeded but no sites fell in the requested distance band.
        """
        sites = self._enumerate_surface_sites_with_pymatgen_subprocess(
            structure=structure,
            adsorbate_indices=adsorbate_indices,
            anchor_pos=anchor_pos,
            min_dist=min_dist,
            max_dist=max_dist,
        )
        if sites is None:
            raise RuntimeError(
                "pymatgen adsorption site subprocess failed after all retries. "
                "Increase CATDT_PYMATGEN_SITE_RETRIES or CATDT_PYMATGEN_SITE_TIMEOUT_SEC, "
                "or inspect the subprocess stderr tail logged above."
            )
        return sites

    @staticmethod
    def _wrap_cartesian_position(structure: Atoms, position: np.ndarray) -> np.ndarray:
        if len(structure) == 0 or (not bool(np.any(structure.pbc))):
            return np.array(position, dtype=float)

        cell = np.array(structure.cell, dtype=float)
        if cell.shape != (3, 3):
            return np.array(position, dtype=float)

        try:
            frac = np.linalg.pinv(cell.T) @ np.array(position, dtype=float)
        except Exception:
            return np.array(position, dtype=float)

        for dim in range(3):
            if bool(structure.pbc[dim]):
                frac[dim] = float(frac[dim] % 1.0)
        return np.array(cell.T @ frac, dtype=float)

    @staticmethod
    def _pbc_cartesian_delta(structure: Atoms, cart_a: np.ndarray, cart_b: np.ndarray) -> np.ndarray:
        delta = np.array(cart_b, dtype=float) - np.array(cart_a, dtype=float)
        if len(structure) == 0 or (not bool(np.any(structure.pbc))):
            return delta

        cell = np.array(structure.cell, dtype=float)
        if cell.shape != (3, 3):
            return delta

        try:
            inv_t = np.linalg.pinv(cell.T)
        except Exception:
            return delta

        frac_delta = inv_t @ delta
        for dim in range(3):
            if bool(structure.pbc[dim]):
                frac_delta[dim] = float(frac_delta[dim] - np.round(frac_delta[dim]))
        return np.array(cell.T @ frac_delta, dtype=float)

    def _pbc_distance(self, structure: Atoms, cart_a: np.ndarray, cart_b: np.ndarray) -> float:
        return float(np.linalg.norm(self._pbc_cartesian_delta(structure, cart_a, cart_b)))

    def _enumerate_surface_sites_with_pymatgen_subprocess(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        anchor_pos: np.ndarray,
        min_dist: float,
        max_dist: float,
    ) -> Optional[List[Dict[str, Any]]]:
        ads_set = set(adsorbate_indices or [])
        slab_indices = [i for i in range(len(structure)) if i not in ads_set]
        if not slab_indices:
            return []

        slab_ase = structure[slab_indices].copy()
        slab_ase.pbc = structure.pbc
        slab_ase.cell = structure.cell

        timeout_s = max(5, int(os.getenv("CATDT_PYMATGEN_SITE_TIMEOUT_SEC", "30")))
        retry_count = max(1, int(os.getenv("CATDT_PYMATGEN_SITE_RETRIES", "6")))
        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "1")
        env.setdefault("OPENBLAS_NUM_THREADS", "1")
        env.setdefault("MKL_NUM_THREADS", "1")
        env.setdefault("NUMEXPR_NUM_THREADS", "1")
        # Prevent child process from touching GPU — avoids CUDA fork SIGSEGV
        env["CUDA_VISIBLE_DEVICES"] = ""

        # Sync + release CUDA context before fork to avoid nondeterministic
        # SIGSEGV in child processes (CUDA driver atexit handlers)
        try:
            import torch
            if torch.cuda.is_available() and torch.cuda.current_device() is not None:
                torch.cuda.synchronize()
        except Exception:
            pass

        with tempfile.TemporaryDirectory(prefix="catdt_pmg_sites_") as tmpdir:
            slab_path = Path(tmpdir) / "slab.traj"
            try:
                from ase.io import write
                write(str(slab_path), slab_ase)
            except Exception as exc:
                logger.warning("Failed to serialize slab for pymatgen site subprocess: %s", exc)
                return None

            payload = {
                "slab_path": str(slab_path),
                "anchor_pos": [float(anchor_pos[0]), float(anchor_pos[1]), float(anchor_pos[2])],
                "min_dist": float(min_dist),
                "max_dist": float(max_dist),
            }
            worker_path = Path(__file__).with_name("pymatgen_site_worker.py")
            proc: Optional[subprocess.CompletedProcess[str]] = None
            for attempt in range(1, retry_count + 1):
                try:
                    proc = subprocess.run(
                        [
                            sys.executable,
                            "-I",
                            "-B",
                            "-X",
                            "faulthandler",
                            str(worker_path),
                            json.dumps(payload),
                        ],
                        capture_output=True,
                        text=True,
                        timeout=timeout_s,
                        env=env,
                        cwd=tmpdir,
                        check=False,
                    )
                except subprocess.TimeoutExpired:
                    logger.warning(
                        "pymatgen adsorption site subprocess timed out after %ds (attempt %d/%d)",
                        timeout_s,
                        attempt,
                        retry_count,
                    )
                    proc = None
                except Exception as exc:
                    logger.warning(
                        "Failed to launch pymatgen adsorption site subprocess (attempt %d/%d): %s",
                        attempt,
                        retry_count,
                        exc,
                    )
                    proc = None

                if proc is not None and proc.returncode == 0:
                    break

                if proc is not None:
                    stderr_tail = (proc.stderr or "").strip().splitlines()[-1] if proc.stderr else ""
                    logger.warning(
                        "pymatgen adsorption site subprocess failed (code=%s, attempt %d/%d). stderr_tail=%s",
                        proc.returncode,
                        attempt,
                        retry_count,
                        stderr_tail,
                    )

            if proc is None or proc.returncode != 0:
                return None

        stdout = (proc.stdout or "").strip()
        if not stdout:
            return None

        try:
            parsed = json.loads(stdout)
        except Exception as exc:
            logger.warning("Failed to parse pymatgen adsorption site subprocess JSON: %s", exc)
            return None

        if not isinstance(parsed, list):
            return None

        sites: List[Dict[str, Any]] = []
        for item in parsed:
            if not isinstance(item, dict):
                continue
            pos = item.get("position")
            if not isinstance(pos, list) or len(pos) != 3:
                continue
            sites.append(
                {
                    "position": [float(pos[0]), float(pos[1]), float(pos[2])],
                    "dist_from_anchor": float(item.get("dist_from_anchor", 0.0)),
                    "site_label": str(item.get("site_label", "unknown")),
                }
            )
        return sites

    def _get_remote_position(self, structure: Atoms) -> List[float]:
        positions = structure.get_positions()
        max_z = np.max(positions[:, 2])
        center_x = float(np.mean(positions[:, 0]))
        center_y = float(np.mean(positions[:, 1]))
        return [center_x + 1.8, center_y + 1.6, float(max_z + 2.8)]

    def _get_local_staging_position(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        core_adsorbate_indices: List[int],
        offset: int = 0,
    ) -> List[float]:
        positions = structure.get_positions()
        anchor_indices = core_adsorbate_indices or adsorbate_indices

        if anchor_indices:
            anchor = np.mean(positions[anchor_indices], axis=0)
            surface_z = self._estimate_surface_top_z(structure, adsorbate_indices)
            radius = 1.35 + 0.15 * (offset % 3)
            theta = (np.pi / 3.0) * (offset % 6)
            x = float(anchor[0] + radius * np.cos(theta))
            y = float(anchor[1] + radius * np.sin(theta))
            z_raw = max(anchor[2] + 0.45, surface_z + 1.0 + 0.05 * offset)
            z = float(np.clip(z_raw, surface_z + 0.9, surface_z + 3.2))
            return [x, y, z]

        fallback = self._get_unbonded_position(structure, adsorbate_indices)
        fallback[0] += 0.2 * offset
        fallback[1] -= 0.2 * offset
        return fallback

    def _stage_outward_from_center(
        self,
        base_position: np.ndarray,
        center_position: np.ndarray,
        surface_z: float,
        radial_distance: float,
        min_height: float = 0.5,
        max_height: float = 3.4,
    ) -> np.ndarray:
        base = np.array(base_position, dtype=float)
        center = np.array(center_position, dtype=float)

        direction = base - center
        direction[2] = max(direction[2], 0.2)
        norm = float(np.linalg.norm(direction))
        if norm < 1e-8:
            direction = np.array([1.0, 0.0, 0.3], dtype=float)
            norm = float(np.linalg.norm(direction))
        direction = direction / norm

        staged = base + direction * float(radial_distance)
        staged[2] = float(np.clip(staged[2], surface_z + min_height, surface_z + max_height))
        return staged

    def _nudge_adsorbate_set_outward(
        self,
        structure: Atoms,
        move_indices: List[int],
        reference_indices: List[int],
        adsorbate_indices: List[int],
        radial_base: float = 0.8,
        min_sep: float = 0.95,
        min_height: float = 0.75,
        max_height: float = 3.4,
    ) -> None:
        if not move_indices:
            return

        valid_ads = [i for i in adsorbate_indices if 0 <= i < len(structure)]
        valid_move = [i for i in move_indices if i in set(valid_ads)]
        if not valid_move:
            return

        ref_pool = [i for i in reference_indices if i in set(valid_ads)]
        if not ref_pool:
            ref_pool = [i for i in valid_ads if i not in set(valid_move)]
        if not ref_pool:
            ref_pool = list(valid_move)

        positions = structure.get_positions()
        center = np.mean(positions[ref_pool], axis=0)
        surface_z = self._estimate_surface_top_z(structure, valid_ads)

        placed = [i for i in valid_ads if i not in set(valid_move)]
        for offset, idx_move in enumerate(valid_move):
            base = positions[idx_move].copy()
            staged = self._stage_outward_from_center(
                base_position=base,
                center_position=center,
                surface_z=surface_z,
                radial_distance=radial_base + 0.35 * offset,
                min_height=min_height,
                max_height=max_height,
            )

            for _ in range(8):
                close_idx = -1
                close_dist = float("inf")
                for idx_other in placed:
                    dist = float(np.linalg.norm(staged - positions[idx_other]))
                    if dist < close_dist:
                        close_dist = dist
                        close_idx = idx_other

                if close_idx < 0 or close_dist >= min_sep:
                    break

                direction = staged - positions[close_idx]
                direction[2] = max(float(direction[2]), 0.12)
                norm = float(np.linalg.norm(direction))
                if norm < 1e-8:
                    direction = np.array([1.0, 0.0, 0.2], dtype=float)
                    norm = float(np.linalg.norm(direction))
                staged = staged + direction / norm * (min_sep - close_dist + 0.08)
                staged[2] = float(np.clip(staged[2], surface_z + min_height, surface_z + max_height))

            structure.positions[idx_move] = staged
            positions[idx_move] = staged
            placed.append(idx_move)

    def _staged_connectivity_metrics(
        self,
        structure: Atoms,
        staged_indices: List[int],
    ) -> Tuple[int, float, float]:
        """Return (connected_components, max_nearest_neighbor, min_pair_distance) for staged atoms."""
        staged = sorted({idx for idx in staged_indices if 0 <= idx < len(structure)})
        if len(staged) < 2:
            return (1 if staged else 0), 0.0, -1.0

        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()

        adjacency: Dict[int, set[int]] = {idx: set() for idx in staged}
        nearest: Dict[int, float] = {idx: float("inf") for idx in staged}
        min_pair = float("inf")

        for i, idx_i in enumerate(staged):
            for idx_j in staged[i + 1 :]:
                dist = float(np.linalg.norm(positions[idx_i] - positions[idx_j]))
                min_pair = min(min_pair, dist)
                nearest[idx_i] = min(nearest[idx_i], dist)
                nearest[idx_j] = min(nearest[idx_j], dist)

                bond_cutoff = min(
                    2.30,
                    1.30 * self._minimum_distance(symbols[idx_i], symbols[idx_j], scale=1.0),
                )
                if dist <= bond_cutoff:
                    adjacency[idx_i].add(idx_j)
                    adjacency[idx_j].add(idx_i)

        seen: set[int] = set()
        components = 0
        for idx in staged:
            if idx in seen:
                continue
            components += 1
            stack = [idx]
            seen.add(idx)
            while stack:
                current = stack.pop()
                for neigh in adjacency[current]:
                    if neigh not in seen:
                        seen.add(neigh)
                        stack.append(neigh)

        max_nearest = max((v for v in nearest.values() if np.isfinite(v)), default=0.0)
        min_pair_dist = min_pair if np.isfinite(min_pair) else -1.0
        return components, float(max_nearest), float(min_pair_dist)

    def _enforce_staged_fragment_cohesion(
        self,
        structure: Atoms,
        staged_indices: List[int],
        adsorbate_indices: List[int],
        core_reference_indices: List[int],
        fixed_count: int,
        allows_gas: bool,
    ) -> None:
        """Keep staged extras geometrically coherent without moving fixed core coordinates."""
        if allows_gas:
            return

        ads_set = {idx for idx in adsorbate_indices if 0 <= idx < len(structure)}
        movable = sorted(
            {
                idx
                for idx in staged_indices
                if idx in ads_set and idx >= max(0, int(fixed_count)) and 0 <= idx < len(structure)
            }
        )
        if not movable:
            return

        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()
        ref_pool = [idx for idx in core_reference_indices if idx in ads_set and idx not in set(movable)]
        if not ref_pool:
            ref_pool = [idx for idx in ads_set if idx not in set(movable)]
        if not ref_pool:
            ref_pool = list(movable)

        center = np.mean(positions[ref_pool], axis=0)
        surface_z = self._estimate_surface_top_z(structure, list(ads_set))

        local_base = np.array(
            self._get_local_staging_position(
                structure=structure,
                adsorbate_indices=list(ads_set),
                core_adsorbate_indices=ref_pool,
                offset=0,
            ),
            dtype=float,
        )
        radial = local_base[:2] - center[:2]
        radial_norm = float(np.linalg.norm(radial))
        if radial_norm < 1e-8:
            radial = np.array([1.0, 0.0], dtype=float)
            radial_norm = 1.0

        if len(movable) == 1:
            idx_move = movable[0]
            cur = np.array(structure.positions[idx_move], dtype=float)
            lateral = float(np.linalg.norm(cur[:2] - center[:2]))
            height = float(cur[2] - surface_z)
            needs_restage = lateral < 1.10 or lateral > 4.20 or height < 0.80 or height > 2.60
            if needs_restage:
                target = np.array(
                    [
                        float(center[0] + radial[0] / radial_norm * 2.05),
                        float(center[1] + radial[1] / radial_norm * 2.05),
                        float(np.clip(local_base[2], surface_z + 1.00, surface_z + 2.25)),
                    ],
                    dtype=float,
                )
                structure.positions[idx_move] = target
                self._relax_close_contacts(structure, [idx_move], max_iter=10)
                structure.positions[idx_move, 2] = float(
                    np.clip(structure.positions[idx_move, 2], surface_z + 0.95, surface_z + 2.40)
                )
            return

        components, max_nearest, _ = self._staged_connectivity_metrics(structure, movable)
        if components <= 1 and max_nearest <= 3.20:
            for idx_move in movable:
                structure.positions[idx_move, 2] = float(
                    np.clip(structure.positions[idx_move, 2], surface_z + 0.95, surface_z + 2.40)
                )
            return

        base_xy_radius = 2.0
        first = movable[0]
        first_pos = np.array(
            [
                float(center[0] + radial[0] / radial_norm * base_xy_radius),
                float(center[1] + radial[1] / radial_norm * base_xy_radius),
                float(np.clip(local_base[2], surface_z + 0.95, surface_z + 2.40)),
            ],
            dtype=float,
        )
        structure.positions[first] = first_pos

        prev_idx = first
        prev_pos = first_pos.copy()
        for order, idx_move in enumerate(movable[1:], start=1):
            bond_len = max(
                0.95,
                self._minimum_distance(symbols[prev_idx], symbols[idx_move], scale=0.95) + 0.10,
            )
            theta = float((np.pi / 2.5) * order)
            direction = np.array([np.cos(theta), np.sin(theta), 0.22], dtype=float)
            direction = direction / float(np.linalg.norm(direction))

            staged = prev_pos + direction * bond_len
            staged[2] = float(np.clip(staged[2], surface_z + 0.95, surface_z + 2.40))

            structure.positions[idx_move] = staged
            prev_idx = idx_move
            prev_pos = staged

        self._relax_close_contacts(structure, movable, max_iter=14)
        for idx_move in movable:
            structure.positions[idx_move, 2] = float(
                np.clip(structure.positions[idx_move, 2], surface_z + 0.95, surface_z + 2.40)
            )

    def _best_interpolation_mapping_metrics(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_ads: List[int],
        product_ads: List[int],
        n_images: int = 7,
    ) -> Tuple[float, int, float, int]:
        candidates: List[Tuple[float, int, float, int]] = []
        for reorder_product_atoms in (False, True):
            try:
                prepared = prepare_endpoint_pair_for_interpolation(
                    reactant=reactant,
                    product=product,
                    reactant_adsorbate_indices=reactant_ads,
                    product_adsorbate_indices=product_ads,
                    reorder_product_atoms=reorder_product_atoms,
                )
                prepared_ads = sorted(
                    {
                        idx
                        for idx in (
                            list(prepared.reactant_adsorbate_indices)
                            + list(prepared.product_adsorbate_indices)
                        )
                        if 0 <= idx < len(prepared.reactant)
                    }
                )
                candidates.append(
                    interpolated_path_min_distances(
                        reactant=prepared.reactant,
                        product=prepared.product,
                        adsorbate_indices=prepared_ads,
                        n_images=n_images,
                    )
                )
            except Exception:
                continue

        if not candidates:
            raw_ads = sorted({i for i in (list(reactant_ads) + list(product_ads)) if 0 <= i < len(reactant)})
            candidates.append(
                self._interpolated_path_min_distances(
                    reactant=reactant.copy(),
                    product=product.copy(),
                    adsorbate_indices=raw_ads,
                    n_images=n_images,
                )
            )

        def _score(item: Tuple[float, int, float, int]) -> Tuple[float, float]:
            min_pair, _, min_ads_surf, _ = item
            pair_score = float(min_pair if min_pair > 0 else -1.0)
            surf_score = float(min_ads_surf if min_ads_surf > 0 else -1.0)
            return pair_score, surf_score

        return max(candidates, key=_score)

    @staticmethod
    def _canonical_element(symbol: str) -> str:
        s = str(symbol or "").strip()
        if not s:
            return ""
        return s[0].upper() + (s[1:].lower() if len(s) > 1 else "")

    def _extract_species_tokens(self, species: str) -> List[str]:
        raw = str(species or "").strip()
        if not raw:
            return []

        token_symbols: List[str] = []
        primary = self._extract_primary_adsorbate_species(raw)
        source = primary or raw

        # direct element
        canon = self._canonical_element(source)
        if canon in atomic_numbers:
            return [canon]

        # parse formula-like token
        formula = source.replace("*", "").replace("(", "").replace(")", "")
        for elem, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
            symbol = self._canonical_element(elem)
            if symbol in atomic_numbers:
                n = int(count) if count else 1
                token_symbols.extend([symbol] * max(n, 1))

        return token_symbols

    def _infer_addition_elements(
        self,
        step: PathwayStepSpec,
        reactant: Atoms,
        product: Atoms,
    ) -> List[str]:
        inferred: List[str] = []

        # First try explicit species labels from LLM output
        for atom_spec in step.atoms_to_add:
            tokens = self._resolve_addition_elements(atom_spec)
            inferred.extend(tokens)

        if inferred:
            return inferred

        # Generic fallback: infer by element-count delta between product and reactant
        reactant_counts = Counter(reactant.get_chemical_symbols())
        product_counts = Counter(product.get_chemical_symbols())
        all_elements = set(reactant_counts.keys()) | set(product_counts.keys())

        for element in sorted(all_elements):
            delta = product_counts.get(element, 0) - reactant_counts.get(element, 0)
            if delta > 0:
                inferred.extend([element] * delta)

        return inferred

    def _select_formula_core_adsorbate_indices(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        formula: str,
    ) -> List[int]:
        if not adsorbate_indices:
            return []

        expected_tokens = self._extract_species_tokens(formula)
        if not expected_tokens:
            return list(adsorbate_indices)

        expected_counts = Counter(expected_tokens)
        symbols = structure.get_chemical_symbols()
        positions = structure.get_positions()
        centroid = np.mean(positions[adsorbate_indices], axis=0)

        selected: List[int] = []
        used: set[int] = set()

        for elem, count in expected_counts.items():
            candidates = [idx for idx in adsorbate_indices if symbols[idx] == elem and idx not in used]
            if not candidates:
                continue
            candidates.sort(key=lambda idx: float(np.linalg.norm(positions[idx] - centroid)))
            keep = candidates[:count]
            selected.extend(keep)
            used.update(keep)

        return selected if selected else list(adsorbate_indices)

    def _ensure_formula_coverage_with_staged_atoms(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        staged_indices: List[int],
        formula: str,
        fixed_count: int,
    ) -> List[int]:
        """Ensure formula-required adsorbate elements exist, using staged atoms only when needed."""
        valid_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(structure)]
        if not valid_ads:
            return []

        expected_tokens = self._extract_species_tokens(formula)
        if self._contains_gas_phase_hint(formula) or (not expected_tokens):
            return self._select_formula_core_adsorbate_indices(structure, valid_ads, formula)

        expected_counts = Counter(expected_tokens)
        core_ads = self._select_formula_core_adsorbate_indices(structure, valid_ads, formula)
        core_set = set(core_ads)
        symbols = structure.get_chemical_symbols()
        observed_counts = Counter(symbols[idx] for idx in core_ads if 0 <= idx < len(structure))

        for elem in sorted(expected_counts.keys()):
            missing = expected_counts[elem] - observed_counts.get(elem, 0)
            if missing <= 0:
                continue

            same_elem_spare = [
                idx for idx in valid_ads
                if symbols[idx] == elem and idx not in core_set
            ]
            for idx in same_elem_spare:
                core_ads.append(idx)
                core_set.add(idx)
                observed_counts[elem] += 1
                missing -= 1
                if missing <= 0:
                    break

            for local_offset in range(max(0, missing)):
                staged = np.array(
                    self._get_local_staging_position(
                        structure=structure,
                        adsorbate_indices=valid_ads,
                        core_adsorbate_indices=core_ads,
                        offset=len(staged_indices) + local_offset,
                    ),
                    dtype=float,
                )
                surface_z = self._estimate_surface_top_z(structure, valid_ads)
                staged[2] = float(np.clip(staged[2], surface_z + 0.95, surface_z + 3.2))
                structure.append(Atom(elem, position=staged.tolist()))
                new_idx = len(structure) - 1
                valid_ads.append(new_idx)
                adsorbate_indices.append(new_idx)
                staged_indices.append(new_idx)
                core_ads.append(new_idx)
                core_set.add(new_idx)
                observed_counts[elem] += 1

        movable = sorted(
            {
                idx
                for idx in staged_indices
                if isinstance(idx, int) and idx >= max(0, fixed_count) and 0 <= idx < len(structure)
            }
        )
        if movable:
            self._nudge_adsorbate_set_outward(
                structure=structure,
                move_indices=movable,
                reference_indices=list(core_ads),
                adsorbate_indices=list(valid_ads),
                radial_base=0.95,
                min_sep=1.00,
                min_height=0.85,
                max_height=3.2,
            )
            self._relax_close_contacts(structure, movable, max_iter=12)

        return self._select_formula_core_adsorbate_indices(structure, list(valid_ads), formula)

    def _agent4_postprocess_step_structures(
        self,
        state: WorkflowState,
        base_structure: Atoms,
        step_structures: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        _ = state, base_structure
        return self.tools.postprocess_step_structures(self, step_structures)

    def _agent5_programmatic_validation(self, step_structures: List[Dict[str, Any]]) -> Dict[str, Any]:
        return self.tools.programmatic_validation(self, step_structures)

    def _suggest_cross_endpoint_position(
        self,
        source: Atoms,
        source_adsorbate_indices: List[int],
        element: str,
        prefer_index_at_least: int,
    ) -> Optional[List[float]]:
        candidates = [
            idx
            for idx in source_adsorbate_indices
            if 0 <= idx < len(source) and source[idx].symbol == element
        ]
        if not candidates:
            return None

        positions = source.get_positions()
        anchor_pool = [idx for idx in source_adsorbate_indices if 0 <= idx < len(source)]
        if anchor_pool:
            anchor = np.mean(positions[anchor_pool], axis=0)
        else:
            anchor = np.mean(positions, axis=0)

        def _rank(idx: int) -> Tuple[int, float]:
            prefer_rank = 0 if idx >= max(0, int(prefer_index_at_least)) else 1
            dist_rank = float(np.linalg.norm(positions[idx] - anchor))
            return (prefer_rank, dist_rank)

        best_idx = sorted(candidates, key=_rank)[0]
        best_pos = positions[best_idx]
        return [float(best_pos[0]), float(best_pos[1]), float(best_pos[2])]

    def _apply_step_ops_on_fixed_product(
        self,
        current_physical: Atoms,
        fixed_product: Atoms,
        step: PathwayStepSpec,
    ) -> Dict[str, Any]:
        return self.tools.apply_step_ops_on_fixed_product(
            workflow=self,
            current_physical=current_physical,
            fixed_product=fixed_product,
            step=step,
        )

    def _build_step_structures(
        self,
        state: WorkflowState,
        base_structure: Atoms,
        steps: List[PathwayStepSpec],
        preplaced_products: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> List[Dict[str, Any]]:
        if preplaced_products is None:
            raise RuntimeError("Tool preplaced products are required for Agent4 endpoint construction.")
        return self.tools.build_step_structures(
            workflow=self,
            state=state,
            base_structure=base_structure,
            steps=steps,
            preplaced_products=preplaced_products,
        )
