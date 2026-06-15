"""Species parsing, index management, structure serialization, and reaction context normalization."""

from __future__ import annotations

import hashlib
import re
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase import Atom, Atoms
from ase.io import read, write

from camel_agents.schemas import (
    AtomAddSpec,
    PathwayDesign,
    PathwayStepSpec,
    ReactionContext,
    WorkflowState,
)

from .common import logger

class WorkflowContextMixin:
    def _set_surface_indices(self, atoms: Atoms, surface_indices: List[int], adsorbate_indices: List[int]) -> None:
        atoms.info["surface_indices"] = list(surface_indices)
        atoms.info["adsorbate_indices"] = list(adsorbate_indices)

    def _make_endpoint(self, atoms: Atoms, surface_indices: List[int], adsorbate_indices: List[int], staged_indices: Optional[List[int]] = None) -> "StructureEndpoint":
        from camel_agents.surface_model import StructureEndpoint
        return StructureEndpoint.from_parts(atoms, surface_indices, adsorbate_indices, staged_indices)

    def _update_adsorbate_indices_from_count(self, atoms: Atoms, adsorbate_count: int) -> List[int]:
        if adsorbate_count <= 0:
            adsorbate_indices: List[int] = []
        else:
            adsorbate_indices = list(range(len(atoms) - adsorbate_count, len(atoms)))
        surface_indices = [i for i in range(len(atoms)) if i not in adsorbate_indices]
        self._set_surface_indices(atoms, surface_indices, adsorbate_indices)
        return adsorbate_indices

    def _infer_adsorbate_indices_after_reconstruction(
        self,
        reference_atoms: Optional[Atoms],
        reference_adsorbate_indices: Optional[List[int]],
        reconstructed_atoms: Atoms,
        fallback_surface_count: Optional[int] = None,
    ) -> Tuple[List[int], List[int]]:
        """Infer adsorbate indices after MC reconstruction using reference mapping first.

        Why:
        - Reconstruction can append/reorder atoms (e.g., added slab atoms),
          so tail-based indexing is unsafe.
        - We first map previous adsorbate atoms to reconstructed atoms by
          element + nearest-distance, then fallback to generic heuristics.
        """
        expected_count = len(reference_adsorbate_indices or [])
        mapped_adsorbate: List[int] = []

        if (
            reference_atoms is not None
            and reference_adsorbate_indices
            and len(reference_atoms) > 0
            and len(reconstructed_atoms) > 0
        ):
            all_reconstructed_indices = list(range(len(reconstructed_atoms)))
            pairs = self._match_adsorbate_atoms_by_element(
                reference_atoms,
                reconstructed_atoms,
                list(reference_adsorbate_indices),
                all_reconstructed_indices,
            )

            pair_by_reference: Dict[int, int] = {int(r_idx): int(p_idx) for r_idx, p_idx in pairs}
            ordered = [pair_by_reference[r_idx] for r_idx in reference_adsorbate_indices if r_idx in pair_by_reference]

            seen: set[int] = set()
            mapped_adsorbate = []
            for idx in ordered:
                if idx not in seen:
                    mapped_adsorbate.append(idx)
                    seen.add(idx)

        if expected_count > 0 and len(mapped_adsorbate) < expected_count:
            fallback_surface, fallback_adsorbate = self._infer_surface_adsorbate_indices(
                reconstructed_atoms,
                fallback_surface_count=fallback_surface_count,
            )
            _ = fallback_surface
            for idx in fallback_adsorbate:
                if idx not in mapped_adsorbate:
                    mapped_adsorbate.append(idx)
                if len(mapped_adsorbate) >= expected_count:
                    break

        if not mapped_adsorbate:
            return self._infer_surface_adsorbate_indices(
                reconstructed_atoms,
                fallback_surface_count=fallback_surface_count,
            )

        if expected_count > 0:
            mapped_adsorbate = mapped_adsorbate[:expected_count]

        surface_indices = [i for i in range(len(reconstructed_atoms)) if i not in set(mapped_adsorbate)]
        return surface_indices, mapped_adsorbate

    def _infer_surface_adsorbate_indices(
        self,
        atoms: Atoms,
        fallback_surface_count: Optional[int] = None,
    ) -> Tuple[List[int], List[int]]:
        tags = atoms.get_tags() if len(atoms) else np.array([])
        if len(tags) == len(atoms) and len(np.unique(tags)) > 1:
            max_tag = int(np.max(tags))
            candidate_ads = [i for i, tag in enumerate(tags.tolist()) if int(tag) == max_tag]
            if 0 < len(candidate_ads) < len(atoms):
                surface = [i for i in range(len(atoms)) if i not in set(candidate_ads)]
                return surface, candidate_ads

        symbols = atoms.get_chemical_symbols()
        if not symbols:
            return [], []

        dominant_symbol = Counter(symbols).most_common(1)[0][0]
        positions = atoms.get_positions()

        dominant_indices = [i for i, s in enumerate(symbols) if s == dominant_symbol]
        if dominant_indices:
            surface_z = float(np.max(positions[dominant_indices, 2]))
        else:
            surface_z = float(np.max(positions[:, 2]))

        adsorbate = [
            i for i, s in enumerate(symbols)
            if s != dominant_symbol and positions[i, 2] >= surface_z - 0.6
        ]

        if not adsorbate:
            adsorbate = [i for i, s in enumerate(symbols) if s != dominant_symbol]

        if (not adsorbate) and fallback_surface_count is not None and 0 < fallback_surface_count < len(atoms):
            adsorbate = list(range(fallback_surface_count, len(atoms)))

        surface = [i for i in range(len(atoms)) if i not in set(adsorbate)]
        return surface, adsorbate

    def _ensure_clean_slab_consistency(
        self,
        state: WorkflowState,
        ads_atoms: Atoms,
        surface_indices: List[int],
        adsorption_dir: Path,
    ) -> None:
        surface = self._get_surface() if hasattr(self, "_get_surface") else None
        clean_path = getattr(surface, "clean_slab_path", None) if surface is not None else None
        if not clean_path:
            return

        try:
            clean_atoms = read(clean_path)
            if len(clean_atoms) == len(surface_indices):
                return
        except Exception:
            pass

        if not surface_indices:
            return

        derived_clean = ads_atoms[surface_indices]
        derived_path = adsorption_dir / "derived_clean_slab.vasp"
        write(derived_path, derived_clean)
        if surface is not None:
            surface.clean_slab_path = str(derived_path)

    def _format_adsorbate_info(self, atoms: Atoms, adsorbate_indices: List[int]) -> str:
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        lines = []
        for idx in adsorbate_indices:
            pos = positions[idx]
            lines.append(f"{idx}: {symbols[idx]} at ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
        return "\n".join(lines) if lines else "(none)"

    @staticmethod
    def _read_text_file(path_text: str) -> str:
        path = Path(str(path_text or "").strip())
        if not path_text or not path.exists() or not path.is_file():
            return "(not provided)"
        try:
            return path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            return path.read_text(errors="ignore")
        except Exception as exc:
            return f"(failed to read file: {exc})"

    def _format_atom_coordinates_table(
        self,
        atoms: Atoms,
        adsorbate_indices: List[int],
        selected_indices: Optional[List[int]] = None,
    ) -> str:
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        ads_set = set(adsorbate_indices or [])
        if selected_indices is None:
            indices = list(range(len(atoms)))
        else:
            indices = [idx for idx in selected_indices if 0 <= idx < len(atoms)]
        lines = ["index symbol x y z role"]
        for idx in indices:
            sym = symbols[idx]
            pos = positions[idx]
            role = "adsorbate" if idx in ads_set else "surface"
            lines.append(
                f"{idx:4d} {sym:>3s} {float(pos[0]): .8f} {float(pos[1]): .8f} {float(pos[2]): .8f} {role}"
            )
        return "\n".join(lines)

    def _fixed_coordinate_signature(
        self,
        atoms: Atoms,
        fixed_count: int,
    ) -> str:
        n = max(0, min(int(fixed_count), len(atoms)))
        if n <= 0:
            return "none"

        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        digest_lines = [
            f"{idx}:{symbols[idx]}:{float(positions[idx, 0]):.10f}:{float(positions[idx, 1]):.10f}:{float(positions[idx, 2]):.10f}"
            for idx in range(n)
        ]
        payload = "\n".join(digest_lines).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()[:20]

    def _format_fixed_coordinate_block(
        self,
        atoms: Atoms,
        adsorbate_indices: List[int],
        fixed_count: int,
        title: str,
    ) -> str:
        n = max(0, min(int(fixed_count), len(atoms)))
        selected = list(range(n))
        signature = self._fixed_coordinate_signature(atoms, n)
        coord_text = self._format_atom_coordinates_table(
            atoms,
            adsorbate_indices,
            selected_indices=selected,
        )
        return (
            f"[{title}]\n"
            f"fixed_atom_count={n}\n"
            f"immutable_signature={signature}\n"
            f"## FIXED_COORDINATES\n"
            f"{coord_text}"
        )

    def _atoms_to_poscar_text(self, atoms: Atoms) -> str:
        buffer = StringIO()
        write(buffer, atoms, format="vasp", direct=True, sort=False, vasp5=True)
        return buffer.getvalue()

    def _atoms_to_structure_text(
        self,
        atoms: Atoms,
        title: str,
        adsorbate_indices: Optional[List[int]] = None,
        selected_indices: Optional[List[int]] = None,
    ) -> str:
        ads_indices = [idx for idx in (adsorbate_indices or []) if 0 <= idx < len(atoms)]
        selected = [idx for idx in (selected_indices or []) if 0 <= idx < len(atoms)]

        serialized_atoms = atoms
        header_lines = [f"[{title}]"]
        if selected:
            serialized_atoms = atoms[selected].copy()
            header_lines.append(f"selected_original_indices={selected}")
            selected_syms = [atoms[idx].symbol for idx in selected]
            header_lines.append(f"selected_symbols={selected_syms}")
        elif ads_indices:
            header_lines.append(f"adsorbate_indices={ads_indices}")

        poscar_text = self._atoms_to_poscar_text(serialized_atoms)
        header_lines.append(f"## POSCAR\n{poscar_text}")
        return "\n".join(header_lines)

    def _build_agent5_structures_payload(self, state: WorkflowState) -> str:
        return self.tools.build_steps_payload(self, state.step_structures)

    def _build_clean_surface_from_known_indices(self, base_structure: Atoms, state: WorkflowState) -> Atoms:
        """Build clean slab using upstream-known index partition (no re-detection)."""
        n = len(base_structure)
        surface = self._get_surface() if hasattr(self, "_get_surface") else None
        surf_src = list(surface.surface_indices) if surface is not None else []
        ads_src = list(surface.adsorbate_indices) if surface is not None else []
        known_surface = [i for i in surf_src if 0 <= i < n]
        known_ads = [i for i in ads_src if 0 <= i < n]

        if known_surface:
            clean = base_structure[known_surface].copy()
        elif known_ads:
            clean = base_structure.copy()
            for idx in sorted(set(known_ads), reverse=True):
                if 0 <= idx < len(clean):
                    del clean[idx]
        else:
            clean = base_structure.copy()

        clean.info["surface_indices"] = list(range(len(clean)))
        clean.info["adsorbate_indices"] = []
        return clean

    def _estimate_anchor_from_known_adsorbate(
        self,
        base_structure: Atoms,
        adsorbate_indices: List[int],
    ) -> Optional[List[float]]:
        valid = [idx for idx in (adsorbate_indices or []) if 0 <= idx < len(base_structure)]
        if not valid:
            return None

        positions = base_structure.get_positions()[valid]
        anchor_local_idx = int(np.argmin(positions[:, 2]))
        anchor = positions[anchor_local_idx]
        return [float(anchor[0]), float(anchor[1]), float(anchor[2])]

    def _species_safe_name(self, species: str) -> str:
        base = self._canonical_species_label(species) or str(species or "").strip()
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", base)
        return safe or "species"

    def _agent4_tool_preplace_products(
        self,
        state: WorkflowState,
        base_structure: Atoms,
        iteration: int,
    ) -> Dict[str, Dict[str, Any]]:
        """Agent4 phase-1: tool-based fixed-site product placement for intermediates."""
        intermediates = self._get_working_intermediate_sequence(state)
        if not intermediates:
            return {}

        placed: Dict[str, Dict[str, Any]] = {}

        surface = self._get_surface() if hasattr(self, "_get_surface") else None
        base_copy = base_structure.copy()
        surf_idx0 = list(surface.surface_indices) if surface is not None else []
        ads_idx0 = list(surface.adsorbate_indices) if surface is not None else []
        self._set_surface_indices(base_copy, surf_idx0, ads_idx0)

        # If first intermediate is a co-adsorption label, add secondary species
        first_label = intermediates[0]
        _, first_secondaries = self._parse_coadsorption_label(first_label)
        if first_secondaries:
            primary_ads_indices = list(base_copy.info.get("adsorbate_indices", []))
            # Fallback: if adsorbate_indices empty, infer from surface atom count
            n_surface_atoms = surface.n_surface_atoms if surface is not None else 0
            if not primary_ads_indices and n_surface_atoms:
                primary_ads_indices = list(range(int(n_surface_atoms), len(base_copy)))
            for sec_formula, sec_count in first_secondaries:
                for i_sec in range(sec_count):
                    site = self._find_nearest_neighbor_site(
                        base_copy, primary_ads_indices, offset=i_sec,
                    )
                    if site is not None:
                        sec_atoms = self._formula_to_atom_list(sec_formula)
                        base_pos = np.array(site, dtype=float)
                        if len(sec_atoms) == 1:
                            base_copy.append(Atom(sec_atoms[0], position=base_pos))
                        else:
                            base_copy.append(Atom(sec_atoms[0], position=base_pos))
                            for j, elem in enumerate(sec_atoms[1:], 1):
                                angle = 2.0 * np.pi * j / len(sec_atoms)
                                bond_len = 0.97 if elem == "H" else 1.20
                                offset_vec = np.array([
                                    bond_len * np.sin(angle),
                                    bond_len * np.cos(angle),
                                    0.0,
                                ])
                                base_copy.append(Atom(elem, position=base_pos + offset_vec))
            # Update indices after adding co-adsorbates
            surface_count = len(surf_idx0)
            surface_indices = list(range(surface_count))
            adsorbate_indices = list(range(surface_count, len(base_copy)))
            self._set_surface_indices(base_copy, surface_indices, adsorbate_indices)

        first_key = self._canonical_species_label(first_label)
        placed[first_key] = {
            "label": first_label,
            "structure": base_copy,
            "source": "agent3_input",
            "path": (getattr(surface, "reconstructed_surface_path", None) or "") if surface is not None else "",
            "surface_indices": list(base_copy.info.get("surface_indices", [])),
            "adsorbate_indices": list(base_copy.info.get("adsorbate_indices", [])),
            "fixed_atom_count": int(len(base_copy)),
            "immutable_signature": self._fixed_coordinate_signature(base_copy, len(base_copy)),
        }

        if len(intermediates) <= 1:
            return placed

        iter_dir = Path(state.output_base_dir) / state.run_id / "04_pathway" / f"iter_{iteration:02d}_tool_preplacement"
        iter_dir.mkdir(parents=True, exist_ok=True)

        # Separate primary species from co-adsorption labels for tool placement.
        # E.g. "*CO+*O" → primary "*CO" goes to adsorption tool; "*O" placed later.
        tool_intermediates: List[str] = []
        coadsorption_map: Dict[str, List[Tuple[str, int]]] = {}  # label → [(formula, count)]
        for species in intermediates:
            if self._contains_gas_phase_hint(species):
                continue
            primary, secondaries = self._parse_coadsorption_label(species)
            if secondaries:
                # Feed only the primary species to the adsorption tool
                primary_label = f"*{primary}"
                tool_intermediates.append(primary_label)
                coadsorption_map[species] = secondaries
            else:
                tool_intermediates.append(species)

        if self.tools is not None and len(tool_intermediates) > 1:
            clean_surface = self._build_clean_surface_from_known_indices(base_structure, state)
            clean_surface_path = iter_dir / "clean_surface_known_indices.vasp"
            write(clean_surface_path, clean_surface)
            surface_atom_count = len(clean_surface)

            # Get UMA predictor for sequential relax
            predictor = None
            try:
                predictor = self.tools.get_shared_fairchem_predictor(
                    cache_key="uma_shared", model_name="uma-s-1p1",
                    use_gpu=True, device="cuda",
                    work_subdir="_shared_fairchem_global",
                    keep_files=False, verbose=False,
                )
                predictor._load_model()
            except Exception as exc:
                logger.warning("Cannot get FairchemPredictor for sequential preplacement: %s", exc)

            # Sequential construction: each intermediate from previous relaxed structure
            prev_structure = base_copy.copy()
            prev_ads_indices = list(base_copy.info.get("adsorbate_indices", []))

            for label in intermediates[1:]:
                if self._contains_gas_phase_hint(label):
                    continue
                key = self._canonical_species_label(label)

                # Parse target adsorbate formula
                primary_label = label
                if label in coadsorption_map:
                    primary, _ = self._parse_coadsorption_label(label)
                    primary_label = f"*{primary}"

                # Determine element delta between current adsorbate and target
                prev_ads_symbols = [prev_structure.get_chemical_symbols()[i]
                                    for i in prev_ads_indices if 0 <= i < len(prev_structure)]
                target_atoms_list = self._formula_to_atom_list(
                    re.sub(r"^\*|\*$", "", self._extract_primary_adsorbate_species(primary_label) or primary_label)
                )
                prev_counts = Counter(prev_ads_symbols)
                target_counts = Counter(target_atoms_list)

                # Build new structure from previous
                new_structure = prev_structure.copy()
                new_ads_indices = list(prev_ads_indices)

                # Compute element delta and modify structure
                all_pos = new_structure.get_positions()
                all_sym = new_structure.get_chemical_symbols()
                ads_positions = all_pos[new_ads_indices] if new_ads_indices else np.zeros((0, 3))
                ads_centroid = np.mean(ads_positions, axis=0) if len(ads_positions) > 0 else all_pos.mean(axis=0)

                # First pass: remove excess atoms (furthest from centroid, matching element)
                for elem in sorted(set(list(prev_counts.keys()) + list(target_counts.keys()))):
                    delta = target_counts.get(elem, 0) - prev_counts.get(elem, 0)
                    if delta < 0:
                        removable = [i for i in new_ads_indices
                                     if 0 <= i < len(new_structure) and new_structure[i].symbol == elem]
                        if removable:
                            removable.sort(key=lambda i: -np.linalg.norm(
                                new_structure.positions[i] - ads_centroid))
                            for i_del in range(min(-delta, len(removable))):
                                idx = removable[i_del]
                                del new_structure[idx]
                                new_ads_indices = [j if j < idx else j - 1
                                                   for j in new_ads_indices if j != idx]

                # Refresh positions after deletions
                all_pos = new_structure.get_positions()
                all_sym = new_structure.get_chemical_symbols()
                ads_positions = all_pos[new_ads_indices] if new_ads_indices else np.zeros((0, 3))
                ads_centroid = np.mean(ads_positions, axis=0) if len(ads_positions) > 0 else all_pos.mean(axis=0)

                # Second pass: add missing atoms near the nearest heavy adsorbate atom
                for elem in sorted(set(list(prev_counts.keys()) + list(target_counts.keys()))):
                    delta = target_counts.get(elem, 0) - prev_counts.get(elem, 0)
                    if delta > 0:
                        # Find the nearest heavy (non-H) adsorbate atom as anchor
                        heavy_ads = [i for i in new_ads_indices
                                     if 0 <= i < len(new_structure) and all_sym[i] != "H"]
                        if not heavy_ads:
                            heavy_ads = new_ads_indices[:1] if new_ads_indices else []
                        anchor_idx = heavy_ads[0] if heavy_ads else None
                        anchor_pos = all_pos[anchor_idx] if anchor_idx is not None else ads_centroid

                        for i_add in range(delta):
                            bond_len = 0.97 if elem == "H" else 1.43 if elem == "O" else 1.54
                            # Place at bond length from anchor, above the surface
                            angle = 2.0 * np.pi * i_add / max(delta, 1)
                            pos = anchor_pos + np.array([
                                bond_len * np.cos(angle) * 0.7,
                                bond_len * np.sin(angle) * 0.7,
                                bond_len * 0.7,
                            ])
                            new_structure.append(Atom(elem, position=pos))
                            new_ads_indices.append(len(new_structure) - 1)

                # Place co-adsorbed secondary species
                if label in coadsorption_map:
                    for sec_formula, sec_count in coadsorption_map[label]:
                        for i_sec in range(sec_count):
                            site = self._find_nearest_neighbor_site(
                                new_structure, new_ads_indices, offset=i_sec,
                            )
                            if site is not None:
                                sec_atoms = self._formula_to_atom_list(sec_formula)
                                base_pos = np.array(site, dtype=float)
                                for j, elem in enumerate(sec_atoms):
                                    if j == 0:
                                        new_structure.append(Atom(elem, position=base_pos))
                                    else:
                                        angle = 2.0 * np.pi * j / len(sec_atoms)
                                        bond_len = 0.97 if elem == "H" else 1.20
                                        new_structure.append(Atom(elem, position=base_pos + np.array([
                                            bond_len * np.sin(angle), bond_len * np.cos(angle), 0.0])))
                                    new_ads_indices.append(len(new_structure) - 1)

                surface_indices = [i for i in range(len(new_structure)) if i not in set(new_ads_indices)]
                self._set_surface_indices(new_structure, surface_indices, new_ads_indices)

                # Relax with fixed surface (only adsorbate + top 2 layers free)
                if predictor is not None:
                    try:
                        from ase.optimize import BFGS
                        from ase.constraints import FixAtoms
                        relax_struct = new_structure.copy()
                        relax_struct.pbc = [True, True, True]
                        relax_struct.calc = predictor._calculator

                        # Fix everything except adsorbate atoms
                        fixed = [i for i in range(len(relax_struct)) if i not in set(new_ads_indices)]
                        relax_struct.set_constraint(FixAtoms(indices=fixed))

                        opt = BFGS(relax_struct, logfile=str(iter_dir / f"{self._species_safe_name(str(label))}_relax.log"))
                        opt.run(fmax=0.05, steps=200)

                        new_structure = relax_struct.copy()
                        try:
                            new_structure.set_constraint()
                        except Exception:
                            pass
                        self._set_surface_indices(new_structure, surface_indices, new_ads_indices)
                        logger.info("Sequential preplacement: relaxed %s (fmax=%.3f)",
                                    label, np.max(np.linalg.norm(relax_struct.get_forces(), axis=1)))
                    except Exception as exc:
                        logger.warning("Sequential preplacement relax failed for %s: %s", label, exc)

                product_path = iter_dir / f"{self._species_safe_name(str(label))}_product.vasp"
                write(product_path, new_structure)

                placed[key] = {
                    "label": label,
                    "structure": new_structure,
                    "source": "sequential_from_reactant",
                    "path": str(product_path),
                    "surface_indices": surface_indices,
                    "adsorbate_indices": new_ads_indices,
                    "fixed_atom_count": int(len(new_structure)),
                    "immutable_signature": self._fixed_coordinate_signature(new_structure, len(new_structure)),
                }

                # Chain: next intermediate starts from this relaxed structure
                prev_structure = new_structure.copy()
                prev_ads_indices = list(new_ads_indices)

        required = [s for s in intermediates if not self._contains_gas_phase_hint(s)]
        missing = [s for s in required if self._canonical_species_label(s) not in placed]
        if missing:
            raise RuntimeError(f"Tool preplacement missing intermediates: {missing}")

        # Align consecutive intermediates: shift each product's adsorbate so that
        # its bottom-most adsorbate atom matches the previous step's bottom-most.
        # This ensures the surface-bonded anchor stays at consistent coordinates
        # across steps, preventing interpolation collisions in NEB.
        self._align_consecutive_adsorbates(intermediates, placed)

        return placed

    def _find_nearest_neighbor_site(
        self,
        structure: Atoms,
        primary_ads_indices: List[int],
        offset: int = 0,
        adsorption_height: float = 1.8,
        min_dist_from_ads: float = 2.0,
    ) -> Optional[List[float]]:
        """Find nearest-neighbor hollow/bridge site relative to primary adsorbate.

        Uses surface metal atoms near the primary adsorbate centroid to compute
        hollow site positions, then picks the closest one that does NOT overlap
        with any existing adsorbate atom (distance >= *min_dist_from_ads*).
        ``offset`` skips the first N valid sites for placing multiple co-adsorbates.
        """
        if not primary_ads_indices or len(structure) == 0:
            return None

        positions = structure.get_positions()
        symbols = structure.get_chemical_symbols()

        # Primary adsorbate centroid
        valid_ads = [i for i in primary_ads_indices if i < len(structure)]
        if not valid_ads:
            return None
        ads_positions = positions[valid_ads]
        centroid_xy = np.mean(ads_positions[:, :2], axis=0)

        # Surface metal atoms (top layer)
        ads_set = set(primary_ads_indices)
        surface_indices = [i for i in range(len(structure)) if i not in ads_set]
        if not surface_indices:
            return None

        metal_counts = Counter(symbols[i] for i in surface_indices)
        dominant = metal_counts.most_common(1)[0][0]
        metal_indices = [i for i in surface_indices if symbols[i] == dominant]

        metal_z = [positions[i][2] for i in metal_indices]
        if not metal_z:
            return None
        top_z = max(metal_z)
        top_layer = [i for i in metal_indices if positions[i][2] > top_z - 1.0]
        if len(top_layer) < 2:
            return None

        # Generate candidate hollow sites from all triplets of nearby top-layer atoms
        from itertools import combinations

        top_positions = np.array([positions[i] for i in top_layer])
        candidate_sites: List[Tuple[float, np.ndarray]] = []

        for combo in combinations(range(len(top_layer)), 3):
            tri = top_positions[list(combo)]
            # Only consider compact triangles (max edge < 3.5 Å, typical fcc nearest-neighbor)
            edges = [float(np.linalg.norm(tri[i] - tri[j])) for i, j in [(0, 1), (0, 2), (1, 2)]]
            if max(edges) > 3.5:
                continue
            site_xy = np.mean(tri[:, :2], axis=0)
            site_z = float(np.mean(tri[:, 2])) + adsorption_height
            site_pos = np.array([site_xy[0], site_xy[1], site_z])

            # Distance from adsorbate centroid (for sorting: prefer near but not under)
            d_centroid = float(np.linalg.norm(site_xy - centroid_xy))

            # Minimum distance from ALL adsorbate atoms (must not overlap)
            d_min_ads = float(min(np.linalg.norm(site_pos - ads_positions[j]) for j in range(len(ads_positions))))

            if d_min_ads >= min_dist_from_ads:
                candidate_sites.append((d_centroid, site_pos))

        if not candidate_sites:
            return None

        # Sort by distance from centroid (nearest neighbor first)
        candidate_sites.sort(key=lambda x: x[0])

        # Skip `offset` sites for placing multiple co-adsorbates at different positions
        idx = min(offset, len(candidate_sites) - 1)
        return [float(candidate_sites[idx][1][k]) for k in range(3)]

    @staticmethod
    def _formula_to_atom_list(formula: str) -> List[str]:
        """Parse a chemical formula into a list of element symbols.

        ``"OH"`` → ``["O", "H"]``, ``"H2"`` → ``["H", "H"]``, ``"O"`` → ``["O"]``.
        """
        atoms: List[str] = []
        for m in re.finditer(r"([A-Z][a-z]?)(\d*)", formula):
            elem = m.group(1)
            count = int(m.group(2)) if m.group(2) else 1
            atoms.extend([elem] * count)
        return atoms

    def _align_consecutive_adsorbates(
        self,
        intermediates: List[str],
        placed: Dict[str, Dict[str, Any]],
    ) -> None:
        """Align adsorbate positions across consecutive intermediates.

        For each pair of consecutive non-gas intermediates, shift the later
        one's adsorbate atoms so that its bottom-most adsorbate atom aligns
        with the previous intermediate's bottom-most adsorbate atom.  All other
        adsorbate atoms keep their relative coordinates (rigid shift).
        Surface atoms are never modified.

        This fixes the independent-relaxation problem where the same anchor atom
        (e.g. C bonded to surface) ends up at different positions in different
        intermediates.
        """
        import numpy as np

        prev_key = None
        for species in intermediates:
            if self._contains_gas_phase_hint(species):
                continue
            key = self._canonical_species_label(species)
            entry = placed.get(key)
            if entry is None or not isinstance(entry.get("structure"), Atoms):
                prev_key = key
                continue

            if prev_key is not None:
                prev_entry = placed.get(prev_key)
                if prev_entry is not None and isinstance(prev_entry.get("structure"), Atoms):
                    self._align_adsorbate_to_reference(prev_entry, entry)

            prev_key = key

    @staticmethod
    def _align_adsorbate_to_reference(
        ref_entry: Dict[str, Any],
        target_entry: Dict[str, Any],
    ) -> None:
        """Shift target adsorbate so its anchor atom aligns with ref's anchor atom.

        Anchor selection priority (general, not reaction-specific):
        1. Heaviest common non-H adsorbate element with lowest z (e.g., C in both CO2 and CO)
        2. Lowest-z adsorbate atom (fallback)
        """
        import numpy as np

        ref_struct: Atoms = ref_entry["structure"]
        target_struct: Atoms = target_entry["structure"]
        ref_ads = list(ref_entry.get("adsorbate_indices", []))
        target_ads = list(target_entry.get("adsorbate_indices", []))

        if not ref_ads or not target_ads:
            return

        ref_positions = ref_struct.get_positions()
        target_positions = target_struct.get_positions()
        ref_symbols = ref_struct.get_chemical_symbols()
        target_symbols = target_struct.get_chemical_symbols()

        ref_valid = [i for i in ref_ads if 0 <= i < len(ref_struct)]
        target_valid = [i for i in target_ads if 0 <= i < len(target_struct)]
        if not ref_valid or not target_valid:
            return

        # Find common non-H elements between ref and target adsorbates
        ref_heavy = {ref_symbols[i] for i in ref_valid if ref_symbols[i] != "H"}
        target_heavy = {target_symbols[i] for i in target_valid if target_symbols[i] != "H"}
        common_heavy = ref_heavy & target_heavy

        if common_heavy:
            # Pick heaviest common element (C > N > O > S by atomic number)
            from ase.data import atomic_numbers
            anchor_elem = max(common_heavy, key=lambda e: atomic_numbers.get(e, 0))
            ref_candidates = [i for i in ref_valid if ref_symbols[i] == anchor_elem]
            target_candidates = [i for i in target_valid if target_symbols[i] == anchor_elem]
            ref_bottom_idx = min(ref_candidates, key=lambda i: ref_positions[i, 2])
            target_bottom_idx = min(target_candidates, key=lambda i: target_positions[i, 2])
        else:
            # Fallback: lowest-z atom
            ref_bottom_idx = min(ref_valid, key=lambda i: ref_positions[i, 2])
            target_bottom_idx = min(target_valid, key=lambda i: target_positions[i, 2])

        ref_anchor = ref_positions[ref_bottom_idx]
        target_anchor = target_positions[target_bottom_idx]

        # Compute shift: move target anchor to ref anchor position
        shift = ref_anchor - target_anchor

        # Only shift if displacement is significant (> 0.1 Å)
        if np.linalg.norm(shift) < 0.1:
            return

        # Apply rigid shift to ALL target adsorbate atoms (surface stays fixed)
        new_positions = target_positions.copy()
        for idx in target_valid:
            new_positions[idx] = target_positions[idx] + shift
        target_struct.set_positions(new_positions)

        # Update signature since positions changed
        from .common import logger
        logger.info(
            "Aligned adsorbate anchor: shifted %d adsorbate atoms by [%.3f, %.3f, %.3f] Å",
            len(target_valid), shift[0], shift[1], shift[2],
        )

    def _summarize_preplaced_products(
        self,
        preplaced_products: Dict[str, Dict[str, Any]],
        intermediates: List[str],
    ) -> str:
        lines: List[str] = []
        for species in intermediates or []:
            key = self._canonical_species_label(species)
            entry = preplaced_products.get(key)
            if not entry:
                lines.append(f"{species}: missing")
                continue

            struct = entry.get("structure")
            ads_indices = list(entry.get("adsorbate_indices", []))
            centroid_text = "(n/a)"
            if isinstance(struct, Atoms) and ads_indices:
                centroid = np.mean(struct.get_positions()[ads_indices], axis=0)
                centroid_text = f"({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})"

            lines.append(
                f"{species}: source={entry.get('source', 'unknown')}; "
                f"path={entry.get('path', '')}; adsorbate_indices={ads_indices}; centroid={centroid_text}; "
                f"fixed_atom_count={entry.get('fixed_atom_count', len(struct) if isinstance(struct, Atoms) else 'n/a')}; "
                f"immutable_signature={entry.get('immutable_signature', self._fixed_coordinate_signature(struct, len(struct)) if isinstance(struct, Atoms) else 'n/a')}"
            )
            if isinstance(struct, Atoms):
                lines.append(
                    self._atoms_to_structure_text(
                        atoms=struct,
                        title=f"{species}_tool_preplaced_adsorbate_only",
                        adsorbate_indices=ads_indices,
                        selected_indices=ads_indices,
                    )
                )

        return "\n".join(lines) if lines else "(none)"

    def _get_preplaced_entry(
        self,
        preplaced_products: Optional[Dict[str, Dict[str, Any]]],
        species: str,
    ) -> Optional[Dict[str, Any]]:
        if not preplaced_products:
            return None
        key = self._canonical_species_label(species)
        if key in preplaced_products:
            return preplaced_products[key]
        raw = str(species or "").strip()
        return preplaced_products.get(raw)

    @staticmethod
    def _parse_coadsorption_label(label: str) -> Tuple[str, List[Tuple[str, int]]]:
        """Parse co-adsorption labels like ``*CO+*O``, ``*CO2+2*H``.

        Returns (primary_formula, [(secondary_formula, count), ...]).
        For plain labels like ``*CO`` returns (``CO``, []).
        """
        raw = str(label).strip().replace(" ", "")
        # Split on '+' but only between starred species
        parts = re.split(r"\+(?=\d*\*)", raw)
        if len(parts) <= 1:
            # No co-adsorption
            cleaned = re.sub(r"^\*+", "", raw)
            cleaned = re.sub(r"\([^\)]*\)$", "", cleaned)
            return cleaned.upper(), []

        primary = re.sub(r"^\*+", "", parts[0])
        primary = re.sub(r"\([^\)]*\)$", "", primary).upper()

        secondaries: List[Tuple[str, int]] = []
        for part in parts[1:]:
            m = re.match(r"(\d*)\*?([A-Za-z0-9]+)", part)
            if m:
                count = int(m.group(1)) if m.group(1) else 1
                formula = m.group(2).upper()
                secondaries.append((formula, count))
        return primary, secondaries

    @staticmethod
    def _is_coadsorption_label(label: str) -> bool:
        """Return True if label contains co-adsorption notation (``+*``)."""
        raw = str(label).strip().replace(" ", "")
        return bool(re.search(r"\+\d*\*", raw))

    @staticmethod
    def _canonical_species_label(label: Optional[str]) -> str:
        if not label:
            return ""
        raw = str(label).strip()
        is_gas_phase = bool(
            re.search(r"\([A-Za-z0-9]*g\)", raw, flags=re.IGNORECASE)
            or re.search(r"\bgas\b", raw, flags=re.IGNORECASE)
        )
        compact = raw.replace(" ", "")

        # Preserve co-adsorption labels: *CO+*O → "CO+O", *CO2+2*H → "CO2+2H"
        if re.search(r"\+\d*\*", compact):
            parts = re.split(r"\+(?=\d*\*)", compact)
            normalized_parts = []
            for part in parts:
                cleaned = re.sub(r"\*+", "", part)  # remove all * (including mid-string like 2*H)
                cleaned = re.sub(r"\([^\)]*\)$", "", cleaned)
                normalized_parts.append(cleaned.upper())
            return "+".join(normalized_parts)

        starred = re.findall(r"\*[A-Za-z0-9]+(?:\*[A-Za-z0-9]+)?", compact)
        if starred:
            compact = starred[0]
        else:
            for sep in ("->", "=", ",", ";", "→", "to"):
                if sep in raw:
                    compact = raw.split(sep, 1)[0].strip().replace(" ", "")
                    break

        cleaned = re.sub(r"^\*+", "", compact)
        cleaned = re.sub(r"\([^\)]*\)$", "", cleaned)
        normalized = cleaned.upper()
        if is_gas_phase and normalized:
            return f"GAS:{normalized}"
        return normalized

    @staticmethod
    def _extract_primary_adsorbate_species(text: Optional[str]) -> str:
        if not text:
            return ""
        raw = str(text).strip()
        compact = raw.replace(" ", "")

        starred = re.findall(r"\*[A-Za-z0-9]+(?:\*[A-Za-z0-9]+)?", compact)
        if starred:
            return starred[0]

        for sep in ("->", "=", "+", ",", ";", "→", "to"):
            if sep in raw:
                head = raw.split(sep, 1)[0].strip()
                if head:
                    return head
        tokens = raw.split()
        return tokens[0] if tokens else raw

    @staticmethod
    def _ensure_star_notation(species: str) -> str:
        label = str(species or "").strip()
        if not label:
            return label
        if re.search(r"\([A-Za-z0-9]*g\)", label, flags=re.IGNORECASE):
            return label
        if re.search(r"\bgas\b", label, flags=re.IGNORECASE):
            return label
        if label == "*" or label.startswith("*"):
            return label
        return f"*{label}"

    def _extract_intermediate_sequence_from_text(self, description: str) -> List[str]:
        raw = str(description or "")
        if not raw:
            return []

        normalized = raw.replace("→", "->").replace("<->", "->")

        _PURE_H_RE = re.compile(r"^\*?\d*H\d*\*?$")
        _BARE_SITE = {"*", "2*", "3*"}

        elementary_steps = re.split(r";\s*", normalized)
        step_species: List[Tuple[List[str], List[str]]] = []
        for step_text in elementary_steps:
            sides = re.split(r"\s*->\s*", step_text.strip())
            if len(sides) != 2:
                continue
            lhs_tokens = [t.strip() for t in re.split(r"\s*\+\s*", sides[0]) if t.strip()]
            rhs_tokens = [t.strip() for t in re.split(r"\s*\+\s*", sides[1]) if t.strip()]

            def _extract_species(tokens: List[str]) -> List[str]:
                result = []
                for t in tokens:
                    if t in _BARE_SITE:
                        continue
                    if _PURE_H_RE.match(t):
                        continue
                    if "*" in t or re.match(r"^[A-Z][A-Za-z0-9]*$", t):
                        result.append(t)
                return result

            lhs_sp = _extract_species(lhs_tokens)
            rhs_sp = _extract_species(rhs_tokens)
            if lhs_sp or rhs_sp:
                step_species.append((lhs_sp, rhs_sp))

        if len(step_species) >= 2:
            sequence: List[str] = []
            for i, (lhs_sp, rhs_sp) in enumerate(step_species):
                main_reactant = lhs_sp[0] if lhs_sp else ""
                if i == 0 and main_reactant and main_reactant not in sequence:
                    sequence.append(main_reactant)
                main_product = rhs_sp[0] if rhs_sp else ""
                if main_product and main_product not in sequence:
                    sequence.append(main_product)
            cleaned: List[str] = []
            for token in sequence:
                species = self._extract_primary_adsorbate_species(token)
                if not species:
                    continue
                cleaned.append(self._ensure_star_notation(species))
            deduped: List[str] = []
            prev_key = ""
            for item in cleaned:
                key = self._canonical_species_label(item)
                if key and key == prev_key:
                    continue
                deduped.append(item)
                prev_key = key
            if len(deduped) >= 2:
                return deduped

        chain_candidates = re.findall(r"([*A-Za-z0-9()\-]+(?:\s*->\s*[*A-Za-z0-9()\-]+)+)", normalized)
        best: List[str] = []

        for candidate in chain_candidates:
            parts = [p.strip(" \t\n,.;") for p in candidate.split("->") if p.strip()]
            if len(parts) < 2:
                continue
            if not any(("*" in p) or ("(g" in p.lower()) for p in parts):
                continue
            if len(parts) > len(best):
                best = parts

        if not best and "->" in normalized:
            parts = [p.strip(" \t\n,.;") for p in normalized.split("->") if p.strip()]
            if len(parts) >= 2:
                best = parts

        cleaned = []
        for token in best:
            species = self._extract_primary_adsorbate_species(token)
            if not species:
                continue
            cleaned.append(self._ensure_star_notation(species))

        deduped = []
        prev_key = ""
        for item in cleaned:
            key = self._canonical_species_label(item)
            if key and key == prev_key:
                continue
            deduped.append(item)
            prev_key = key

        return deduped

    def _fallback_reaction_context(self, reaction_description: str) -> ReactionContext:
        _ = reaction_description
        raise RuntimeError("Fallback reaction-context generation is disabled by policy.")

    def _fallback_pathway_design(self, state: WorkflowState) -> PathwayDesign:
        _ = state
        raise RuntimeError("Fallback pathway design is disabled by policy.")

    def _normalize_reaction_context(self, state: WorkflowState) -> None:
        if not state.reaction_context:
            return

        context = state.reaction_context
        initial = self._extract_primary_adsorbate_species(context.initial_adsorbate)

        normalized_intermediates: List[str] = []
        for item in context.intermediates:
            if not isinstance(item, str):
                continue
            species = self._extract_primary_adsorbate_species(item)
            if species:
                normalized_intermediates.append(species)

        if initial:
            initial_key = self._canonical_species_label(initial)
            if not normalized_intermediates:
                normalized_intermediates = [initial]
            elif self._canonical_species_label(normalized_intermediates[0]) != initial_key:
                normalized_intermediates.insert(0, initial)

        deduped: List[str] = []
        prev_key = ""
        for item in normalized_intermediates:
            key = self._canonical_species_label(item)
            if key and key == prev_key:
                continue
            deduped.append(item)
            prev_key = key

        starred_intermediates = [self._ensure_star_notation(item) for item in deduped]
        parsed_from_text = self._extract_intermediate_sequence_from_text(state.reaction_description)

        if len(parsed_from_text) >= 2:
            starred_intermediates = parsed_from_text
            initial = parsed_from_text[0]

        context.initial_adsorbate = self._ensure_star_notation(
            initial or (starred_intermediates[0] if starred_intermediates else context.initial_adsorbate)
        )
        context.intermediates = starred_intermediates
        state.intermediates = starred_intermediates

    def _derive_intermediate_sequence_from_design(self, state: WorkflowState) -> List[str]:
        design = state.pathway_design
        if not design or not design.steps:
            return []

        sequence: List[str] = []
        for spec in design.steps:
            reactant = str(spec.reactant_formula or "").strip()
            product = str(spec.product_formula or "").strip()

            for species in (reactant, product):
                if not species:
                    continue
                if not sequence:
                    sequence.append(species)
                    continue

                if self._canonical_species_label(sequence[-1]) == self._canonical_species_label(species):
                    continue
                sequence.append(species)

        return sequence

    def _canonical_step_pair(self, reactant: str, product: str) -> Tuple[str, str]:
        reactant_primary = self._extract_primary_adsorbate_species(reactant)
        product_primary = self._extract_primary_adsorbate_species(product)
        return (
            self._canonical_species_label(reactant_primary),
            self._canonical_species_label(product_primary),
        )

    def _formula_element_key(self, label: str) -> Tuple[str, ...]:
        """Convert a species label to a sorted element tuple for fuzzy matching.

        Both 'CHOH' and 'CH2O' produce ('C', 'H', 'H', 'O').
        """
        tokens = self._extract_species_tokens(label)
        return tuple(sorted(tokens))

    def _align_pathway_design_to_context(self, state: WorkflowState, design: PathwayDesign) -> PathwayDesign:
        if not design:
            raise RuntimeError("Agent4 returned empty pathway design.")

        context_seq = list(state.reaction_context.intermediates) if state.reaction_context else list(state.intermediates or [])
        if len(context_seq) < 2:
            return design

        required_pairs = [(context_seq[i], context_seq[i + 1]) for i in range(len(context_seq) - 1)]
        if not required_pairs:
            return design

        pair_to_step: Dict[Tuple[str, str], PathwayStepSpec] = {}
        # Also build element-count-based fallback index
        elem_pair_to_step: Dict[Tuple[Tuple[str, ...], Tuple[str, ...]], PathwayStepSpec] = {}
        for spec in design.steps or []:
            pair = self._canonical_step_pair(spec.reactant_formula, spec.product_formula)
            if pair[0] and pair[1] and pair not in pair_to_step:
                pair_to_step[pair] = spec
            elem_pair = (
                self._formula_element_key(spec.reactant_formula),
                self._formula_element_key(spec.product_formula),
            )
            if elem_pair[0] and elem_pair[1] and elem_pair not in elem_pair_to_step:
                elem_pair_to_step[elem_pair] = spec

        aligned_steps: List[PathwayStepSpec] = []
        design_steps = list(design.steps or [])

        # Try pair-based matching first
        all_matched = True
        for idx, (reactant_label, product_label) in enumerate(required_pairs):
            expected_pair = self._canonical_step_pair(reactant_label, product_label)
            candidate = pair_to_step.get(expected_pair)

            # Fallback: match by element composition (CHOH == CH2O)
            if candidate is None:
                expected_elem_pair = (
                    self._formula_element_key(reactant_label),
                    self._formula_element_key(product_label),
                )
                candidate = elem_pair_to_step.get(expected_elem_pair)

            if candidate is None:
                all_matched = False
                break

            aligned_steps.append(
                PathwayStepSpec(
                    step_name=f"{reactant_label}_to_{product_label}",
                    reactant_formula=str(reactant_label),
                    product_formula=str(product_label),
                    reaction_type=str(candidate.reaction_type or "elementary_transition"),
                    atoms_to_add=list(candidate.atoms_to_add or []),
                    atoms_to_modify=list(candidate.atoms_to_modify or []),
                    atoms_to_remove=list(candidate.atoms_to_remove or []),
                    confidence=float(candidate.confidence or 0.0),
                )
            )

        # Fallback: sequential matching when LLM uses NEB-endpoint formulas (X->X)
        if not all_matched and len(design_steps) == len(required_pairs):
            import logging
            _log = logging.getLogger("catdt.species")
            _log.info(
                "Formula pair matching failed; falling back to sequential alignment "
                "(%d design steps == %d required pairs)",
                len(design_steps), len(required_pairs),
            )
            aligned_steps = []
            for idx, (reactant_label, product_label) in enumerate(required_pairs):
                candidate = design_steps[idx]
                aligned_steps.append(
                    PathwayStepSpec(
                        step_name=f"{reactant_label}_to_{product_label}",
                        reactant_formula=str(reactant_label),
                        product_formula=str(product_label),
                        reaction_type=str(candidate.reaction_type or "elementary_transition"),
                        atoms_to_add=list(candidate.atoms_to_add or []),
                        atoms_to_modify=list(candidate.atoms_to_modify or []),
                        atoms_to_remove=list(candidate.atoms_to_remove or []),
                        confidence=float(candidate.confidence or 0.0),
                    )
                )
        elif not all_matched:
            import logging
            _log = logging.getLogger("catdt.species")
            _log.error(
                "Alignment FAIL: available_pairs=%s; available_elem_pairs=%s; "
                "design_steps=[%s]; required=%d vs design=%d",
                list(pair_to_step.keys()),
                list(elem_pair_to_step.keys()),
                "; ".join(
                    f"{s.reactant_formula}->{s.product_formula}" for s in (design.steps or [])
                ),
                len(required_pairs), len(design_steps),
            )
            raise RuntimeError(
                f"Agent4 pathway design has {len(design_steps)} steps but "
                f"{len(required_pairs)} are required, and formula matching failed"
            )

        aligned = PathwayDesign(
            pathway_name=design.pathway_name or "aligned_pathway",
            description=(str(design.description or "").strip() + " | aligned_to_context").strip(" |"),
            overall_reaction=" -> ".join(context_seq),
            steps=aligned_steps,
            confidence=float(design.confidence or 0.0),
        )
        return self._normalize_pathway_step_operations(aligned)

    def _infer_formula_element_deltas(self, reactant_formula: str, product_formula: str) -> Tuple[List[str], List[str]]:
        reactant_tokens = Counter(self._extract_species_tokens(reactant_formula))
        product_tokens = Counter(self._extract_species_tokens(product_formula))

        to_add: List[str] = []
        to_remove: List[str] = []
        for element in sorted(set(reactant_tokens.keys()) | set(product_tokens.keys())):
            delta = product_tokens.get(element, 0) - reactant_tokens.get(element, 0)
            if delta > 0:
                to_add.extend([element] * delta)
            elif delta < 0:
                to_remove.extend([element] * (-delta))
        return to_add, to_remove

    def _count_explicit_remove_elements(self, atoms_to_remove: List[Any]) -> Counter:
        counts: Counter = Counter()
        for item in atoms_to_remove or []:
            if not isinstance(item, dict):
                continue

            explicit_symbol = item.get("element") or item.get("symbol") or item.get("species")
            tokens = self._extract_species_tokens(str(explicit_symbol or ""))
            if not tokens:
                continue

            count_raw = item.get("count", 1)
            try:
                n = max(int(count_raw), 1)
            except Exception:
                n = 1

            for token in tokens:
                counts[token] += n
        return counts

    def _normalize_pathway_step_operations(self, design: PathwayDesign) -> PathwayDesign:
        if not design or not design.steps:
            return design

        normalized_steps: List[PathwayStepSpec] = []
        for step in design.steps:
            step_copy = PathwayStepSpec(
                step_name=str(step.step_name),
                reactant_formula=str(step.reactant_formula),
                product_formula=str(step.product_formula),
                reaction_type=str(step.reaction_type or "elementary_transition"),
                atoms_to_add=list(step.atoms_to_add or []),
                atoms_to_modify=[],
                atoms_to_remove=list(step.atoms_to_remove or []),
                confidence=float(step.confidence or 0.0),
            )

            normalized_steps.append(step_copy)

        return PathwayDesign(
            pathway_name=design.pathway_name,
            description=design.description,
            overall_reaction=design.overall_reaction,
            steps=normalized_steps,
            confidence=float(design.confidence or 0.0),
        )

    @staticmethod
    def _is_adsorbate_species(label: Optional[str]) -> bool:
        if not label:
            return False
        return str(label).strip().startswith("*")

    def _get_working_intermediate_sequence(self, state: WorkflowState) -> List[str]:
        context_seq = list(state.reaction_context.intermediates) if state.reaction_context else []
        design_seq = self._derive_intermediate_sequence_from_design(state)

        if len(context_seq) >= 2:
            chosen = context_seq
        elif len(design_seq) >= 2:
            chosen = design_seq
        else:
            chosen = context_seq or design_seq

        if state.reaction_context is not None:
            state.reaction_context.intermediates = chosen
            if chosen and not state.reaction_context.initial_adsorbate:
                state.reaction_context.initial_adsorbate = chosen[0]

        state.intermediates = chosen
        return chosen

    def _resolve_removal_indices(self, step: PathwayStepSpec, adsorbate_indices: List[int]) -> List[int]:
        resolved: List[int] = []
        for item in step.atoms_to_remove:
            idx_candidates: List[Any] = []
            if isinstance(item, dict):
                if "index" in item:
                    idx_candidates = [item.get("index")]
                elif isinstance(item.get("indices"), list):
                    idx_candidates = list(item.get("indices", []))
            else:
                idx_candidates = [item]

            for idx_raw in idx_candidates:
                resolved_idx = self._resolve_atom_index(
                    idx_raw,
                    adsorbate_indices,
                    total_atoms=max(adsorbate_indices) + 1 if adsorbate_indices else 0,
                )
                if resolved_idx >= 0:
                    resolved.append(resolved_idx)
        return sorted(set(i for i in resolved if 0 <= i))

    def _resolve_formula_guided_removal_indices(
        self,
        step: PathwayStepSpec,
        structure: Atoms,
        adsorbate_indices: List[int],
    ) -> List[int]:
        explicit_indices = self._resolve_removal_indices(step, adsorbate_indices)
        chosen: List[int] = list(explicit_indices)
        used = set(chosen)

        remove_targets = self._count_explicit_remove_elements(step.atoms_to_remove)
        _, expected_remove = self._infer_formula_element_deltas(step.reactant_formula, step.product_formula)
        for element, count in Counter(expected_remove).items():
            remove_targets[element] = max(remove_targets.get(element, 0), count)

        if not remove_targets:
            return sorted(set(chosen))

        positions = structure.get_positions()
        anchor = np.mean(positions[adsorbate_indices], axis=0) if adsorbate_indices else np.mean(positions, axis=0)

        for element, target_count in remove_targets.items():
            existing = sum(1 for idx in chosen if 0 <= idx < len(structure) and structure[idx].symbol == element)
            need = max(0, int(target_count) - existing)
            if need <= 0:
                continue

            candidates = [
                idx for idx in adsorbate_indices
                if idx not in used and 0 <= idx < len(structure) and structure[idx].symbol == element
            ]
            if not candidates:
                continue

            candidates.sort(key=lambda idx: float(np.linalg.norm(positions[idx] - anchor)))
            for idx in candidates[:need]:
                chosen.append(idx)
                used.add(idx)

        return sorted(set(i for i in chosen if 0 <= i < len(structure)))
