"""Multi-pathway bundle builder and exporter.

Converts candidate pathways from the search engine into structured
output directories compatible with existing Agent4/5 consumption.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class PathwayBundleBuilder:
    """Build and export multi-pathway bundles to disk.

    Output layout::

        <output_base>/<run_id>/04_pathway/
            search_manifest.json
            iter_01_step_visualizations/
                path_00/
                    state_00_CO.vasp
                    step01_CO_to_CHO_reactant.vasp
                    step01_CO_to_CHO_product.vasp
                    ...
                path_01/
                    ...
    """

    def __init__(
        self,
        output_base_dir: str,
        run_id: str,
        iteration: int = 1,
    ):
        self.output_base_dir = output_base_dir
        self.run_id = run_id
        self.iteration = iteration
        self._pathway_dir = (
            Path(output_base_dir)
            / run_id
            / "04_pathway"
            / f"iter_{iteration:02d}_step_visualizations"
        )

    def export_pathway(
        self,
        pathway_id: str,
        pathway_index: int,
        states: List[Dict[str, Any]],
        steps: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Export a single pathway's structures and manifest entry.

        Args:
            pathway_id: unique pathway identifier
            pathway_index: numeric index (for path_XX directory)
            states: list of state dicts (species_label, structure_path, energy, etc.)
            steps: list of step dicts (reactant/product structure paths, step info)
            metadata: extra info to include in manifest

        Returns:
            Manifest entry dict for this pathway.
        """
        path_dir = self._pathway_dir / f"path_{pathway_index:02d}"
        path_dir.mkdir(parents=True, exist_ok=True)

        exported_states: List[Dict[str, Any]] = []
        exported_steps: List[Dict[str, Any]] = []

        # Export state structures
        for idx, state in enumerate(states):
            label = state.get("species_label", f"state_{idx}")
            safe_label = _safe_name(label)
            structure_path = state.get("structure_path")

            state_entry: Dict[str, Any] = {
                "index": idx,
                "label": label,
                "energy_eV": state.get("free_energy_eV"),
            }

            if structure_path and Path(structure_path).exists():
                try:
                    from ase.io import read, write
                    atoms = read(structure_path)
                    dest = path_dir / f"state_{idx:02d}_{safe_label}.vasp"
                    write(str(dest), atoms)
                    state_entry["structure_file"] = str(dest)
                except Exception as exc:
                    logger.debug("Failed to export state structure %s: %s", label, exc)

            exported_states.append(state_entry)

        # Export step structures (reactant/product pairs)
        for idx, step in enumerate(steps):
            step_name = step.get("step_name", f"step_{idx}")
            safe_name = _safe_name(step_name)

            step_entry: Dict[str, Any] = {
                "index": idx,
                "step_name": step_name,
                "operation_type": step.get("operation_type", ""),
            }

            for endpoint in ("reactant", "product"):
                src_path = step.get(f"{endpoint}_structure_path")
                if src_path and Path(src_path).exists():
                    try:
                        from ase.io import read, write
                        atoms = read(src_path)
                        dest = path_dir / f"step{idx + 1:02d}_{safe_name}_{endpoint}.vasp"
                        write(str(dest), atoms)
                        step_entry[f"{endpoint}_file"] = str(dest)
                    except Exception as exc:
                        logger.debug("Failed to export %s for %s: %s", endpoint, step_name, exc)

            exported_steps.append(step_entry)

        manifest_entry = {
            "pathway_id": pathway_id,
            "pathway_index": pathway_index,
            "directory": str(path_dir),
            "states": exported_states,
            "steps": exported_steps,
            "metadata": metadata or {},
        }

        # Write per-pathway manifest
        manifest_file = path_dir / "pathway_manifest.json"
        try:
            manifest_file.write_text(
                json.dumps(manifest_entry, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8",
            )
        except Exception as exc:
            logger.debug("Failed to write pathway manifest: %s", exc)

        return manifest_entry

    def export_multi_pathway_bundle(
        self,
        pathways: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Export multiple pathways and write a top-level search manifest.

        Args:
            pathways: list of dicts, each with:
                - pathway_id
                - is_retained
                - states (list)
                - steps (list)
                - metadata (optional)

        Returns:
            Complete search manifest dict.
        """
        manifest_entries: List[Dict[str, Any]] = []

        for idx, pw in enumerate(pathways):
            entry = self.export_pathway(
                pathway_id=pw.get("pathway_id", f"path_{idx}"),
                pathway_index=idx,
                states=pw.get("states", []),
                steps=pw.get("steps", []),
                metadata={
                    "is_retained": pw.get("is_retained", True),
                    "prune_reason": pw.get("prune_reason", ""),
                    "total_energy_change_eV": pw.get("total_free_energy_change_eV"),
                    **(pw.get("metadata", {})),
                },
            )
            manifest_entries.append(entry)

        search_manifest = {
            "run_id": self.run_id,
            "iteration": self.iteration,
            "n_pathways": len(pathways),
            "n_retained": sum(1 for p in pathways if p.get("is_retained", True)),
            "pathways": manifest_entries,
        }

        # Write top-level manifest
        manifest_path = self._pathway_dir.parent / "search_manifest.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            manifest_path.write_text(
                json.dumps(search_manifest, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8",
            )
        except Exception as exc:
            logger.debug("Failed to write search manifest: %s", exc)

        return search_manifest


def _safe_name(name: str) -> str:
    """Convert a label to filesystem-safe name."""
    return re.sub(r"[^\w\-]+", "_", name).strip("_")


__all__ = ["PathwayBundleBuilder"]
