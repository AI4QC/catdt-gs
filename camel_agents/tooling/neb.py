
"""NEB and endpoint validation tool mixin."""

import json
import os
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms, Atom
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import FIRE

from core.pathway.enhanced_neb_validator import EnhancedNEBValidator, validate_neb_endpoints, ValidationStatus, ReactionType

from camel_agents.schemas import PathwayStepSpec, WorkflowState

from core.fairchem_config import DEFAULT_FAIRCHEM_MODEL

from .common import DEPS_BASE_PATH, logger


class NEBToolsMixin:
    def _reorder_product_to_match_reactant(self, *args: Any, **kwargs: Any):
        """Delegate to the single MIC-aware implementation in geometry.WorkflowGeometryMixin."""
        from .geometry import WorkflowGeometryMixin
        return WorkflowGeometryMixin._reorder_product_to_match_reactant(self, *args, **kwargs)

    @staticmethod
    def _valid_atom_mapping_constraints(
        reactant: Atoms,
        product: Atoms,
        atom_mapping_constraints: Optional[List[Dict[str, Any]]],
    ) -> List[Tuple[int, int]]:
        """Delegate to the single implementation in geometry.WorkflowGeometryMixin."""
        from .geometry import WorkflowGeometryMixin
        return WorkflowGeometryMixin._valid_atom_mapping_constraints(
            reactant=reactant,
            product=product,
            atom_mapping_constraints=atom_mapping_constraints,
        )

    def run_neb_for_steps(
        self,
        steps: List[Dict[str, Any]],
        output_dir: str,
        n_frames: int = 10,
        fmax: float = 0.1,
        max_steps: int = 300,
    ) -> Dict[str, Any]:
        from core.pathway.barrier_predictor import BarrierPredictor

        def _build_barrier():
            predictor = self.get_shared_uma_predictor()
            return BarrierPredictor(
                fairchem_predictor=predictor,
                work_dir=output_dir,
                use_llm_controller=False,
            )

        results: Dict[str, Any] = {}
        barrier = None

        try:
            barrier = _build_barrier()
        except Exception as exc_gpu:
            logger.warning("Failed to initialize NEB predictor on GPU: %s", exc_gpu)
            for step in steps:
                results[step["name"]] = {"error": f"predictor_init_failed: {exc_gpu}"}
            return results

        for step_idx, step in enumerate(steps):
            if step_idx > 0:
                try:
                    import gc
                    gc.collect()
                    import torch
                    torch.cuda.empty_cache()
                except Exception:
                    pass
            name = step["name"]
            reactant = step["reactant"].copy()
            product = step["product"].copy()
            reactant_ads = list(step.get("reactant_adsorbate_indices", []))
            product_ads = list(step.get("product_adsorbate_indices", []))
            reactant_staged = list(step.get("reactant_staged_indices", []))
            product_staged = list(step.get("product_staged_indices", []))

            product, product_ads, product_staged, reordered = self._reorder_product_to_match_reactant(
                reactant=reactant,
                product=product,
                product_adsorbate_indices=product_ads,
                product_staged_indices=product_staged,
                atom_mapping_constraints=list(step.get("atom_mapping_constraints", [])),
            )
            if reordered:
                logger.info("Reordered product atom order for NEB step %s", name)

            reactant_ads_set = set(idx for idx in reactant_ads if 0 <= idx < len(reactant))
            product_ads_set = set(idx for idx in product_ads if 0 <= idx < len(product))
            reactant_fixed_surface = [idx for idx in range(len(reactant)) if idx not in reactant_ads_set]
            product_fixed_surface = [idx for idx in range(len(product)) if idx not in product_ads_set]

            try:
                if reactant_fixed_surface:
                    reactant.set_constraint(FixAtoms(indices=reactant_fixed_surface))
                if product_fixed_surface:
                    product.set_constraint(FixAtoms(indices=product_fixed_surface))
            except Exception as exc:
                logger.warning("Failed to set surface constraints for step %s: %s", name, exc)

            # Pre-NEB geometry screening: max per-atom displacement (MIC/PBC-aware).
            # If any atom must travel >4 Å between endpoints, NEB linear interpolation
            # forces atoms through physically unreasonable paths (atoms overlap, bonds
            # stretch wildly), producing 5-20 eV pseudo-barriers that drown the real
            # chemistry. This usually means staging placed atoms in inconsistent sites
            # between R and P — the fix belongs upstream (Agent4/5), so abort NEB here.
            try:
                if len(reactant) == len(product):
                    r_pos = reactant.get_positions()
                    p_pos = product.get_positions()
                    cell = np.array(reactant.cell.array, dtype=float)
                    pbc = tuple(bool(x) for x in reactant.pbc)
                    delta = p_pos - r_pos
                    if cell.shape == (3, 3) and any(pbc):
                        # Minimum-image convention: wrap fractional delta to [-0.5, 0.5]
                        inv_t = np.linalg.pinv(cell.T)
                        frac = (inv_t @ delta.T).T
                        for dim in range(3):
                            if pbc[dim]:
                                frac[:, dim] -= np.round(frac[:, dim])
                        delta = frac @ cell
                    disp = np.linalg.norm(delta, axis=1)
                    d_max = float(disp.max()) if len(disp) else 0.0
                    d_max_idx = int(np.argmax(disp)) if len(disp) else -1
                    if d_max > 4.0:
                        sym = reactant.get_chemical_symbols()[d_max_idx] if d_max_idx >= 0 else "?"
                        logger.warning(
                            "Step %s: max per-atom PBC displacement %.2f Å (atom %d '%s') "
                            "exceeds 4.0 Å — staging inconsistency between endpoints. "
                            "Skipping NEB (would produce spurious high barrier).",
                            name, d_max, d_max_idx, sym,
                        )
                        results[name] = {
                            "error": "endpoint_displacement_too_large",
                            "max_disp": d_max,
                            "worst_atom_index": d_max_idx,
                            "worst_atom_symbol": sym,
                        }
                        continue
            except Exception as exc_disp:
                logger.warning("Pre-NEB displacement screen failed for step %s: %s (proceeding)", name, exc_disp)

            # Pre-NEB energy screening: skip obviously bad endpoints
            try:
                barrier.fairchem._load_model()
                _r_tmp = reactant.copy()
                _p_tmp = product.copy()
                _r_tmp.calc = barrier.fairchem._calculator
                _p_tmp.calc = barrier.fairchem._calculator
                e_r = float(_r_tmp.get_potential_energy())
                e_p = float(_p_tmp.get_potential_energy())
                dE = abs(e_r - e_p)
                if dE > 8.0:
                    logger.warning(
                        "Step %s: endpoint energy gap %.2f eV exceeds 8.0 eV threshold "
                        "(E_reactant=%.3f, E_product=%.3f). Skipping NEB.",
                        name, dE, e_r, e_p,
                    )
                    results[name] = {"error": "endpoint_energy_gap_too_large", "dE": dE}
                    continue
            except Exception as exc_screen:
                logger.warning("Pre-NEB energy screening failed for step %s: %s (proceeding with NEB)", name, exc_screen)

            has_staged = bool(reactant_staged or product_staged)
            try:
                result = barrier.predict_from_structures(
                    reactant=reactant,
                    product=product,
                    n_frames=n_frames,
                    fmax=fmax,
                    max_steps=max_steps,
                    relax_endpoints=has_staged,
                    reaction_name=name,
                    reactant_adsorbate_indices=reactant_ads,
                    product_adsorbate_indices=product_ads,
                    reactant_staged_indices=reactant_staged,
                    product_staged_indices=product_staged,
                )
                results[name] = result
            except Exception as exc:
                logger.warning("NEB failed for step %s on GPU: %s", name, exc)
                results[name] = {"error": str(exc)}
            finally:
                try:
                    import gc
                    import torch
                    # Clear any result references that hold CUDA tensors
                    if name in results and hasattr(results[name], 'neb_frames'):
                        for fr in (results[name].neb_frames or []):
                            if hasattr(fr, 'calc') and fr.calc is not None:
                                fr.calc = None
                    gc.collect()
                    torch.cuda.empty_cache()
                    if torch.cuda.is_available():
                        alloc = torch.cuda.memory_allocated() / 1e6
                        resrv = torch.cuda.memory_reserved() / 1e6
                        if alloc > 2000:
                            logger.warning(
                                "GPU memory high after step %d (%s): %.0f/%.0f MB",
                                step_idx, name, alloc, resrv,
                            )
                except Exception:
                    pass

        return results

    def validate_neb_endpoints(
        self,
        initial_structure_path: str,
        final_structure_path: str,
        expected_reaction_type: Optional[str] = None,
        auto_fix: bool = True,
        run_id: str = "default_run",
        workflow_step: str = "04_neb_validation"
    ) -> str:
        """
        验证 NEB 初末态结构，并提供智能修复建议

        Parameters
        ----------
        initial_structure_path : str
            初态结构文件路径（VASP格式）
        final_structure_path : str
            末态结构文件路径（VASP格式）
        expected_reaction_type : str, optional
            期望的反应类型（"hydrogenation", "dehydrogenation", "dissociation", "coupling"）
        auto_fix : bool
            是否尝试自动修复问题
        run_id : str
            运行ID
        workflow_step : str
            工作流步骤标识

        Returns
        -------
        str
            验证报告JSON文件路径
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Validating NEB endpoints: {initial_structure_path} -> {final_structure_path}")

        # 读取结构
        initial = read(initial_structure_path)
        final = read(final_structure_path)

        # 检查必要的元数据
        if 'adsorbate_indices' not in initial.info:
            logger.warning("Initial structure missing 'adsorbate_indices', attempting to auto-detect...")
            initial.info['adsorbate_indices'] = self._detect_adsorbate_indices(initial)

        if 'adsorbate_indices' not in final.info:
            logger.warning("Final structure missing 'adsorbate_indices', attempting to auto-detect...")
            final.info['adsorbate_indices'] = self._detect_adsorbate_indices(final)

        # 执行验证
        validator = EnhancedNEBValidator(auto_fix=auto_fix, fix_attempts=3)
        report = validator.validate_endpoint_pair(initial, final, expected_reaction_type)

        # EnhancedNEBValidator does NOT expose the fixed Atoms objects, so an
        # AUTO_FIXED status only means fixes were applied internally during
        # re-validation — the ORIGINAL (unfixed) structures are what downstream
        # NEB will use. Reflect that honestly by downgrading to WARNING.
        if report.status == ValidationStatus.AUTO_FIXED and auto_fix:
            logger.warning(
                "Validator reported AUTO_FIXED, but fixed structures are not "
                "exposed by EnhancedNEBValidator; the original (unfixed) "
                "structures will be used. Downgrading status to WARNING."
            )
            report.status = ValidationStatus.WARNING

        # 保存验证报告
        report_path = output_dir / "validation_report.json"
        report.save(str(report_path))

        # 生成可视化（如果验证失败）
        if report.status != ValidationStatus.PASS:
            self._generate_validation_visualization(initial, final, report, output_dir)

        logger.info(f"Validation status: {report.status.value}")
        logger.info(f"Validation report saved: {report_path}")

        return str(report_path)

    def _detect_adsorbate_indices(self, atoms: Atoms) -> List[int]:
        """自动检测吸附物原子索引。

        Heuristic order:
        1. ``atoms.info['adsorbate_indices']`` — trusted if present and valid.
        2. ASE tags — atoms carrying the maximum tag are treated as adsorbate
           when tags partition the structure non-trivially.
        3. z-position fallback — slab elements are identified as elements with
           at least one atom in the LOWER half of the z-range (slab bottom);
           this works on oxides where O dominates the element count. Adsorbate
           atoms are non-slab-element atoms above the slab top (+0.5 Å).
        """
        n = len(atoms)
        if n == 0:
            return []

        # 1. Known indices from atoms.info
        known = atoms.info.get("adsorbate_indices")
        if isinstance(known, (list, tuple)):
            valid_known = [int(i) for i in known if 0 <= int(i) < n]
            if valid_known:
                logger.info(f"Adsorbate indices from atoms.info: {valid_known}")
                return sorted(valid_known)

        # 2. ASE tags
        tags = atoms.get_tags()
        if len(tags) == n and len(np.unique(tags)) > 1:
            max_tag = int(np.max(tags))
            tagged = [i for i, tag in enumerate(tags.tolist()) if int(tag) == max_tag]
            if 0 < len(tagged) < n:
                logger.info(f"Adsorbate indices from ASE tags (tag={max_tag}): {tagged}")
                return sorted(tagged)

        # 3. z-position fallback (element-count heuristic, slab-bottom aware)
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        z = positions[:, 2]
        z_min, z_max = float(np.min(z)), float(np.max(z))
        z_cut = z_min + 0.5 * (z_max - z_min)
        slab_elements = {symbols[i] for i in range(n) if z[i] <= z_cut}
        if not slab_elements:
            # Degenerate (flat) structure: fall back to dominant element only.
            slab_elements = {Counter(symbols).most_common(1)[0][0]}

        slab_z = [z[i] for i in range(n) if symbols[i] in slab_elements]
        surface_z = max(slab_z) if slab_z else z_max

        adsorbate_indices = [
            i for i in range(n)
            if symbols[i] not in slab_elements and z[i] > surface_z + 0.5
        ]

        logger.info(f"Auto-detected adsorbate indices: {adsorbate_indices}")
        return adsorbate_indices

    def _generate_validation_visualization(
        self,
        initial: Atoms,
        final: Atoms,
        report: Any,
        output_dir: Path
    ):
        """生成验证问题的可视化"""
        try:
            viz_manager = self._get_viz_manager_instance(output_dir)

            # 高亮问题的结构图
            initial_viz = output_dir / "initial_with_issues.png"
            final_viz = output_dir / "final_with_issues.png"

            viz_manager.structure_visualizer.visualize_structure(
                initial,
                output_file=str(initial_viz),
                title=f"Initial - Status: {report.status.value}",
                show=False,
            )
            viz_manager.structure_visualizer.visualize_structure(
                final,
                output_file=str(final_viz),
                title=f"Final - Status: {report.status.value}",
                show=False,
            )

            logger.info(f"Validation visualizations saved to {output_dir}")
        except Exception as e:
            logger.warning(f"Failed to generate validation visualization: {e}")

    def run_neb_with_validation(
        self,
        initial_structure_path: str,
        final_structure_path: str,
        expected_reaction_type: Optional[str] = None,
        n_images: int = 7,
        run_id: str = "default_run",
        workflow_step: str = "05_neb_calculation"
    ) -> str:
        """
        运行 NEB 计算，带自动预验证和修复

        这是增强版的 NEB 计算，会自动：
        1. 验证初末态结构
        2. 尝试自动修复问题
        3. 使用修复后的结构运行 NEB
        4. 如果失败，提供详细的诊断报告

        Parameters
        ----------
        initial_structure_path : str
            初态结构路径
        final_structure_path : str
            末态结构路径
        expected_reaction_type : str, optional
            期望的反应类型
        n_images : int
            NEB 图像数量
        run_id : str
            运行ID
        workflow_step : str
            工作流步骤

        Returns
        -------
        str
            NEB 结果 pickle 文件路径
        """
        output_dir = self.output_base_dir / run_id / workflow_step
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info("=" * 60)
        logger.info("Running NEB with Pre-validation")
        logger.info("=" * 60)

        # 步骤 1: 验证结构
        logger.info("Step 1: Validating endpoint structures...")
        validation_report_path = self.validate_neb_endpoints(
            initial_structure_path,
            final_structure_path,
            expected_reaction_type,
            auto_fix=True,
            run_id=run_id,
            workflow_step=workflow_step + "_validation"
        )

        # 加载验证报告
        with open(validation_report_path, 'r') as f:
            validation_result = json.load(f)

        # 步骤 2: 根据验证结果决定下一步
        status = validation_result.get('status')

        if status == 'fail':
            logger.error("Structure validation FAILED. Cannot proceed with NEB.")
            logger.error("Issues found:")
            for issue in validation_result.get('issues', []):
                if issue['severity'] == 'critical':
                    logger.error(f"  - {issue['description']}")
            logger.error("Please fix the issues manually or check the validation report.")
            raise RuntimeError(f"NEB validation failed. See {validation_report_path}")

        elif status == 'warning':
            logger.warning("Structure validation has warnings, but proceeding with NEB...")
            for rec in validation_result.get('recommendations', []):
                logger.warning(f"  - {rec}")

        elif status == 'auto_fixed':
            logger.warning(
                "Validation reported auto_fixed, but auto-fix is unavailable at "
                "the NEB level (fixed structures are not exposed); proceeding "
                "with the ORIGINAL structures."
            )

        else:  # pass
            logger.info("Structure validation PASSED. Proceeding with NEB.")

        # 步骤 3: 运行 NEB（使用现有的 run_neb_for_steps 逻辑）
        logger.info("Step 2: Running NEB calculation...")

        # Build the steps payload expected by run_neb_for_steps from the
        # structure files (each step needs 'name', 'reactant', 'product').
        initial_atoms = read(initial_structure_path)
        final_atoms = read(final_structure_path)
        if 'adsorbate_indices' not in initial_atoms.info:
            initial_atoms.info['adsorbate_indices'] = self._detect_adsorbate_indices(initial_atoms)
        if 'adsorbate_indices' not in final_atoms.info:
            final_atoms.info['adsorbate_indices'] = self._detect_adsorbate_indices(final_atoms)

        step_name = f"{Path(initial_structure_path).stem}_to_{Path(final_structure_path).stem}"
        steps_payload = [{
            "name": step_name,
            "reactant": initial_atoms,
            "product": final_atoms,
            "reactant_adsorbate_indices": list(initial_atoms.info.get("adsorbate_indices", [])),
            "product_adsorbate_indices": list(final_atoms.info.get("adsorbate_indices", [])),
        }]

        neb_results = self.run_neb_for_steps(
            steps=steps_payload,
            output_dir=str(output_dir),
            n_frames=int(n_images),
        )

        neb_result_path = self._save_result_to_pickle(neb_results, output_dir / "neb_results.pkl")
        logger.info(f"NEB calculation completed: {neb_result_path}")
        return neb_result_path


# ---------------------------------------------------------------------------
# Agent4/5 generic geometry staging and validation tools
# ---------------------------------------------------------------------------


class Agent45WorkflowToolsMixin:
    def _agent45_geometry(self, workflow: Any) -> "Agent45GeometryTools":
        return Agent45GeometryTools(workflow)

    def build_steps_payload(self, workflow: Any, step_structures: List[Dict[str, Any]]) -> str:
        return self._agent45_geometry(workflow).build_steps_payload(step_structures)

    def build_tool_baseline_steps(
        self,
        workflow: Any,
        state: WorkflowState,
        base_structure: Atoms,
        preplaced_products: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        return self._agent45_geometry(workflow).build_tool_baseline_steps(
            state=state,
            base_structure=base_structure,
            preplaced_products=preplaced_products,
        )

    def build_single_step_baseline(
        self,
        workflow: Any,
        state: WorkflowState,
        current_structure: Atoms,
        preplaced_products: Dict[str, Dict[str, Any]],
        step_idx: int,
    ) -> Dict[str, Any]:
        return self._agent45_geometry(workflow).build_single_step_baseline(
            state=state,
            current_structure=current_structure,
            preplaced_products=preplaced_products,
            step_idx=step_idx,
        )

    def apply_step_ops_on_fixed_product(
        self,
        workflow: Any,
        current_physical: Atoms,
        fixed_product: Atoms,
        step: PathwayStepSpec,
    ) -> Dict[str, Any]:
        return self._agent45_geometry(workflow).apply_step_ops_on_fixed_product(
            current_physical=current_physical,
            fixed_product=fixed_product,
            step=step,
        )

    def build_step_structures(
        self,
        workflow: Any,
        state: WorkflowState,
        base_structure: Atoms,
        steps: List[PathwayStepSpec],
        preplaced_products: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        return self._agent45_geometry(workflow).build_step_structures(
            state=state,
            base_structure=base_structure,
            steps=steps,
            preplaced_products=preplaced_products,
        )

    def postprocess_step_structures(
        self,
        workflow: Any,
        step_structures: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        return self._agent45_geometry(workflow).postprocess_step_structures(step_structures)

    def programmatic_validation(self, workflow: Any, step_structures: List[Dict[str, Any]]) -> Dict[str, Any]:
        return self._agent45_geometry(workflow).programmatic_validation(step_structures)

    def get_step_element_deltas(self, workflow: Any) -> List[Dict[str, Any]]:
        return self._agent45_geometry(workflow).get_step_element_deltas()

    def suggest_staged_positions(
        self, workflow: Any, requests: List[Dict[str, str]],
    ) -> Dict[str, Any]:
        return self._agent45_geometry(workflow).suggest_staged_positions(requests=requests)

    def validate_proposed_positions(
        self,
        workflow: Any,
        steps: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        return self._agent45_geometry(workflow).validate_proposed_positions(steps=steps)

    def run_agent45_energy_gate(
        self,
        workflow: Any,
        step_structures: List[Dict[str, Any]],
        output_dir: str = "",
        min_pair_fatal_threshold: float = 0.80,
        force_fatal_threshold: float = 1.0,
        force_warning_threshold: float = 0.5,
    ) -> Dict[str, Any]:
        return self._agent45_geometry(workflow).run_agent45_energy_gate(
            step_structures=step_structures,
            output_dir=output_dir,
            min_pair_fatal_threshold=min_pair_fatal_threshold,
            force_fatal_threshold=force_fatal_threshold,
            force_warning_threshold=force_warning_threshold,
        )

    def reset_site_retry_state(self, workflow: Any = None) -> None:
        """Clear retained Agent4/5 site retry ranks (call between independent runs)."""
        state = getattr(self, "_agent45_site_retry_state", None)
        if isinstance(state, dict):
            state.clear()

    def find_bond_change_sites(
        self,
        workflow: Any,
        clean_reactant: Atoms,
        clean_reactant_ads: List[int],
        clean_product: Atoms,
        clean_product_ads: List[int],
        min_site_dist: float = 2.5,
        max_site_dist: float = 4.5,
        max_sites_per_hint: int = 5,
    ) -> List[Dict[str, Any]]:
        return self._agent45_geometry(workflow).find_bond_change_sites(
            clean_reactant=clean_reactant,
            clean_reactant_ads=clean_reactant_ads,
            clean_product=clean_product,
            clean_product_ads=clean_product_ads,
            min_site_dist=min_site_dist,
            max_site_dist=max_site_dist,
            max_sites_per_hint=max_sites_per_hint,
        )



class Agent45GeometryTools:
    def __init__(self, workflow: Any):
        self.workflow = workflow

    # ── Phase-2 Stage B: geometry-based staging hint generator ──────────────

    def find_bond_change_sites(
        self,
        clean_reactant: Atoms,
        clean_reactant_ads: List[int],
        clean_product: Atoms,
        clean_product_ads: List[int],
        min_site_dist: float = 2.5,
        max_site_dist: float = 4.5,
        max_sites_per_hint: int = 5,
    ) -> List[Dict[str, Any]]:
        """Generate staging hints by RMSD-matching adsorbate atoms across endpoints.

        Algorithm (Q3 option A — geometry-based, element-agnostic beyond matching):
        1. Greedy same-element RMSD match between reactant and product adsorbates.
           Pairs are consumed in ascending distance order.
        2. Unmatched atoms on one side are atoms that must be added (as staging)
           to the OTHER side so that NEB sees matching element multisets.
        3. For each unmatched atom, the anchor is its own position in the side
           that contains it; that anchor position is then used to enumerate
           nearby surface adsorption sites on the OPPOSITE (staging) side.

        Returns
        -------
        list of dict, each with::

            {
              "element": str,                # element of the staging atom
              "side":    "reactant"|"product",  # which side gets the staging atom
              "anchor_position":  [x,y,z],   # position in the source (has-atom) side
              "source_side":      "reactant"|"product",
              "bond_partner_element": str,   # nearest heavy-atom bonding partner (if any)
              "bond_partner_position": [x,y,z] | None,
              "candidate_sites": [
                  {"position": [x,y,z], "site_type": str,
                   "distance_to_anchor": float},
                  ...
              ],
            }
        """
        r_pos = clean_reactant.get_positions()
        p_pos = clean_product.get_positions()
        r_sym = clean_reactant.get_chemical_symbols()
        p_sym = clean_product.get_chemical_symbols()

        # Normalize to valid index lists
        r_valid = [i for i in clean_reactant_ads if 0 <= i < len(clean_reactant)]
        p_valid = [i for i in clean_product_ads if 0 <= i < len(clean_product)]

        # Greedy same-element matching, shortest distance first (under PBC).
        pairs: List[Tuple[float, int, int]] = []
        for ri in r_valid:
            for pi in p_valid:
                if r_sym[ri] != p_sym[pi]:
                    continue
                dist = float(
                    self.workflow._pbc_distance(clean_reactant, r_pos[ri], p_pos[pi])
                )
                pairs.append((dist, ri, pi))
        pairs.sort(key=lambda t: t[0])

        matched_r: set[int] = set()
        matched_p: set[int] = set()
        for dist, ri, pi in pairs:
            if ri in matched_r or pi in matched_p:
                continue
            matched_r.add(ri)
            matched_p.add(pi)

        unmatched_r = [i for i in r_valid if i not in matched_r]
        unmatched_p = [i for i in p_valid if i not in matched_p]

        hints: List[Dict[str, Any]] = []

        # Atoms present in reactant but not product → staging copy goes on PRODUCT side.
        for ri in unmatched_r:
            hint = self._build_bond_change_hint(
                source_side="reactant",
                target_side="product",
                source_atoms=clean_reactant,
                source_ads=r_valid,
                source_atom_idx=ri,
                target_atoms=clean_product,
                target_ads=p_valid,
                element=r_sym[ri],
                min_dist=min_site_dist,
                max_dist=max_site_dist,
                max_sites=max_sites_per_hint,
            )
            if hint is not None:
                hints.append(hint)

        # Atoms present in product but not reactant → staging copy goes on REACTANT side.
        for pi in unmatched_p:
            hint = self._build_bond_change_hint(
                source_side="product",
                target_side="reactant",
                source_atoms=clean_product,
                source_ads=p_valid,
                source_atom_idx=pi,
                target_atoms=clean_reactant,
                target_ads=r_valid,
                element=p_sym[pi],
                min_dist=min_site_dist,
                max_dist=max_site_dist,
                max_sites=max_sites_per_hint,
            )
            if hint is not None:
                hints.append(hint)

        return hints

    def _build_bond_change_hint(
        self,
        source_side: str,
        target_side: str,
        source_atoms: Atoms,
        source_ads: List[int],
        source_atom_idx: int,
        target_atoms: Atoms,
        target_ads: List[int],
        element: str,
        min_dist: float,
        max_dist: float,
        max_sites: int,
    ) -> Optional[Dict[str, Any]]:
        """Assemble a single bond-change hint (helper for ``find_bond_change_sites``)."""
        source_pos = source_atoms.get_positions()[source_atom_idx]
        source_syms = source_atoms.get_chemical_symbols()

        # Find the closest heavy-atom bonding partner in the SOURCE side.
        partner_idx: Optional[int] = None
        partner_dist = float("inf")
        for j in source_ads:
            if j == source_atom_idx:
                continue
            if source_syms[j] == "H":
                continue
            d = float(
                self.workflow._pbc_distance(source_atoms, source_pos, source_atoms.get_positions()[j])
            )
            if d < partner_dist:
                partner_dist = d
                partner_idx = j

        partner_element: Optional[str] = None
        partner_position: Optional[List[float]] = None
        if partner_idx is not None and partner_dist < 2.5:
            partner_element = source_syms[partner_idx]
            partner_position = [
                float(x) for x in source_atoms.get_positions()[partner_idx]
            ]

        # Enumerate nearby surface sites on the TARGET side near the anchor position.
        # The anchor position is taken directly from the source side (same cell).
        source_anchor = np.array(source_pos, dtype=float)
        search_windows = [(float(min_dist), float(max_dist))]
        widened = (max(0.8, float(min_dist) - 1.5), float(max_dist) + 2.0)
        if widened != search_windows[0]:
            search_windows.append(widened)

        sites: List[Dict[str, Any]] = []
        for lo, hi in search_windows:
            sites = self.workflow.enumerate_surface_sites_near_anchor(
                target_atoms,
                target_ads,
                source_anchor,
                element,
                lo,
                hi,
            )
            if sites:
                break

        if str(element) == "H" and self._has_physical_surface_context(target_atoms, target_ads):
            local_sites = self._enumerate_local_h_surface_sites(
                target_atoms,
                target_ads,
                source_anchor,
                max_sites=max(max_sites, 3),
            )
            if local_sites:
                by_key: Dict[Tuple[float, float, float], Dict[str, Any]] = {}
                for site in list(local_sites) + list(sites or []):
                    pos = site.get("position", [])
                    if not isinstance(pos, list) or len(pos) != 3:
                        continue
                    key = (round(float(pos[0]), 3), round(float(pos[1]), 3), round(float(pos[2]), 3))
                    by_key.setdefault(key, site)
                sites = sorted(
                    by_key.values(),
                    key=lambda site: (
                        0 if bool(site.get("active_center_site", False)) else 1,
                        0 if str(site.get("site_label", "")).startswith("local_H_") else 1,
                        float(site.get("dist_from_anchor", float("inf"))),
                        float(site.get("position", [0.0, 0.0, 0.0])[2]),
                    ),
                )

        candidate_sites: List[Dict[str, Any]] = []
        if sites:
            for s in sites[:max_sites]:
                candidate_sites.append(
                    {
                        "position": [float(x) for x in s.get("position", [])],
                        "site_type": str(s.get("site_label", "unknown")),
                        "distance_to_anchor": float(s.get("dist_from_anchor", 0.0)),
                    }
                )

        if not candidate_sites:
            return None

        return {
            "element": element,
            "side": target_side,
            "source_side": source_side,
            "source_atom_index": int(source_atom_idx),
            "anchor_position": [float(x) for x in source_pos],
            "bond_partner_element": partner_element,
            "bond_partner_index": int(partner_idx) if partner_idx is not None else None,
            "bond_partner_position": partner_position,
            "candidate_sites": candidate_sites,
        }

    def _nearest_reference_position(
        self,
        structure: Atoms,
        position: np.ndarray,
        reference_positions: List[np.ndarray] | np.ndarray | List[List[float]],
    ) -> Tuple[Optional[int], Optional[np.ndarray], float]:
        """Return nearest reference position using the workflow's canonical PBC distance."""
        if structure is None:
            return None, None, float("inf")

        try:
            refs_iter = list(reference_positions)
        except TypeError:
            return None, None, float("inf")

        if not refs_iter:
            return None, None, float("inf")

        pos_arr = np.array(position, dtype=float)
        best_idx: Optional[int] = None
        best_pos: Optional[np.ndarray] = None
        best_dist = float("inf")

        for idx, ref_pos in enumerate(refs_iter):
            ref_arr = np.array(ref_pos, dtype=float)
            if ref_arr.shape != (3,):
                continue
            dist = float(self.workflow._pbc_distance(structure, pos_arr, ref_arr))
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
                best_pos = ref_arr

        return best_idx, best_pos, best_dist

    def _validate_staged_position(
        self,
        resolved_pos: List[float],
        fallback_pos: List[float],
        structure: Atoms,
        adsorbate_indices: List[int],
        element: str,
        step_name: str,
        endpoint: str,
        max_distance_from_ads: float = 5.0,
        min_distance_from_any: float = 0.8,
    ) -> List[float]:
        """Validate Agent4's resolved position; fall back to staging plan if bad.

        Checks:
        - Position must be within *max_distance_from_ads* of nearest adsorbate atom.
        - Position must not overlap with any existing atom (< *min_distance_from_any*).
        If either check fails, the staging plan position is used instead.
        """
        pos_arr = np.array(resolved_pos, dtype=float)
        positions = structure.get_positions()

        # Check distance to nearest adsorbate atom
        valid_ads = [i for i in adsorbate_indices if 0 <= i < len(structure)]
        if valid_ads:
            _, _, min_ads_dist = self._nearest_reference_position(
                structure,
                pos_arr,
                [positions[idx] for idx in valid_ads],
            )
            if min_ads_dist > max_distance_from_ads:
                logger.warning(
                    "Step '%s' %s: Agent4 position for %s is %.2f Å from nearest adsorbate "
                    "(threshold %.1f Å); replacing with staging plan position.",
                    step_name, endpoint, element, min_ads_dist, max_distance_from_ads,
                )
                return list(fallback_pos)

        # Check overlap with any existing atom
        if len(positions) > 0:
            _, _, min_any_dist = self._nearest_reference_position(structure, pos_arr, positions)
            if min_any_dist < min_distance_from_any:
                logger.warning(
                    "Step '%s' %s: Agent4 position for %s overlaps with existing atom "
                    "(distance %.3f Å < %.1f Å threshold); replacing with staging plan position.",
                    step_name, endpoint, element, min_any_dist, min_distance_from_any,
                )
                return list(fallback_pos)

        return list(resolved_pos)

    # ------------------------------------------------------------------
    # Tool-callable helpers for Agent4/5 (fully generic, no reaction-specific logic)
    # ------------------------------------------------------------------

    def get_step_element_deltas(self) -> List[Dict[str, Any]]:
        """Return per-step element deltas: what atoms must be added to which side.

        Agent4 calls this first to know exactly what staging atoms are needed.
        Pure arithmetic from reactant/product element counts — no hardcoded rules.

        Returns
        -------
        list of dict, each with:
            step_name, reactant_formula, product_formula,
            to_add_to_reactant: list of element symbols,
            to_add_to_product: list of element symbols,
            total_staged: int
        """
        owner = getattr(self.workflow, "tools", self.workflow)
        baseline_steps = list(getattr(owner, "_agent45_baseline_steps", []) or [])

        result: List[Dict[str, Any]] = []
        for entry in baseline_steps:
            name = str(entry.get("name", ""))
            reactant = entry.get("reactant")
            product = entry.get("product")
            if not isinstance(reactant, Atoms) or not isinstance(product, Atoms):
                result.append({"step_name": name, "error": "missing structure"})
                continue

            r_cnt = Counter(reactant.get_chemical_symbols())
            p_cnt = Counter(product.get_chemical_symbols())
            all_elems = sorted(set(list(r_cnt) + list(p_cnt)))

            to_add_r: List[str] = []
            to_add_p: List[str] = []
            for elem in all_elems:
                diff = p_cnt.get(elem, 0) - r_cnt.get(elem, 0)
                if diff > 0:
                    to_add_r.extend([elem] * diff)
                elif diff < 0:
                    to_add_p.extend([elem] * (-diff))

            result.append({
                "step_name": name,
                "reactant_formula": str(entry.get("reactant_formula", "")),
                "product_formula": str(entry.get("product_formula", "")),
                "to_add_to_reactant": to_add_r,
                "to_add_to_product": to_add_p,
                "total_staged": len(to_add_r) + len(to_add_p),
            })
        return result

    def _site_retry_owner(self) -> Any:
        owner = getattr(self.workflow, "tools", None)
        return owner if owner is not None else self.workflow

    def _site_retry_state(self) -> Dict[Tuple[str, ...], int]:
        owner = self._site_retry_owner()
        state = getattr(owner, "_agent45_site_retry_state", None)
        if isinstance(state, dict):
            return state
        state = {}
        setattr(owner, "_agent45_site_retry_state", state)
        return state

    def reset_site_retry_state(self) -> None:
        """Clear all retained site retry ranks (call between independent runs)."""
        self._site_retry_state().clear()

    def _site_retry_run_scope(self) -> str:
        """Best-effort run identifier so retry ranks don't bleed across runs."""
        for owner in (self.workflow, getattr(self.workflow, "tools", None)):
            if owner is None:
                continue
            for attr in ("run_id", "_current_run_id"):
                value = getattr(owner, attr, None)
                if value:
                    return str(value)
            state = getattr(owner, "state", None)
            state_run_id = getattr(state, "run_id", None)
            if state_run_id:
                return str(state_run_id)
        return ""

    def _site_retry_key(self, step_name: str, side: str, fragment_label: str) -> Tuple[str, ...]:
        base = (str(step_name), str(side), str(fragment_label))
        scope = self._site_retry_run_scope()
        # Prefix the run id when available so retry ranks never bleed across
        # reactions/runs sharing one tools object.
        return (scope,) + base if scope else base

    def _get_allowed_site_rank(
        self,
        step_name: str,
        side: str,
        fragment_label: str,
        candidate_sites: List[Dict[str, Any]],
    ) -> int:
        if not candidate_sites:
            return 1
        state = self._site_retry_state()
        # Exact fragment-label equality only — no substring fallback (a bare
        # "H" must not inherit retry ranks recorded for "OH"/"CH" fragments).
        raw_rank = state.get(self._site_retry_key(step_name, side, fragment_label), None)
        if raw_rank is None:
            raw_rank = 1
        try:
            rank = int(raw_rank)
        except Exception:
            rank = 1
        return max(1, min(rank, len(candidate_sites)))

    def _get_allowed_site_rank_for_hint(
        self,
        step_entry: Dict[str, Any],
        hint: Dict[str, Any],
    ) -> int:
        return self._get_allowed_site_rank(
            step_name=str(step_entry.get("name", "")),
            side=str(hint.get("side", "reactant")),
            fragment_label=str(hint.get("fragment_label", "")),
            candidate_sites=list(hint.get("candidate_sites", [])),
        )

    def _promote_allowed_site_rank_for_hint(
        self,
        step_entry: Dict[str, Any],
        hint: Dict[str, Any],
        next_rank: int,
    ) -> None:
        candidate_sites = list(hint.get("candidate_sites", []))
        if not candidate_sites:
            return
        state = self._site_retry_state()
        key = self._site_retry_key(
            str(step_entry.get("name", "")),
            str(hint.get("side", "reactant")),
            str(hint.get("fragment_label", "")),
        )
        state[key] = max(1, min(int(next_rank), len(candidate_sites)))

    @staticmethod
    def _candidate_site_distance(site: Dict[str, Any]) -> float:
        return float(site.get("dist_from_anchor", site.get("distance_to_anchor", float("inf"))))

    @staticmethod
    def _candidate_site_label(site: Dict[str, Any]) -> str:
        return str(site.get("site_label", site.get("site_type", "")))

    @classmethod
    def _candidate_site_priority(cls, site: Dict[str, Any]) -> int:
        label = cls._candidate_site_label(site)
        if label.startswith(("local_H_on_O", "local_H_on_N", "local_H_on_S")):
            return 0
        if label.startswith("local_H_"):
            return 1
        return 2

    @classmethod
    def _candidate_site_is_h_acceptor(cls, site: Dict[str, Any]) -> bool:
        return cls._candidate_site_label(site).startswith(
            ("local_H_on_O", "local_H_on_N", "local_H_on_S")
        )

    @classmethod
    def _nearest_nonactive_local_h_site_distance(cls, candidate_sites: List[Dict[str, Any]]) -> float:
        return cls._nearest_nonactive_local_h_site_metrics(candidate_sites)[0]

    @classmethod
    def _nearest_nonactive_local_h_site_metrics(cls, candidate_sites: List[Dict[str, Any]]) -> Tuple[float, float]:
        best_dist = float("inf")
        best_z = float("inf")
        for site in candidate_sites:
            if bool(site.get("active_center_site", False)):
                continue
            if not cls._candidate_site_is_h_acceptor(site):
                continue
            dist = cls._candidate_site_distance(site)
            if dist >= best_dist:
                continue
            pos = site.get("position", [0.0, 0.0, float("inf")])
            try:
                z = float(pos[2])
            except Exception:
                z = float("inf")
            best_dist = dist
            best_z = z
        return best_dist, best_z

    @classmethod
    def _candidate_site_is_clearly_far_active_acceptor(
        cls,
        site: Dict[str, Any],
        *,
        nearest_nonactive_local_dist: float,
        nearest_nonactive_local_z: float,
    ) -> bool:
        if not bool(site.get("active_center_site", False)):
            return False
        if not cls._candidate_site_is_h_acceptor(site):
            return False
        if nearest_nonactive_local_dist == float("inf"):
            return False
        dist = cls._candidate_site_distance(site)
        try:
            site_z = float(site.get("position", [0.0, 0.0, float("inf")])[2])
        except Exception:
            site_z = float("inf")
        nonactive_is_same_local_shell = (
            nearest_nonactive_local_z < float("inf")
            and site_z < float("inf")
            and (nearest_nonactive_local_z - site_z) <= 2.0
        )
        return bool(
            nonactive_is_same_local_shell
            and (
                (dist >= 3.50 and (dist - nearest_nonactive_local_dist) >= 0.75)
                or (
                    nearest_nonactive_local_z <= site_z + 2.50
                    and (dist - nearest_nonactive_local_dist) >= 0.25
                )
            )
        )

    @classmethod
    def _candidate_site_active_bucket(
        cls,
        site: Dict[str, Any],
        *,
        nearest_nonactive_local_dist: float,
        nearest_nonactive_local_z: float = float("inf"),
        active_acceptor_demoted: bool = False,
    ) -> int:
        """Prefer active-center H sites unless they are clearly not local.

        Curved SMSI pockets can contain a chemically active Ti/O pair that is
        several Angstrom away from the outgoing C-H atom after a carried Stage B
        spectator.  In that case forcing the active O as site 1 makes the staged
        H jump across the pocket and breaks NEB atom mapping.  A much closer
        local H-on-O/Ti site should be offered first; the active-center sites
        remain available as later candidates.
        """
        is_active = bool(site.get("active_center_site", False))
        priority = cls._candidate_site_priority(site)
        clearly_far_active = cls._candidate_site_is_clearly_far_active_acceptor(
            site,
            nearest_nonactive_local_dist=nearest_nonactive_local_dist,
            nearest_nonactive_local_z=nearest_nonactive_local_z,
        )
        if is_active and not clearly_far_active:
            if active_acceptor_demoted and priority == 1:
                return 2
            return 0
        if not is_active and priority <= 1:
            return 1
        if is_active:
            return 2
        return 3

    @staticmethod
    def _rank_candidate_sites_by_direction(
        candidate_sites: List[Dict[str, Any]],
        anchor_pos: Optional[np.ndarray],
        preferred_reference_pos: Optional[np.ndarray],
    ) -> List[Dict[str, Any]]:
        if not candidate_sites:
            return []

        ranked = [dict(site) for site in candidate_sites]
        nearest_nonactive_local, nearest_nonactive_local_z = (
            Agent45GeometryTools._nearest_nonactive_local_h_site_metrics(ranked)
        )
        active_acceptor_demoted = any(
            Agent45GeometryTools._candidate_site_is_clearly_far_active_acceptor(
                site,
                nearest_nonactive_local_dist=nearest_nonactive_local,
                nearest_nonactive_local_z=nearest_nonactive_local_z,
            )
            for site in ranked
        )
        if anchor_pos is None or preferred_reference_pos is None:
            ranked.sort(
                key=lambda site: (
                    Agent45GeometryTools._candidate_site_active_bucket(
                        site,
                        nearest_nonactive_local_dist=nearest_nonactive_local,
                        nearest_nonactive_local_z=nearest_nonactive_local_z,
                        active_acceptor_demoted=active_acceptor_demoted,
                    ),
                    Agent45GeometryTools._candidate_site_priority(site),
                    Agent45GeometryTools._candidate_site_distance(site),
                    float(site.get("position", [0.0, 0.0, 0.0])[0]),
                    float(site.get("position", [0.0, 0.0, 0.0])[1]),
                )
            )
        else:
            anchor = np.array(anchor_pos, dtype=float)
            preferred = np.array(preferred_reference_pos, dtype=float) - anchor
            preferred[2] = 0.0
            preferred_norm = float(np.linalg.norm(preferred))
            if preferred_norm <= 1e-8:
                ranked.sort(
                    key=lambda site: (
                        Agent45GeometryTools._candidate_site_active_bucket(
                            site,
                            nearest_nonactive_local_dist=nearest_nonactive_local,
                            nearest_nonactive_local_z=nearest_nonactive_local_z,
                            active_acceptor_demoted=active_acceptor_demoted,
                        ),
                        Agent45GeometryTools._candidate_site_priority(site),
                        Agent45GeometryTools._candidate_site_distance(site),
                        float(site.get("position", [0.0, 0.0, 0.0])[0]),
                        float(site.get("position", [0.0, 0.0, 0.0])[1]),
                    )
                )
            else:
                preferred = preferred / preferred_norm

                def _alignment(site: Dict[str, Any]) -> float:
                    pos = np.array(site.get("position", []), dtype=float)
                    if pos.shape != (3,):
                        return -1.0
                    direction = pos - anchor
                    direction[2] = 0.0
                    norm = float(np.linalg.norm(direction))
                    if norm <= 1e-8:
                        return -1.0
                    return float(np.dot(direction / norm, preferred))

                ranked.sort(
                    key=lambda site: (
                        Agent45GeometryTools._candidate_site_active_bucket(
                            site,
                            nearest_nonactive_local_dist=nearest_nonactive_local,
                            nearest_nonactive_local_z=nearest_nonactive_local_z,
                            active_acceptor_demoted=active_acceptor_demoted,
                        ),
                        Agent45GeometryTools._candidate_site_priority(site),
                        round(Agent45GeometryTools._candidate_site_distance(site), 2),
                        -_alignment(site),
                        float(site.get("position", [0.0, 0.0, 0.0])[0]),
                        float(site.get("position", [0.0, 0.0, 0.0])[1]),
                    )
                )

        for rank, site in enumerate(ranked, start=1):
            site["rank"] = rank
        return ranked

    def _local_surface_top_z(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
        position: np.ndarray,
        *,
        radius: float = 3.0,
        fallback: Optional[float] = None,
    ) -> float:
        positions = structure.get_positions()
        ads_set = set(adsorbate_indices or [])
        surf_indices = [i for i in range(len(structure)) if i not in ads_set]
        if not surf_indices:
            return float(fallback if fallback is not None else 0.0)

        local: List[int] = []
        for idx in surf_indices:
            try:
                dist = float(self.workflow._pbc_distance(structure, position, positions[idx]))
            except Exception:
                dist = float(np.linalg.norm(position - positions[idx]))
            if dist <= radius:
                local.append(idx)

        if local:
            return float(np.max(positions[local, 2]))
        if fallback is not None:
            return float(fallback)
        return float(np.max(positions[surf_indices, 2]))

    def _enumerate_local_h_surface_sites(
        self,
        target_atoms: Atoms,
        target_ads: List[int],
        anchor_pos: np.ndarray,
        *,
        max_sites: int,
        active_center_indices: Optional[List[int]] = None,
    ) -> List[Dict[str, Any]]:
        """Generate local H-on-surface candidates around a bond-change anchor.

        Pymatgen adsorption sites are generic molecular adsorption positions and
        can put a single H too high above an oxide overlayer.  For staged H, add
        local atop candidates on nearby exposed Ti/O/metal atoms, then let the
        existing Agent4/5 validation and energy gate decide whether they pass.
        """
        positions = target_atoms.get_positions()
        symbols = target_atoms.get_chemical_symbols()
        ads_set = set(target_ads or [])
        surf_indices = [i for i in range(len(target_atoms)) if i not in ads_set]
        if not surf_indices:
            return []

        surf_top_z = float(np.max(positions[surf_indices, 2]))
        active_centers = [
            int(idx)
            for idx in (active_center_indices or [])
            if 0 <= int(idx) < len(target_atoms) and int(idx) in surf_indices and symbols[int(idx)] not in {"H", "C"}
        ]
        exposed: List[int] = []
        if active_centers:
            local_pool: set[int] = set(active_centers)
            for center_idx in active_centers:
                center_pos = positions[center_idx]
                for idx in surf_indices:
                    if symbols[idx] in {"H", "C"}:
                        continue
                    try:
                        dist = float(self.workflow._pbc_distance(target_atoms, center_pos, positions[idx]))
                    except Exception:
                        dist = float(np.linalg.norm(center_pos - positions[idx]))
                    if dist <= 3.2:
                        local_pool.add(idx)
            exposed = sorted(local_pool, key=lambda idx: (idx not in active_centers, symbols[idx] != "O", idx))
        if not exposed:
            exposed = [
                idx for idx in surf_indices
                if positions[idx, 2] >= surf_top_z - 2.5 and symbols[idx] not in {"H", "C"}
            ]
        if not exposed:
            exposed = [idx for idx in surf_indices if symbols[idx] not in {"H", "C"}]

        height_by_symbol = {
            "O": [0.98, 1.20, 1.60],
            "N": [1.02, 1.25],
            "S": [1.35, 1.60],
            "Ti": [1.80, 2.10, 2.40],
        }
        default_heights = [1.60, 1.90, 2.20]

        ads_positions = [positions[i] for i in target_ads if 0 <= i < len(target_atoms)]
        ads_h_positions = [
            positions[i]
            for i in target_ads
            if 0 <= i < len(target_atoms) and symbols[i] == "H"
        ]
        candidates: List[Dict[str, Any]] = []
        for idx in exposed:
            sym = symbols[idx]
            for height in height_by_symbol.get(sym, default_heights):
                site_pos = np.array(
                    [positions[idx, 0], positions[idx, 1], positions[idx, 2] + float(height)],
                    dtype=float,
                )
                local_top_z = self._local_surface_top_z(
                    target_atoms,
                    target_ads,
                    site_pos,
                    fallback=float(np.max(positions[exposed, 2])) if active_centers and exposed else surf_top_z,
                )
                is_active_center_site = bool(idx in active_centers)
                is_proton_acceptor_site = sym in {"O", "N", "S"}
                if site_pos[2] < local_top_z + 0.5 and not is_proton_acceptor_site:
                    continue
                try:
                    dist_anchor = float(self.workflow._pbc_distance(target_atoms, anchor_pos, site_pos))
                except Exception:
                    dist_anchor = float(np.linalg.norm(anchor_pos - site_pos))
                if not (1.2 <= dist_anchor <= 6.5):
                    continue

                min_any = float("inf")
                for atom_pos in positions:
                    try:
                        d_any = float(self.workflow._pbc_distance(target_atoms, site_pos, atom_pos))
                    except Exception:
                        d_any = float(np.linalg.norm(site_pos - atom_pos))
                    min_any = min(min_any, d_any)
                if min_any < 0.8:
                    continue

                too_close_to_existing_h = False
                for h_pos in ads_h_positions:
                    try:
                        h_dist = float(self.workflow._pbc_distance(target_atoms, site_pos, h_pos))
                    except Exception:
                        h_dist = float(np.linalg.norm(site_pos - h_pos))
                    if h_dist < 1.70:
                        too_close_to_existing_h = True
                        break
                if too_close_to_existing_h:
                    continue

                min_ads = float("inf")
                for atom_pos in ads_positions:
                    try:
                        d_ads = float(self.workflow._pbc_distance(target_atoms, site_pos, atom_pos))
                    except Exception:
                        d_ads = float(np.linalg.norm(site_pos - atom_pos))
                    min_ads = min(min_ads, d_ads)
                if ads_positions and min_ads > 5.0:
                    continue

                candidates.append(
                    {
                        "position": [float(v) for v in site_pos],
                        "site_label": f"local_H_on_{sym}",
                        "dist_from_anchor": dist_anchor,
                        "active_center_site": is_active_center_site,
                        "_host_symbol": sym,
                        "_host_index": idx,
                        "_min_ads": min_ads,
                    }
                )

        nearest_nonactive_local, nearest_nonactive_local_z = self._nearest_nonactive_local_h_site_metrics(candidates)
        active_acceptor_demoted = any(
            self._candidate_site_is_clearly_far_active_acceptor(
                site,
                nearest_nonactive_local_dist=nearest_nonactive_local,
                nearest_nonactive_local_z=nearest_nonactive_local_z,
            )
            for site in candidates
        )
        candidates.sort(
            key=lambda site: (
                self._candidate_site_active_bucket(
                    site,
                    nearest_nonactive_local_dist=nearest_nonactive_local,
                    nearest_nonactive_local_z=nearest_nonactive_local_z,
                    active_acceptor_demoted=active_acceptor_demoted,
                ),
                self._candidate_site_priority(site),
                float(site.get("dist_from_anchor", float("inf"))),
                float(site.get("_min_ads", float("inf"))),
                float(site.get("position", [0.0, 0.0, 0.0])[2]),
            )
        )
        for site in candidates:
            site.pop("_min_ads", None)
            site.pop("_host_symbol", None)
            site.pop("_host_index", None)
        return candidates[:max_sites]

    def _has_physical_surface_context(
        self,
        structure: Atoms,
        adsorbate_indices: List[int],
    ) -> bool:
        positions = structure.get_positions()
        ads_set = set(adsorbate_indices or [])
        surf_indices = [i for i in range(len(structure)) if i not in ads_set]
        valid_ads = [i for i in adsorbate_indices if 0 <= i < len(structure)]
        if not surf_indices or not valid_ads:
            return False
        surf_top_z = float(np.max(positions[surf_indices, 2]))
        ads_min_z = float(np.min(positions[valid_ads, 2]))
        return bool(surf_top_z >= ads_min_z - 1.5)

    def _build_step_staging_hints(self, step_entry: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Build staging hints for one step using geometry-based RMSD matching.

        Internally delegates to ``find_bond_change_sites`` (Phase-2 Q3=A,
        element-agnostic RMSD match) and reshapes the output to the legacy hint
        schema so downstream consumers (Agent4 prompts, ``suggest_staged_positions_tool``,
        ``_find_matching_staging_hint``) keep working unchanged.

        The existing pymatgen site-enumeration infrastructure
        (``enumerate_surface_sites_near_anchor``) is preserved and reused by
        ``find_bond_change_sites``.
        """
        reactant = step_entry.get("reactant")
        product = step_entry.get("product")
        if not isinstance(reactant, Atoms) or not isinstance(product, Atoms):
            return []

        reactant_ads = list(step_entry.get("reactant_adsorbate_indices", []))
        product_ads = list(step_entry.get("product_adsorbate_indices", []))
        active_center_indices = list(step_entry.get("active_center_indices", []))

        raw_hints = self.find_bond_change_sites(
            clean_reactant=reactant,
            clean_reactant_ads=reactant_ads,
            clean_product=product,
            clean_product_ads=product_ads,
            min_site_dist=2.5,
            max_site_dist=4.5,
            max_sites_per_hint=3,
        )
        if not raw_hints and active_center_indices:
            raw_hints = self._build_active_center_fallback_h_hints(
                reactant=reactant,
                reactant_ads=reactant_ads,
                product=product,
                product_ads=product_ads,
                active_center_indices=active_center_indices,
                max_sites_per_hint=3,
            )

        hints: List[Dict[str, Any]] = []
        for rh in raw_hints:
            element = str(rh.get("element", ""))
            if not element:
                continue
            side = str(rh.get("side", ""))
            source_side = str(rh.get("source_side", ""))
            target_atoms = reactant if side == "reactant" else product
            target_ads = reactant_ads if side == "reactant" else product_ads
            legacy_sites: List[Dict[str, Any]] = []
            raw_sites = list(rh.get("candidate_sites", []))
            injected_local_h_sites = False
            if element == "H" and self._has_physical_surface_context(target_atoms, target_ads):
                anchor_position = np.array(rh.get("anchor_position", []), dtype=float)
                if anchor_position.shape == (3,):
                    local_sites = self._enumerate_local_h_surface_sites(
                        target_atoms,
                        target_ads,
                        anchor_position,
                        max_sites=3,
                        active_center_indices=active_center_indices,
                    )
                    merged: Dict[Tuple[float, float, float], Dict[str, Any]] = {}
                    for site in local_sites:
                        pos = site.get("position", [])
                        if isinstance(pos, list) and len(pos) == 3:
                            merged[(round(float(pos[0]), 3), round(float(pos[1]), 3), round(float(pos[2]), 3))] = {
                                "position": [float(v) for v in pos],
                                "site_type": str(site.get("site_label", "?")),
                                "distance_to_anchor": float(site.get("dist_from_anchor", 0.0)),
                                "active_center_site": bool(site.get("active_center_site", False)),
                            }
                    for site in raw_sites:
                        pos = site.get("position", [])
                        if isinstance(pos, list) and len(pos) == 3:
                            merged.setdefault(
                                (round(float(pos[0]), 3), round(float(pos[1]), 3), round(float(pos[2]), 3)),
                                dict(site),
                            )
                    injected_local_h_sites = bool(local_sites)
                    if injected_local_h_sites:
                        raw_sites = sorted(
                            merged.values(),
                            key=lambda site: (
                                0 if bool(site.get("active_center_site", False)) else 1,
                                0 if str(site.get("site_type", "")).startswith("local_H_") else 1,
                                float(site.get("distance_to_anchor", site.get("dist_from_anchor", float("inf")))),
                                float(site.get("position", [0.0, 0.0, 0.0])[2]),
                            ),
                        )

            anchor_position = np.array(rh.get("anchor_position", []), dtype=float)
            bond_partner_position = rh.get("bond_partner_position")
            preferred_reference_pos = None
            rank_anchor_pos = anchor_position if anchor_position.shape == (3,) else None
            if isinstance(bond_partner_position, list) and len(bond_partner_position) == 3:
                bond_partner_arr = np.array(bond_partner_position, dtype=float)
                if source_side == "product":
                    rank_anchor_pos = bond_partner_arr
                    preferred_reference_pos = anchor_position if anchor_position.shape == (3,) else None
                else:
                    preferred_reference_pos = bond_partner_arr

            ranked_sites = self._rank_candidate_sites_by_direction(
                candidate_sites=raw_sites,
                anchor_pos=rank_anchor_pos,
                preferred_reference_pos=preferred_reference_pos,
            )
            for rank, site in enumerate(ranked_sites, start=1):
                legacy_sites.append(
                    {
                        "rank": rank,
                        "position": [round(float(v), 2) for v in site.get("position", [])],
                        "site_label": str(site.get("site_type", site.get("site_label", "?"))),
                        "dist_from_anchor": round(
                            float(site.get("distance_to_anchor", site.get("dist_from_anchor", 0.0))),
                            2,
                        ),
                        "active_center_site": bool(site.get("active_center_site", False)),
                    }
                )
            hints.append(
                {
                    "side": side,
                    "fragment_elements": [element],
                    "fragment_label": element,
                    "source_side": source_side,
                    "source_atom_index": rh.get("source_atom_index"),
                    "bond_partner_index": rh.get("bond_partner_index"),
                    "anchor_symbol": str(rh.get("bond_partner_element") or "?"),
                    "candidate_sites": legacy_sites,
                }
            )
        return hints

    def _build_active_center_fallback_h_hints(
        self,
        *,
        reactant: Atoms,
        reactant_ads: List[int],
        product: Atoms,
        product_ads: List[int],
        active_center_indices: List[int],
        max_sites_per_hint: int,
    ) -> List[Dict[str, Any]]:
        """Build H staging hints from active-center sites when generic hints are empty."""
        r_syms = reactant.get_chemical_symbols()
        p_syms = product.get_chemical_symbols()
        r_h = [int(i) for i in reactant_ads if 0 <= int(i) < len(reactant) and r_syms[int(i)] == "H"]
        p_h = [int(i) for i in product_ads if 0 <= int(i) < len(product) and p_syms[int(i)] == "H"]
        if len(r_h) == len(p_h):
            return []

        if len(r_h) > len(p_h):
            source_atoms = reactant
            source_ads = reactant_ads
            target_atoms = product
            target_ads = product_ads
            source_side = "reactant"
            target_side = "product"
            source_h = r_h
            target_h = p_h
        else:
            source_atoms = product
            source_ads = product_ads
            target_atoms = reactant
            target_ads = reactant_ads
            source_side = "product"
            target_side = "reactant"
            source_h = p_h
            target_h = r_h

        target_positions = target_atoms.get_positions()
        source_positions = source_atoms.get_positions()
        unmatched: List[Tuple[float, int]] = []
        for h_idx in source_h:
            if target_h:
                nearest = min(
                    float(self.workflow._pbc_distance(source_atoms, source_positions[h_idx], target_positions[j]))
                    for j in target_h
                )
            else:
                nearest = float("inf")
            unmatched.append((nearest, h_idx))
        unmatched.sort(reverse=True)
        if not unmatched:
            return []

        hints: List[Dict[str, Any]] = []
        for _, h_idx in unmatched[: abs(len(r_h) - len(p_h))]:
            anchor_pos = np.array(source_positions[h_idx], dtype=float)
            sites = self._enumerate_local_h_surface_sites(
                target_atoms,
                target_ads,
                anchor_pos,
                max_sites=max_sites_per_hint,
                active_center_indices=active_center_indices,
            )
            if not sites:
                continue
            candidate_sites = [
                {
                    "position": [float(v) for v in site.get("position", [])],
                    "site_type": str(site.get("site_label", "unknown")),
                    "distance_to_anchor": float(site.get("dist_from_anchor", 0.0)),
                    "active_center_site": bool(site.get("active_center_site", False)),
                }
                for site in sites[:max_sites_per_hint]
            ]
            if not candidate_sites:
                continue
            hints.append(
                {
                    "element": "H",
                    "side": target_side,
                    "source_side": source_side,
                    "source_atom_index": int(h_idx),
                    "anchor_position": [float(v) for v in anchor_pos],
                    "bond_partner_element": "H",
                    "bond_partner_index": None,
                    "bond_partner_position": [float(v) for v in anchor_pos],
                    "candidate_sites": candidate_sites,
                }
            )
        return hints

    def _ensure_step_staging_hints(self, step_entry: Dict[str, Any]) -> List[Dict[str, Any]]:
        hints = step_entry.get("staging_hints")
        if isinstance(hints, list):
            return hints
        hints = self._build_step_staging_hints(step_entry)
        step_entry["staging_hints"] = hints
        return hints

    def _find_matching_staging_hint(
        self,
        step_entry: Dict[str, Any],
        side: str,
        element: str,
    ) -> Optional[Dict[str, Any]]:
        for hint in self._ensure_step_staging_hints(step_entry):
            if str(hint.get("side", "")) != side:
                continue
            fragment_elements = list(hint.get("fragment_elements", []))
            if fragment_elements and str(fragment_elements[0]) == element:
                return hint
        return None

    def _nearest_candidate_site(
        self,
        structure: Atoms,
        candidate_sites: List[Dict[str, Any]],
        position: np.ndarray,
    ) -> Tuple[Optional[Dict[str, Any]], float]:
        if not candidate_sites:
            return None, float("inf")

        valid_sites: List[Dict[str, Any]] = []
        valid_positions: List[np.ndarray] = []
        for site in candidate_sites:
            site_pos = np.array(site.get("position", []), dtype=float)
            if site_pos.shape != (3,):
                continue
            valid_sites.append(site)
            valid_positions.append(site_pos)

        if not valid_sites:
            return None, float("inf")

        best_idx, _, best_dist = self._nearest_reference_position(structure, position, valid_positions)
        if best_idx is None:
            return None, float("inf")
        return valid_sites[best_idx], best_dist

    @staticmethod
    def _candidate_site_rank(
        candidate_sites: List[Dict[str, Any]],
        site: Optional[Dict[str, Any]],
    ) -> int:
        if site is None or not candidate_sites:
            return 1
        try:
            rank = int(site.get("rank"))
            if 1 <= rank <= len(candidate_sites):
                return rank
        except Exception:
            pass

        site_pos = list(site.get("position", []))
        for idx, candidate in enumerate(candidate_sites, start=1):
            if candidate is site:
                return idx
            if list(candidate.get("position", [])) == site_pos:
                return idx
        return 1

    def suggest_staged_positions(
        self,
        requests: List[Dict[str, str]],
    ) -> Dict[str, Any]:
        """Suggest positions for ALL staged atoms across ALL steps in one call.

        Parameters
        ----------
        requests : list of dict
            Each dict: ``{"step_name": "...", "element": "H", "side": "reactant"}``.

        Returns
        -------
        dict with "suggestions": list of per-request results, each with
        position, anchor_atom, bond_distance, reason.
        """
        # Try multiple sources for baseline steps
        baseline_steps: List[Dict[str, Any]] = []
        for src in [self.workflow, getattr(self.workflow, "tools", None)]:
            if src is not None:
                candidate = getattr(src, "_agent45_baseline_steps", None)
                if candidate:
                    baseline_steps = list(candidate)
                    break
        step_lookup: Dict[str, Any] = {}
        for e in baseline_steps:
            name = str(e.get("name", ""))
            step_lookup[name] = e
            step_lookup[name.replace("*", "")] = e

        suggestions: List[Dict[str, Any]] = []
        for req in (requests or []):
            step_name = str(req.get("step_name", "")).replace("*", "")
            element = str(req.get("element", ""))
            side = str(req.get("side", "reactant"))

            step_entry = step_lookup.get(step_name)
            if step_entry is None:
                suggestions.append({"step_name": step_name, "error": f"Step not found"})
                continue

            atoms = step_entry.get(side)
            if not isinstance(atoms, Atoms):
                suggestions.append({"step_name": step_name, "error": f"No {side} structure"})
                continue

            ads_indices = list(step_entry.get(f"{side}_adsorbate_indices", []))

            hint = self._find_matching_staging_hint(step_entry, side, element)
            if not hint:
                suggestions.append(
                    {
                        "step_name": step_name,
                        "element": element,
                        "side": side,
                        "error": (
                            f"No staging hint found for {element} on {side} side. "
                            "find_bond_change_sites did not identify this element as a "
                            "bond-change atom; check the clean reactant/product pair."
                        ),
                    }
                )
                continue
            candidate_sites = list(hint.get("candidate_sites", []))
            if not candidate_sites:
                suggestions.append(
                    {
                        "step_name": step_name,
                        "element": element,
                        "side": side,
                        "error": (
                            f"Hint for {element} on {side} has no candidate surface sites. "
                            "pymatgen site enumerator returned empty for the bond-change anchor."
                        ),
                    }
                )
                continue

            allowed_rank = self._get_allowed_site_rank_for_hint(step_entry, hint)
            pos = list(candidate_sites[allowed_rank - 1]["position"])
            suggestions.append({
                "step_name": step_name,
                "element": element,
                "side": side,
                "position": [round(float(v), 4) for v in pos],
                "anchor": f"candidate site {allowed_rank} near {hint.get('anchor_symbol', '?')}",
                "candidates": candidate_sites,
            })

        return {"suggestions": suggestions}

    def validate_proposed_positions(
        self,
        steps: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Validate proposed atom positions for ALL steps in one call.

        Parameters
        ----------
        steps : list of dict
            Each dict: ``{"step_name": "...", "atoms": [{"species": "H", "side": "reactant", "position": [x,y,z]}]}``.

        Returns
        -------
        dict with: ok (bool), per_step (list of per-step results), suggestion (str)
        """
        baseline_steps: List[Dict[str, Any]] = []
        for src in [self.workflow, getattr(self.workflow, "tools", None)]:
            if src is not None:
                candidate = getattr(src, "_agent45_baseline_steps", None)
                if candidate:
                    baseline_steps = list(candidate)
                    break
        step_lookup: Dict[str, Any] = {}
        for e in baseline_steps:
            name = str(e.get("name", ""))
            step_lookup[name] = e
            step_lookup[name.replace("*", "")] = e

        per_step: List[Dict[str, Any]] = []
        all_ok = True

        for step_req in (steps or []):
            step_name = str(step_req.get("step_name", "")).replace("*", "")
            proposed_atoms = list(step_req.get("atoms", []))
            step_entry = step_lookup.get(step_name)

            if step_entry is None:
                per_step.append({"step_name": step_name, "ok": False, "issue": "Step not found"})
                all_ok = False
                continue

            step_ok = True
            atom_results: List[Dict[str, Any]] = []

            for atom_spec in proposed_atoms:
                element = str(atom_spec.get("species", atom_spec.get("element", "?")))
                side = str(atom_spec.get("side", "reactant"))
                position = atom_spec.get("position", [])

                if not (isinstance(position, list) and len(position) == 3):
                    atom_results.append({"element": element, "side": side, "ok": False, "issue": "Invalid [x,y,z]"})
                    step_ok = False
                    continue

                try:
                    pos = np.array([float(position[0]), float(position[1]), float(position[2])])
                except (TypeError, ValueError):
                    atom_results.append({"element": element, "side": side, "ok": False, "issue": "Non-numeric position"})
                    step_ok = False
                    continue

                atoms: Optional[Atoms] = step_entry.get(side)
                if not isinstance(atoms, Atoms):
                    atom_results.append({"element": element, "side": side, "ok": False, "issue": f"No {side} structure"})
                    step_ok = False
                    continue

                ads_indices = list(step_entry.get(f"{side}_adsorbate_indices", []))
                positions = atoms.get_positions()
                symbols = atoms.get_chemical_symbols()

                # Distance to nearest adsorbate
                valid_ads = [i for i in ads_indices if 0 <= i < len(atoms)]
                min_ads_dist = float("inf")
                nearest_ads = "?"
                if valid_ads:
                    idx_rel, p, min_ads_dist = self._nearest_reference_position(
                        atoms,
                        pos,
                        [positions[idx] for idx in valid_ads],
                    )
                    idx = valid_ads[int(idx_rel)] if idx_rel is not None else valid_ads[0]
                    if p is None:
                        p = positions[idx]
                    nearest_ads = f"{symbols[idx]} at [{p[0]:.2f},{p[1]:.2f},{p[2]:.2f}]"

                # Overlap with any atom
                min_any_dist = float("inf")
                nearest_any = "?"
                if len(positions) > 0:
                    idx, p, min_any_dist = self._nearest_reference_position(atoms, pos, positions)
                    if idx is None:
                        idx = 0
                    if p is None:
                        p = positions[idx]
                    nearest_any = f"{symbols[idx]} at [{p[0]:.2f},{p[1]:.2f},{p[2]:.2f}]"

                hint = self._find_matching_staging_hint(step_entry, side, element)

                # Surface z
                surf_indices = [i for i in range(len(atoms)) if i not in set(ads_indices)]
                surf_top_z = float(np.max(positions[surf_indices, 2])) if surf_indices else 0.0
                local_top_z = surf_top_z
                is_active_recessed_acceptor_site = False
                if hint:
                    candidate_sites = list(hint.get("candidate_sites", []))
                    nearest_site, site_dist = self._nearest_candidate_site(atoms, candidate_sites, pos)
                    if nearest_site is not None and site_dist <= 1.5:
                        site_label = str(nearest_site.get("site_label", nearest_site.get("label", "")))
                        is_active_recessed_acceptor_site = bool(
                            nearest_site.get("active_center_site", False)
                            and site_label.startswith("local_H_on_O")
                        )
                        local_top_z = self._local_surface_top_z(
                            atoms,
                            ads_indices,
                            pos,
                            fallback=surf_top_z,
                        )

                atom_ok = True
                issues: List[str] = []
                if min_ads_dist > 5.0:
                    atom_ok = False
                    issues.append(f"Too far ({min_ads_dist:.2f} Å) from nearest adsorbate {nearest_ads}")
                if min_any_dist < 0.8:
                    atom_ok = False
                    issues.append(f"Overlaps ({min_any_dist:.3f} Å) with {nearest_any}")
                if pos[2] < local_top_z + 0.5 and not is_active_recessed_acceptor_site:
                    atom_ok = False
                    issues.append(f"Below local surface (z={pos[2]:.2f}, local surface top={local_top_z:.2f})")

                if hint:
                    candidate_sites = list(hint.get("candidate_sites", []))
                    allowed_rank = self._get_allowed_site_rank_for_hint(step_entry, hint)
                    allowed_site = candidate_sites[allowed_rank - 1] if candidate_sites else None
                    nearest_site, site_dist = self._nearest_candidate_site(atoms, candidate_sites, pos)
                    nearest_rank = self._candidate_site_rank(candidate_sites, nearest_site)
                    if nearest_site is not None and site_dist > 1.5:
                        atom_ok = False
                        site_pos = nearest_site["position"]
                        if allowed_rank <= 1:
                            issues.append(
                                f"Not at a suggested candidate site: nearest candidate site {nearest_rank} "
                                f"is [{site_pos[0]:.2f},{site_pos[1]:.2f},{site_pos[2]:.2f}] "
                                f"({site_dist:.2f} Å away). Use site 1 first; only switch to site 2 or 3 "
                                f"after Agent5/energy-gate feedback."
                            )
                        else:
                            issues.append(
                                f"Not at the currently allowed candidate site: nearest candidate site {nearest_rank} "
                                f"is [{site_pos[0]:.2f},{site_pos[1]:.2f},{site_pos[2]:.2f}] "
                                f"({site_dist:.2f} Å away). This step is currently assigned to site {allowed_rank}."
                            )
                    elif nearest_site is not None and allowed_site is not None and nearest_rank != allowed_rank:
                        atom_ok = False
                        allowed_pos = allowed_site["position"]
                        if allowed_rank == 1:
                            issues.append(
                                f"Use suggested candidate site 1 [{allowed_pos[0]:.2f},{allowed_pos[1]:.2f},{allowed_pos[2]:.2f}] "
                                f"on the first attempt. site {nearest_rank} is only allowed after explicit "
                                f"Agent5/energy-gate failure for this step."
                            )
                        else:
                            issues.append(
                                f"This step has already escalated to candidate site {allowed_rank}. "
                                f"Use [{allowed_pos[0]:.2f},{allowed_pos[1]:.2f},{allowed_pos[2]:.2f}] now."
                            )

                if not atom_ok:
                    step_ok = False

                atom_results.append({
                    "element": element, "side": side, "ok": atom_ok,
                    "dist_ads": round(min_ads_dist, 2), "dist_any": round(min_any_dist, 2),
                    "issues": issues,
                })

            # --- Interpolation collision check ---
            # Build temporary structures with staged atoms appended, then check
            # if linear interpolation between reactant and product causes collisions.
            interp_issue = ""
            if step_ok and proposed_atoms:
                try:
                    reactant_base = step_entry.get("reactant")
                    product_base = step_entry.get("product")
                    if isinstance(reactant_base, Atoms) and isinstance(product_base, Atoms):
                        tmp_r = reactant_base.copy()
                        tmp_p = product_base.copy()
                        r_ads = list(step_entry.get("reactant_adsorbate_indices", []))
                        p_ads = list(step_entry.get("product_adsorbate_indices", []))

                        for a in proposed_atoms:
                            elem = str(a.get("species", a.get("element", "")))
                            side = str(a.get("side", ""))
                            apos = a.get("position", [])
                            if not (isinstance(apos, list) and len(apos) == 3):
                                continue
                            apos_f = [float(apos[0]), float(apos[1]), float(apos[2])]
                            if side == "reactant":
                                tmp_r.append(Atom(elem, position=apos_f))
                                r_ads.append(len(tmp_r) - 1)
                            elif side == "product":
                                tmp_p.append(Atom(elem, position=apos_f))
                                p_ads.append(len(tmp_p) - 1)

                        # Only check if both endpoints have same element count
                        if (len(tmp_r) == len(tmp_p)
                                and sorted(tmp_r.get_chemical_symbols()) == sorted(tmp_p.get_chemical_symbols())):
                            interp_mins = self.workflow._best_interpolation_mapping_metrics(
                                tmp_r, tmp_p, r_ads, p_ads, n_images=5,
                            )
                            finite = [float(v) for v in interp_mins if np.isfinite(v) and v > 0]
                            if finite:
                                min_interp = min(finite)
                                if min_interp < 0.50:
                                    step_ok = False
                                    interp_issue = (
                                        f"Interpolation collision: min separation {min_interp:.2f} Å along "
                                        f"NEB path (< 0.50 Å). Staged atoms cross existing atoms during "
                                        f"linear interpolation. Reposition staged atoms so they do not "
                                        f"pass through other atoms between reactant and product positions."
                                    )
                                elif min_interp < 0.80:
                                    interp_issue = (
                                        f"Interpolation warning: min separation {min_interp:.2f} Å "
                                        f"(< 0.80 Å). NEB may still work but positions could be improved."
                                    )
                except Exception as exc:
                    interp_issue = f"Interpolation check failed: {exc}"

            if interp_issue:
                for ar in atom_results:
                    ar.setdefault("issues", []).append(interp_issue)

            if not step_ok:
                all_ok = False
            per_step.append({
                "step_name": step_name, "ok": step_ok, "atoms": atom_results,
                "interpolation": interp_issue if interp_issue else "OK",
            })

        suggestion = "All positions OK." if all_ok else "Some positions need adjustment — see per_step details."
        return {"ok": all_ok, "per_step": per_step, "suggestion": suggestion}

    def _build_gas_phase_product_from_current(self, current: Atoms) -> Atoms:
        """Create a synthetic desorbed product for gas-phase endpoints."""
        product = current.copy()
        adsorbate_indices = list(product.info.get("adsorbate_indices", []))
        if not adsorbate_indices:
            return product

        positions = product.get_positions()
        surface_top_z = self.workflow._estimate_surface_top_z(product, adsorbate_indices)
        adsorbate_com_z = float(np.mean(positions[adsorbate_indices, 2]))
        cell = np.array(product.cell, dtype=float)
        c_len = float(np.linalg.norm(cell[2])) if cell.shape == (3, 3) else 0.0

        target_z = max(
            surface_top_z + 4.5,
            float(np.max(positions[:, 2]) + 2.0),
        )
        if c_len > 0:
            # Keep gas-phase endpoint inside the simulation cell vacuum region.
            upper_bound = c_len - 1.2
            target_z = min(target_z, upper_bound)
            if target_z <= surface_top_z + 2.5:
                target_z = min(upper_bound, surface_top_z + 3.2)
        dz = float(target_z - adsorbate_com_z)
        positions[adsorbate_indices, 2] = positions[adsorbate_indices, 2] + dz
        product.set_positions(positions)
        return product

    def _resolve_fixed_product_for_step(
        self,
        current: Atoms,
        preplaced_products: Dict[str, Dict[str, Any]],
        product_label: str,
    ) -> Tuple[Atoms, str, str, List[int], List[int]]:
        fixed_entry = self.workflow._get_preplaced_entry(preplaced_products, product_label)
        if fixed_entry and isinstance(fixed_entry.get("structure"), Atoms):
            fixed_product = fixed_entry["structure"].copy()
            product_source = str(fixed_entry.get("source", "tool_preplacement"))
            product_path = str(fixed_entry.get("path", ""))
            surface_indices = list(fixed_entry.get("surface_indices", []))
            adsorbate_indices = list(fixed_entry.get("adsorbate_indices", []))
            if not surface_indices and not adsorbate_indices:
                raise RuntimeError(
                    f"Tool-preplaced structure for '{product_label}' lacks surface/adsorbate indices"
                )
            return fixed_product, product_source, product_path, surface_indices, adsorbate_indices

        if self.workflow._contains_gas_phase_hint(product_label):
            fallback_product = self._build_gas_phase_product_from_current(current)
            surface_indices = list(fallback_product.info.get("surface_indices", []))
            adsorbate_indices = list(fallback_product.info.get("adsorbate_indices", []))
            if not surface_indices and not adsorbate_indices:
                surface_indices = [i for i in range(len(fallback_product)) if i not in set(adsorbate_indices)]
            return fallback_product, "synthetic_gas_phase", "", surface_indices, adsorbate_indices

        raise RuntimeError(f"Missing tool-preplaced product structure for step product '{product_label}'")

    @staticmethod
    def _pymatgen_sort(
        atoms: Atoms,
        ads_indices: List[int],
        staged_indices: List[int],
    ) -> Tuple[Atoms, List[int], List[int]]:
        """Sort structure by element (atomic number) using pymatgen for consistent NEB ordering.

        After sorting, atoms are grouped by element (H, C, O, Cu, …) with a
        deterministic sub-order (fractional z, y, x).  Adsorbate / staged roles
        are tracked via ``site.properties`` so they survive the sort.
        """
        from pymatgen.io.ase import AseAtomsAdaptor

        adaptor = AseAtomsAdaptor()
        struct = adaptor.get_structure(atoms)

        # Tag each site with its original index so we can rebuild ads/staged lists
        ads_set = set(ads_indices)
        staged_set = set(staged_indices)
        for i, site in enumerate(struct):
            if i in staged_set:
                site.properties["_role"] = 2   # staged (also adsorbate)
            elif i in ads_set:
                site.properties["_role"] = 1   # adsorbate
            else:
                site.properties["_role"] = 0   # surface

        struct.sort()  # sorts by species (atomic number), then frac coords

        sorted_atoms = adaptor.get_atoms(struct)
        # Copy non-index metadata only: index-bearing keys recorded before the
        # sort would point at stale pre-sort indices. adsorbate_indices /
        # surface_indices are rewritten below from the tracked roles.
        carried_info = dict(atoms.info)
        for stale_key in (
            "adsorbate_indices",
            "surface_indices",
            "staged_indices",
            "_reorder_sequence",
            "atom_mapping_constraints",
        ):
            carried_info.pop(stale_key, None)
        sorted_atoms.info = carried_info

        new_ads: List[int] = []
        new_staged: List[int] = []
        for new_idx, site in enumerate(struct):
            role = site.properties.get("_role", 0)
            if role >= 1:
                new_ads.append(new_idx)
            if role == 2:
                new_staged.append(new_idx)

        new_ads.sort()
        new_staged.sort()
        ads_set_new = set(new_ads)
        sorted_atoms.info["adsorbate_indices"] = new_ads
        sorted_atoms.info["surface_indices"] = [
            i for i in range(len(sorted_atoms)) if i not in ads_set_new
        ]
        return sorted_atoms, new_ads, new_staged

    def _sort_endpoints_by_element(
        self,
        reactant: Atoms,
        product: Atoms,
        reactant_ads: List[int],
        product_ads: List[int],
        reactant_staged: List[int],
        product_staged: List[int],
    ) -> Tuple[Atoms, Atoms, List[int], List[int], List[int], List[int]]:
        """Sort both NEB endpoints by element (pymatgen) for consistent ordering."""
        r, ra, rs = self._pymatgen_sort(reactant, reactant_ads, reactant_staged)
        p, pa, ps = self._pymatgen_sort(product, product_ads, product_staged)
        return r, p, ra, pa, rs, ps

    def build_steps_payload(self, step_structures: List[Dict[str, Any]]) -> str:
        """Build a compact per-step summary with only adsorbate coordinates.

        Redundancy reduction:
        - Step i product = Step i+1 reactant, so only product adsorbate is shown per step
          (step 1 also includes reactant adsorbate as the starting point)
        - Internal metadata (signatures, fixed atom counts) omitted
        - Compact xyz format instead of POSCAR
        """
        blocks: List[str] = []
        for idx, step in enumerate(step_structures or [], start=1):
            step_name = str(step.get("name", f"step_{idx}"))
            reactant = step.get("reactant")
            product = step.get("product")
            r_formula = step.get("reactant_formula", "")
            p_formula = step.get("product_formula", "")

            # Compute element delta (what Agent4 needs to add/remove)
            delta_str = "none"
            if isinstance(reactant, Atoms) and isinstance(product, Atoms):
                r_cnt = Counter(reactant.get_chemical_symbols())
                p_cnt = Counter(product.get_chemical_symbols())
                parts = []
                for elem in sorted(set(list(r_cnt) + list(p_cnt))):
                    diff = p_cnt.get(elem, 0) - r_cnt.get(elem, 0)
                    if diff != 0:
                        parts.append(f"{'+' if diff > 0 else ''}{diff}{elem}")
                delta_str = " ".join(parts) if parts else "none"

            header = f"=== STEP {idx}: {step_name} ({r_formula} -> {p_formula}), element_delta={delta_str} ==="
            blocks.append(header)

            # Show BOTH reactant and product adsorbate for every step.
            # Agent4 needs both endpoints to decide where to place new atoms.
            if isinstance(reactant, Atoms):
                r_ads = list(step.get("reactant_adsorbate_indices", []))
                blocks.append(self._compact_adsorbate_xyz(reactant, r_ads, "reactant"))

            if isinstance(product, Atoms):
                p_ads = list(step.get("product_adsorbate_indices", []))
                blocks.append(self._compact_adsorbate_xyz(product, p_ads, "product"))

            if delta_str != "none" and isinstance(reactant, Atoms) and isinstance(product, Atoms):
                hints: List[str] = []
                for hint in self._ensure_step_staging_hints(step):
                    side = str(hint.get("side", "reactant"))
                    frag_label = str(hint.get("fragment_label", ""))
                    anchor_sym = str(hint.get("anchor_symbol", "?"))
                    candidate_sites = list(hint.get("candidate_sites", []))
                    allowed_rank = self._get_allowed_site_rank_for_hint(step, hint)

                    if candidate_sites:
                        if allowed_rank <= 1:
                            hints.append(
                                f"  [hint] Add {frag_label} to {side} at nearest-neighbor surface site "
                                f"in bonding direction from {anchor_sym}. Use site 1 first; "
                                f"retry site 2/3 only after failure. "
                                f"Candidate sites (nearest to {anchor_sym}):"
                            )
                        else:
                            hints.append(
                                f"  [hint] Add {frag_label} to {side} at nearest-neighbor surface site "
                                f"from {anchor_sym}. Previous failure → use site {allowed_rank}. "
                                f"Candidate sites:"
                            )
                        for site in candidate_sites[:3]:
                            p = site["position"]
                            hints.append(
                                f"    site {site['rank']}: [{p[0]:.2f},{p[1]:.2f},{p[2]:.2f}] "
                                f"{site['site_label']}, {site['dist_from_anchor']:.1f} Å from {anchor_sym}"
                            )
                    else:
                        hints.append(
                            f"  [hint] Add {frag_label} to {side}: place near {anchor_sym} anchor at the nearest "
                            "real co-adsorption site within 4.5 Å."
                        )

                if hints:
                    blocks.extend(hints)

        return "\n".join(blocks) if blocks else "(none)"

    def _compact_adsorbate_xyz(self, atoms: Atoms, ads_indices: List[int], label: str) -> str:
        """Compact adsorbate coordinate listing: one line per atom."""
        lines = [f"  [{label}_adsorbate] indices={ads_indices}"]
        for i in ads_indices:
            if 0 <= i < len(atoms):
                sym = atoms[i].symbol
                x, y, z = atoms.positions[i]
                lines.append(f"    {i}: {sym}  {x:.4f}  {y:.4f}  {z:.4f}")
        return "\n".join(lines)

    def build_agent5_structures_payload(self, state: WorkflowState) -> str:
        return self.build_steps_payload(state.step_structures)

    def _energy_gate_owner(self) -> Any:
        owner = getattr(self.workflow, "tools", None)
        if owner is not None:
            return owner
        return self.workflow

    def _get_uma_predictor_for_energy_gate(self, output_dir: str = "") -> Any:
        owner = self._energy_gate_owner()
        shared_getter = getattr(owner, "get_shared_uma_predictor", None)
        if callable(shared_getter):
            return shared_getter()

        predictor = getattr(owner, "_agent45_uma_predictor", None)
        if predictor is not None:
            return predictor

        from core.pathway.fairchem_predictor import FairchemPredictor

        predictor_work_dir = Path(output_dir).resolve() if output_dir else Path(getattr(owner, "output_base_dir", ".")).resolve()
        predictor_work_dir = predictor_work_dir / "_agent45_energy_gate"
        predictor_work_dir.mkdir(parents=True, exist_ok=True)

        predictor = FairchemPredictor(
            fairchem_root=str(DEPS_BASE_PATH / "fairchem"),
            model_name=DEFAULT_FAIRCHEM_MODEL,
            use_gpu=True,
            device="cuda",
            work_dir=str(predictor_work_dir),
            keep_files=False,
            verbose=False,
        )
        setattr(owner, "_agent45_uma_predictor", predictor)
        return predictor

    @staticmethod
    def _select_relaxable_indices(
        adsorbate_indices: List[int],
        staged_indices: List[int],
        fixed_atom_count: int,
        total_atoms: int,
    ) -> List[int]:
        """Return ONLY staged atom indices for relaxation.

        Baseline adsorbate core atoms are NOT relaxed — their positions are
        determined by tool preplacement and must remain immutable so that
        NEB endpoints reflect the intended reaction path.
        Only staged atoms (added by Agent4 for element balancing) are relaxed.
        """
        staged_set = {idx for idx in staged_indices if isinstance(idx, int) and 0 <= idx < total_atoms}
        return sorted(staged_set)

    def _constrained_subset_relax(
        self,
        atoms: Atoms,
        movable_indices: List[int],
        predictor: Any,
        fmax: float,
        max_steps: int,
    ) -> Tuple[Atoms, Dict[str, Any]]:
        if not movable_indices:
            return atoms.copy(), {
                "ran": False,
                "steps": 0,
                "reason": "no_movable_indices",
            }

        relaxed = atoms.copy()
        predictor._load_model()
        relaxed.calc = predictor._calculator

        movable_set = set(movable_indices)
        fixed = [idx for idx in range(len(relaxed)) if idx not in movable_set]
        relaxed.set_constraint(FixAtoms(indices=fixed))

        optimizer = FIRE(relaxed, logfile=None)
        ran = True
        try:
            optimizer.run(fmax=float(fmax), steps=max(1, int(max_steps)))
        except Exception as exc:
            return atoms.copy(), {
                "ran": ran,
                "steps": int(optimizer.get_number_of_steps()),
                "error": str(exc),
            }
        finally:
            try:
                relaxed.set_constraint()
            except Exception:
                pass
            relaxed.calc = None

        return relaxed, {
            "ran": ran,
            "steps": int(optimizer.get_number_of_steps()),
            "error": "",
        }

    def _evaluate_endpoint_with_uma(
        self,
        predictor: Any,
        atoms: Atoms,
        adsorbate_indices: List[int],
        focus_indices: Optional[List[int]] = None,
    ) -> Dict[str, float]:
        result = predictor.predict_energy(
            structure=atoms.copy(),
            relax=False,
        )
        forces = np.array(result.forces, dtype=float)
        force_norm = np.linalg.norm(forces, axis=1) if forces.size else np.zeros((len(atoms),), dtype=float)

        valid_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(atoms)]
        if valid_ads:
            ads_forces = force_norm[valid_ads]
            max_force_ads = float(np.max(ads_forces))
            mean_force_ads = float(np.mean(ads_forces))
        else:
            max_force_ads = float(np.max(force_norm)) if len(force_norm) else 0.0
            mean_force_ads = float(np.mean(force_norm)) if len(force_norm) else 0.0

        valid_focus = [
            idx for idx in (focus_indices or [])
            if 0 <= idx < len(atoms)
        ]
        if valid_focus:
            focus_forces = force_norm[valid_focus]
            max_force_focus = float(np.max(focus_forces))
            mean_force_focus = float(np.mean(focus_forces))
        else:
            max_force_focus = max_force_ads
            mean_force_focus = mean_force_ads

        return {
            "energy": float(result.energy),
            "max_force_all": float(np.max(force_norm)) if len(force_norm) else 0.0,
            "max_force_adsorbate": max_force_ads,
            "mean_force_adsorbate": mean_force_ads,
            "max_force_focus": max_force_focus,
            "mean_force_focus": mean_force_focus,
            "min_pair_distance": float(self.workflow._min_distance_any_pair(atoms)),
            "min_adsorbate_surface_distance": float(
                self.workflow._min_adsorbate_surface_distance(atoms, valid_ads)
            ),
        }

    def run_agent45_energy_gate(
        self,
        step_structures: List[Dict[str, Any]],
        output_dir: str = "",
        min_pair_fatal_threshold: float = 0.80,
        force_fatal_threshold: float = 1.0,
        force_warning_threshold: float = 0.5,
    ) -> Dict[str, Any]:
        """Agent4/5 endpoint UMA energy gate: evaluate energies and forces.

        Staged atoms (added by Agent4 for element balancing) are first pre-relaxed with
        every other atom held fixed (FIRE, 40 steps, fmax 0.3 eV/Å; SI Table S3), and the
        relaxed positions are written back to the step. Adsorbate core and slab stay exactly
        as Agent4 designed them. Energies, forces and overlaps are then evaluated.
        Set CATDT_AGENT4_STAGED_PRERELAX=0 to evaluate the designs strictly as-is.

        Parameters
        ----------
        min_pair_fatal_threshold : float
            Minimum allowed pair distance (Å); below this is FATAL overlap.
        force_fatal_threshold : float
            Max adsorbate force (eV/Å) at or above which the endpoint is FATAL.
        force_warning_threshold : float
            Max adsorbate force (eV/Å) at or above which a warning is issued.
        """
        if not step_structures:
            return {"status": "PASS", "summary": "(no step structures)",
                    "fatal_issues": [], "warning_issues": [], "steps": []}

        predictor = self._get_uma_predictor_for_energy_gate(output_dir=output_dir)
        predictor._load_model()

        fatal_issues: List[str] = []
        warning_issues: List[str] = []
        summaries: List[str] = []
        step_reports: List[Dict[str, Any]] = []

        for step_idx, step in enumerate(step_structures, start=1):
            step_name = str(step.get("name", f"step_{step_idx}"))
            step_report: Dict[str, Any] = {"step_index": step_idx, "step_name": step_name, "endpoints": {}}

            for endpoint in ("reactant", "product"):
                atoms = step.get(endpoint) if isinstance(step, dict) else None
                if not isinstance(atoms, Atoms):
                    fatal_issues.append(
                        f"{step_name}:{endpoint} malformed payload — missing ASE Atoms structure"
                    )
                    step_report["endpoints"][endpoint] = {"error": "missing_structure"}
                    continue
                adsorbate_indices = list(step.get(f"{endpoint}_adsorbate_indices", []))

                staged = [int(i) for i in step.get(f"{endpoint}_staged_indices", []) or []
                          if isinstance(i, (int, np.integer)) and 0 <= int(i) < len(atoms)]
                prerelax_info: Dict[str, Any] = {"ran": False, "steps": 0}
                if staged and os.getenv("CATDT_AGENT4_STAGED_PRERELAX", "1") != "0":
                    relaxed, prerelax_info = self._constrained_subset_relax(
                        atoms=atoms, movable_indices=staged, predictor=predictor,
                        fmax=float(os.getenv("CATDT_AGENT4_STAGED_PRERELAX_FMAX", "0.3")),
                        max_steps=int(os.getenv("CATDT_AGENT4_STAGED_PRERELAX_STEPS", "40")),
                    )
                    if not prerelax_info.get("error"):
                        # keep Agent4's design: the staging-site rule judges the design,
                        # the physical checks judge the relaxed structure
                        step[f"{endpoint}_design_positions"] = np.asarray(atoms.get_positions()).copy()
                        # Only staged atoms may have moved. After product reordering a staged atom
                        # can sit inside the fixed-count range covered by the immutable signature,
                        # so re-sign that range, but only after checking every other atom is unmoved.
                        others = [i for i in range(len(atoms)) if i not in set(staged)]
                        unmoved = bool(np.allclose(relaxed.positions[others], atoms.positions[others], atol=1e-8))
                        prerelax_info["non_staged_unmoved"] = unmoved
                        sig_key = f"{endpoint}_fixed_reference_signature"
                        if unmoved and step.get(sig_key):
                            fixed_count = int(step.get(f"{endpoint}_fixed_atom_count", len(relaxed)))
                            step[sig_key] = self.workflow._fixed_coordinate_signature(relaxed, fixed_count)
                            prerelax_info["signature_updated"] = True
                        relaxed.info.update(atoms.info)
                        step[endpoint] = relaxed
                        atoms = relaxed

                result = self._evaluate_endpoint_with_uma(
                    predictor=predictor, atoms=atoms,
                    adsorbate_indices=adsorbate_indices,
                    focus_indices=adsorbate_indices,
                )

                # Overlap check
                min_pair = result["min_pair_distance"]
                if 0 < min_pair < float(min_pair_fatal_threshold):
                    fatal_issues.append(
                        f"{step_name}:{endpoint} min pair distance {min_pair:.2f} Å "
                        f"< {float(min_pair_fatal_threshold):.2f} Å"
                    )

                max_force = max(
                    float(result.get("max_force_adsorbate", 0.0)),
                    float(result.get("max_force_focus", 0.0)),
                )
                if max_force >= float(force_fatal_threshold):
                    fatal_issues.append(
                        f"{step_name}:{endpoint} max adsorbate force {max_force:.2f} eV/Å "
                        f">= {float(force_fatal_threshold):.2f} eV/Å"
                    )
                elif max_force >= float(force_warning_threshold):
                    warning_issues.append(
                        f"{step_name}:{endpoint} max adsorbate force {max_force:.2f} eV/Å"
                    )

                summaries.append(
                    f"{step_name}:{endpoint} E={result['energy']:.3f} eV, "
                    f"fmax={result['max_force_adsorbate']:.2f} eV/Å, "
                    f"min_pair={min_pair:.2f} Å"
                )
                result["staged_prerelax"] = {"staged_indices": staged, **prerelax_info}
                step_report["endpoints"][endpoint] = result

            step_reports.append(step_report)

        status = "FAIL" if fatal_issues else "PASS"
        report = {
            "status": status,
            "summary": "\n".join(summaries),
            "fatal_issues": fatal_issues,
            "warning_issues": warning_issues,
            "steps": step_reports,
        }

        if output_dir:
            try:
                gate_path = Path(output_dir) / "agent45_energy_gate.json"
                gate_path.parent.mkdir(parents=True, exist_ok=True)
                with open(gate_path, "w", encoding="utf-8") as fh:
                    json.dump(report, fh, ensure_ascii=False, indent=2)
            except Exception as exc:
                logger.debug("Failed to persist agent45 energy gate report: %s", exc)

        return report

    def build_tool_baseline_steps(
        self,
        state: WorkflowState,
        base_structure: Atoms,
        preplaced_products: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        intermediates = self.workflow._get_working_intermediate_sequence(state)
        if len(intermediates) < 2:
            return []

        baseline_steps: List[Dict[str, Any]] = []
        current = base_structure.copy()
        surface = self.workflow._get_surface()
        self.workflow._set_surface_indices(
            current,
            list(surface.surface_indices) if surface is not None else [],
            list(surface.adsorbate_indices) if surface is not None else [],
        )

        for idx in range(len(intermediates) - 1):
            reactant_label = intermediates[idx]
            product_label = intermediates[idx + 1]
            fixed_product, _, _, surface_indices, adsorbate_indices = self._resolve_fixed_product_for_step(
                current=current,
                preplaced_products=preplaced_products,
                product_label=product_label,
            )
            self.workflow._set_surface_indices(fixed_product, surface_indices, adsorbate_indices)

            reactant = current.copy()
            product = fixed_product.copy()

            reactant_ads = list(reactant.info.get("adsorbate_indices", []))
            product_ads = list(product.info.get("adsorbate_indices", []))

            step_name = f"{reactant_label}_to_{product_label}"
            entry = {
                "name": step_name,
                "reactant": reactant,
                "product": product,
                "reactant_formula": str(reactant_label),
                "product_formula": str(product_label),
                "reaction_type": "tool_baseline",
                "reactant_adsorbate_indices": reactant_ads,
                "product_adsorbate_indices": product_ads,
                "reactant_fixed_atom_count": int(len(reactant)),
                "product_fixed_atom_count": int(len(product)),
                "reactant_staged_indices": [],
                "product_staged_indices": [],
                "reactant_fixed_reference_signature": self.workflow._fixed_coordinate_signature(reactant, len(reactant)),
                "product_fixed_reference_signature": self.workflow._fixed_coordinate_signature(product, len(product)),
            }
            entry["staging_hints"] = self._build_step_staging_hints(entry)
            baseline_steps.append(entry)
            current = fixed_product.copy()

        return baseline_steps

    def build_single_step_baseline(
        self,
        state: WorkflowState,
        current_structure: Atoms,
        preplaced_products: Dict[str, Dict[str, Any]],
        step_idx: int,
    ) -> Dict[str, Any]:
        """Build baseline for a single step. Used by sequential step-by-step flow.

        Product is built by transplanting the preplaced adsorbate onto the
        reactant's surface, ensuring both endpoints share the same surface.
        """
        intermediates = self.workflow._get_working_intermediate_sequence(state)
        if step_idx < 0 or step_idx >= len(intermediates) - 1:
            raise IndexError(f"step_idx {step_idx} out of range for {len(intermediates)} intermediates")

        reactant_label = intermediates[step_idx]
        product_label = intermediates[step_idx + 1]

        fixed_product, _, _, prod_surface_indices, prod_adsorbate_indices = self._resolve_fixed_product_for_step(
            current=current_structure,
            preplaced_products=preplaced_products,
            product_label=product_label,
        )

        reactant = current_structure.copy()
        reactant_ads = list(reactant.info.get("adsorbate_indices", []))
        reactant_surf = list(reactant.info.get("surface_indices", []))

        # Build product by: reactant's surface + preplaced adsorbate (aligned to reactant)
        product = reactant.copy()
        for idx in sorted(reactant_ads, reverse=True):
            if 0 <= idx < len(product):
                del product[idx]

        # Collect preplaced adsorbate atoms
        prod_ads_positions = []
        prod_ads_symbols = []
        for idx in prod_adsorbate_indices:
            if 0 <= idx < len(fixed_product):
                prod_ads_positions.append(fixed_product.positions[idx].copy())
                prod_ads_symbols.append(fixed_product[idx].symbol)

        # Align product adsorbate to reactant adsorbate anchor (lowest-z heavy atom)
        if prod_ads_positions and reactant_ads:
            r_pos = reactant.get_positions()
            r_sym = reactant.get_chemical_symbols()
            r_heavy = [i for i in reactant_ads if 0 <= i < len(reactant) and r_sym[i] != "H"]
            p_heavy = [j for j, s in enumerate(prod_ads_symbols) if s != "H"]

            if r_heavy and p_heavy:
                from ase.data import atomic_numbers as _an
                common_elems = set(r_sym[i] for i in r_heavy) & set(prod_ads_symbols[j] for j in p_heavy)
                if common_elems:
                    anchor_elem = max(common_elems, key=lambda e: _an.get(e, 0))
                    r_anchor_idx = min(
                        [i for i in r_heavy if r_sym[i] == anchor_elem],
                        key=lambda i: r_pos[i, 2],
                    )
                    p_anchor_j = min(
                        [j for j in p_heavy if prod_ads_symbols[j] == anchor_elem],
                        key=lambda j: prod_ads_positions[j][2],
                    )
                    shift = r_pos[r_anchor_idx] - prod_ads_positions[p_anchor_j]
                    if np.linalg.norm(shift) > 0.1:
                        for j in range(len(prod_ads_positions)):
                            prod_ads_positions[j] = prod_ads_positions[j] + shift

        new_ads_indices = []
        for sym, pos in zip(prod_ads_symbols, prod_ads_positions):
            product.append(Atom(sym, position=pos))
            new_ads_indices.append(len(product) - 1)
        product_ads = new_ads_indices
        product_surf = list(range(len(reactant_surf)))
        self.workflow._set_surface_indices(product, product_surf, product_ads)

        step_name = f"{reactant_label}_to_{product_label}"
        entry = {
            "name": step_name,
            "reactant": reactant,
            "product": product,
            "reactant_formula": str(reactant_label),
            "product_formula": str(product_label),
            "reaction_type": "tool_baseline",
            "reactant_adsorbate_indices": reactant_ads,
            "product_adsorbate_indices": product_ads,
            "reactant_fixed_atom_count": int(len(reactant)),
            "product_fixed_atom_count": int(len(product)),
            "reactant_staged_indices": [],
            "product_staged_indices": [],
            "reactant_fixed_reference_signature": self.workflow._fixed_coordinate_signature(reactant, len(reactant)),
            "product_fixed_reference_signature": self.workflow._fixed_coordinate_signature(product, len(product)),
        }
        entry["staging_hints"] = self._build_step_staging_hints(entry)
        return entry

    def apply_step_ops_on_fixed_product(
        self,
        current_physical: Atoms,
        fixed_product: Atoms,
        step: PathwayStepSpec,
    ) -> Dict[str, Any]:
        """Build one step endpoints while preserving immutable core coordinates.

        Agent4 must provide explicit atoms_to_add with positions. No hardcoded fallback.
        """
        reactant_neb = current_physical.copy()
        product_neb = fixed_product.copy()

        reactant_ads = list(reactant_neb.info.get("adsorbate_indices", []))
        product_ads = list(product_neb.info.get("adsorbate_indices", []))

        reactant_fixed_len = len(reactant_neb)
        product_fixed_len = len(product_neb)
        reactant_fixed_snapshot = np.array(reactant_neb.positions[:reactant_fixed_len], dtype=float)
        product_fixed_snapshot = np.array(product_neb.positions[:product_fixed_len], dtype=float)

        reactant_staged_indices: List[int] = []
        product_staged_indices: List[int] = []
        atom_mapping_constraints: List[Dict[str, Any]] = []

        has_llm_additions = bool(step.atoms_to_add)

        # Compute element delta from ADSORBATE atoms only — not the full
        # structure which includes slab. Using full-structure counts can
        # give wrong staging direction when the slab transplant in
        # build_single_step_baseline produces slightly different slab
        # compositions between reactant and product.
        r_syms = reactant_neb.get_chemical_symbols()
        p_syms = product_neb.get_chemical_symbols()
        base_reactant_counts = Counter(r_syms[i] for i in reactant_ads if 0 <= i < len(reactant_neb))
        base_product_counts = Counter(p_syms[i] for i in product_ads if 0 <= i < len(product_neb))
        pending_staging: List[Tuple[str, int]] = []
        for element in sorted(set(base_reactant_counts.keys()) | set(base_product_counts.keys())):
            delta = int(base_product_counts.get(element, 0) - base_reactant_counts.get(element, 0))
            if delta != 0:
                pending_staging.append((element, delta))
        if pending_staging:
            logger.info(
                "Step '%s': adsorbate delta=%s (R_ads=%s, P_ads=%s)",
                step.step_name, pending_staging,
                dict(base_reactant_counts), dict(base_product_counts),
            )

        if pending_staging and not has_llm_additions:
            raise RuntimeError(
                f"Step '{step.step_name}' requires atoms_to_add for element deltas {pending_staging}, "
                "but Agent4 returned no additions."
            )

        llm_positions: Dict[str, Dict[str, List[List[float]]]] = {"reactant": {}, "product": {}}

        def _append_llm_position(side: str, element: str, value: Any) -> None:
            if not (isinstance(value, list) and len(value) == 3):
                return
            try:
                pos = [float(value[0]), float(value[1]), float(value[2])]
            except Exception:
                return
            if not all(np.isfinite(v) for v in pos):
                return
            llm_positions[side].setdefault(element, []).append(pos)

        if has_llm_additions:
            for atom_spec in step.atoms_to_add:
                elements = self.workflow._resolve_addition_elements(atom_spec)
                if not elements:
                    continue

                base_candidate = atom_spec.position
                reactant_candidate = atom_spec.reactant_position if atom_spec.reactant_position is not None else base_candidate
                product_candidate = atom_spec.product_position if atom_spec.product_position is not None else base_candidate

                for element in elements:
                    _append_llm_position("reactant", element, reactant_candidate)
                    _append_llm_position("product", element, product_candidate)

        incoming_offsets: Counter = Counter()
        outgoing_offsets: Counter = Counter()
        step_entry = {
            "name": step.step_name,
            "reactant": reactant_neb,
            "product": product_neb,
            "reactant_adsorbate_indices": reactant_ads,
            "product_adsorbate_indices": product_ads,
        }
        for src in [self.workflow, getattr(self.workflow, "tools", None)]:
            if src is None:
                continue
            for baseline in list(getattr(src, "_agent45_baseline_steps", []) or []):
                if not isinstance(baseline, dict):
                    continue
                if str(baseline.get("name", "")) != str(step.step_name):
                    continue
                hints = baseline.get("staging_hints")
                if isinstance(hints, list):
                    step_entry["staging_hints"] = [
                        {
                            **dict(hint),
                            "fragment_elements": list(hint.get("fragment_elements", [])),
                            "candidate_sites": [dict(site) for site in list(hint.get("candidate_sites", []))],
                        }
                        for hint in hints
                        if isinstance(hint, dict)
                    ]
                break
            if "staging_hints" in step_entry:
                break

        def _pop_position(source: Dict[str, Dict[str, List[List[float]]]], side: str, element: str) -> Optional[List[float]]:
            slots = source.get(side, {}).get(element, [])
            if slots:
                return slots.pop(0)
            return None

        for element, delta in pending_staging:
            if delta > 0:
                for _ in range(delta):
                    incoming_offsets[element] += 1

                    reactant_candidate = (
                        _pop_position(llm_positions, "reactant", element)
                        or _pop_position(llm_positions, "product", element)
                    )
                    if reactant_candidate is None:
                        raise RuntimeError(
                            f"Step '{step.step_name}': Agent4 must provide position for "
                            f"adding {element} to reactant (element delta +{delta}). "
                            f"Use suggest_staged_positions_tool if unsure."
                        )

                    reactant_pos = self.workflow._safe_position(
                        reactant_candidate,
                        reactant_candidate,
                        structure=reactant_neb,
                    )

                    reactant_neb.append(Atom(element, position=reactant_pos))
                    reactant_ads.append(len(reactant_neb) - 1)
                    reactant_staged_indices.append(len(reactant_neb) - 1)
                    hint = self._find_matching_staging_hint(step_entry, "reactant", element)
                    if hint and str(hint.get("source_side", "")) == "product":
                        source_idx = hint.get("source_atom_index")
                        if isinstance(source_idx, int):
                            atom_mapping_constraints.append(
                                {
                                    "reactant_index": len(reactant_neb) - 1,
                                    "product_index": int(source_idx),
                                    "element": element,
                                    "reason": "Stage B staged atom maps to the unmatched Stage A atom from the clean product.",
                                }
                            )
            else:
                for _ in range(-delta):
                    outgoing_offsets[element] += 1

                    product_candidate = (
                        _pop_position(llm_positions, "product", element)
                        or _pop_position(llm_positions, "reactant", element)
                    )
                    if product_candidate is None:
                        raise RuntimeError(
                            f"Step '{step.step_name}': Agent4 must provide position for "
                            f"adding {element} to product (element delta {delta}). "
                            f"Use suggest_staged_positions_tool if unsure."
                        )

                    product_pos = self.workflow._safe_position(
                        product_candidate,
                        product_candidate,
                        structure=product_neb,
                    )

                    product_neb.append(Atom(element, position=product_pos))
                    product_ads.append(len(product_neb) - 1)
                    product_staged_indices.append(len(product_neb) - 1)
                    hint = self._find_matching_staging_hint(step_entry, "product", element)
                    if hint and str(hint.get("source_side", "")) == "reactant":
                        source_idx = hint.get("source_atom_index")
                        if isinstance(source_idx, int):
                            atom_mapping_constraints.append(
                                {
                                    "reactant_index": int(source_idx),
                                    "product_index": len(product_neb) - 1,
                                    "element": element,
                                    "reason": "Stage B staged atom maps to the unmatched Stage A atom from the clean reactant.",
                                }
                            )

        ignored_additions = sum(
            len(pos_list)
            for side_positions in llm_positions.values()
            for pos_list in side_positions.values()
        )
        if ignored_additions > 0:
            logger.warning(
                "Step '%s': ignored %d redundant atoms_to_add entries beyond required element deltas.",
                step.step_name,
                ignored_additions,
            )

        remove_indices = self.workflow._resolve_removal_indices(step, reactant_ads)
        if remove_indices:
            logger.warning(
                "Step '%s': atoms_to_remove should be empty per prompt constraints, "
                "but Agent4 specified removal of %d atoms. Proceeding with caution.",
                step.step_name, len(remove_indices),
            )
            illegal = [idx for idx in remove_indices if idx < reactant_fixed_len]
            if illegal:
                raise RuntimeError(
                    f"Step '{step.step_name}' attempts to remove immutable core atoms: {illegal}"
                )

            removed_records: List[Dict[str, Any]] = []
            for idx in sorted(set(remove_indices)):
                if 0 <= idx < len(reactant_neb):
                    removed_records.append(
                        {
                            "symbol": str(reactant_neb[idx].symbol),
                            "position": [
                                float(reactant_neb.positions[idx][0]),
                                float(reactant_neb.positions[idx][1]),
                                float(reactant_neb.positions[idx][2]),
                            ],
                        }
                    )

            reactant_ads, reactant_staged_indices = self.workflow._delete_atoms_with_index_tracking(
                structure=reactant_neb,
                adsorbate_indices=reactant_ads,
                staged_indices=reactant_staged_indices,
                delete_indices=remove_indices,
            )

            for offset, removed in enumerate(removed_records):
                fallback = self.workflow._get_unbonded_position(product_neb, product_ads)
                fallback[0] += 0.18 * offset
                fallback[1] -= 0.18 * offset
                position = self.workflow._safe_position(
                    removed.get("position"),
                    fallback,
                    structure=product_neb,
                )
                product_neb.append(Atom(str(removed["symbol"]), position=position))
                product_ads.append(len(product_neb) - 1)
                product_staged_indices.append(len(product_neb) - 1)

        # Agent4 positions are used as-is. No alignment or geometric relaxation.
        # Agent5 will validate and provide feedback if positions are problematic.

        reactant_surface = [i for i in range(len(reactant_neb)) if i not in set(reactant_ads)]
        product_surface = [i for i in range(len(product_neb)) if i not in set(product_ads)]
        self.workflow._set_surface_indices(reactant_neb, reactant_surface, reactant_ads)
        self.workflow._set_surface_indices(product_neb, product_surface, product_ads)

        if len(reactant_neb) < reactant_fixed_len or len(product_neb) < product_fixed_len:
            raise RuntimeError(f"Step '{step.step_name}' corrupted immutable core atom counts")

        if not np.allclose(reactant_neb.positions[:reactant_fixed_len], reactant_fixed_snapshot, atol=1e-10):
            raise RuntimeError(f"Step '{step.step_name}' changed immutable reactant core coordinates")
        if not np.allclose(product_neb.positions[:product_fixed_len], product_fixed_snapshot, atol=1e-10):
            raise RuntimeError(f"Step '{step.step_name}' changed immutable product core coordinates")

        if sorted(reactant_neb.get_chemical_symbols()) != sorted(product_neb.get_chemical_symbols()):
            raise RuntimeError(
                f"Step '{step.step_name}' endpoint element multisets differ after staged add/remove"
            )

        # Bounds-check + diagnostic before the return. Any out-of-range index
        # indicates a stale atoms.info carry from upstream — we log the full
        # state so the actual bug location is visible in the log.
        r_valid = [i for i in reactant_ads if 0 <= i < len(reactant_neb)]
        p_valid = [i for i in product_ads if 0 <= i < len(product_neb)]
        if len(r_valid) != len(reactant_ads) or len(p_valid) != len(product_ads):
            logger.error(
                "Step '%s': stale adsorbate indices. reactant_neb has %d atoms, "
                "reactant_ads=%s (valid: %s); product_neb has %d atoms, "
                "product_ads=%s (valid: %s). reactant_staged=%s product_staged=%s "
                "reactant_fixed_len=%d product_fixed_len=%d",
                step.step_name,
                len(reactant_neb), list(reactant_ads), r_valid,
                len(product_neb), list(product_ads), p_valid,
                list(reactant_staged_indices), list(product_staged_indices),
                reactant_fixed_len, product_fixed_len,
            )
            raise RuntimeError(
                f"Step '{step.step_name}': adsorbate indices out of range "
                f"after staging. See log above for details."
            )

        return {
            "reactant": reactant_neb,
            "product": product_neb,
            "reactant_adsorbate_indices": sorted(set(reactant_ads)),
            "product_adsorbate_indices": sorted(set(product_ads)),
            "reactant_fixed_atom_count": int(reactant_fixed_len),
            "product_fixed_atom_count": int(product_fixed_len),
            "reactant_staged_indices": sorted(set(reactant_staged_indices)),
            "product_staged_indices": sorted(set(product_staged_indices)),
            "atom_mapping_constraints": atom_mapping_constraints,
            "reactant_min_distance": self.workflow._min_pair_distance(reactant_neb, reactant_ads),
            "product_min_distance": self.workflow._min_pair_distance(product_neb, product_ads),
        }

    @staticmethod
    def _mic_xy(atoms: Atoms, vec_xy: Any) -> np.ndarray:
        """Shortest periodic image of a lateral vector (surface cell vectors a, b)."""
        a, b = np.asarray(atoms.cell[0][:2], float), np.asarray(atoms.cell[1][:2], float)
        v = np.asarray(vec_xy, float)
        best = v
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                cand = v + i * a + j * b
                if np.linalg.norm(cand) < np.linalg.norm(best):
                    best = cand
        return best

    def _equivalent_surface_translation(
        self, atoms: Atoms, surface_indices: List[int], wanted_xy: Any, tol: float = 0.25,
    ) -> Optional[np.ndarray]:
        """Lateral translation closest to ``wanted_xy`` that maps the slab onto itself.

        Candidates are the vectors between same-element atoms of the top slab layer; each is
        accepted only if every slab atom lands on a slab atom of the same element (minimum
        image, within ``tol``). Returns None when no candidate is a symmetry of the slab.
        """
        idx = [i for i in surface_indices if 0 <= i < len(atoms)]
        if not idx:
            return None
        pos = atoms.positions[idx]
        sym = [atoms[i].symbol for i in idx]
        top = [k for k in range(len(idx)) if pos[k, 2] > pos[:, 2].max() - 0.6]
        wanted = np.asarray(wanted_xy, float)
        cands = {}
        for k0 in top[:1]:
            for k in top:
                if sym[k] != sym[k0]:
                    continue
                v = self._mic_xy(atoms, pos[k, :2] - pos[k0, :2])
                cands[(round(v[0], 3), round(v[1], 3))] = v
        best, best_d = None, float("inf")
        for v in sorted(cands.values(), key=lambda v: float(np.linalg.norm(self._mic_xy(atoms, wanted - v)))):
            d = float(np.linalg.norm(self._mic_xy(atoms, wanted - v)))
            if d >= best_d:
                continue
            ok = True
            for k in range(len(idx)):
                target = pos[k].copy()
                target[:2] += v
                dists = [np.linalg.norm(np.r_[self._mic_xy(atoms, target[:2] - pos[m, :2]), target[2] - pos[m, 2]])
                         for m in range(len(idx)) if sym[m] == sym[k]]
                if not dists or min(dists) > tol:
                    ok = False
                    break
            if ok:
                best, best_d = v, d
                break
        return best

    def build_step_structures(
        self,
        state: WorkflowState,
        base_structure: Atoms,
        steps: List[PathwayStepSpec],
        preplaced_products: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        step_structures: List[Dict[str, Any]] = []
        baseline_lookup = {
            str(entry.get("name", "")): entry
            for entry in list(state.tool_baseline_steps or [])
            if isinstance(entry, dict)
        }

        current = base_structure.copy()
        surface = self.workflow._get_surface()
        self.workflow._set_surface_indices(
            current,
            list(surface.surface_indices) if surface is not None else [],
            list(surface.adsorbate_indices) if surface is not None else [],
        )

        for step_idx, step in enumerate(steps):
            fixed_product, product_source, product_path, surface_indices, adsorbate_indices = (
                self._resolve_fixed_product_for_step(
                    current=current,
                    preplaced_products=preplaced_products,
                    product_label=step.product_formula,
                )
            )
            self.workflow._set_surface_indices(fixed_product, surface_indices, adsorbate_indices)

            built = self.apply_step_ops_on_fixed_product(
                current_physical=current,
                fixed_product=fixed_product,
                step=step,
            )

            reactant_neb = built["reactant"]
            product_neb = built["product"]

            # Fix adsorbate drift: if product adsorbate centroid xy is far from
            # reactant, shift product adsorbate xy only (preserve z to avoid
            # embedding atoms into the surface for species with different
            # adsorption heights, e.g. CH3 on surface vs CH4 in gas phase).
            r_ads_idx = list(built["reactant_adsorbate_indices"])
            p_ads_idx = list(built["product_adsorbate_indices"])
            if r_ads_idx and p_ads_idx:
                r_centroid = np.mean([reactant_neb.positions[i] for i in r_ads_idx if i < len(reactant_neb)], axis=0)
                p_centroid = np.mean([product_neb.positions[i] for i in p_ads_idx if i < len(product_neb)], axis=0)
                # Minimum-image lateral offset: a product sitting next to the reactant across the
                # periodic boundary has not drifted.
                wanted_xy = self._mic_xy(product_neb, r_centroid[:2] - p_centroid[:2])
                xy_drift = float(np.linalg.norm(wanted_xy))
                if xy_drift > 2.0:
                    # Move the product adsorbate only by a translation that maps the slab onto
                    # itself, so it stays on an equivalent site (an arbitrary shift can push a
                    # relaxed adsorbate from a hollow onto a top site). No such translation
                    # (e.g. a defective surface): leave the product where it is.
                    surf_idx = [i for i in range(len(product_neb)) if i not in set(p_ads_idx)]
                    shift_xy = self._equivalent_surface_translation(product_neb, surf_idx, wanted_xy)
                    if shift_xy is not None and float(np.linalg.norm(shift_xy)) > 1e-6:
                        for idx in p_ads_idx:
                            if idx < len(product_neb):
                                product_neb.positions[idx][:2] += shift_xy
                    logger.info(
                        "Step '%s': product adsorbate xy drifted %.2f A (minimum image) from reactant; "
                        "%s",
                        step.step_name, xy_drift,
                        "shifted by the surface-equivalent translation (%.2f, %.2f)" % (shift_xy[0], shift_xy[1])
                        if shift_xy is not None else "no surface-equivalent translation, product left in place",
                    )

            # Reorder product atoms to match reactant ordering so that
            # interpolation and NEB see consistent atom indices.
            product_ads_built = list(built["product_adsorbate_indices"])
            product_staged_built = list(built["product_staged_indices"])
            mapping_constraints_built = list(built.get("atom_mapping_constraints", []))
            reordered_prod, reordered_ads, reordered_staged, did_reorder = self.workflow._reorder_product_to_match_reactant(
                reactant=reactant_neb,
                product=product_neb,
                product_adsorbate_indices=product_ads_built,
                product_staged_indices=product_staged_built,
                atom_mapping_constraints=mapping_constraints_built,
            )
            if did_reorder:
                product_neb = reordered_prod
                built["product_adsorbate_indices"] = reordered_ads
                built["product_staged_indices"] = reordered_staged
                old_to_new = {old: new for new, old in enumerate(reordered_prod.info.get("_reorder_sequence", []))}
                if old_to_new:
                    for constraint in mapping_constraints_built:
                        try:
                            old_p = int(constraint.get("product_index"))
                        except Exception:
                            continue
                        if old_p in old_to_new:
                            constraint["product_index"] = int(old_to_new[old_p])

            baseline_entry = baseline_lookup.get(step.step_name, {})
            baseline_hints = list(baseline_entry.get("staging_hints", [])) if isinstance(baseline_entry, dict) else []
            copied_hints = [
                {
                    **dict(hint),
                    "fragment_elements": list(hint.get("fragment_elements", [])),
                    "candidate_sites": [dict(site) for site in list(hint.get("candidate_sites", []))],
                }
                for hint in baseline_hints
            ]

            step_structures.append(
                {
                    "name": step.step_name,
                    "reactant": reactant_neb,
                    "product": product_neb,
                    "reactant_formula": step.reactant_formula,
                    "product_formula": step.product_formula,
                    "reaction_type": step.reaction_type,
                    "atoms_to_add_count": len(step.atoms_to_add),
                    "atoms_to_remove_count": len(step.atoms_to_remove),
                    "reactant_adsorbate_indices": list(built["reactant_adsorbate_indices"]),
                    "product_adsorbate_indices": list(built["product_adsorbate_indices"]),
                    "reactant_fixed_atom_count": int(built["reactant_fixed_atom_count"]),
                    "product_fixed_atom_count": int(built["product_fixed_atom_count"]),
                    "reactant_staged_indices": list(built["reactant_staged_indices"]),
                    "product_staged_indices": list(built["product_staged_indices"]),
                    "atom_mapping_constraints": list(built.get("atom_mapping_constraints", [])),
                    "reactant_min_distance": float(built["reactant_min_distance"]),
                    "product_min_distance": float(built["product_min_distance"]),
                    "product_source": product_source,
                    "product_structure_path": product_path,
                    "reactant_fixed_reference_signature": self.workflow._fixed_coordinate_signature(
                        reactant_neb,
                        int(built["reactant_fixed_atom_count"]),
                    ),
                    "product_fixed_reference_signature": self.workflow._fixed_coordinate_signature(
                        product_neb,
                        int(built["product_fixed_atom_count"]),
                    ),
                    "staging_hints": copied_hints,
                }
            )

            current = fixed_product.copy()

        return step_structures

    def postprocess_step_structures(
        self,
        step_structures: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Strict no-op postprocess: verify immutable signatures, do not alter coordinates."""
        processed: List[Dict[str, Any]] = []
        for entry in step_structures:
            reactant: Atoms = entry["reactant"]
            product: Atoms = entry["product"]

            r_fixed = int(entry.get("reactant_fixed_atom_count", len(reactant)))
            p_fixed = int(entry.get("product_fixed_atom_count", len(product)))

            r_ref = str(entry.get("reactant_fixed_reference_signature", "")).strip()
            p_ref = str(entry.get("product_fixed_reference_signature", "")).strip()

            r_now = self.workflow._fixed_coordinate_signature(reactant, r_fixed)
            p_now = self.workflow._fixed_coordinate_signature(product, p_fixed)

            if r_ref and r_ref != r_now:
                raise RuntimeError(
                    f"Immutable signature mismatch on reactant for step '{entry.get('name', 'unknown')}'"
                )
            if p_ref and p_ref != p_now:
                raise RuntimeError(
                    f"Immutable signature mismatch on product for step '{entry.get('name', 'unknown')}'"
                )

            processed_entry = dict(entry)
            processed_entry["agent4_postprocessed"] = True
            processed_entry["reactant_min_distance"] = self.workflow._min_pair_distance(
                reactant,
                list(entry.get("reactant_adsorbate_indices", [])),
            )
            processed_entry["product_min_distance"] = self.workflow._min_pair_distance(
                product,
                list(entry.get("product_adsorbate_indices", [])),
            )
            processed.append(processed_entry)
        return processed

    def _overlap_suggestion(
        self,
        reactant: Atoms,
        product: Atoms,
        step_name: str,
        threshold: float,
    ) -> str:
        """Generate actionable suggestion for the closest overlapping atom pair.

        Uses ``get_all_distances(mic=True)`` so the closest pair is found
        consistently with the caller's MIC-based overlap detection.
        """
        for label, atoms in [("reactant", reactant), ("product", product)]:
            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()
            n = len(atoms)
            if n < 2:
                continue
            dists = atoms.get_all_distances(mic=True)
            np.fill_diagonal(dists, np.inf)
            best_i, best_j = np.unravel_index(int(np.argmin(dists)), dists.shape)
            best_i, best_j = int(best_i), int(best_j)
            min_d = float(dists[best_i, best_j])
            if min_d < threshold:
                pi = positions[best_i]
                pj = positions[best_j]
                dz_needed = max(0.5, threshold - min_d + 0.3)
                higher = best_j if pj[2] >= pi[2] else best_i
                elem_h = symbols[higher]
                ph = positions[higher]
                return (
                    f"Move atom {elem_h} at [{ph[0]:.2f},{ph[1]:.2f},{ph[2]:.2f}] "
                    f"by +{dz_needed:.1f} Å in z to resolve overlap with atom "
                    f"{symbols[best_i if higher == best_j else best_j]} "
                    f"(current distance: {min_d:.2f} Å)"
                )
        return ""

    def programmatic_validation(self, step_structures: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Generic geometry gate for Agent5."""
        strict_interpolation_fatal = str(
            os.getenv("CATDT_STRICT_INTERPOLATION_FATAL", "0")
        ).strip().lower() in {"1", "true", "yes", "on"}
        if step_structures:
            first = step_structures[0]
            if not isinstance(first, dict) or not isinstance(first.get("reactant"), Atoms) or not isinstance(first.get("product"), Atoms):
                return {
                    "status": "FAIL",
                    "fatal_issues": ["Invalid step_structures payload: expected ASE Atoms under 'reactant' and 'product'"],
                    "warning_issues": [],
                    "summary": "(invalid payload)",
                    "feedback": "Provide canonical step_structures objects generated by tooling before calling validation.",
                }

        def _core_signature_order_invariant(atoms: Atoms, fixed_count: int, decimals: int = 8) -> tuple:
            n = max(0, min(int(fixed_count), len(atoms)))
            if n <= 0:
                return tuple()
            symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            rows = [
                (
                    str(symbols[idx]),
                    round(float(positions[idx][0]), decimals),
                    round(float(positions[idx][1]), decimals),
                    round(float(positions[idx][2]), decimals),
                )
                for idx in range(n)
            ]
            rows.sort()
            return tuple(rows)

        fatal_issues: List[str] = []
        warning_issues: List[str] = []
        summaries: List[str] = []

        previous_product: Atoms | None = None
        previous_product_fixed = 0

        for idx, step in enumerate(step_structures, start=1):
            step_name = str(step.get("name", f"step_{idx}"))
            reactant: Atoms = step["reactant"]
            product: Atoms = step["product"]

            reactant_ads = list(step.get("reactant_adsorbate_indices", []))
            product_ads = list(step.get("product_adsorbate_indices", []))
            reactant_fixed = int(step.get("reactant_fixed_atom_count", len(reactant)))
            product_fixed = int(step.get("product_fixed_atom_count", len(product)))

            reactant_ref_sig = str(step.get("reactant_fixed_reference_signature", "")).strip()
            product_ref_sig = str(step.get("product_fixed_reference_signature", "")).strip()
            reactant_sig = self.workflow._fixed_coordinate_signature(reactant, reactant_fixed)
            product_sig = self.workflow._fixed_coordinate_signature(product, product_fixed)
            if reactant_ref_sig and reactant_ref_sig != reactant_sig:
                fatal_issues.append(f"{step_name}: reactant immutable signature mismatch")
            if product_ref_sig and product_ref_sig != product_sig:
                fatal_issues.append(f"{step_name}: product immutable signature mismatch")

            if len(reactant) != len(product):
                fatal_issues.append(
                    f"{step_name}: endpoint atom counts differ ({len(reactant)} vs {len(product)})"
                )
            elif sorted(reactant.get_chemical_symbols()) != sorted(product.get_chemical_symbols()):
                fatal_issues.append(f"{step_name}: endpoint element multiset mismatch")

            reactant_label = self.workflow._canonical_species_label(str(step.get("reactant_formula", "")))
            product_label = self.workflow._canonical_species_label(str(step.get("product_formula", "")))
            add_count = int(step.get("atoms_to_add_count", 0) or 0)
            remove_count = int(step.get("atoms_to_remove_count", 0) or 0)
            raw_reactant_formula = str(step.get("reactant_formula", "")).strip()
            raw_product_formula = str(step.get("product_formula", "")).strip()
            phase_change_only = (
                self.workflow._contains_gas_phase_hint(raw_reactant_formula)
                or self.workflow._contains_gas_phase_hint(raw_product_formula)
            )
            if reactant_label != product_label and add_count == 0 and remove_count == 0 and not phase_change_only:
                fatal_issues.append(
                    f"{step_name}: labels differ ({reactant_label}->{product_label}) but no add/remove operations are provided"
                )

            endpoints_identical = (
                len(reactant) == len(product)
                and reactant.get_chemical_symbols() == product.get_chemical_symbols()
                and np.allclose(reactant.positions, product.positions, atol=1e-10)
            )
            if endpoints_identical and not phase_change_only:
                if raw_reactant_formula != raw_product_formula:
                    fatal_issues.append(
                        f"{step_name}: reactant/product coordinates are identical while endpoint formulas differ "
                        f"({raw_reactant_formula}->{raw_product_formula})"
                    )
                elif reactant_label != product_label:
                    fatal_issues.append(
                        f"{step_name}: reactant/product coordinates are identical while species labels differ"
                    )
            elif endpoints_identical and phase_change_only:
                # Gas-phase adsorption/desorption: identical coordinates are expected
                # at this stage; the actual gas-phase displacement is applied later
                # by _build_gas_phase_product_from_current.
                warning_issues.append(
                    f"{step_name}: adsorption/desorption step — gas-phase displacement applied at NEB stage"
                )
            elif endpoints_identical:
                warning_issues.append(f"{step_name}: reactant and product endpoints are numerically identical")

            # Compute each endpoint's distance matrix once and reuse.
            reactant_min = self.workflow._min_distance_any_pair(reactant)
            product_min = self.workflow._min_distance_any_pair(product)
            positive_mins = [v for v in (reactant_min, product_min) if v > 0]
            endpoint_min = min(positive_mins) if positive_mins else -1.0

            # H-X bonds are shorter (O-H ~0.97, C-H ~1.09, N-H ~1.01)
            # so use a lower threshold when H is involved in the closest pair
            h_involved = False
            if 0 < endpoint_min < 0.80:
                for struct in (reactant, product):
                    dists = struct.get_all_distances(mic=True)
                    np.fill_diagonal(dists, 999.0)
                    i_min, j_min = np.unravel_index(dists.argmin(), dists.shape)
                    sym = struct.get_chemical_symbols()
                    if dists[i_min, j_min] < 0.80 and (sym[i_min] == "H" or sym[j_min] == "H"):
                        h_involved = True
                        break
            fatal_threshold = 0.50 if h_involved else 0.65
            warning_threshold = 0.70 if h_involved else 0.80

            if 0 < endpoint_min < fatal_threshold:
                overlap_suggestion = self._overlap_suggestion(reactant, product, step_name, fatal_threshold)
                fatal_issues.append(
                    f"{step_name}: severe endpoint overlap (min pair {endpoint_min:.2f} Å). {overlap_suggestion}"
                )
            elif 0 < endpoint_min < warning_threshold:
                overlap_suggestion = self._overlap_suggestion(reactant, product, step_name, warning_threshold)
                warning_issues.append(
                    f"{step_name}: tight endpoint contacts (min pair {endpoint_min:.2f} Å). {overlap_suggestion}"
                )

            # Check staged atoms distance from adsorbate
            # Fatal if >5.0 Å (disconnected), warning if >4.0 Å
            # Staged atoms at 2-4 Å are expected (approaching/departing species)
            for endpoint_label, atoms_ep, ads_ep, staged_key in [
                ("reactant", reactant, reactant_ads, "reactant_staged_indices"),
                ("product", product, product_ads, "product_staged_indices"),
            ]:
                staged_ep = list(step.get(staged_key, []))
                if not staged_ep or not ads_ep:
                    continue
                core_ads = [i for i in ads_ep if i not in set(staged_ep) and 0 <= i < len(atoms_ep)]
                if not core_ads:
                    continue
                positions_ep = atoms_ep.get_positions()
                symbols_ep = atoms_ep.get_chemical_symbols()
                for si in staged_ep:
                    if si < 0 or si >= len(atoms_ep):
                        continue
                    sp = positions_ep[si]
                    idx_rel, nearest_pos, min_d = self._nearest_reference_position(
                        atoms_ep,
                        sp,
                        [positions_ep[idx] for idx in core_ads],
                    )
                    if idx_rel is None:
                        continue
                    nearest_idx = core_ads[int(idx_rel)]
                    nearest_sym = symbols_ep[nearest_idx]
                    if nearest_pos is None:
                        nearest_pos = positions_ep[nearest_idx]
                    suggest_pos = nearest_pos.copy()
                    suggest_pos[2] += 1.5
                    suggest_str = (
                        f"Move {symbols_ep[si]} closer to {nearest_sym} at "
                        f"[{nearest_pos[0]:.2f},{nearest_pos[1]:.2f},{nearest_pos[2]:.2f}]; "
                        f"suggested position: [{suggest_pos[0]:.2f},{suggest_pos[1]:.2f},{suggest_pos[2]:.2f}]"
                    )
                    if min_d > 5.0:
                        fatal_issues.append(
                            f"{step_name} {endpoint_label}: Staged {symbols_ep[si]} at "
                            f"[{sp[0]:.2f},{sp[1]:.2f},{sp[2]:.2f}] is {min_d:.1f} Å from nearest "
                            f"adsorbate atom — disconnected (max 5.0 Å). {suggest_str}"
                        )
                    elif min_d > 4.0:
                        warning_issues.append(
                            f"{step_name} {endpoint_label}: Staged {symbols_ep[si]} at "
                            f"[{sp[0]:.2f},{sp[1]:.2f},{sp[2]:.2f}] is {min_d:.1f} Å from nearest "
                            f"adsorbate — consider moving closer. {suggest_str}"
                        )

            for hint in list(step.get("staging_hints", [])):
                side = str(hint.get("side", "reactant"))
                fragment_elements = list(hint.get("fragment_elements", []))
                candidate_sites = list(hint.get("candidate_sites", []))
                if not fragment_elements or not candidate_sites:
                    continue

                endpoint_atoms = reactant if side == "reactant" else product
                staged_indices = list(step.get(f"{side}_staged_indices", []))
                if not staged_indices:
                    continue

                anchor_element = str(fragment_elements[0])
                matching_staged = [
                    idx_ep for idx_ep in staged_indices
                    if 0 <= idx_ep < len(endpoint_atoms) and endpoint_atoms[idx_ep].symbol == anchor_element
                ]
                if not matching_staged:
                    continue

                # The site rule constrains Agent4's DESIGN. If the energy gate has since
                # pre-relaxed the staged atoms, judge the positions Agent4 proposed.
                endpoint_positions = endpoint_atoms.get_positions()
                design_positions = step.get(f"{side}_design_positions")
                if design_positions is not None and np.shape(design_positions) == endpoint_positions.shape:
                    endpoint_positions = np.asarray(design_positions, float)
                best_idx = min(
                    matching_staged,
                    key=lambda idx_ep: self._nearest_candidate_site(endpoint_atoms, candidate_sites, endpoint_positions[idx_ep])[1],
                )
                nearest_site, site_dist = self._nearest_candidate_site(endpoint_atoms, candidate_sites, endpoint_positions[best_idx])
                nearest_rank = self._candidate_site_rank(candidate_sites, nearest_site)
                if nearest_site is None:
                    continue
                if site_dist > 1.5:
                    sp = endpoint_positions[best_idx]
                    site_pos = nearest_site["position"]
                    fatal_issues.append(
                        f"{step_name} {side}: staged {anchor_element} at "
                        f"[{sp[0]:.2f},{sp[1]:.2f},{sp[2]:.2f}] is not on a suggested co-adsorption site. "
                        f"Nearest candidate site {nearest_rank} is "
                        f"[{site_pos[0]:.2f},{site_pos[1]:.2f},{site_pos[2]:.2f}] ({site_dist:.2f} Å away)."
                    )
                else:
                    allowed_rank = self._get_allowed_site_rank_for_hint(step, hint)
                    if nearest_rank != allowed_rank:
                        sp = endpoint_positions[best_idx]
                        allowed_site = candidate_sites[allowed_rank - 1]
                        fatal_issues.append(
                            f"{step_name} {side}: staged {anchor_element} at "
                            f"[{sp[0]:.2f},{sp[1]:.2f},{sp[2]:.2f}] is on candidate site {nearest_rank}, "
                            f"but this step must use site {allowed_rank} at "
                            f"[{allowed_site['position'][0]:.2f},{allowed_site['position'][1]:.2f},{allowed_site['position'][2]:.2f}] "
                            f"before moving to other sites."
                        )

            interp_mins = self.workflow._best_interpolation_mapping_metrics(
                reactant,
                product,
                reactant_ads,
                product_ads,
                n_images=7,
            )
            finite_interp = [float(v) for v in interp_mins if np.isfinite(v) and v > 0]
            if finite_interp:
                min_interp = min(finite_interp)
                if min_interp < 0.50:
                    message = (
                        f"{step_name}: interpolated path has severe overlap in the geometric precheck "
                        f"(min {min_interp:.2f} Å < 0.50 Å). Rebuild staged endpoints before NEB; "
                        "fragment-aware alignment cannot make atom-crossing endpoints physically meaningful."
                    )
                    if strict_interpolation_fatal:
                        fatal_issues.append(message)
                    else:
                        warning_issues.append(message)
                elif min_interp < 0.80:
                    warning_issues.append(
                        f"{step_name}: interpolated path has close approaches (min {min_interp:.2f} Å). "
                        "Treat as advisory only; formal NEB frame generation will attempt to resolve this."
                    )

            if previous_product is not None:
                compare_n = min(previous_product_fixed, reactant_fixed, len(previous_product), len(reactant))
                if compare_n > 0:
                    previous_sig = _core_signature_order_invariant(previous_product, compare_n)
                    current_sig = _core_signature_order_invariant(reactant, compare_n)
                    if previous_sig != current_sig:
                        warning_issues.append(
                            f"{step_name}: previous step product core does not match current step reactant core"
                        )

            summaries.append(
                f"{step_name}: endpoint_min={endpoint_min:.3f}; "
                f"reactant_atoms={len(reactant)}; product_atoms={len(product)}"
            )

            previous_product = product.copy()
            previous_product_fixed = product_fixed

        status = "FAIL" if fatal_issues else "PASS"
        feedback_parts: List[str] = []
        if fatal_issues:
            feedback_parts.append("Fix all fatal geometry/mapping issues before NEB.")
        if warning_issues:
            feedback_parts.append("Advisory issues remain; improve staged atom placement.")

        fatal_text = "\n".join(fatal_issues)
        retry_site_lines: List[str] = []
        if fatal_text:
            for step in step_structures:
                step_name = str(step.get("name", ""))
                if not step_name or step_name not in fatal_text:
                    continue
                for hint in list(step.get("staging_hints", [])):
                    candidate_sites = list(hint.get("candidate_sites", []))
                    fragment_elements = list(hint.get("fragment_elements", []))
                    side = str(hint.get("side", "reactant"))
                    if len(candidate_sites) < 2 or not fragment_elements:
                        continue
                    endpoint_atoms = step.get(side)
                    if not isinstance(endpoint_atoms, Atoms):
                        continue
                    staged_indices = list(step.get(f"{side}_staged_indices", []))
                    anchor_element = str(fragment_elements[0])
                    matching_staged = [
                        idx_ep for idx_ep in staged_indices
                        if 0 <= idx_ep < len(endpoint_atoms) and endpoint_atoms[idx_ep].symbol == anchor_element
                    ]
                    if not matching_staged:
                        continue
                    endpoint_positions = endpoint_atoms.get_positions()
                    best_idx = min(
                        matching_staged,
                        key=lambda idx_ep: self._nearest_candidate_site(endpoint_atoms, candidate_sites, endpoint_positions[idx_ep])[1],
                    )
                    nearest_site, _ = self._nearest_candidate_site(endpoint_atoms, candidate_sites, endpoint_positions[best_idx])
                    if nearest_site is None:
                        continue
                    allowed_rank = self._get_allowed_site_rank_for_hint(step, hint)
                    current_rank = self._candidate_site_rank(candidate_sites, nearest_site)
                    if current_rank != allowed_rank:
                        # Agent4 keeps proposing a different candidate site, so the "try site N
                        # first" rank is never reached and can never be promoted: the design is
                        # rejected forever. Advance the rank after two such rejections.
                        retry_state = self._site_retry_state()
                        miss_key = self._site_retry_key(
                            step_name, side, str(hint.get("fragment_label", ""))
                        ) + ("rank_miss",)  # _site_retry_key returns a tuple
                        misses = int(retry_state.get(miss_key, 0)) + 1
                        retry_state[miss_key] = misses
                        if misses >= 2 and allowed_rank < len(candidate_sites):
                            retry_state[miss_key] = 0
                            self._promote_allowed_site_rank_for_hint(step, hint, allowed_rank + 1)
                        continue
                    next_rank = min(current_rank + 1, len(candidate_sites))
                    if next_rank <= current_rank:
                        continue
                    self._promote_allowed_site_rank_for_hint(step, hint, next_rank)
                    next_site = candidate_sites[next_rank - 1]
                    line = (
                        f"{step_name} {side}: after fixing any geometry errors, if this endpoint still fails "
                        f"Agent5 or the energy gate at site {current_rank}, retry site {next_rank} at "
                        f"[{next_site['position'][0]:.2f},{next_site['position'][1]:.2f},{next_site['position'][2]:.2f}]"
                    )
                    if len(candidate_sites) >= next_rank + 1:
                        site3 = candidate_sites[next_rank]
                        line += (
                            f"; if that also fails, use site {next_rank + 1} at "
                            f"[{site3['position'][0]:.2f},{site3['position'][1]:.2f},{site3['position'][2]:.2f}]"
                        )
                    if line not in retry_site_lines:
                        retry_site_lines.append(line)
        if retry_site_lines:
            feedback_parts.extend(retry_site_lines)

        return {
            "status": status,
            "fatal_issues": fatal_issues,
            "warning_issues": warning_issues,
            "summary": "\n".join(summaries) if summaries else "(none)",
            "feedback": " ".join(feedback_parts).strip(),
        }
