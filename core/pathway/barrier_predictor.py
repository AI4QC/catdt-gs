"""
Barrier Predictor - 反应能垒预测封装类

这个模块提供了一个完整的封装类，用于使用 Fairchem 的 CatTSunami 框架预测反应能垒。

CatTSunami 使用 Nudged Elastic Band (NEB) 方法来寻找反应过渡态和计算活化能。

使用方式:
    from barrier_predictor import BarrierPredictor
    from fairchem_predictor import FairchemPredictor

    # 初始化预测器（共享 Fairchem 模型）
    fairchem = FairchemPredictor(
        fairchem_root="/path/to/fairchem",
        model_name="uma-s-1p1",
    )

    barrier_predictor = BarrierPredictor(
        fairchem_predictor=fairchem,
    )

    # 方式1: 从反应式预测能垒（需要 CatTSunami 数据库）
    result = barrier_predictor.predict_from_reaction(
        reaction_str="*CH -> *C + *H",
        surface="Pt_111.vasp",
        n_frames=10,
    )
    print(f"Activation energy: {result.activation_energy:.3f} eV")

    # 方式2: 从已知的反应物和产物结构预测能垒
    result = barrier_predictor.predict_from_structures(
        reactant="reactant.vasp",
        product="product.vasp",
        n_frames=10,
    )
    print(f"Activation energy: {result.activation_energy:.3f} eV")

    # 方式3: 手动提供初始 NEB 路径
    result = barrier_predictor.run_neb(
        initial_frames=[atoms1, atoms2, ..., atoms_n],
        fmax=0.1,
    )

作者: Claude
"""

import os
import sys
import tempfile
import logging
from typing import Union, Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
import warnings

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.mep import DyNEB
from ase.optimize import BFGS, FIRE

# Import LLM-controlled NEB structure preparation
from .llm_neb_controller import LLMNEBController, NEBStructurePlan
from .neb_frame_prep import prepare_endpoint_pair_for_interpolation


@dataclass
class NEBResult:
    """NEB 计算结果"""
    reaction_name: str
    reactant: Atoms
    product: Atoms
    neb_frames: List[Atoms]
    n_frames: int
    energies: List[float]  # eV
    activation_energy_forward: float  # eV (reactant -> product)
    activation_energy_reverse: float  # eV (product -> reactant)
    reaction_energy: float  # eV (E_product - E_reactant)
    transition_state_index: int
    transition_state: Atoms
    converged: bool
    fmax_final: float
    optimization_steps: int
    trajectory_file: Optional[str] = None

    def __repr__(self):
        return (f"NEBResult(reaction={self.reaction_name}, "
                f"E_act_forward={self.activation_energy_forward:.3f} eV, "
                f"E_act_reverse={self.activation_energy_reverse:.3f} eV, "
                f"converged={self.converged})")

    def summary(self) -> str:
        """生成结果摘要"""
        lines = [
            "NEB Barrier Prediction Results",
            "=" * 60,
            f"Reaction: {self.reaction_name}",
            f"Number of frames: {self.n_frames}",
            f"Converged: {self.converged}",
            "",
            "Energetics:",
            f"  Reactant energy: {self.energies[0]:.4f} eV",
            f"  Product energy: {self.energies[-1]:.4f} eV",
            f"  Transition state energy: {self.energies[self.transition_state_index]:.4f} eV",
            "",
            f"  Reaction energy (ΔE): {self.reaction_energy:.4f} eV",
            f"  Forward barrier (E_act): {self.activation_energy_forward:.4f} eV",
            f"  Reverse barrier: {self.activation_energy_reverse:.4f} eV",
            "",
            f"Transition state at frame {self.transition_state_index}/{self.n_frames}",
            f"Final fmax: {self.fmax_final:.4f} eV/Å",
            f"Optimization steps: {self.optimization_steps}",
        ]

        if self.trajectory_file:
            lines.append(f"Trajectory saved to: {self.trajectory_file}")

        return "\n".join(lines)


@dataclass
class ReactionBarriersResult:
    """多个反应能垒的结果"""
    surface_formula: str
    reactions: List[str]
    barriers: Dict[str, NEBResult]
    output_dir: str

    def summary(self) -> str:
        """生成结果摘要"""
        lines = [
            "Reaction Barriers Prediction Results",
            "=" * 60,
            f"Surface: {self.surface_formula}",
            f"Number of reactions: {len(self.reactions)}",
            "",
            "Activation energies:",
        ]

        # 按正向能垒排序
        sorted_reactions = sorted(
            self.barriers.items(),
            key=lambda x: x[1].activation_energy_forward
        )

        for reaction, result in sorted_reactions:
            lines.append(
                f"  {reaction:30s}: E_act = {result.activation_energy_forward:6.3f} eV "
                f"(ΔE = {result.reaction_energy:6.3f} eV)"
            )

        lines.append("")
        lines.append(f"Output directory: {self.output_dir}")

        return "\n".join(lines)


class BarrierPredictor:
    """
    反应能垒预测器

    这个类封装了 Fairchem CatTSunami 的功能，提供简洁的 API 用于预测反应能垒。

    Parameters
    ----------
    fairchem_predictor : FairchemPredictor
        Fairchem 预测器实例（用于能量计算）
    fairchem_root : str, optional
        Fairchem 项目根目录（如果不提供 fairchem_predictor）
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """

    def __init__(
        self,
        fairchem_predictor=None,
        fairchem_root: Optional[str] = None,
        work_dir: Optional[str] = None,
        keep_files: bool = True,
        verbose: bool = True,
        use_llm_controller: bool = False,  # New: enable LLM-controlled NEB preparation
        llm_model: str = None,  # New: LLM model for structure preparation
    ):
        # 使用共享的 Fairchem 预测器或创建新的
        if fairchem_predictor is not None:
            self.fairchem = fairchem_predictor
            self.fairchem_root = fairchem_predictor.fairchem_root
        elif fairchem_root is not None:
            # 延迟导入避免循环依赖
            from .fairchem_predictor import FairchemPredictor
            self.fairchem = FairchemPredictor(
                fairchem_root=fairchem_root,
                verbose=verbose
            )
            self.fairchem_root = fairchem_root
        else:
            raise ValueError("Must provide either fairchem_predictor or fairchem_root")

        self.keep_files = keep_files
        self.verbose = verbose

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="barrier_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 添加 Fairchem 路径
        fairchem_src = os.path.join(self.fairchem_root, "src")
        if fairchem_src not in sys.path:
            sys.path.insert(0, fairchem_src)

        # 初始化日志
        self.logger = logging.getLogger(__name__)
        if self.verbose:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.WARNING)

        # Initialize LLM controller for NEB structure preparation (after logger is set up)
        self.use_llm_controller = use_llm_controller
        if self.use_llm_controller:
            self.llm_controller = LLMNEBController(
                model=llm_model,
                max_retries=3,
                temperature=0.1,
            )
            self.logger.info(f"LLM controller enabled (model: {llm_model or 'default'})")
        else:
            self.llm_controller = None

        # LLM配置（用于智能原子调整）
        self._llm_client = None
        self._llm_model = None
        self._llm_enabled = True

    def _init_llm(self):
        """初始化LLM用于智能原子调整 - 直接使用OpenAI客户端"""
        if self._llm_client is not None:
            return

        try:
            import os
            from dotenv import load_dotenv
            load_dotenv()

            from openai import OpenAI

            api_key = os.getenv("OPENAI_API_KEY")
            base_url = os.getenv("OPENAI_BASE_URL")
            model = os.getenv("OPENAI_MODEL", "claude-opus-4-5-20251101")

            if not api_key:
                self.logger.warning("OPENAI_API_KEY not set, LLM atom adjustment disabled")
                self._llm_enabled = False
                return

            self._llm_client = OpenAI(api_key=api_key, base_url=base_url)
            self._llm_model = model
            self.logger.info(f"LLM initialized for intelligent atom adjustment (model: {model})")

        except ImportError:
            self.logger.warning("openai package not available, using rule-based atom adjustment")
            self._llm_enabled = False
        except Exception as e:
            self.logger.warning(f"Failed to initialize LLM: {e}")
            self._llm_enabled = False

    def _analyze_atom_difference(
        self,
        reactant: Atoms,
        product: Atoms
    ) -> Dict[str, Any]:
        """
        分析反应物和产物之间的原子差异

        Returns
        -------
        Dict with:
            - 'atoms_to_add': List of (element, count) to add to reactant
            - 'atoms_to_remove': List of (element, count) to remove from reactant
            - 'net_change': Dict of element -> count change
        """
        from collections import Counter

        reactant_symbols = Counter(reactant.get_chemical_symbols())
        product_symbols = Counter(product.get_chemical_symbols())

        net_change = {}
        atoms_to_add = []
        atoms_to_remove = []

        all_elements = set(reactant_symbols.keys()) | set(product_symbols.keys())

        for elem in all_elements:
            r_count = reactant_symbols.get(elem, 0)
            p_count = product_symbols.get(elem, 0)
            diff = p_count - r_count

            if diff != 0:
                net_change[elem] = diff
                if diff > 0:
                    atoms_to_add.append((elem, diff))
                else:
                    atoms_to_remove.append((elem, -diff))

        return {
            'atoms_to_add': atoms_to_add,
            'atoms_to_remove': atoms_to_remove,
            'net_change': net_change,
            'reactant_formula': reactant.get_chemical_formula(),
            'product_formula': product.get_chemical_formula(),
        }

    def _llm_determine_atom_positions(
        self,
        atoms: Atoms,
        adsorbate_indices: List[int],
        atoms_to_add: List[Tuple[str, int]],
        reaction_name: str,
        reference_atoms: Atoms = None,
        reference_adsorbate_indices: List[int] = None,
    ) -> List[Tuple[str, np.ndarray]]:
        """
        使用LLM决定要添加的原子的位置 - 提供完整结构信息

        Parameters
        ----------
        atoms : Atoms
            当前结构（反应物）
        adsorbate_indices : List[int]
            吸附分子原子的索引
        atoms_to_add : List[Tuple[str, int]]
            要添加的原子列表 [(element, count), ...]
        reaction_name : str
            反应名称
        reference_atoms : Atoms, optional
            参考结构（产物）
        reference_adsorbate_indices : List[int], optional
            参考结构的吸附分子原子索引

        Returns
        -------
        List[Tuple[str, np.ndarray]]
            要添加的原子列表 [(element, position), ...]
        """
        self._init_llm()

        # 生成POSCAR格式的结构信息
        def atoms_to_poscar_string(structure: Atoms, name: str = "Structure") -> str:
            """将ASE Atoms对象转换为POSCAR格式字符串"""
            lines = [name]
            lines.append("1.0")

            # 晶格向量
            cell = structure.cell
            for vec in cell:
                lines.append(f"  {vec[0]:15.10f}  {vec[1]:15.10f}  {vec[2]:15.10f}")

            # 元素和数量
            symbols = structure.get_chemical_symbols()
            unique_elements = []
            counts = []
            for s in symbols:
                if s not in unique_elements:
                    unique_elements.append(s)
                    counts.append(1)
                else:
                    counts[unique_elements.index(s)] += 1

            lines.append("  " + "  ".join(unique_elements))
            lines.append("  " + "  ".join(map(str, counts)))
            lines.append("Cartesian")

            # 原子坐标
            positions = structure.get_positions()
            for i, (pos, sym) in enumerate(zip(positions, symbols)):
                marker = " <-- ADSORBATE" if i in adsorbate_indices else ""
                lines.append(f"  {pos[0]:15.10f}  {pos[1]:15.10f}  {pos[2]:15.10f}  # {i}: {sym}{marker}")

            return "\n".join(lines)

        # 获取结构信息
        reactant_poscar = atoms_to_poscar_string(atoms, "REACTANT")

        product_poscar = ""
        product_adsorbate_info = ""
        if reference_atoms is not None and reference_adsorbate_indices is not None:
            # 为产物生成POSCAR时使用产物的吸附原子索引
            def atoms_to_poscar_product(structure: Atoms, ads_indices: List[int], name: str) -> str:
                lines = [name]
                lines.append("1.0")
                cell = structure.cell
                for vec in cell:
                    lines.append(f"  {vec[0]:15.10f}  {vec[1]:15.10f}  {vec[2]:15.10f}")
                symbols = structure.get_chemical_symbols()
                unique_elements = []
                counts = []
                for s in symbols:
                    if s not in unique_elements:
                        unique_elements.append(s)
                        counts.append(1)
                    else:
                        counts[unique_elements.index(s)] += 1
                lines.append("  " + "  ".join(unique_elements))
                lines.append("  " + "  ".join(map(str, counts)))
                lines.append("Cartesian")
                positions = structure.get_positions()
                for i, (pos, sym) in enumerate(zip(positions, symbols)):
                    marker = " <-- ADSORBATE" if i in ads_indices else ""
                    lines.append(f"  {pos[0]:15.10f}  {pos[1]:15.10f}  {pos[2]:15.10f}  # {i}: {sym}{marker}")
                return "\n".join(lines)

            product_poscar = atoms_to_poscar_product(reference_atoms, reference_adsorbate_indices, "PRODUCT")

            # 提取产物中吸附分子的详细信息
            ref_positions = reference_atoms.get_positions()
            ref_symbols = reference_atoms.get_chemical_symbols()
            product_adsorbate_info = "Product adsorbate atoms (target positions):\n"
            for idx in reference_adsorbate_indices:
                pos = ref_positions[idx]
                product_adsorbate_info += f"  Index {idx}: {ref_symbols[idx]} at ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})\n"

        # 计算表面信息
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        # 找到表面原子（非吸附分子的最高z坐标）
        surface_atom_z = []
        for i, pos in enumerate(positions):
            if i not in adsorbate_indices:
                surface_atom_z.append(pos[2])
        surface_z_max = max(surface_atom_z) if surface_atom_z else positions[:, 2].max()

        # 吸附分子信息
        adsorbate_info = "Reactant adsorbate atoms:\n"
        for idx in adsorbate_indices:
            pos = positions[idx]
            adsorbate_info += f"  Index {idx}: {symbols[idx]} at ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})\n"

        # 构建详细的LLM提示
        prompt = f"""You are a computational chemistry expert helping to set up NEB (Nudged Elastic Band) calculations for surface catalysis.

TASK: Determine the INITIAL position for atoms that need to be added to the REACTANT structure for NEB calculation.

REACTION: {reaction_name}

ATOMS TO ADD TO REACTANT: {atoms_to_add}

=== REACTANT STRUCTURE (POSCAR format) ===
Atoms marked with "<-- ADSORBATE" are the adsorbate molecule on the surface.
{reactant_poscar}

{adsorbate_info}

=== PRODUCT STRUCTURE (POSCAR format) ===
This shows where the atoms end up AFTER the reaction.
{product_poscar if product_poscar else "Not provided"}

{product_adsorbate_info}

=== SURFACE INFORMATION ===
- Surface top z-coordinate (highest non-adsorbate atom): {surface_z_max:.3f} Å
- Cell dimensions: a={atoms.cell[0,0]:.2f}, b={atoms.cell[1,1]:.2f}, c={atoms.cell[2,2]:.2f} Å

=== CRITICAL GUIDELINES ===

**IMPORTANT**: The added atom position determines the NEB initial path. Poor positioning can create artificial barriers.

1. **General placement for any added atom/group**:
   - Infer where each added atom should appear from the PRODUCT structure.
   - Keep reactant x,y close to the corresponding product location when approaching bond formation.
   - Keep reactant z only moderately offset (typically ~0.5-1.2 Å), not unrealistically far away.

2. **General handling for atoms/groups being removed**:
   - Use a "just detached" geometry near the detachment region.
   - Keep x,y near the local site and apply only moderate z offset.

3. **Interpolation-aware positioning**:
   - NEB linearly interpolates between endpoints.
   - Large displacements can create nonphysical high-energy frames.
   - Prefer minimal, chemically plausible endpoint differences.

4. **Site continuity**:
   - Preserve local adsorption site continuity unless migration/desorption is explicitly intended.
   - Avoid random relocation of the whole adsorbate.

Return a JSON array with positions:
[
  {{"element": "<element_symbol>", "position": [x, y, z], "reasoning": "brief explanation"}}
]

Return ONLY the JSON array, no other text."""

        atoms_with_positions = []

        if self._llm_enabled and self._llm_client is not None:
            try:
                # 调用LLM
                response = self._llm_client.chat.completions.create(
                    model=self._llm_model,
                    messages=[{"role": "user", "content": prompt}],
                    temperature=0.1,
                )
                response_text = response.choices[0].message.content

                # 解析JSON
                import json
                import re

                # 尝试提取JSON（处理可能的```json代码块）
                cleaned_text = response_text.strip()
                if cleaned_text.startswith("```"):
                    cleaned_text = re.sub(r'^```(?:json)?\s*', '', cleaned_text)
                    cleaned_text = re.sub(r'\s*```$', '', cleaned_text)

                json_match = re.search(r'\[.*\]', cleaned_text, re.DOTALL)
                if json_match:
                    result = json.loads(json_match.group())
                    for item in result:
                        elem = item['element']
                        pos = np.array(item['position'])

                        # 验证位置合理性
                        # 使用统一且保守的最小距离阈值，避免对特定元素写死规则
                        min_dist_threshold = 0.9
                        is_valid, validation_msg = self._validate_atom_position(
                            pos, atoms, surface_z_max, min_distance=min_dist_threshold
                        )

                        if is_valid:
                            atoms_with_positions.append((elem, pos))
                            self.logger.info(f"LLM positioned {elem} at {pos}: {item.get('reasoning', '')}")
                        else:
                            self.logger.warning(f"LLM position for {elem} invalid: {validation_msg}")
                            # 尝试修正位置
                            corrected_pos = self._correct_atom_position(
                                pos, elem, atoms, adsorbate_indices, reference_atoms, reference_adsorbate_indices
                            )
                            atoms_with_positions.append((elem, corrected_pos))
                            self.logger.info(f"Corrected {elem} position to {corrected_pos}")

                    if atoms_with_positions:
                        return atoms_with_positions

            except Exception as e:
                self.logger.warning(f"LLM atom positioning failed: {e}")

        # 回退：使用参考结构中的位置
        self.logger.info("Using reference structure positions for atom placement")
        return self._fallback_atom_positions(
            atoms, adsorbate_indices, atoms_to_add,
            reference_atoms, reference_adsorbate_indices
        )

    def _validate_atom_position(
        self,
        pos: np.ndarray,
        atoms: Atoms,
        surface_z_max: float,
        min_distance: float = 1.0,
    ) -> Tuple[bool, str]:
        """
        验证原子位置是否合理

        Returns
        -------
        Tuple[bool, str]
            (是否有效, 原因说明)
        """
        # 检查是否在表面上方
        if pos[2] < surface_z_max - 0.5:
            return False, f"z-coordinate {pos[2]:.2f} is below surface ({surface_z_max:.2f})"

        # 检查与现有原子的距离
        positions = atoms.get_positions()
        distances = np.linalg.norm(positions - pos, axis=1)
        min_dist = distances.min()

        if min_dist < min_distance:
            return False, f"Too close to existing atom (min distance: {min_dist:.2f} Å)"

        # 检查是否在合理范围内（不要太远）
        if min_dist > 10.0:
            return False, f"Too far from structure (min distance: {min_dist:.2f} Å)"

        return True, "Position is valid"

    def _correct_atom_position(
        self,
        pos: np.ndarray,
        element: str,
        atoms: Atoms,
        adsorbate_indices: List[int],
        reference_atoms: Atoms = None,
        reference_adsorbate_indices: List[int] = None,
    ) -> np.ndarray:
        """
        修正不合理的原子位置

        策略：
        1. 如果有参考结构，使用参考位置但偏移2-3 Å
        2. 否则在吸附物质心上方添加
        """
        positions = atoms.get_positions()

        # 尝试从参考结构获取目标位置
        target_pos = None
        if reference_atoms is not None and reference_adsorbate_indices is not None:
            ref_symbols = reference_atoms.get_chemical_symbols()
            ref_positions = reference_atoms.get_positions()

            for idx in reference_adsorbate_indices:
                if ref_symbols[idx] == element:
                    target_pos = ref_positions[idx]
                    break

        if target_pos is not None:
            # 从目标位置向上偏移0.5-0.8 Å（代表预反应状态）
            # 关键：偏移量要小，0.5-0.8 Å 足够让NEB找到真正的过渡态
            # 太大的偏移（如2.5 Å）会导致NEB找不到正确的反应路径
            corrected = target_pos.copy()
            # H原子用更小的偏移，因为H-C/H-O键长约1.1 Å
            if element == 'H':
                corrected[2] += 0.6  # 小偏移，接近产物位置
            else:
                corrected[2] += 0.8
            return corrected

        # 回退：在吸附物质心上方
        if adsorbate_indices:
            ads_center = positions[adsorbate_indices].mean(axis=0)
            corrected = ads_center.copy()
            # 使用较小的偏移
            corrected[2] += 1.0 if element != 'H' else 0.8
            return corrected

        # 最后回退：在原位置上方
        corrected = pos.copy()
        corrected[2] = positions[:, 2].max() + 1.0
        return corrected

    def _fallback_atom_positions(
        self,
        atoms: Atoms,
        adsorbate_indices: List[int],
        atoms_to_add: List[Tuple[str, int]],
        reference_atoms: Atoms = None,
        reference_adsorbate_indices: List[int] = None,
    ) -> List[Tuple[str, np.ndarray]]:
        """
        回退方法：从参考结构中获取原子位置

        如果没有参考结构，则在吸附分子质心附近添加
        """
        atoms_with_positions = []
        positions = atoms.get_positions()

        # 计算吸附分子质心
        if adsorbate_indices:
            ads_center = positions[adsorbate_indices].mean(axis=0)
        else:
            # 没有吸附分子索引时，使用表面顶部中心
            z_max = positions[:, 2].max()
            cell = atoms.cell
            ads_center = np.array([cell[0, 0] / 2, cell[1, 1] / 2, z_max + 1.5])

        for elem, count in atoms_to_add:
            for i in range(count):
                new_pos = None

                # 尝试从参考结构获取位置
                if reference_atoms is not None:
                    ref_symbols = reference_atoms.get_chemical_symbols()
                    ref_positions = reference_atoms.get_positions()

                    # 找到参考结构中该元素在吸附分子区域的位置
                    if reference_adsorbate_indices:
                        for idx in reference_adsorbate_indices:
                            if ref_symbols[idx] == elem:
                                candidate_pos = ref_positions[idx]
                                # 检查是否已经使用
                                already_used = any(
                                    np.linalg.norm(candidate_pos - p) < 0.3
                                    for _, p in atoms_with_positions
                                )
                                if not already_used:
                                    new_pos = candidate_pos
                                    break

                # 如果没找到，使用更智能的定位策略
                if new_pos is None:
                    # 策略1：如果有参考结构，找到参考结构中该元素的位置并添加小偏移
                    if reference_atoms is not None:
                        ref_symbols = reference_atoms.get_chemical_symbols()
                        ref_positions = reference_atoms.get_positions()

                        # 遍历所有参考原子，找到该元素类型
                        for ref_idx, (ref_sym, ref_pos) in enumerate(zip(ref_symbols, ref_positions)):
                            if ref_sym == elem:
                                candidate_pos = ref_pos.copy()
                                # 添加小的z偏移（0.5-0.8 Å），使其略高于目标位置
                                # 这样NEB可以找到真正的过渡态
                                candidate_pos[2] += 0.6 if elem == 'H' else 0.8

                                already_used = any(
                                    np.linalg.norm(candidate_pos - p) < 0.3
                                    for _, p in atoms_with_positions
                                )
                                if not already_used:
                                    new_pos = candidate_pos
                                    self.logger.info(f"Using reference position for {elem}: {new_pos}")
                                    break

                    # 策略2：如果还是没找到，在吸附分子上方的合理位置添加
                    if new_pos is None:
                        # 对于H原子，应该靠近C或O原子（典型键长1.0-1.1 Å）
                        if elem == 'H' and adsorbate_indices:
                            symbols = atoms.get_chemical_symbols()
                            # 找到C或O原子
                            for idx in adsorbate_indices:
                                if symbols[idx] in ['C', 'O']:
                                    target_atom_pos = positions[idx]
                                    # 在该原子上方1.0 Å处放置H
                                    new_pos = target_atom_pos.copy()
                                    new_pos[2] += 1.0
                                    self.logger.info(f"Placing H above {symbols[idx]} at {new_pos}")
                                    break

                        # 如果仍然没找到，在质心上方
                        if new_pos is None:
                            new_pos = ads_center.copy()
                            new_pos[2] += 1.0
                            self.logger.info(f"Fallback: placing {elem} above adsorbate center at {new_pos}")

                atoms_with_positions.append((elem, new_pos))
                self.logger.info(f"Positioned {elem} at {new_pos}")

        return atoms_with_positions

    def _adjust_structures_for_neb(
        self,
        reactant: Atoms,
        product: Atoms,
        reaction_name: str = "reaction",
        reactant_adsorbate_indices: List[int] = None,
        product_adsorbate_indices: List[int] = None,
    ) -> Tuple[Atoms, Atoms, str]:
        """
        使用LLM智能调整反应物和产物结构，使原子数匹配以进行NEB计算

        完全通用的方法，适用于任何涉及原子数变化的反应

        Parameters
        ----------
        reactant : Atoms
            反应物结构
        product : Atoms
            产物结构
        reaction_name : str
            反应名称
        reactant_adsorbate_indices : List[int], optional
            反应物中吸附分子的原子索引（从工作流传递）
        product_adsorbate_indices : List[int], optional
            产物中吸附分子的原子索引（从工作流传递）

        Returns
        -------
        Tuple[Atoms, Atoms, str]
            (调整后的反应物, 调整后的产物, 调整说明)
        """
        diff = self._analyze_atom_difference(reactant, product)

        if not diff['net_change']:
            return reactant.copy(), product.copy(), "No adjustment needed"

        self.logger.info(f"Atom count adjustment needed for NEB:")
        self.logger.info(f"  Reactant: {diff['reactant_formula']}")
        self.logger.info(f"  Product: {diff['product_formula']}")
        self.logger.info(f"  Net change: {diff['net_change']}")

        adjusted_reactant = reactant.copy()
        adjusted_product = product.copy()
        adjustments = []

        # 如果产物原子更多（如加氢反应），在反应物中添加原子
        if diff['atoms_to_add']:
            self.logger.info(f"  Need to add atoms to reactant: {diff['atoms_to_add']}")

            # 使用LLM确定原子位置
            atoms_to_place = self._llm_determine_atom_positions(
                atoms=adjusted_reactant,
                adsorbate_indices=reactant_adsorbate_indices or [],
                atoms_to_add=diff['atoms_to_add'],
                reaction_name=reaction_name,
                reference_atoms=product,
                reference_adsorbate_indices=product_adsorbate_indices,
            )

            # 添加原子到结构中
            from ase import Atom
            for elem, pos in atoms_to_place:
                new_atom = Atom(elem, position=pos)
                adjusted_reactant.append(new_atom)
                adjustments.append(f"+{elem} to reactant")

            # 为新添加的原子设置tag（吸附物原子应该是tag=2）
            # 获取或创建tags数组
            n_atoms = len(adjusted_reactant)
            if reactant.has('tags'):
                # 扩展原有tags
                old_tags = list(reactant.get_tags())
                new_tags = old_tags + [2] * (n_atoms - len(reactant))  # 新原子tag=2
            else:
                # 创建新tags数组
                new_tags = [0] * len(reactant) + [2] * (n_atoms - len(reactant))
            adjusted_reactant.set_tags(new_tags)

        # 如果反应物原子更多（如脱附反应），在产物中添加原子
        if diff['atoms_to_remove']:
            self.logger.info(f"  Need to add atoms to product: {diff['atoms_to_remove']}")

            # 使用LLM确定原子位置
            atoms_to_place = self._llm_determine_atom_positions(
                atoms=adjusted_product,
                adsorbate_indices=product_adsorbate_indices or [],
                atoms_to_add=diff['atoms_to_remove'],
                reaction_name=reaction_name,
                reference_atoms=reactant,
                reference_adsorbate_indices=reactant_adsorbate_indices,
            )

            # 添加原子到结构中
            from ase import Atom
            for elem, pos in atoms_to_place:
                new_atom = Atom(elem, position=pos)
                adjusted_product.append(new_atom)
                adjustments.append(f"+{elem} to product")

            # 为新添加的原子设置tag（吸附物原子应该是tag=2）
            # 获取或创建tags数组
            n_atoms = len(adjusted_product)
            if product.has('tags'):
                # 扩展原有tags
                old_tags = list(product.get_tags())
                new_tags = old_tags + [2] * (n_atoms - len(product))  # 新原子tag=2
            else:
                # 创建新tags数组
                new_tags = [0] * len(product) + [2] * (n_atoms - len(product))
            adjusted_product.set_tags(new_tags)

        adjustment_desc = ", ".join(adjustments) if adjustments else "No adjustment"

        # 验证调整后原子数匹配
        if len(adjusted_reactant) != len(adjusted_product):
            self.logger.warning(f"Adjustment failed: {len(adjusted_reactant)} vs {len(adjusted_product)} atoms")
            return reactant.copy(), product.copy(), "Adjustment failed"

        # 对齐原子顺序，确保NEB可以正确运行
        # NEB要求两个结构的原子顺序完全一致
        adjusted_reactant, adjusted_product = self._align_atoms_for_neb(
            adjusted_reactant, adjusted_product
        )

        self.logger.info(f"  Adjusted: {len(adjusted_reactant)} atoms each")
        return adjusted_reactant, adjusted_product, adjustment_desc

    def _align_atoms_for_neb(
        self,
        reactant: Atoms,
        product: Atoms,
    ) -> Tuple[Atoms, Atoms]:
        """
        对齐两个结构的原子顺序，使其适合NEB计算

        通过元素类型和位置匹配原子，确保两个结构的原子索引对应
        """
        from scipy.optimize import linear_sum_assignment

        n_atoms = len(reactant)
        assert n_atoms == len(product), "Structures must have same number of atoms"

        r_symbols = reactant.get_chemical_symbols()
        p_symbols = product.get_chemical_symbols()
        r_pos = reactant.get_positions()
        p_pos = product.get_positions()

        # 构建成本矩阵：只有相同元素的原子才能匹配
        cost_matrix = np.full((n_atoms, n_atoms), 1e10)

        for i in range(n_atoms):
            for j in range(n_atoms):
                if r_symbols[i] == p_symbols[j]:
                    # 成本 = 位置距离
                    dist = np.linalg.norm(r_pos[i] - p_pos[j])
                    cost_matrix[i, j] = dist

        # 使用匈牙利算法找到最优匹配
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # 重排产物结构以匹配反应物顺序
        new_product = product[col_ind]
        new_product.cell = product.cell
        new_product.pbc = product.pbc

        # 重新计算tags（重要：fairchem需要正确的tags）
        # 不能简单重排tags，因为匹配是基于位置而非原子类型
        # 需要根据原子类型和z坐标重新确定tags
        if reactant.has('tags'):
            # 从reactant获取表面元素信息
            reactant_tags = reactant.get_tags()
            reactant_symbols = reactant.get_chemical_symbols()

            # 识别表面元素（tag=0或1的元素）
            surface_elements = set()
            for i, tag in enumerate(reactant_tags):
                if tag in [0, 1]:
                    surface_elements.add(reactant_symbols[i])

            # 为new_product设置tags
            new_tags = []
            new_symbols = new_product.get_chemical_symbols()
            new_positions = new_product.positions

            # 找到所有表面原子和吸附原子
            surface_atom_indices = []
            adsorbate_atom_indices = []
            for i, symbol in enumerate(new_symbols):
                if symbol in surface_elements:
                    surface_atom_indices.append(i)
                else:
                    adsorbate_atom_indices.append(i)

            # 对表面原子：根据z坐标分bulk和surface
            if surface_atom_indices:
                surface_z = new_positions[surface_atom_indices, 2]
                z_threshold = np.percentile(surface_z, 70)  # 顶部30%为surface

                for i in range(len(new_product)):
                    if i in adsorbate_atom_indices:
                        new_tags.append(2)  # adsorbate
                    elif i in surface_atom_indices:
                        if new_positions[i, 2] > z_threshold:
                            new_tags.append(1)  # surface layer
                        else:
                            new_tags.append(0)  # bulk layer
            else:
                # 全部是吸附原子（不太可能）
                new_tags = [2] * len(new_product)

            new_product.set_tags(new_tags)
            self.logger.info(f"  Recalculated tags for aligned product: {np.unique(new_tags)} (counts: {np.bincount(new_tags)})")

        # 保留constraints
        if product.constraints:
            new_product.set_constraint(product.constraints)

        self.logger.info(f"  Aligned atoms for NEB (max cost: {cost_matrix[row_ind, col_ind].max():.2f} Å)")

        return reactant, new_product

    def _prerelax_staged_atoms(
        self,
        reactant: Atoms,
        r_staged: List[int],
        product: Atoms,
        p_staged: List[int],
        fmax: float = 0.05,
        max_steps: int = 100,
        reactant_adsorbate_indices: Optional[List[int]] = None,
        product_adsorbate_indices: Optional[List[int]] = None,
    ) -> None:
        """Constrained relaxation of adsorbate + staged atoms to find local minima.

        All adsorbate atoms (including staged) are free; only slab atoms are
        fixed. This ensures both the adsorbate core and staged atoms reach a
        true energy minimum, which is critical for meaningful NEB barriers.

        Modifies reactant and product IN PLACE.
        """
        from ase.constraints import FixAtoms

        self.fairchem._load_model()
        self.logger.info("Pre-relaxing adsorbate + staged atoms (surface fixed)...")

        for atoms, staged, ads_indices, label in [
            (reactant, r_staged, reactant_adsorbate_indices, "reactant"),
            (product, p_staged, product_adsorbate_indices, "product"),
        ]:
            # Determine free atoms: all adsorbate indices (includes staged)
            free_set = set(ads_indices or []) | set(staged or [])
            if not free_set:
                continue

            # Use gentler relaxation for endpoints with staged atoms to avoid
            # the staged atom collapsing into the bonded position (which would
            # eliminate the reactant-product difference needed for NEB).
            has_staged = bool(staged)
            effective_fmax = max(fmax, 0.3) if has_staged else fmax
            effective_steps = min(max_steps, 40) if has_staged else max_steps

            fixed = [i for i in range(len(atoms)) if i not in free_set]
            atoms.calc = self.fairchem._calculator
            atoms.set_constraint(FixAtoms(indices=fixed))

            try:
                e_before = atoms.get_potential_energy()
            except Exception:
                e_before = None

            opt = FIRE(atoms, logfile=None)
            try:
                opt.run(fmax=effective_fmax, steps=effective_steps)
                n_steps = opt.get_number_of_steps()
                try:
                    e_after = atoms.get_potential_energy()
                except Exception:
                    e_after = None

                if e_before is not None and e_after is not None:
                    self.logger.info(
                        f"  Pre-relaxed {label} adsorbate ({len(free_set)} free atoms, "
                        f"{len(staged)} staged): "
                        f"{n_steps} steps, dE = {e_after - e_before:+.3f} eV "
                        f"({e_before:.3f} -> {e_after:.3f} eV)"
                    )
                else:
                    self.logger.info(
                        f"  Pre-relaxed {label} adsorbate: {n_steps} steps"
                    )
            except Exception as e:
                self.logger.warning(f"  Pre-relaxation failed for {label}: {e}")

            atoms.set_constraint()
            atoms.calc = None

    def _interpolate_with_overlap_correction(
        self,
        initial_frame: Atoms,
        final_frame: Atoms,
        num_frames: int,
        n_iterations: int = 100,
        rate: float = 0.1,
    ) -> list:
        """CatTSunami-style interpolation: linear + iterative overlap correction.

        Based on fairchem/applications/cattsunami/core/autoframe.py::interpolate().
        Uses linear interpolation as initial guess, then iteratively pushes
        apart atoms that are closer than their target distances.

        Args:
            initial_frame: reactant Atoms
            final_frame: product Atoms
            num_frames: total number of frames (including endpoints)
            n_iterations: number of overlap correction iterations
            rate: correction rate per iteration

        Returns:
            list of Atoms frames
        """
        import torch

        start_pos = torch.from_numpy(initial_frame.get_positions()).double()
        end_pos = torch.from_numpy(final_frame.get_positions()).double()
        num_atoms = len(start_pos)

        # Target distance matrices (linearly interpolated)
        start_dist = torch.from_numpy(
            initial_frame.get_all_distances(mic=True)
        ).double().reshape(-1)
        end_dist = torch.from_numpy(
            final_frame.get_all_distances(mic=True)
        ).double().reshape(-1)

        alpha = torch.linspace(0, 1, num_frames, dtype=torch.float64)

        # Linear interpolation: frame_i = initial * (1-alpha_i) + final * alpha_i
        frames = start_pos.unsqueeze(0) * (1.0 - alpha.view(-1, 1, 1)) + \
                 end_pos.unsqueeze(0) * alpha.view(-1, 1, 1)

        target_dist = start_dist.unsqueeze(0) * (1.0 - alpha.view(-1, 1)) + \
                      end_dist.unsqueeze(0) * alpha.view(-1, 1)

        # Build initial Atoms list
        def frames_to_atoms(frames_tensor):
            atoms_list = []
            for i in range(num_frames):
                a = initial_frame.copy()
                a.set_positions(frames_tensor[i].numpy(), apply_constraint=False)
                atoms_list.append(a)
            return atoms_list

        # Iterative overlap correction
        for _it in range(n_iterations):
            atoms_list = frames_to_atoms(frames)

            frame_dist = []
            frame_vec = []
            for a in atoms_list:
                d = a.get_all_distances(mic=True)
                v = a.get_all_distances(mic=True, vector=True)
                frame_dist.append(d.ravel())
                frame_vec.append(v)

            frame_dist = torch.from_numpy(np.array(frame_dist)).double()
            frame_vec = torch.from_numpy(np.array(frame_vec)).double()

            # Normalize direction vectors
            norms = torch.norm(frame_vec, dim=3, keepdim=True) + 1e-8
            frame_vec = frame_vec / norms

            # Only correct atoms that are too close (delta < 0)
            delta = torch.clamp(frame_dist - target_dist, max=0)
            weight = torch.exp(-(target_dist ** 2) / (0.5 * 2.0 * 2.0))
            weight = weight + torch.exp(-(frame_dist ** 2) / (0.5 * 2.0 * 2.0))
            delta = delta * weight

            delta_3d = rate * frame_vec * delta.view(num_frames, num_atoms, num_atoms, 1)
            delta_3d = torch.sum(delta_3d, dim=1)  # sum over atom pairs

            # Freeze endpoints
            delta_3d[0] = 0.0
            delta_3d[num_frames - 1] = 0.0

            frames = frames - delta_3d

            # Smooth adjacent frames to avoid large jumps
            if num_frames > 2:
                mean_pos = (frames[:-2] + frames[2:]) / 2.0
                frames[1:-1] = 0.9 * frames[1:-1] + 0.1 * mean_pos

        result = frames_to_atoms(frames)

        # Log min distances for each frame
        for i, a in enumerate(result):
            if i == 0 or i == num_frames - 1:
                continue
            dists = a.get_all_distances(mic=True)
            np.fill_diagonal(dists, 999.0)
            min_d = dists.min()
            if min_d < 0.8:
                self.logger.warning(
                    f"  Frame {i}: min interatomic distance = {min_d:.3f} Å (< 0.8 Å!)"
                )

        for atoms_now in result:
            try:
                atoms_now.wrap()
            except Exception:
                pass

        return result

    def _build_neb_image_calculators(self, n_images: int) -> tuple[list[object], bool]:
        """Create one calculator per NEB image when the backend supports it.

        ASE's NEB implementation assumes each image owns its own calculator
        instance. Sharing a single FAIRChemCalculator across all images means
        mutable calculator state (`self.atoms`, `self.results`) is reused across
        different endpoint/image objects, which is explicitly discouraged by ASE
        and has shown up as nondeterministic native crashes during GPU NEB runs.

        To keep GPU memory usage reasonable, we do not duplicate the underlying
        ML model. Instead we create lightweight calculator wrappers that all
        point at the same loaded predict unit.
        """
        self.fairchem._load_model()

        base_calc = self.fairchem._calculator
        if n_images <= 0:
            return [], False

        calc_cls = type(base_calc)
        shared_predictor = getattr(base_calc, "predictor", None)
        task_name = getattr(base_calc, "task_name", None)

        if shared_predictor is not None:
            calculators: list[object] = []
            try:
                for _ in range(int(n_images)):
                    calculators.append(calc_cls(shared_predictor, task_name=task_name))
            except Exception as exc:
                self.logger.warning(
                    "Failed to create isolated NEB calculators from shared predict unit; "
                    "falling back to shared calculator: %s",
                    exc,
                )
            else:
                return calculators, False

        return [base_calc for _ in range(int(n_images))], True

    def run_neb(
        self,
        initial_frames: List[Atoms],
        fmax: float = 0.1,
        max_steps: int = 300,
        delta_fmax_climb: float = 0.1,
        k: float = 1.0,
        climb: bool = True,
        use_two_phase: bool = True,
        trajectory_file: Optional[str] = None,
        reaction_name: str = "reaction",
    ) -> NEBResult:
        """
        运行 NEB 计算

        Uses DyNEB + FIRE with two-phase approach:
          Phase 1: Plain NEB (no climbing) to fmax + delta_fmax_climb (0.2 eV/Å)
          Phase 2: CI-NEB (climbing image) to fmax (0.1 eV/Å)

        Parameters
        ----------
        initial_frames : list of Atoms
            初始 NEB 路径（包括反应物和产物）
        fmax : float, default=0.1
            CI-NEB 力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        delta_fmax_climb: float, default=0.1
            Phase 1 的 fmax 增量，Phase 1 目标为 fmax + delta_fmax_climb
        k : float, default=1.0
            弹簧常数 (eV/Å²)
        climb : bool, default=True
            是否使用爬升像方法
        use_two_phase : bool, default=True
            是否使用两阶段优化（Phase 1无climbing + Phase 2有climbing）
            True: 先以 fmax+delta 稳定路径，再启用 climbing 优化到 fmax（推荐）
            False: 直接使用 CI-NEB 优化到 fmax
        trajectory_file : str, optional
            轨迹文件路径
        reaction_name : str, default="reaction"
            反应名称

        Returns
        -------
        NEBResult
            NEB 计算结果
        """
        self.logger.info(f"Running {'Two-Phase ' if use_two_phase else ''}CI-NEB for {reaction_name}")
        self.logger.info(f"Number of frames: {len(initial_frames)}")

        # 复制帧并为每个 image 分配计算器
        frames = [atoms.copy() for atoms in initial_frames]
        image_calculators, using_shared_calculator = self._build_neb_image_calculators(len(frames))
        for frame, calculator in zip(frames, image_calculators):
            frame.calc = calculator

        # 使用 DyNEB（动态 NEB）
        neb = DyNEB(
            frames,
            k=k,
            climb=climb if not use_two_phase else False,  # 如果不用两阶段，直接启用climbing
            parallel=False,           # 共享 GPU predictor 不支持并行（会数据竞争）
            dynamic_relaxation=True,  # 启用动态松弛
            scale_fmax=0.0,           # 禁用 fmax 缩放
            allow_shared_calculator=using_shared_calculator,
        )

        # 设置轨迹文件
        if trajectory_file is None:
            trajectory_file = os.path.join(
                self.work_dir,
                f"{reaction_name}_neb.traj"
            )

        # 优化器（FIRE 比 BFGS 对 NEB 更稳定，不易振荡）
        optimizer = FIRE(neb, trajectory=trajectory_file)

        def _run_and_check_convergence(target_fmax: float, step_limit: int, phase_name: str) -> bool:
            step_limit = max(0, int(step_limit))
            if step_limit <= 0:
                return False

            monitor_interval = 2
            min_steps_for_stagnation = min(500, max(100, step_limit // 3))
            stagnant_checks_limit = 50
            improvement_tol = max(0.01, 0.15 * float(target_fmax))
            severe_multiplier = 10.0

            monitor_state = {
                "best_fmax": float("inf"),
                "stagnant_checks": 0,
                "early_stopped": False,
                "latest_fmax": float("inf"),
            }

            def _stagnation_monitor() -> None:
                try:
                    forces = np.array(neb.get_forces(), dtype=float)
                    if forces.size == 0:
                        return
                    current_fmax = float(np.max(np.linalg.norm(forces, axis=1)))
                except Exception:
                    return

                monitor_state["latest_fmax"] = current_fmax
                best_fmax = float(monitor_state["best_fmax"])
                if current_fmax < (best_fmax - improvement_tol):
                    monitor_state["best_fmax"] = current_fmax
                    monitor_state["stagnant_checks"] = 0
                else:
                    monitor_state["stagnant_checks"] = int(monitor_state["stagnant_checks"]) + 1

                current_steps = int(optimizer.get_number_of_steps())
                if (
                    current_steps >= min_steps_for_stagnation
                    and int(monitor_state["stagnant_checks"]) >= stagnant_checks_limit
                    and current_fmax > (severe_multiplier * float(target_fmax))
                ):
                    monitor_state["early_stopped"] = True
                    raise StopIteration(
                        f"{phase_name}: stagnated at fmax={current_fmax:.3f} eV/Å for {int(monitor_state['stagnant_checks'])} checks"
                    )

            optimizer.attach(_stagnation_monitor, interval=monitor_interval)
            try:
                run_result = optimizer.run(fmax=target_fmax, steps=step_limit)
            except Exception as exc:
                msg = str(exc)
                if isinstance(exc, StopIteration) or "StopIteration" in msg:
                    self.logger.warning("Early stop in %s due to stagnation: %s", phase_name, msg)
                    return False
                raise
            finally:
                # ASE has no detach() method — manually remove from observers list
                try:
                    optimizer.observers = [
                        obs for obs in optimizer.observers
                        if obs[0] is not _stagnation_monitor
                    ]
                except Exception:
                    pass

            if isinstance(run_result, bool):
                return run_result
            # Fallback for optimizer implementations that don't return bool
            return optimizer.get_number_of_steps() < max(1, int(step_limit))

        if use_two_phase:
            # 两阶段优化：先不用climbing，再用climbing
            # Phase 1 gets at most 60% of total budget; Phase 2 always gets at least 40%
            # This ensures the climbing image (critical for saddle point) always runs
            phase1_budget = int(max_steps * 0.6)
            min_phase2_budget = max_steps - phase1_budget  # at least 40%

            # 第一阶段：不使用爬升像
            self.logger.info(f"Phase 1: Optimizing without climbing image (budget: {phase1_budget} steps)...")
            try:
                converged_phase1 = _run_and_check_convergence(
                    target_fmax=fmax + delta_fmax_climb,
                    step_limit=phase1_budget,
                    phase_name="phase1",
                )
            except Exception as e:
                self.logger.warning(f"Phase 1 optimization failed: {e}")
                converged_phase1 = False

            optimization_steps = optimizer.get_number_of_steps()

            # 第二阶段：使用爬升像
            # Phase 2 ALWAYS runs — climbing image is essential for locating saddle points
            converged = converged_phase1
            remaining_steps = max(max_steps - optimization_steps, min_phase2_budget)
            if climb:
                if not converged_phase1:
                    self.logger.info(f"Phase 1 did not fully converge, proceeding to Phase 2 (CI-NEB needs climbing image)")
                self.logger.info(f"Phase 2: Optimizing with climbing image ({remaining_steps} steps remaining)...")
                neb.climb = True
                try:
                    converged_phase2 = _run_and_check_convergence(
                        target_fmax=fmax,
                        step_limit=remaining_steps,
                        phase_name="phase2",
                    )
                    converged = converged_phase2
                except Exception as e:
                    self.logger.warning(f"Phase 2 optimization failed: {e}")
                    converged = False

                optimization_steps = optimizer.get_number_of_steps()
        else:
            # 单阶段CI-NEB优化（推荐）
            self.logger.info(f"Optimizing with CI-NEB (climb={climb})...")
            try:
                converged = _run_and_check_convergence(
                    target_fmax=fmax,
                    step_limit=max_steps,
                    phase_name="single_phase",
                )
            except Exception as e:
                self.logger.warning(f"CI-NEB optimization failed: {e}")
                converged = False

            optimization_steps = optimizer.get_number_of_steps()

        self.logger.info(f"NEB optimizer convergence flag: {converged}")
        self.logger.info(f"NEB optimizer steps used: {optimization_steps}/{max_steps}")

        # 提取结果
        energies = []
        for frame in frames:
            try:
                e = frame.get_potential_energy()
                energies.append(e)
            except Exception as e:
                self.logger.error(f"Failed to get energy for frame: {e}")
                energies.append(np.nan)

        energies = np.array(energies)

        # 打印NEB路径能量用于诊断
        self.logger.info(f"NEB path energies: {energies}")
        self.logger.info(f"  Energy range: {energies.min():.3f} to {energies.max():.3f} eV")

        # 找到过渡态 — 仅在 interior images 中寻找（排除端点）
        # 与 ASE 的 NEBState.imax 一致: 1 + argmax(energies[1:-1])
        # 端点是固定的初末态，不应作为过渡态候选
        reactant_energy = energies[0]
        product_energy = energies[-1]

        if len(energies) > 2:
            ts_index = 1 + np.argmax(energies[1:-1])
        else:
            ts_index = np.argmax(energies)
        ts_energy = energies[ts_index]

        self.logger.info(f"  Reactant (frame 0) energy = {reactant_energy:.3f} eV")
        self.logger.info(f"  Product (frame {len(energies)-1}) energy = {product_energy:.3f} eV")
        self.logger.info(f"  Transition state at frame {ts_index} (interior), energy = {ts_energy:.3f} eV")

        # 计算能垒和反应能
        activation_energy_forward = max(0.0, ts_energy - reactant_energy)
        activation_energy_reverse = max(0.0, ts_energy - product_energy)
        reaction_energy = product_energy - reactant_energy

        self.logger.info(f"  Ea_fwd = {activation_energy_forward:.3f} eV, Ea_rev = {activation_energy_reverse:.3f} eV")

        # 警告：端点能量高于所有 interior images → 端点可能不是局部极小
        if reactant_energy > ts_energy:
            self.logger.warning(
                f"  WARNING: Reactant energy ({reactant_energy:.3f}) > TS energy ({ts_energy:.3f}). "
                f"Reactant may not be a local minimum — staged atoms may be in gas phase!"
            )
        if product_energy > ts_energy:
            self.logger.warning(
                f"  WARNING: Product energy ({product_energy:.3f}) > TS energy ({ts_energy:.3f}). "
                f"Product may not be a local minimum — staged atoms may be in gas phase!"
            )

        # 计算最终 fmax
        forces = frames[ts_index].get_forces()
        fmax_final = np.max(np.linalg.norm(forces, axis=1))

        self.logger.info(f"DyNEB completed. E_act = {activation_energy_forward:.3f} eV")

        # Release CUDA memory held by per-frame calculators.
        # NEBResult keeps the Atoms objects but we clear their calc references
        # (energies/forces are already extracted above).
        for frame in frames:
            if hasattr(frame, 'calc') and frame.calc is not None:
                frame.calc.results = {}
                frame.calc = None

        return NEBResult(
            reaction_name=reaction_name,
            reactant=frames[0],
            product=frames[-1],
            neb_frames=frames,
            n_frames=len(frames),
            energies=energies.tolist(),
            activation_energy_forward=activation_energy_forward,
            activation_energy_reverse=activation_energy_reverse,
            reaction_energy=reaction_energy,
            transition_state_index=ts_index,
            transition_state=frames[ts_index],
            converged=converged,
            fmax_final=fmax_final,
            optimization_steps=optimization_steps,
            trajectory_file=trajectory_file,
        )

    def predict_from_structures(
        self,
        reactant: Union[str, Atoms],
        product: Union[str, Atoms],
        n_frames: int = 10,
        fmax: float = 0.1,
        max_steps: int = 300,
        relax_endpoints: bool = False,
        reaction_name: Optional[str] = None,
        reactant_adsorbate_indices: List[int] = None,
        product_adsorbate_indices: List[int] = None,
        reactant_staged_indices: List[int] = None,
        product_staged_indices: List[int] = None,
        previous_step_product: Atoms = None,
        previous_step_product_adsorbate_indices: List[int] = None,
        step_index: int = None,
        use_two_phase: bool = True,
        delta_fmax_climb: float = 0.1,
        k: float = 1.0,
    ) -> NEBResult:
        """
        从反应物和产物结构预测能垒

        Parameters
        ----------
        reactant : str or Atoms
            反应物结构
        product : str or Atoms
            产物结构
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            CI-NEB 力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        relax_endpoints : bool, default=True
            是否先优化反应物和产物
        reaction_name : str, optional
            反应名称
        reactant_adsorbate_indices : List[int], optional
            反应物中吸附分子的原子索引（从工作流传递，用于LLM智能原子调整）
        product_adsorbate_indices : List[int], optional
            产物中吸附分子的原子索引（从工作流传递，用于LLM智能原子调整）
        reactant_staged_indices : List[int], optional
            反应物中暂驻原子的索引（仅这些原子在端点松弛时可移动）
        product_staged_indices : List[int], optional
            产物中暂驻原子的索引（仅这些原子在端点松弛时可移动）
        previous_step_product : Atoms, optional
            前一步反应的产物结构（用于确保路径连续性）
        previous_step_product_adsorbate_indices : List[int], optional
            前一步产物的吸附质原子索引
        step_index : int, optional
            当前步骤在反应路径中的索引
        use_two_phase : bool, default=True
            是否使用两阶段优化（Phase 1 fmax=0.2 → Phase 2 CI-NEB fmax=0.1）
        delta_fmax_climb : float, default=0.1
            Phase 1 的 fmax 增量
        k : float, default=1.0
            弹簧常数 (eV/Å²)

        Returns
        -------
        NEBResult
            NEB 计算结果
        """
        # 读取结构
        if isinstance(reactant, str):
            reactant_atoms = read(reactant)
            if reaction_name is None:
                reaction_name = os.path.basename(reactant).replace(".vasp", "")
        else:
            reactant_atoms = reactant.copy()
            if reaction_name is None:
                reaction_name = "reaction"

        if isinstance(product, str):
            product_atoms = read(product)
        else:
            product_atoms = product.copy()

        # 保存原始结构用于能量校正
        reactant_atoms_original = reactant_atoms.copy()
        product_atoms_original = product_atoms.copy()
        adjusted_reactant = None
        adjusted_product = None
        reactant_result = None
        product_result = None

        # 优化端点
        if relax_endpoints:
            r_staged = set(reactant_staged_indices or [])
            p_staged = set(product_staged_indices or [])
            has_staged = bool(r_staged or p_staged)

            if has_staged:
                # Post-staging relaxation is handled UPSTREAM by
                # workflow._relax_neb_endpoints (which frees staged atoms +
                # top-2 surface layers). By the time we get here the endpoints
                # are already at local energy minima that include the staging
                # atoms. Doing a second relaxation here with different
                # constraints would either be redundant or counterproductive.
                self.logger.info(
                    "Staged atoms present (reactant_staged=%d, product_staged=%d). "
                    "Skipping internal relax — endpoints already pre-relaxed by workflow.",
                    len(r_staged), len(p_staged),
                )
            else:
                # Full unconstrained relaxation (original behavior)
                self.logger.info("Relaxing reactant (unconstrained)...")
                reactant_result = self.fairchem.predict_energy(
                    reactant_atoms,
                    relax=True,
                    fmax=fmax,
                    max_steps=max_steps,
                )
                reactant_atoms = reactant_result.final_structure

                self.logger.info("Relaxing product (unconstrained)...")
                product_result = self.fairchem.predict_energy(
                    product_atoms,
                    relax=True,
                    fmax=fmax,
                    max_steps=max_steps,
                )
                product_atoms = product_result.final_structure

        # 检查原子数是否匹配
        if len(reactant_atoms) != len(product_atoms):
            self.logger.info(
                f"Reactant ({len(reactant_atoms)} atoms) and product ({len(product_atoms)} atoms) "
                f"have different number of atoms. Attempting intelligent atom adjustment..."
            )

            # Use LLM controller if enabled, otherwise use legacy method
            if self.use_llm_controller and self.llm_controller is not None:
                self.logger.info("Using LLM-controlled structure preparation...")

                try:
                    # Use LLM controller to prepare structures
                    adjusted_reactant, adjusted_product, new_r_idx, new_p_idx, plan = \
                        self.llm_controller.prepare_neb_structures(
                            reactant=reactant_atoms,
                            product=product_atoms,
                            reactant_adsorbate_indices=reactant_adsorbate_indices or [],
                            product_adsorbate_indices=product_adsorbate_indices or [],
                            reaction_name=reaction_name,
                            previous_step_product=previous_step_product,
                            previous_step_product_adsorbate_indices=previous_step_product_adsorbate_indices,
                            step_index=step_index,
                        )

                    # 保存LLM调整后的初末态结构（用于可视化）
                    from ase.io import write as ase_write
                    try:
                        adjusted_reactant_file = os.path.join(
                            self.work_dir,
                            f"{reaction_name}_neb_initial.vasp"
                        )
                        adjusted_product_file = os.path.join(
                            self.work_dir,
                            f"{reaction_name}_neb_final.vasp"
                        )

                        ase_write(adjusted_reactant_file, adjusted_reactant)
                        ase_write(adjusted_product_file, adjusted_product)

                        self.logger.info(f"Saved NEB initial state: {adjusted_reactant_file}")
                        self.logger.info(f"Saved NEB final state: {adjusted_product_file}")
                    except Exception as e:
                        self.logger.warning(f"Failed to save adjusted structures: {e}")

                    # ⭐ FIXED: apply_plan already handles atom reordering
                    # No need to call align_structures_for_neb which uses Hungarian algorithm
                    # that would undo the element-type based reordering from apply_plan
                    self.logger.info("Structures already aligned by LLM controller apply_plan()")

                    adjustment_desc = f"LLM-controlled: {plan.reaction_type}"
                    self.logger.info(f"LLM plan: {plan.reasoning}")

                except Exception as e:
                    self.logger.error(f"LLM-controlled preparation failed: {e}")
                    self.logger.info("Falling back to legacy method...")
                    # Fallback to legacy method
                    adjusted_reactant, adjusted_product, adjustment_desc = self._adjust_structures_for_neb(
                        reactant_atoms, product_atoms, reaction_name,
                        reactant_adsorbate_indices=reactant_adsorbate_indices,
                        product_adsorbate_indices=product_adsorbate_indices,
                    )
            else:
                # Legacy method: rule-based atom adjustment
                self.logger.info("Using legacy structure adjustment method...")
                adjusted_reactant, adjusted_product, adjustment_desc = self._adjust_structures_for_neb(
                    reactant_atoms, product_atoms, reaction_name,
                    reactant_adsorbate_indices=reactant_adsorbate_indices,
                    product_adsorbate_indices=product_adsorbate_indices,
                )

            if len(adjusted_reactant) == len(adjusted_product):
                self.logger.info(f"Atom adjustment successful: {adjustment_desc}")

                # 对调整后的结构进行优化
                self.logger.info("Relaxing adjusted reactant...")
                # 调试：打印tags
                if adjusted_reactant.has('tags'):
                    self.logger.info(f"  Adjusted reactant tags: {adjusted_reactant.get_tags()}")
                else:
                    self.logger.warning("  Adjusted reactant has NO tags!")

                adj_reactant_result = self.fairchem.predict_energy(
                    adjusted_reactant,
                    relax=True,
                    fmax=fmax,
                    max_steps=max_steps,
                )
                reactant_atoms = adj_reactant_result.final_structure

                self.logger.info("Relaxing adjusted product...")
                # 调试：打印tags
                if adjusted_product.has('tags'):
                    self.logger.info(f"  Adjusted product tags: {adjusted_product.get_tags()}")
                else:
                    self.logger.warning("  Adjusted product has NO tags!")

                adj_product_result = self.fairchem.predict_energy(
                    adjusted_product,
                    relax=True,
                    fmax=fmax,
                    max_steps=max_steps,
                )
                product_atoms = adj_product_result.final_structure

                self.logger.info(f"Adjusted structures optimized. Proceeding with NEB calculation...")
            else:
                # 调整失败，回退到仅计算反应能的模式
                self.logger.warning(
                    f"Atom adjustment failed: {adjustment_desc}. "
                    f"Returning reaction energy only (no barrier calculation)."
                )

                # 计算反应能但不计算能垒
                if relax_endpoints:
                    reactant_energy = reactant_result.energy
                    product_energy = product_result.energy
                else:
                    reactant_pred = self.fairchem.predict_energy(reactant_atoms, relax=False)
                    product_pred = self.fairchem.predict_energy(product_atoms, relax=False)
                    reactant_energy = reactant_pred.energy
                    product_energy = product_pred.energy

                reaction_energy = product_energy - reactant_energy

                return NEBResult(
                    reaction_name=reaction_name,
                    reactant=reactant_atoms,
                    product=product_atoms,
                    neb_frames=[reactant_atoms, product_atoms],
                    n_frames=2,
                    energies=[reactant_energy, product_energy],
                    activation_energy_forward=None,
                    activation_energy_reverse=None,
                    reaction_energy=reaction_energy,
                    transition_state_index=None,
                    transition_state=None,
                    converged=False,
                    fmax_final=np.nan,
                    optimization_steps=0,
                    trajectory_file=None,
                )

        # NOTE: Pre-relaxation of endpoints is intentionally DISABLED.
        # Testing showed that pre-relaxing adsorbate/staged atoms causes them
        # to collapse into bonded positions, eliminating the energy difference
        # between endpoints that NEB needs to find transition states.
        # Without pre-relax, NEB finds more non-zero barriers (5/8 vs 2/8)
        # with more physically reasonable values (0.07-0.72 eV vs 0-2.1 eV).

        try:
            prepared_pair = prepare_endpoint_pair_for_interpolation(
                reactant=reactant_atoms,
                product=product_atoms,
                reactant_adsorbate_indices=reactant_adsorbate_indices or reactant_atoms.info.get("adsorbate_indices", []),
                product_adsorbate_indices=product_adsorbate_indices or product_atoms.info.get("adsorbate_indices", []),
                reorder_product_atoms=True,
            )
            reactant_atoms = prepared_pair.reactant
            product_atoms = prepared_pair.product
            if prepared_pair.reordered:
                self.logger.info("Reordered product endpoint to match reactant atom order before interpolation")
        except Exception as exc:
            self.logger.warning("Shared NEB endpoint preparation failed, using raw endpoints: %s", exc)

        # 生成插值帧（使用 CatTSunami 风格的线性插值 + 迭代重叠校正）
        frames = self._interpolate_with_overlap_correction(
            reactant_atoms, product_atoms, n_frames
        )

        # 可视化插值路径
        self.logger.info(f"Visualizing interpolated NEB path before optimization...")
        try:
            from core.viz.visualization_manager import CatalysisVisualizationManager

            # 创建可视化管理器（使用work_dir作为输出目录）
            viz_manager = CatalysisVisualizationManager(
                output_dir=self.work_dir,
                quality='medium',
                renderer='tachyon',  # 使用Tachyon高质量渲染
            )

            # 可视化插值路径
            interpolated_gif = viz_manager.visualize_neb_trajectory(
                trajectory=frames,
                reaction_name=reaction_name,
                fps=2,
                is_optimized=False,  # 插值路径，非优化后的
            )

            if interpolated_gif:
                self.logger.info(f"✓ Saved interpolated path: {interpolated_gif}")
                self.logger.info(f"  Please review this GIF to check if initial structures are reasonable")
            else:
                self.logger.warning(f"  Failed to generate interpolated path visualization")

        except Exception as e:
            self.logger.warning(f"Failed to visualize interpolated path: {e}")
            self.logger.warning("Continuing with NEB calculation anyway...")

        # 运行 NEB
        neb_result = self.run_neb(
            initial_frames=frames,
            fmax=fmax,
            max_steps=max_steps,
            use_two_phase=use_two_phase,
            delta_fmax_climb=delta_fmax_climb,
            k=k,
            reaction_name=reaction_name,
        )

        # 能量校正：处理虚拟原子导致的能量偏移
        # 当使用LLM控制器添加虚拟原子时，NEB的反应物能量会改变
        # 需要应用能量偏移以使E_act与吸附能量系统一致
        if (
            adjusted_reactant is not None
            and adjusted_product is not None
            and len(adjusted_reactant) != len(reactant_atoms_original)
        ):
            n_virtual_reactant = len(adjusted_reactant) - len(reactant_atoms_original)
            n_virtual_product = (
                len(adjusted_product) - len(product_atoms_original)
                if len(adjusted_product) != len(product_atoms_original)
                else 0
            )

            if n_virtual_reactant > 0:
                # 虚拟原子添加到反应物
                # 从NEB能量减去虚拟原子的贡献（大约每个H原子-2到-4eV）
                # 更精确的做法：使用原始反应物能量作为参考

                # 获取原始反应物的DFT能量（从relax_endpoints的结果）
                if hasattr(reactant_result, 'energy'):
                    original_reactant_energy = reactant_result.energy
                else:
                    original_reactant_energy = self.fairchem.predict_energy(reactant_atoms_original, relax=False).energy

                # NEB Frame 0 的虚拟原子能量贡献
                virtual_energy_contribution = neb_result.energies[0] - (original_reactant_energy + (product_result.energy if hasattr(product_result, 'energy') else 0))

                # 实际的修正：使用原始反应物能量替代Frame 0
                # 然后重新计算能垒
                self.logger.info(f"Applying energy correction for {n_virtual_reactant} virtual atoms in reactant")

                # 从NEB能量中移除虚拟原子的贡献
                # Frame 0应该对应原始反应物的能量
                energies_corrected = np.array(neb_result.energies, dtype=float)

                # 简单的校正：保持所有帧相对关系，但重新基准化Frame 0
                energy_shift = original_reactant_energy - energies_corrected[0]
                energies_corrected = energies_corrected + energy_shift

                self.logger.info(f"  Energy shift applied: {energy_shift:+.3f} eV")
                self.logger.info(f"  Corrected reactant energy: {original_reactant_energy:.3f} eV (was {neb_result.energies[0]:.3f})")

                # 重新计算能垒
                ts_index = np.argmax(energies_corrected)
                ts_energy = energies_corrected[ts_index]

                activation_energy_forward = ts_energy - energies_corrected[0]
                activation_energy_reverse = ts_energy - energies_corrected[-1]
                reaction_energy = energies_corrected[-1] - energies_corrected[0]

                self.logger.info(f"  Corrected E_act: {activation_energy_forward:+.3f} eV (was {neb_result.activation_energy_forward:.3f})")
                self.logger.info(f"  Corrected ΔE: {reaction_energy:+.3f} eV (was {neb_result.reaction_energy:.3f})")

                # 返回校正后的结果
                return NEBResult(
                    reaction_name=reaction_name,
                    reactant=neb_result.reactant,
                    product=neb_result.product,
                    neb_frames=neb_result.neb_frames,
                    n_frames=neb_result.n_frames,
                    energies=energies_corrected.tolist(),
                    activation_energy_forward=activation_energy_forward,
                    activation_energy_reverse=activation_energy_reverse,
                    reaction_energy=reaction_energy,
                    transition_state_index=ts_index,
                    transition_state=neb_result.neb_frames[ts_index],
                    converged=neb_result.converged,
                    fmax_final=neb_result.fmax_final,
                    optimization_steps=neb_result.optimization_steps,
                    trajectory_file=neb_result.trajectory_file,
                )

        return neb_result

    def predict_from_reaction(
        self,
        reaction_str: str,
        surface: Union[str, Atoms],
        n_frames: int = 10,
        n_pdt1_sites: int = 4,
        n_pdt2_sites: int = 4,
        fmax: float = 0.1,
        max_steps: int = 300,
    ) -> NEBResult:
        """
        从反应式预测能垒（使用 CatTSunami AutoFrame）

        Parameters
        ----------
        reaction_str : str
            反应式，如 "*CH -> *C + *H"
        surface : str or Atoms
            表面结构
        n_frames : int, default=10
            NEB 帧数
        n_pdt1_sites : int, default=4
            产物1的位点数
        n_pdt2_sites : int, default=4
            产物2的位点数
        fmax : float, default=0.1
            力收敛标准
        max_steps : int, default=300
            最大优化步数

        Returns
        -------
        NEBResult
            NEB 计算结果
        """
        try:
            from fairchem.applications.cattsunami.core import Reaction
            from fairchem.applications.cattsunami.core.autoframe import (
                AutoFrameDissociation
            )
            from fairchem.applications.cattsunami.databases import (
                DISSOCIATION_REACTION_DB_PATH
            )
            from fairchem.data.oc.core import (
                Adsorbate, Bulk, Slab, AdsorbateSlabConfig
            )
            from fairchem.data.oc.databases.pkls import (
                ADSORBATE_PKL_PATH, BULK_PKL_PATH
            )
        except ImportError as e:
            raise ImportError(
                f"Failed to import cattsunami modules. Error: {e}"
            )

        self.logger.info(f"Setting up reaction: {reaction_str}")

        # 创建反应对象
        reaction = Reaction(
            reaction_str_from_db=reaction_str,
            reaction_db_path=DISSOCIATION_REACTION_DB_PATH,
            adsorbate_db_path=ADSORBATE_PKL_PATH
        )

        # 读取表面
        if isinstance(surface, str):
            surface_atoms = read(surface)
        else:
            surface_atoms = surface.copy()

        # 创建 Slab 对象
        slab = Slab(atoms=surface_atoms)

        # 1. 生成并优化反应物配置
        self.logger.info("Generating reactant configurations...")
        reactant = Adsorbate(
            adsorbate_id_from_db=reaction.reactant1_idx,
            adsorbate_db_path=ADSORBATE_PKL_PATH
        )

        reactant_configs = AdsorbateSlabConfig(
            slab=slab,
            adsorbate=reactant,
            mode="random_site_heuristic_placement",
            num_sites=10
        ).atoms_list

        # 优化反应物
        reactant_energies = []
        for config in reactant_configs:
            result = self.fairchem.predict_energy(
                config, relax=True, fmax=fmax, max_steps=max_steps
            )
            reactant_energies.append(result.energy)

        best_reactant_idx = np.argmin(reactant_energies)
        reactant_system = reactant_configs[best_reactant_idx]

        # 2. 生成并优化产物配置
        self.logger.info("Generating product configurations...")

        product1 = Adsorbate(
            adsorbate_id_from_db=reaction.product1_idx,
            adsorbate_db_path=ADSORBATE_PKL_PATH
        )
        product2 = Adsorbate(
            adsorbate_id_from_db=reaction.product2_idx,
            adsorbate_db_path=ADSORBATE_PKL_PATH
        )

        product1_configs = AdsorbateSlabConfig(
            slab=slab,
            adsorbate=product1,
            mode="random_site_heuristic_placement",
            num_sites=n_pdt1_sites
        ).atoms_list

        product2_configs = AdsorbateSlabConfig(
            slab=slab,
            adsorbate=product2,
            mode="random_site_heuristic_placement",
            num_sites=n_pdt2_sites
        ).atoms_list

        # 优化产物
        product1_energies = []
        for config in product1_configs:
            result = self.fairchem.predict_energy(
                config, relax=True, fmax=fmax, max_steps=max_steps
            )
            product1_energies.append(result.energy)

        product2_energies = []
        for config in product2_configs:
            result = self.fairchem.predict_energy(
                config, relax=True, fmax=fmax, max_steps=max_steps
            )
            product2_energies.append(result.energy)

        # 3. 使用 AutoFrame 生成 NEB 帧
        self.logger.info("Generating NEB frames with AutoFrame...")

        af = AutoFrameDissociation(
            reaction=reaction,
            reactant_system=reactant_system,
            product1_systems=product1_configs,
            product1_energies=product1_energies,
            product2_systems=product2_configs,
            product2_energies=product2_energies,
        )

        frame_sets, mapping_idxs = af.get_neb_frames(
            self.fairchem._calculator,
            n_frames=n_frames,
            n_pdt1_sites=n_pdt1_sites,
            n_pdt2_sites=n_pdt2_sites,
        )

        if not frame_sets:
            raise RuntimeError("Failed to generate NEB frames")

        # 4. 运行 NEB（使用第一组帧）
        self.logger.info(f"Running NEB with {len(frame_sets[0])} frames...")

        return self.run_neb(
            initial_frames=frame_sets[0],
            fmax=fmax,
            max_steps=max_steps,
            reaction_name=reaction_str,
        )

    def predict_multiple_barriers(
        self,
        surface: Union[str, Atoms],
        reactions: List[Tuple[Union[str, Atoms], Union[str, Atoms], str]],
        n_frames: int = 10,
        fmax: float = 0.1,
        max_steps: int = 300,
        output_dir: Optional[str] = None,
    ) -> ReactionBarriersResult:
        """
        预测多个反应的能垒

        Parameters
        ----------
        surface : str or Atoms
            表面结构
        reactions : list of tuple
            反应列表，每个元组为 (reactant, product, name)
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            力收敛标准
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录

        Returns
        -------
        ReactionBarriersResult
            多个反应能垒的结果
        """
        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "barriers")
        os.makedirs(output_dir, exist_ok=True)

        # 读取表面
        if isinstance(surface, str):
            surface_atoms = read(surface)
        else:
            surface_atoms = surface.copy()

        surface_formula = surface_atoms.get_chemical_formula()

        # 预测每个反应
        barriers = {}
        reaction_names = []

        for reactant, product, name in reactions:
            self.logger.info("=" * 60)
            self.logger.info(f"Processing reaction: {name}")

            # 计算吸附分子的原子索引
            surface_atom_count = len(surface_atoms)
            reactant_atoms = read(reactant) if isinstance(reactant, str) else reactant
            product_atoms = read(product) if isinstance(product, str) else product
            reactant_adsorbate_indices = list(range(surface_atom_count, len(reactant_atoms)))
            product_adsorbate_indices = list(range(surface_atom_count, len(product_atoms)))

            try:
                result = self.predict_from_structures(
                    reactant=reactant,
                    product=product,
                    n_frames=n_frames,
                    fmax=fmax,
                    max_steps=max_steps,
                    reaction_name=name,
                    reactant_adsorbate_indices=reactant_adsorbate_indices,
                    product_adsorbate_indices=product_adsorbate_indices,
                )

                barriers[name] = result
                reaction_names.append(name)

                # 保存轨迹
                traj_file = os.path.join(output_dir, f"{name}_neb.traj")
                if result.trajectory_file and result.trajectory_file != traj_file:
                    import shutil
                    shutil.copy(result.trajectory_file, traj_file)

                self.logger.info(
                    f"  E_act = {result.activation_energy_forward:.3f} eV"
                )

            except Exception as e:
                self.logger.error(f"Failed to process {name}: {e}")
                continue

        return ReactionBarriersResult(
            surface_formula=surface_formula,
            reactions=reaction_names,
            barriers=barriers,
            output_dir=output_dir,
        )

    def __del__(self):
        """清理临时文件"""
        if not self.keep_files and self._temp_dir and os.path.exists(self._temp_dir):
            import shutil
            shutil.rmtree(self._temp_dir)


if __name__ == "__main__":
    print("Barrier Predictor module loaded successfully")
    print("See documentation for usage examples")
