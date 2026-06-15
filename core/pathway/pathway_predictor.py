"""
Pathway Predictor - 反应路径完整分析模块

这个模块整合了自由能预测和能垒预测，提供完整的反应路径分析功能。

对于给定的催化剂表面和反应路径，这个模块可以：
1. 预测每个反应中间体的吸附能（自由能）
2. 预测相邻中间体之间的反应能垒
3. 生成完整的能量剖面图
4. 识别速率决定步骤

使用方式:
    from pathway_predictor import PathwayPredictor

    # 初始化预测器
    predictor = PathwayPredictor(
        fairchem_root="/path/to/fairchem",
        model_name="uma-s-1p1",
        use_gpu=True,
    )

    # 方式1: 从吸附质列表预测完整反应路径
    result = predictor.predict_pathway(
        surface="Pt_111.vasp",
        adsorbates=["*O", "*OH", "*OOH", "*"],  # 按反应顺序
        calculate_barriers=True,  # 是否计算能垒
    )

    # 打印结果
    print(result.summary())

    # 绘制能量剖面图
    predictor.plot_energy_profile(result, output="energy_profile.png")

    # 方式2: 手动指定反应步骤
    result = predictor.predict_pathway_from_reactions(
        surface="Pt_111.vasp",
        reactions=[
            ("*O", "*OH", "O_hydrogenation"),
            ("*OH", "*OOH", "OH_hydrogenation"),
            ("*OOH", "*", "OOH_desorption"),
        ]
    )

作者: Claude
"""

import os
import sys
import tempfile
import logging
import re
from typing import Union, Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from collections import Counter
import warnings

import numpy as np
from ase import Atoms
from ase.io import read, write

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False


@dataclass
class PathwayStep:
    """反应路径中的单个步骤"""
    step_index: int
    name: str
    reactant_adsorbate: str
    product_adsorbate: str
    reactant_energy: float  # eV
    product_energy: float  # eV
    reaction_energy: float  # eV — gas-phase-corrected ΔE (pathway convention)
    activation_energy: Optional[float] = None  # eV
    transition_state_energy: Optional[float] = None  # eV
    is_rate_determining: bool = False
    # Raw NEB reactant->product energy difference (same atom count on both
    # endpoints because NEB stages the transferred atoms). This is the correct
    # ΔE for CatMAP's split elementary reactions like A_s -> B_s + H_s and
    # must be used in place of `reaction_energy` when propagating consistent
    # formation energies across a chain with atomic intermediates.
    neb_reaction_energy: Optional[float] = None

    def __repr__(self):
        barrier_str = f", E_act={self.activation_energy:.3f} eV" if self.activation_energy else ""
        return (f"PathwayStep({self.name}: {self.reactant_adsorbate} -> {self.product_adsorbate}, "
                f"ΔE={self.reaction_energy:.3f} eV{barrier_str})")


@dataclass
class CompletePathwayResult:
    """完整反应路径分析结果"""
    surface_formula: str
    adsorbates: List[str]
    adsorbate_energies: Dict[str, float]  # adsorbate -> total energy (eV)
    adsorption_energies: Dict[str, float]  # adsorbate -> E_ads (eV)
    steps: List[PathwayStep]
    rate_determining_step: Optional[PathwayStep]
    max_barrier: Optional[float]  # eV
    overall_reaction_energy: float  # eV
    output_dir: str
    # Vacuum-phase vibrational frequencies (cm^-1) per adsorbate label.
    # Fed to CatMAP as the `frequencies` column so ideal_gas thermo can
    # compute translational/rotational/vibrational entropy at T.
    gas_frequencies: Dict[str, List[float]] = field(default_factory=dict)

    def summary(self) -> str:
        """生成结果摘要"""
        lines = [
            "Complete Reaction Pathway Analysis",
            "=" * 80,
            f"Surface: {self.surface_formula}",
            f"Number of intermediates: {len(self.adsorbates)}",
            f"Number of steps: {len(self.steps)}",
            f"Overall reaction energy: {self.overall_reaction_energy:.4f} eV",
            "",
            "=" * 80,
            "Adsorption Energies:",
            "-" * 80,
        ]

        for ads in self.adsorbates:
            e_ads = self.adsorption_energies.get(ads, 0.0)
            e_total = self.adsorbate_energies.get(ads, 0.0)
            lines.append(f"  {ads:15s}: E_ads = {e_ads:7.3f} eV, E_total = {e_total:10.3f} eV")

        lines.append("")
        lines.append("=" * 80)
        lines.append("Reaction Steps:")
        lines.append("-" * 80)

        for step in self.steps:
            rds_marker = " [RDS]" if step.is_rate_determining else ""
            lines.append(f"\nStep {step.step_index + 1}: {step.name}{rds_marker}")
            lines.append(f"  {step.reactant_adsorbate} -> {step.product_adsorbate}")
            lines.append(f"  Reaction energy (ΔE): {step.reaction_energy:+.4f} eV")

            if step.activation_energy is not None:
                lines.append(f"  Activation energy (E_act): {step.activation_energy:.4f} eV")
                if step.transition_state_energy is not None:
                    lines.append(f"  Transition state energy: {step.transition_state_energy:.4f} eV")

        if self.rate_determining_step:
            lines.append("")
            lines.append("=" * 80)
            lines.append("Rate-Determining Step:")
            lines.append(f"  Step {self.rate_determining_step.step_index + 1}: {self.rate_determining_step.name}")
            lines.append(f"  Barrier: {self.rate_determining_step.activation_energy:.4f} eV")

        if self.max_barrier is not None:
            lines.append(f"\nMaximum barrier: {self.max_barrier:.4f} eV")

        lines.append("")
        lines.append("=" * 80)
        lines.append(f"Output directory: {self.output_dir}")

        return "\n".join(lines)


class PathwayPredictor:
    """
    反应路径完整分析预测器

    这个类整合了自由能预测和能垒预测，提供完整的反应路径分析。

    Parameters
    ----------
    fairchem_root : str
        Fairchem 项目根目录
    model_name : str, default="uma-s-1p1"
        使用的模型名称
    model_path : str, optional
        本地模型文件路径
    use_gpu : bool, default=True
        是否使用 GPU
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """

    # 原子参考能量 (eV) - 与 FairchemPredictor 保持一致
    # 当原子从表面离开或加入时，使用这些能量进行校正
    ATOMIC_REFERENCE_ENERGIES = {
        "H": -3.477,
        "C": -7.282,
        "N": -8.083,
        "O": -7.204,
        "F": -4.891,
        "S": -4.659,
    }

    @staticmethod
    def parse_adsorbate_formula(formula: str) -> Counter:
        """
        解析吸附物化学式，返回元素组成

        支持格式: *CO, *CHOH, *CH3, *, CO, CHOH 等

        Parameters
        ----------
        formula : str
            吸附物化学式

        Returns
        -------
        Counter
            元素组成计数，如 {'C': 1, 'H': 1, 'O': 1} for *CHOH
        """
        # 移除 * 前缀
        formula = formula.strip().lstrip('*').strip()

        # 空字符串（纯表面）返回空计数
        if not formula:
            return Counter()

        # 解析化学式
        # 匹配模式: 元素符号(大写+可选小写) + 可选数字
        pattern = r'([A-Z][a-z]?)(\d*)'
        composition = Counter()

        for match in re.finditer(pattern, formula):
            element = match.group(1)
            count = int(match.group(2)) if match.group(2) else 1
            if element:  # 忽略空匹配
                composition[element] += count

        return composition

    def calculate_element_change(
        self,
        reactant_formula: str,
        product_formula: str
    ) -> Dict[str, int]:
        """
        计算从反应物到产物的元素变化

        Parameters
        ----------
        reactant_formula : str
            反应物化学式 (如 *CHOH)
        product_formula : str
            产物化学式 (如 *CH)

        Returns
        -------
        Dict[str, int]
            元素变化量，正值表示产物多出的元素（需要从气相加入），
            负值表示反应物多出的元素（离开表面进入气相）
        """
        reactant_comp = self.parse_adsorbate_formula(reactant_formula)
        product_comp = self.parse_adsorbate_formula(product_formula)

        # 计算差异
        change = {}
        all_elements = set(reactant_comp.keys()) | set(product_comp.keys())

        for elem in all_elements:
            diff = product_comp[elem] - reactant_comp[elem]
            if diff != 0:
                change[elem] = diff

        return change

    def calculate_gas_phase_correction(
        self,
        element_change: Dict[str, int]
    ) -> float:
        """
        计算气相能量校正

        当元素从表面离开时，需要加上该元素的参考能量（因为产物包含气相分子）
        当元素加入表面时，需要减去该元素的参考能量（因为反应物包含气相分子）

        对于反应 *CHOH → *CH:
        - 元素变化: O: -1, H: -1 (OH离开表面)
        - 校正 = -(-1) * E_O + -(-1) * E_H = E_O + E_H
        - 即产物侧需要加上离开的原子的能量

        Parameters
        ----------
        element_change : Dict[str, int]
            元素变化量 (product - reactant)

        Returns
        -------
        float
            气相能量校正 (eV)
        """
        correction = 0.0

        for element, delta in element_change.items():
            if element in self.ATOMIC_REFERENCE_ENERGIES:
                # delta < 0 表示原子离开表面（产物有气相分子）
                # delta > 0 表示原子加入表面（反应物有气相分子）
                # 校正 = -delta * E_ref
                # 这样：离开的原子会增加产物能量，加入的原子会增加反应物能量
                correction -= delta * self.ATOMIC_REFERENCE_ENERGIES[element]

                if self.verbose:
                    if delta < 0:
                        self.logger.info(
                            f"    Gas correction: {-delta} {element} atom(s) leaving surface, "
                            f"adding {-delta * self.ATOMIC_REFERENCE_ENERGIES[element]:.3f} eV"
                        )
                    else:
                        self.logger.info(
                            f"    Gas correction: {delta} {element} atom(s) joining surface, "
                            f"subtracting {delta * self.ATOMIC_REFERENCE_ENERGIES[element]:.3f} eV"
                        )
            else:
                if self.verbose:
                    self.logger.warning(
                        f"    No reference energy for element {element}, skipping correction"
                    )

        return correction

    def __init__(
        self,
        fairchem_root: str,
        model_name: str = "uma-s-1p1",
        model_path: Optional[str] = None,
        use_gpu: bool = True,
        work_dir: Optional[str] = None,
        keep_files: bool = True,
        verbose: bool = True,
        use_llm_controller: bool = False,
        llm_model: Optional[str] = None,
    ):
        self.fairchem_root = os.path.abspath(fairchem_root)
        self.model_name = model_name
        self.model_path = model_path
        self.use_gpu = use_gpu
        self.keep_files = keep_files
        self.verbose = verbose
        self.use_llm_controller = use_llm_controller
        self.llm_model = llm_model

        # 设置工作目录
        if work_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="pathway_")
            self.work_dir = self._temp_dir
        else:
            self._temp_dir = None
            self.work_dir = os.path.abspath(work_dir)
            os.makedirs(self.work_dir, exist_ok=True)

        # 初始化日志
        self.logger = logging.getLogger(__name__)
        if self.verbose:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.WARNING)

        # 延迟初始化子模块
        self._fairchem_predictor = None
        self._barrier_predictor = None

    def _init_predictors(self):
        """初始化子预测器"""
        if self._fairchem_predictor is not None:
            return

        from .fairchem_predictor import FairchemPredictor
        from .barrier_predictor import BarrierPredictor

        self.logger.info("Initializing predictors...")

        self._fairchem_predictor = FairchemPredictor(
            fairchem_root=self.fairchem_root,
            model_name=self.model_name,
            model_path=self.model_path,  # 传递本地模型路径
            use_gpu=self.use_gpu,
            work_dir=os.path.join(self.work_dir, "fairchem"),
            keep_files=self.keep_files,
            verbose=self.verbose,
        )

        self._barrier_predictor = BarrierPredictor(
            fairchem_predictor=self._fairchem_predictor,
            work_dir=os.path.join(self.work_dir, "barriers"),
            keep_files=self.keep_files,
            verbose=self.verbose,
            use_llm_controller=self.use_llm_controller,
            llm_model=self.llm_model,
        )

        self.logger.info("Predictors initialized")

    def predict_pathway(
        self,
        surface: Union[str, Atoms],
        adsorbates: List[Union[str, Atoms]],
        calculate_barriers: bool = True,
        num_sites: int = 5,
        n_frames: int = 10,
        fmax: float = 0.1,
        max_steps: int = 300,
        output_dir: Optional[str] = None,
        current_adsorbate_indices: Optional[List[int]] = None,
    ) -> CompletePathwayResult:
        """
        预测完整反应路径

        Parameters
        ----------
        surface : str or Atoms
            催化剂表面结构（可能包含吸附物）
        adsorbates : list of str or Atoms
            反应路径上的吸附质列表（按顺序）
        calculate_barriers : bool, default=True
            是否计算反应能垒
        num_sites : int, default=5
            每个吸附质尝试的位点数
        n_frames : int, default=10
            NEB 计算的帧数
        fmax : float, default=0.1
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录
        current_adsorbate_indices : list of int, optional
            当前吸附物的原子索引（如果表面已有吸附物）

        Returns
        -------
        CompletePathwayResult
            完整反应路径分析结果
        """
        self._init_predictors()

        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "pathway_complete")
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info("=" * 80)
        self.logger.info("Starting complete pathway analysis")
        self.logger.info("=" * 80)

        # 读取表面
        if isinstance(surface, str):
            surface_atoms = read(surface)
        else:
            surface_atoms = surface.copy()

        surface_formula = surface_atoms.get_chemical_formula()

        # 如果表面已有吸附物，找到吸附位点并移除吸附物
        adsorption_site = None
        clean_surface = surface_atoms
        if current_adsorbate_indices is not None and len(current_adsorbate_indices) > 0:
            self.logger.info(f"Surface contains {len(current_adsorbate_indices)} adsorbate atoms")

            # 找到吸附物最下面的原子（z坐标最小）作为吸附位点
            adsorbate_positions = surface_atoms.positions[current_adsorbate_indices]
            lowest_atom_idx = np.argmin(adsorbate_positions[:, 2])
            adsorption_site = adsorbate_positions[lowest_atom_idx]

            self.logger.info(f"  Adsorption site determined from adsorbate: {adsorption_site}")

            # 移除吸附物得到clean表面
            clean_surface = surface_atoms.copy()
            del clean_surface[current_adsorbate_indices]

            self.logger.info(f"  Clean surface after removing adsorbate: {clean_surface.get_chemical_formula()}")

            # 重要：重新设置clean表面的tags，确保有正确的0/1分布
            # 重构后的表面可能所有原子tags都是0（因为都被fixed了）
            # 我们需要确保表面原子是tag=1，体相原子是tag=0
            # 同时，FAIRChem要求原子顺序：bulk(0) -> surface(1) -> adsorbate(2)
            if clean_surface.has('tags'):
                from ase import Atoms as ASEAtoms

                z_coords = clean_surface.positions[:, 2]
                z_max = z_coords.max()
                # 顶部3层原子作为表面（tag=1），其余为体相（tag=0）
                is_surface = z_coords > z_max - 3.0

                # 关键修复：重新排序原子，使bulk原子在前，surface原子在后
                # 这样tags才能正确对应到原子
                bulk_indices = np.where(~is_surface)[0]
                surface_indices = np.where(is_surface)[0]

                # 创建新的排序
                new_order = np.concatenate([bulk_indices, surface_indices])

                # 重新排序：创建一个新的plain ASE Atoms对象
                # 这避免了AtomsBatch等子类的内部状态问题
                old_atoms = clean_surface
                new_positions = old_atoms.positions[new_order]
                new_symbols = [old_atoms.get_chemical_symbols()[i] for i in new_order]

                clean_surface = ASEAtoms(
                    symbols=new_symbols,
                    positions=new_positions,
                    cell=old_atoms.cell,
                    pbc=old_atoms.pbc,
                )

                # 如果有constraints，也需要重新映射
                if old_atoms.constraints:
                    # 创建索引映射：old_index -> new_index
                    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(new_order)}
                    new_constraints = []
                    for constraint in old_atoms.constraints:
                        if hasattr(constraint, 'index'):
                            # 重新映射constraint的索引
                            from ase.constraints import FixAtoms
                            old_indices = constraint.index if hasattr(constraint.index, '__iter__') else [constraint.index]
                            new_indices = [index_map[i] for i in old_indices if i in index_map]
                            if new_indices:
                                new_constraints.append(FixAtoms(indices=new_indices))
                    if new_constraints:
                        clean_surface.set_constraint(new_constraints)

                # 设置tags：前N个是bulk(0)，后M个是surface(1)
                n_bulk = len(bulk_indices)
                n_surface = len(surface_indices)
                new_tags = [0] * n_bulk + [1] * n_surface
                clean_surface.set_tags(new_tags)

                self.logger.info(f"  Reset tags for clean surface: {n_bulk} bulk, {n_surface} surface atoms")

        # 1. 预测所有吸附质的能量
        self.logger.info("\n" + "=" * 80)
        self.logger.info("PHASE 1: Predicting adsorption energies")
        self.logger.info("=" * 80)

        pathway_result = self._fairchem_predictor.predict_pathway_energies(
            surface=clean_surface,
            adsorbates=adsorbates,
            num_sites=num_sites,
            fmax=fmax,
            max_steps=max_steps,
            output_dir=os.path.join(output_dir, "adsorption"),
            fixed_adsorption_site=adsorption_site,  # 使用固定的吸附位点
        )

        adsorbate_energies = pathway_result.adsorbate_energies
        adsorption_energies = pathway_result.adsorption_energies
        best_configs = pathway_result.best_configurations

        # 2. 计算相邻步骤之间的能垒
        steps = []
        previous_product_struct = None  # Track previous step's product for continuity
        previous_product_adsorbate_indices = None

        if calculate_barriers and len(adsorbates) > 1:
            self.logger.info("\n" + "=" * 80)
            self.logger.info("PHASE 2: Calculating reaction barriers")
            self.logger.info("=" * 80)

            for i in range(len(adsorbates) - 1):
                reactant_ads = str(adsorbates[i])
                product_ads = str(adsorbates[i + 1])
                step_name = f"{reactant_ads}_to_{product_ads}"

                self.logger.info(f"\nStep {i + 1}: {reactant_ads} -> {product_ads}")

                reactant_energy = adsorbate_energies[reactant_ads]
                product_energy = adsorbate_energies[product_ads]

                # 计算元素变化并应用气相能量校正
                element_change = self.calculate_element_change(reactant_ads, product_ads)
                gas_correction = 0.0

                if element_change:
                    self.logger.info(f"  Element change detected: {element_change}")
                    gas_correction = self.calculate_gas_phase_correction(element_change)
                    self.logger.info(f"  Total gas-phase correction: {gas_correction:+.3f} eV")

                # 反应能 = 产物能量 - 反应物能量 + 气相校正
                # 气相校正已经包含了离开/加入表面的原子能量
                reaction_energy = product_energy - reactant_energy + gas_correction

                # 获取反应物和产物结构
                reactant_struct = best_configs[reactant_ads]
                product_struct = best_configs[product_ads]

                # 计算吸附分子的原子索引（表面原子之后的所有原子）
                # 使用clean_surface的长度，因为best_configs中的结构也是基于clean_surface的
                surface_atom_count = len(clean_surface)
                reactant_adsorbate_indices = list(range(surface_atom_count, len(reactant_struct)))
                product_adsorbate_indices = list(range(surface_atom_count, len(product_struct)))

                try:
                    # 计算能垒 - 传递前一步的产物信息以确保路径连续性
                    neb_result = self._barrier_predictor.predict_from_structures(
                        reactant=reactant_struct,
                        product=product_struct,
                        n_frames=n_frames,
                        fmax=fmax,
                        max_steps=max_steps,
                        relax_endpoints=False,  # 已经优化过了
                        reaction_name=step_name,
                        reactant_adsorbate_indices=reactant_adsorbate_indices,
                        product_adsorbate_indices=product_adsorbate_indices,
                        previous_step_product=previous_product_struct,
                        previous_step_product_adsorbate_indices=previous_product_adsorbate_indices,
                        step_index=i+1,
                    )

                    # Update previous product for next iteration
                    previous_product_struct = product_struct
                    previous_product_adsorbate_indices = product_adsorbate_indices

                    activation_energy = neb_result.activation_energy_forward
                    # 只有当 transition_state_index 不为 None 时才获取 ts_energy
                    if neb_result.transition_state_index is not None:
                        ts_energy = neb_result.energies[neb_result.transition_state_index]
                    else:
                        ts_energy = None

                    if activation_energy is not None:
                        self.logger.info(f"  ΔE = {reaction_energy:+.3f} eV, E_act = {activation_energy:.3f} eV")
                    else:
                        self.logger.info(f"  ΔE = {reaction_energy:+.3f} eV (barrier calculation not available)")

                except Exception as e:
                    self.logger.error(f"Failed to calculate barrier for step {i + 1}: {e}")
                    activation_energy = None
                    ts_energy = None

                # 创建步骤对象
                step = PathwayStep(
                    step_index=i,
                    name=step_name,
                    reactant_adsorbate=reactant_ads,
                    product_adsorbate=product_ads,
                    reactant_energy=reactant_energy,
                    product_energy=product_energy,
                    reaction_energy=reaction_energy,
                    activation_energy=activation_energy,
                    transition_state_energy=ts_energy,
                )
                steps.append(step)

        else:
            # 不计算能垒，只记录反应能
            for i in range(len(adsorbates) - 1):
                reactant_ads = str(adsorbates[i])
                product_ads = str(adsorbates[i + 1])

                reactant_energy = adsorbate_energies[reactant_ads]
                product_energy = adsorbate_energies[product_ads]

                # 计算元素变化并应用气相能量校正
                element_change = self.calculate_element_change(reactant_ads, product_ads)
                gas_correction = 0.0

                if element_change:
                    self.logger.info(f"  Step {i+1}: {reactant_ads} -> {product_ads}")
                    self.logger.info(f"    Element change: {element_change}")
                    gas_correction = self.calculate_gas_phase_correction(element_change)
                    self.logger.info(f"    Gas-phase correction: {gas_correction:+.3f} eV")

                reaction_energy = product_energy - reactant_energy + gas_correction

                step = PathwayStep(
                    step_index=i,
                    name=f"{reactant_ads}_to_{product_ads}",
                    reactant_adsorbate=reactant_ads,
                    product_adsorbate=product_ads,
                    reactant_energy=reactant_energy,
                    product_energy=product_energy,
                    reaction_energy=reaction_energy,
                )
                steps.append(step)

        # 3. 识别速率决定步骤（最大能垒）
        rate_determining_step = None
        max_barrier = None

        if calculate_barriers:
            valid_barriers = [
                (step, step.activation_energy)
                for step in steps
                if step.activation_energy is not None
            ]

            if valid_barriers:
                rate_determining_step = max(valid_barriers, key=lambda x: x[1])[0]
                rate_determining_step.is_rate_determining = True
                max_barrier = rate_determining_step.activation_energy

                self.logger.info("\n" + "=" * 80)
                self.logger.info("Rate-determining step identified:")
                self.logger.info(f"  {rate_determining_step.name}")
                self.logger.info(f"  Barrier: {max_barrier:.4f} eV")

        # 4. 计算总反应能
        overall_reaction_energy = (
            adsorbate_energies[str(adsorbates[-1])] -
            adsorbate_energies[str(adsorbates[0])]
        )

        self.logger.info("\n" + "=" * 80)
        self.logger.info("Pathway analysis complete")
        self.logger.info("=" * 80)

        return CompletePathwayResult(
            surface_formula=surface_formula,
            adsorbates=[str(ads) for ads in adsorbates],
            adsorbate_energies=adsorbate_energies,
            adsorption_energies=adsorption_energies,
            steps=steps,
            rate_determining_step=rate_determining_step,
            max_barrier=max_barrier,
            overall_reaction_energy=overall_reaction_energy,
            output_dir=output_dir,
        )

    def predict_barriers_from_structures(
        self,
        steps: List[Dict[str, Any]],
        n_frames: int = 10,
        fmax: float = 0.1,
        max_steps: int = 300,
        relax_endpoints: bool = False,
        output_dir: Optional[str] = None,
    ) -> CompletePathwayResult:
        """
        从预计算的反应物/产物结构计算能垒

        适用于已有 Agent4/5 生成的端点结构，跳过吸附能计算阶段，
        直接对每步运行 NEB 并生成完整路径结果。

        Parameters
        ----------
        steps : list of dict
            每步包含:
              - name: str, 步骤名称
              - reactant: Atoms, 反应物结构
              - product: Atoms, 产物结构
              - reactant_formula: str, 反应物标签 (如 "*CO")
              - product_formula: str, 产物标签 (如 "*CHO")
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            CI-NEB 力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        relax_endpoints : bool, default=False
            是否先优化端点结构
        output_dir : str, optional
            输出目录

        Returns
        -------
        CompletePathwayResult
            完整反应路径分析结果
        """
        self._init_predictors()

        if output_dir is None:
            output_dir = os.path.join(self.work_dir, "pathway_from_structures")
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info("=" * 80)
        self.logger.info("Calculating barriers from pre-computed structures")
        self.logger.info(f"  Steps: {len(steps)}, fmax={fmax}, relax_endpoints={relax_endpoints}")
        self.logger.info("=" * 80)

        # Collect adsorbate labels
        adsorbates = [steps[0]["reactant_formula"]]
        for s in steps:
            adsorbates.append(s["product_formula"])

        # Surface formula from first reactant
        surface_formula = steps[0]["reactant"].get_chemical_formula()

        # Placeholder energies (NEB gives relative, not absolute)
        adsorbate_energies = {}
        adsorption_energies = {}

        pathway_steps = []
        previous_product_struct = None
        previous_product_ads_indices = None

        for i, step in enumerate(steps):
            name = step["name"]
            reactant = step["reactant"]
            product = step["product"]
            reactant_label = step["reactant_formula"]
            product_label = step["product_formula"]

            self.logger.info(f"\nStep {i+1}: {reactant_label} -> {product_label}")

            # Detect adsorbate indices (non-surface atoms)
            reactant_ads = list(step.get("reactant_adsorbate_indices", []))
            product_ads = list(step.get("product_adsorbate_indices", []))

            # Staged (parked) atom indices for constrained endpoint relaxation
            r_staged = list(step.get("reactant_staged_indices", []))
            p_staged = list(step.get("product_staged_indices", []))

            # Gas-phase correction for element changes
            element_change = self.calculate_element_change(reactant_label, product_label)
            gas_correction = 0.0
            if element_change:
                gas_correction = self.calculate_gas_phase_correction(element_change)
                self.logger.info(f"  Gas-phase correction: {gas_correction:+.3f} eV")

            activation_energy = None
            ts_energy = None
            reaction_energy = 0.0

            try:
                neb_result = self._barrier_predictor.predict_from_structures(
                    reactant=reactant,
                    product=product,
                    n_frames=n_frames,
                    fmax=fmax,
                    max_steps=max_steps,
                    relax_endpoints=relax_endpoints,
                    reaction_name=name,
                    reactant_adsorbate_indices=reactant_ads,
                    product_adsorbate_indices=product_ads,
                    reactant_staged_indices=r_staged if r_staged else None,
                    product_staged_indices=p_staged if p_staged else None,
                    previous_step_product=previous_product_struct,
                    previous_step_product_adsorbate_indices=previous_product_ads_indices,
                    step_index=i + 1,
                )

                previous_product_struct = product
                previous_product_ads_indices = product_ads

                activation_energy = neb_result.activation_energy_forward
                reaction_energy = neb_result.reaction_energy
                if neb_result.transition_state_index is not None:
                    ts_energy = neb_result.energies[neb_result.transition_state_index]

                self.logger.info(
                    f"  Ea_fwd={activation_energy:.3f} eV, "
                    f"Ea_rev={neb_result.activation_energy_reverse:.3f} eV, "
                    f"dE={reaction_energy:.3f} eV, "
                    f"converged={neb_result.converged}"
                )

            except Exception as e:
                self.logger.error(f"  NEB failed for step {i+1}: {e}")

            # Store energies for profile
            reactant_e = neb_result.energies[0] if activation_energy is not None else 0.0
            product_e = neb_result.energies[-1] if activation_energy is not None else 0.0
            if reactant_label not in adsorbate_energies:
                adsorbate_energies[reactant_label] = reactant_e
            adsorbate_energies[product_label] = product_e

            pathway_steps.append(PathwayStep(
                step_index=i,
                name=name,
                reactant_adsorbate=reactant_label,
                product_adsorbate=product_label,
                reactant_energy=reactant_e,
                product_energy=product_e,
                reaction_energy=reaction_energy,
                activation_energy=activation_energy,
                transition_state_energy=ts_energy,
            ))

        # Identify rate-determining step
        rate_determining_step = None
        max_barrier = None
        valid_barriers = [
            (s, s.activation_energy)
            for s in pathway_steps
            if s.activation_energy is not None and s.activation_energy > 0
        ]
        if valid_barriers:
            rate_determining_step = max(valid_barriers, key=lambda x: x[1])[0]
            rate_determining_step.is_rate_determining = True
            max_barrier = rate_determining_step.activation_energy
            self.logger.info(f"\nRate-determining step: {rate_determining_step.name}")
            self.logger.info(f"  Barrier: {max_barrier:.4f} eV")

        # Overall reaction energy
        first_e = adsorbate_energies.get(adsorbates[0], 0.0)
        last_e = adsorbate_energies.get(adsorbates[-1], 0.0)
        overall_reaction_energy = last_e - first_e

        return CompletePathwayResult(
            surface_formula=surface_formula,
            adsorbates=adsorbates,
            adsorbate_energies=adsorbate_energies,
            adsorption_energies=adsorption_energies,
            steps=pathway_steps,
            rate_determining_step=rate_determining_step,
            max_barrier=max_barrier,
            overall_reaction_energy=overall_reaction_energy,
            output_dir=output_dir,
        )

    def predict_pathway_from_reactions(
        self,
        surface: Union[str, Atoms],
        reactions: List[Tuple[Union[str, Atoms], Union[str, Atoms], str]],
        num_sites: int = 5,
        n_frames: int = 10,
        fmax: float = 0.1,
        max_steps: int = 300,
        output_dir: Optional[str] = None,
    ) -> CompletePathwayResult:
        """
        从反应列表预测路径

        Parameters
        ----------
        surface : str or Atoms
            催化剂表面
        reactions : list of tuple
            反应列表，每个元组为 (reactant, product, name)
        num_sites : int, default=5
            位点数
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
        CompletePathwayResult
            完整反应路径分析结果
        """
        # 提取所有唯一的吸附质
        unique_adsorbates = set()
        for reactant, product, _ in reactions:
            unique_adsorbates.add(reactant)
            unique_adsorbates.add(product)

        # 先按反应顺序排列吸附质
        adsorbates_ordered = []
        for reactant, product, _ in reactions:
            if reactant not in adsorbates_ordered:
                adsorbates_ordered.append(reactant)
            if product not in adsorbates_ordered:
                adsorbates_ordered.append(product)

        # 使用标准路径预测
        return self.predict_pathway(
            surface=surface,
            adsorbates=adsorbates_ordered,
            calculate_barriers=True,
            num_sites=num_sites,
            n_frames=n_frames,
            fmax=fmax,
            max_steps=max_steps,
            output_dir=output_dir,
        )

    def plot_energy_profile(
        self,
        result: CompletePathwayResult,
        output: Optional[str] = None,
        figsize: Tuple[float, float] = (12, 6),
        show_barriers: bool = True,
    ):
        """
        绘制能量剖面图（使用 EnergyDiagramPlotter）

        从每步的 activation_energy 和 reaction_energy 构建累积能量剖面，
        而非使用绝对能量（不同步骤原子数不同，绝对能量不可比）。

        Parameters
        ----------
        result : CompletePathwayResult
            路径分析结果
        output : str, optional
            输出文件路径
        figsize : tuple, default=(12, 6)
            图像大小
        show_barriers : bool, default=True
            是否显示能垒
        """
        from core.viz.energy_diagram_plotter import EnergyDiagramPlotter, DiagramStyle

        # Build cumulative energy profile from relative barriers
        energies = [0.0]
        labels = [result.adsorbates[0]]
        is_ts = [False]

        cumulative_energy = 0.0

        for step in result.steps:
            ea_fwd = step.activation_energy if step.activation_energy is not None else 0.0
            reaction_e = step.reaction_energy if step.reaction_energy is not None else 0.0

            # TS energy = current level + forward barrier
            ts_energy = cumulative_energy + ea_fwd
            energies.append(ts_energy)
            labels.append(f"TS{step.step_index + 1}")
            is_ts.append(True)

            # Product energy = current level + reaction energy
            cumulative_energy += reaction_e
            energies.append(cumulative_energy)
            labels.append(step.product_adsorbate)
            is_ts.append(False)

        style = DiagramStyle()
        style.figsize = (max(figsize[0], len(energies) * 0.8), figsize[1])
        style.color_scheme = "nature"
        style.ylabel_units = "eV"

        plotter = EnergyDiagramPlotter(style=style)
        plotter.plot(energies, labels, is_ts)

        if output:
            plotter.save(output, dpi=300)
            self.logger.info(f"Energy profile saved to {output}")
        else:
            import matplotlib.pyplot as plt
            plt.show()

    def __del__(self):
        """清理临时文件"""
        if not self.keep_files and self._temp_dir and os.path.exists(self._temp_dir):
            import shutil
            shutil.rmtree(self._temp_dir)


if __name__ == "__main__":
    print("Pathway Predictor module loaded successfully")
    print("See documentation for usage examples")
