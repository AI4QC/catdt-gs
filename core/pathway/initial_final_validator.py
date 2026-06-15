"""
初末态结构验证器 (Agent 4.5)

自动验证和修复NEB计算的初始态和最终态结构，确保：
1. 键长合理
2. 无原子重叠
3. 配位数正确
4. 吸附位点合理
5. 氢供体位置合适

支持任意反应、任意体系、任意吸附分子的通用性验证。
"""

import numpy as np
from ase import Atoms
from ase.geometry import get_distances
from scipy.spatial.distance import cdist
from typing import Tuple, Dict, List, Optional
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ReactionRuleDatabase:
    """
    反应类型规则数据库
    """
    
    # 标准共价键长（Å）
    STANDARD_BOND_LENGTHS = {
        'C-H': 1.09,
        'C-C': 1.54,
        'C-O': 1.43,
        'C=O': 1.20,
        'O-H': 0.96,
        'C-N': 1.47,
        'C-surface': 2.0,  # 典型C吸附在金属表面
        'H-surface': 1.8,  # 典型H吸附在金属表面
        'O-surface': 2.0,  # 典型O吸附在金属表面
    }
    
    # 允许的键长范围（倍数）
    BOND_LENGTH_TOLERANCE = 0.3  # ±30%
    
    # 反应类型模板
    REACTION_TEMPLATES = {
        "hydrogenation": {
            "pattern": "*CHx + H → *CHx+1",
            "bonds_forming": ["C-H"],
            "critical_distances": {
                "C-H_forming": (1.0, 2.5),  # 形成中的C-H键
                "H-surface": (1.5, 2.2),    # H供体到表面距离
                "C-surface": (1.8, 2.5),    # C到表面距离
            },
            "coordination_changes": {
                "C": {"initial": (3, 4), "final": (4, 4)},  # sp2 → sp3
            },
            "requires_h_donor": True,
        },
        
        "dehydrogenation": {
            "pattern": "*CHx → *CHx-1 + H",
            "bonds_breaking": ["C-H"],
            "critical_distances": {
                "C-H_breaking": (1.0, 1.3),  # 初态的C-H键
                "C-surface": (1.8, 2.5),
            },
            "coordination_changes": {
                "C": {"initial": (4, 4), "final": (3, 4)},  # sp3 → sp2
            },
        },
        
        "c_o_coupling": {
            "pattern": "*CO + *O → *CO2",
            "bonds_forming": ["C-O"],
            "critical_distances": {
                "C-O_forming": (1.2, 2.8),  # 形成中的C-O键
                "C-surface": (1.8, 2.5),
                "O-surface": (1.8, 2.3),
            },
        },
        
        "generic": {
            "min_interatomic_distance": 0.8,  # 任意两原子最小距离
            "max_bond_length_factor": 1.5,   # 共价键最大长度倍数
        }
    }
    
    @classmethod
    def get_rule(cls, reaction_type: str) -> Dict:
        """获取反应类型的规则"""
        return cls.REACTION_TEMPLATES.get(reaction_type, cls.REACTION_TEMPLATES["generic"])
    
    @classmethod
    def infer_reaction_type(cls, initial_formula: str, final_formula: str) -> str:
        """
        从化学式推断反应类型
        
        Examples:
            "*CH2", "*CH3" → "hydrogenation"
            "*CO + *O", "*CO2" → "c_o_coupling"
        """
        # 简单的启发式推断
        initial_h = initial_formula.count('H')
        final_h = final_formula.count('H')
        
        if final_h > initial_h:
            return "hydrogenation"
        elif final_h < initial_h:
            return "dehydrogenation"
        elif "*CO" in initial_formula and "*O" in initial_formula and "*CO2" in final_formula:
            return "c_o_coupling"
        else:
            return "generic"


class StructureValidator:
    """
    结构合理性验证器
    """
    
    def __init__(self, reaction_type: str = "generic"):
        self.reaction_type = reaction_type
        self.rules = ReactionRuleDatabase.get_rule(reaction_type)
        
    def check_bond_distances(
        self, 
        atoms: Atoms, 
        bond_pairs: List[Tuple[int, int]], 
        bond_type: str
    ) -> Tuple[bool, str]:
        """
        检查键长是否在合理范围
        
        Args:
            atoms: ASE Atoms对象
            bond_pairs: 键的原子索引对列表 [(i, j), ...]
            bond_type: 键类型（如 "C-H_forming"）
        
        Returns:
            (is_valid, message)
        """
        if "critical_distances" not in self.rules:
            return True, "N/A"
        
        if bond_type not in self.rules["critical_distances"]:
            return True, f"No rule for {bond_type}"
        
        min_d, max_d = self.rules["critical_distances"][bond_type]
        
        for i, j in bond_pairs:
            distance = atoms.get_distance(i, j)
            if not (min_d <= distance <= max_d):
                return False, f"{bond_type} distance {distance:.2f} Å out of range [{min_d}, {max_d}]"
        
        return True, f"{bond_type} distances OK"
    
    def check_atom_overlap(self, atoms: Atoms, min_distance: float = 0.8) -> Tuple[bool, str]:
        """
        检查原子是否过于接近（重叠）
        
        Args:
            atoms: ASE Atoms对象
            min_distance: 最小允许距离（Å）
        
        Returns:
            (is_valid, message)
        """
        positions = atoms.get_positions()
        
        # 计算所有原子对的距离
        distances = cdist(positions, positions)
        np.fill_diagonal(distances, np.inf)  # 忽略自身
        
        min_dist = distances.min()
        overlap_indices = np.where(distances == min_dist)
        
        if min_dist < min_distance:
            i, j = overlap_indices[0][0], overlap_indices[1][0]
            return False, f"Atoms {i}({atoms[i].symbol}) and {j}({atoms[j].symbol}) overlap: {min_dist:.2f} Å < {min_distance} Å"
        
        return True, f"No overlap (min distance: {min_dist:.2f} Å)"
    
    def check_coordination(
        self, 
        atoms: Atoms, 
        atom_indices: List[int], 
        expected_range: Tuple[int, int],
        cutoff: float = 2.0
    ) -> Tuple[bool, str]:
        """
        检查原子的配位数
        
        Args:
            atoms: ASE Atoms对象
            atom_indices: 要检查的原子索引列表
            expected_range: 期望的配位数范围 (min, max)
            cutoff: 配位数统计的距离阈值（Å）
        
        Returns:
            (is_valid, message)
        """
        min_coord, max_coord = expected_range
        
        for idx in atom_indices:
            element = atoms[idx].symbol
            
            # 计算近邻数
            positions = atoms.get_positions()
            distances = np.linalg.norm(positions - positions[idx], axis=1)
            neighbors = np.sum((distances > 0.1) & (distances < cutoff))
            
            if not (min_coord <= neighbors <= max_coord):
                return False, f"{element} atom {idx} has {neighbors} neighbors (expected {min_coord}-{max_coord})"
        
        return True, "Coordination numbers OK"
    
    def check_hydrogen_donor_distance(
        self, 
        atoms: Atoms, 
        h_donor_idx: int, 
        c_acceptor_idx: int
    ) -> Tuple[bool, str]:
        """
        检查氢供体到碳受体的距离（加氢反应专用）
        
        Args:
            atoms: ASE Atoms对象
            h_donor_idx: 氢供体原子索引
            c_acceptor_idx: 碳受体原子索引
        
        Returns:
            (is_valid, message)
        """
        if not self.rules.get("requires_h_donor", False):
            return True, "N/A (no H transfer)"
        
        distance = atoms.get_distance(h_donor_idx, c_acceptor_idx)
        
        if distance > 3.5:
            return False, f"H donor too far from C acceptor: {distance:.2f} Å > 3.5 Å"
        elif distance < 1.2:
            return False, f"H donor too close to C acceptor: {distance:.2f} Å < 1.2 Å"
        
        return True, f"H donor distance OK: {distance:.2f} Å"
    
    def check_surface_proximity(
        self, 
        atoms: Atoms, 
        adsorbate_indices: List[int],
        surface_indices: List[int],
        min_distance: float = 1.5,
        max_distance: float = 3.5
    ) -> Tuple[bool, str]:
        """
        检查吸附物到表面的距离
        
        Args:
            atoms: ASE Atoms对象
            adsorbate_indices: 吸附物原子索引
            surface_indices: 表面原子索引
            min_distance: 最小距离（太近可能重叠）
            max_distance: 最大距离（太远无法吸附）
        
        Returns:
            (is_valid, message)
        """
        adsorbate_pos = atoms.get_positions()[adsorbate_indices]
        surface_pos = atoms.get_positions()[surface_indices]
        
        # 计算吸附物到表面的最小距离
        distances = cdist(adsorbate_pos, surface_pos)
        min_dist = distances.min()
        
        if min_dist < min_distance:
            return False, f"Adsorbate too close to surface: {min_dist:.2f} Å < {min_distance} Å"
        elif min_dist > max_distance:
            return False, f"Adsorbate too far from surface: {min_dist:.2f} Å > {max_distance} Å"
        
        return True, f"Surface proximity OK: {min_dist:.2f} Å"


class StructureFixer:
    """
    结构自动修复工具
    """
    
    def __init__(self):
        pass
    
    def fix_bond_distance(
        self, 
        atoms: Atoms, 
        atom_i: int, 
        atom_j: int, 
        target_distance: float,
        freeze_indices: List[int] = []
    ) -> Atoms:
        """
        调整两个原子间的距离到目标值
        
        Args:
            atoms: ASE Atoms对象
            atom_i, atom_j: 两个原子的索引
            target_distance: 目标距离（Å）
            freeze_indices: 冻结的原子索引（通常是表面原子）
        
        Returns:
            修复后的Atoms对象
        """
        atoms = atoms.copy()
        
        # 确定哪个原子可移动
        if atom_i in freeze_indices and atom_j in freeze_indices:
            logger.warning(f"Both atoms {atom_i} and {atom_j} are frozen, cannot fix distance")
            return atoms
        
        if atom_i in freeze_indices:
            movable_idx, anchor_idx = atom_j, atom_i
        elif atom_j in freeze_indices:
            movable_idx, anchor_idx = atom_i, atom_j
        else:
            # 都可移动，移动第二个
            movable_idx, anchor_idx = atom_j, atom_i
        
        # 计算方向向量
        direction = atoms[movable_idx].position - atoms[anchor_idx].position
        current_distance = np.linalg.norm(direction)
        direction /= current_distance
        
        # 设置新位置
        new_position = atoms[anchor_idx].position + direction * target_distance
        atoms[movable_idx].position = new_position
        
        logger.info(f"Fixed distance {atom_i}-{atom_j}: {current_distance:.2f} → {target_distance:.2f} Å")
        
        return atoms
    
    def separate_overlapping_atoms(
        self, 
        atoms: Atoms, 
        atom_i: int, 
        atom_j: int,
        min_distance: float = 1.0,
        freeze_indices: List[int] = []
    ) -> Atoms:
        """
        分离重叠的原子
        
        Args:
            atoms: ASE Atoms对象
            atom_i, atom_j: 重叠的两个原子
            min_distance: 最小安全距离
            freeze_indices: 冻结的原子索引
        
        Returns:
            修复后的Atoms对象
        """
        return self.fix_bond_distance(atoms, atom_i, atom_j, min_distance, freeze_indices)
    
    def adjust_adsorbate_height(
        self, 
        atoms: Atoms, 
        adsorbate_indices: List[int],
        target_height: float = 2.0
    ) -> Atoms:
        """
        调整吸附物到表面的高度（Z方向）
        
        Args:
            atoms: ASE Atoms对象
            adsorbate_indices: 吸附物原子索引
            target_height: 目标高度（Å，相对于表面顶层）
        
        Returns:
            修复后的Atoms对象
        """
        atoms = atoms.copy()
        
        # 找到表面顶层Z坐标
        all_z = atoms.get_positions()[:, 2]
        adsorbate_z = all_z[adsorbate_indices]
        surface_z = np.delete(all_z, adsorbate_indices)
        surface_top = surface_z.max()
        
        # 计算吸附物质心当前高度
        current_height = adsorbate_z.mean() - surface_top
        
        # 整体移动吸附物
        shift = target_height - current_height
        positions = atoms.get_positions()
        positions[adsorbate_indices, 2] += shift
        atoms.set_positions(positions)
        
        logger.info(f"Adjusted adsorbate height: {current_height:.2f} → {target_height:.2f} Å")
        
        return atoms


class InitialFinalStateValidator:
    """
    Agent 4.5: 初末态结构验证与修复
    
    工作流程:
    1. 接收反应初末态结构
    2. 运行一系列验证检查
    3. 对于失败的检查，尝试自动修复
    4. 生成验证报告
    5. 返回修复后的结构或标记为需要人工审查
    """
    
    def __init__(self, reaction_type: str = "generic", auto_fix: bool = True):
        """
        Args:
            reaction_type: 反应类型（hydrogenation, dehydrogenation, c_o_coupling, generic）
            auto_fix: 是否自动修复问题
        """
        self.reaction_type = reaction_type
        self.auto_fix = auto_fix
        self.validator = StructureValidator(reaction_type)
        self.fixer = StructureFixer()
    
    def validate_and_fix(
        self, 
        initial_state: Atoms, 
        final_state: Atoms,
        adsorbate_indices: List[int],
        surface_indices: List[int],
        reaction_info: Optional[Dict] = None
    ) -> Tuple[Atoms, Atoms, Dict]:
        """
        验证并修复初末态结构
        
        Args:
            initial_state: 初始态结构
            final_state: 最终态结构
            adsorbate_indices: 吸附物原子索引
            surface_indices: 表面原子索引（冻结）
            reaction_info: 反应的额外信息（如氢供体位置等）
        
        Returns:
            (fixed_initial, fixed_final, validation_report)
        """
        logger.info(f"Validating initial and final states for {self.reaction_type} reaction...")
        
        # 初始化报告
        report = {
            "reaction_type": self.reaction_type,
            "initial_state": {"errors": [], "warnings": [], "passed": []},
            "final_state": {"errors": [], "warnings": [], "passed": []},
            "fixes_applied": [],
            "status": "UNKNOWN",
        }
        
        # 验证初始态
        initial_state, initial_report = self._validate_state(
            initial_state, 
            adsorbate_indices, 
            surface_indices,
            reaction_info,
            state_name="initial"
        )
        report["initial_state"] = initial_report
        
        # 验证最终态
        final_state, final_report = self._validate_state(
            final_state, 
            adsorbate_indices, 
            surface_indices,
            reaction_info,
            state_name="final"
        )
        report["final_state"] = final_report
        
        # 判断整体状态
        has_errors = (len(initial_report["errors"]) > 0) or (len(final_report["errors"]) > 0)
        
        if has_errors:
            report["status"] = "FAIL - Critical errors remain"
            report["recommendation"] = "Manual review required. Consider alternative pathway or initial placement."
        else:
            report["status"] = "PASS"
            report["recommendation"] = "Structures are ready for NEB calculation."
        
        # 打印摘要
        self._print_report(report)
        
        return initial_state, final_state, report
    
    def _validate_state(
        self, 
        atoms: Atoms, 
        adsorbate_indices: List[int],
        surface_indices: List[int],
        reaction_info: Optional[Dict],
        state_name: str
    ) -> Tuple[Atoms, Dict]:
        """
        验证单个状态（初态或末态）
        
        Returns:
            (fixed_atoms, report_dict)
        """
        report = {
            "errors": [],
            "warnings": [],
            "passed": [],
        }
        
        # 检查1: 原子重叠
        passed, msg = self.validator.check_atom_overlap(atoms)
        if not passed:
            report["errors"].append(f"Atom overlap: {msg}")
            if self.auto_fix:
                # 提取重叠的原子索引
                # 这里简化处理，实际需要解析msg
                logger.warning(f"{state_name}: {msg}, attempting to fix...")
                # TODO: 实现智能分离
        else:
            report["passed"].append(msg)
        
        # 检查2: 表面距离
        passed, msg = self.validator.check_surface_proximity(
            atoms, adsorbate_indices, surface_indices
        )
        if not passed:
            report["errors"].append(f"Surface proximity: {msg}")
            if self.auto_fix:
                logger.info(f"{state_name}: Adjusting adsorbate height...")
                atoms = self.fixer.adjust_adsorbate_height(atoms, adsorbate_indices)
                # 重新检查
                passed, msg = self.validator.check_surface_proximity(
                    atoms, adsorbate_indices, surface_indices
                )
                if passed:
                    report["warnings"].append("Fixed: " + msg)
                    report["errors"] = [e for e in report["errors"] if "Surface proximity" not in e]
        else:
            report["passed"].append(msg)
        
        # 检查3: 氢供体距离（如果是加氢反应）
        if reaction_info and "h_donor_idx" in reaction_info and "c_acceptor_idx" in reaction_info:
            passed, msg = self.validator.check_hydrogen_donor_distance(
                atoms, 
                reaction_info["h_donor_idx"], 
                reaction_info["c_acceptor_idx"]
            )
            if not passed:
                report["errors"].append(f"H donor distance: {msg}")
                if self.auto_fix:
                    logger.info(f"{state_name}: Adjusting H-C distance...")
                    target_distance = 2.0  # 合理的H-C距离
                    atoms = self.fixer.fix_bond_distance(
                        atoms, 
                        reaction_info["h_donor_idx"],
                        reaction_info["c_acceptor_idx"],
                        target_distance,
                        freeze_indices=surface_indices
                    )
                    # 重新检查
                    passed, msg = self.validator.check_hydrogen_donor_distance(
                        atoms, 
                        reaction_info["h_donor_idx"],
                        reaction_info["c_acceptor_idx"]
                    )
                    if passed:
                        report["warnings"].append("Fixed: " + msg)
                        report["errors"] = [e for e in report["errors"] if "H donor" not in e]
            else:
                report["passed"].append(msg)
        
        # 检查4: 配位数（根据反应类型）
        if "coordination_changes" in self.validator.rules:
            coord_rules = self.validator.rules["coordination_changes"].get("C", None)
            if coord_rules:
                # 找到碳原子索引（简化：假设第一个吸附原子是C）
                c_indices = [i for i in adsorbate_indices if atoms[i].symbol == 'C']
                if c_indices:
                    if state_name == "initial":
                        expected = coord_rules.get("initial", (3, 4))
                    else:
                        expected = coord_rules.get("final", (3, 4))
                    
                    passed, msg = self.validator.check_coordination(atoms, c_indices, expected)
                    if not passed:
                        report["warnings"].append(f"Coordination: {msg}")
                    else:
                        report["passed"].append(msg)
        
        return atoms, report
    
    def _print_report(self, report: Dict):
        """打印验证报告"""
        logger.info("\n" + "="*60)
        logger.info(f"VALIDATION REPORT - {report['reaction_type'].upper()} REACTION")
        logger.info("="*60)
        
        for state in ["initial_state", "final_state"]:
            state_name = state.replace("_state", "").upper()
            logger.info(f"\n{state_name} STATE:")
            
            if report[state]["errors"]:
                logger.error(f"  ❌ Errors ({len(report[state]['errors'])}):")
                for err in report[state]["errors"]:
                    logger.error(f"     - {err}")
            
            if report[state]["warnings"]:
                logger.warning(f"  ⚠️  Warnings ({len(report[state]['warnings'])}):")
                for warn in report[state]["warnings"]:
                    logger.warning(f"     - {warn}")
            
            if report[state]["passed"]:
                logger.info(f"  ✅ Passed ({len(report[state]['passed'])}):")
                for p in report[state]["passed"][:3]:  # 只显示前3个
                    logger.info(f"     - {p}")
        
        logger.info(f"\nFINAL STATUS: {report['status']}")
        logger.info(f"RECOMMENDATION: {report['recommendation']}")
        logger.info("="*60 + "\n")


def quick_test():
    """快速测试"""
    from ase.build import fcc111, add_adsorbate
    from ase import Atom
    
    # 创建测试表面
    slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
    
    # 添加*CH2
    add_adsorbate(slab, 'C', height=2.0, position='fcc')
    c_idx = len(slab) - 1
    
    # 添加2个H（模拟CH2）
    c_pos = slab[c_idx].position
    slab.append(Atom('H', position=c_pos + [0.5, 0.5, 1.0]))
    slab.append(Atom('H', position=c_pos + [-0.5, 0.5, 1.0]))
    
    initial_state = slab.copy()
    
    # 创建最终态（添加第3个H，模拟*CH3）
    # 故意让H原子距离很远，测试验证器
    final_state = slab.copy()
    final_state.append(Atom('H', position=c_pos + [0, 0, 4.5]))  # 太远！
    
    # 定义吸附物和表面索引
    n_surface = len(fcc111('Pt', size=(3, 3, 4), vacuum=10.0))
    adsorbate_indices = list(range(n_surface, len(initial_state)))
    surface_indices = list(range(n_surface))
    
    # 运行验证
    validator = InitialFinalStateValidator(reaction_type="hydrogenation", auto_fix=True)
    
    fixed_initial, fixed_final, report = validator.validate_and_fix(
        initial_state,
        final_state,
        adsorbate_indices,
        surface_indices,
        reaction_info={
            "h_donor_idx": len(final_state) - 1,  # 最后一个H
            "c_acceptor_idx": c_idx,
        }
    )
    
    print("\n✅ Test completed!")
    print(f"Initial state: {len(fixed_initial)} atoms")
    print(f"Final state: {len(fixed_final)} atoms")
    print(f"Validation status: {report['status']}")
    
    return fixed_initial, fixed_final, report


if __name__ == "__main__":
    quick_test()
