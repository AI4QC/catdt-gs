"""
增强型 NEB 结构验证与修复系统 (Enhanced NEB Validator)

整合并改进现有的验证功能，提供：
1. 基于 atoms.info 标签的通用验证（无需手动输入）
2. 智能结构修复（自动调整不合理结构）
3. 反应类型感知的验证规则
4. 详细的验证报告和修复建议
"""

import numpy as np
from ase import Atoms, Atom
from ase.geometry import get_distances
from scipy.spatial.distance import cdist
from typing import Tuple, Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from enum import Enum
import json
import logging
from pathlib import Path
import pickle

logger = logging.getLogger(__name__)


class ValidationStatus(Enum):
    """验证状态"""
    PASS = "pass"
    WARNING = "warning"
    FAIL = "fail"
    AUTO_FIXED = "auto_fixed"


class ReactionType(Enum):
    """反应类型"""
    HYDROGENATION = "hydrogenation"  # 加氢
    DEHYDROGENATION = "dehydrogenation"  # 脱氢
    DISSOCIATION = "dissociation"  # 解离
    COUPLING = "coupling"  # 偶联
    ISOMERIZATION = "isomerization"  # 异构化
    UNKNOWN = "unknown"


@dataclass
class GeometryIssue:
    """几何问题描述"""
    issue_type: str  # "too_far", "too_close", "unreasonable_bond", etc.
    atom_indices: Tuple[int, ...]  # 相关原子索引
    description: str
    severity: str  # "critical", "warning", "info"
    suggestion: str
    auto_fixable: bool = False


@dataclass
class ValidationReport:
    """验证报告"""
    status: ValidationStatus
    reaction_type: ReactionType
    initial_valid: bool
    final_valid: bool
    issues: List[GeometryIssue] = field(default_factory=list)
    fixes_applied: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        return {
            "status": self.status.value,
            "reaction_type": self.reaction_type.value,
            "initial_valid": self.initial_valid,
            "final_valid": self.final_valid,
            "issues": [
                {
                    "type": i.issue_type,
                    "atoms": i.atom_indices,
                    "description": i.description,
                    "severity": i.severity,
                    "suggestion": i.suggestion,
                    "auto_fixable": i.auto_fixable
                } for i in self.issues
            ],
            "fixes_applied": self.fixes_applied,
            "recommendations": self.recommendations,
            "metadata": self.metadata
        }
    
    def save(self, path: str):
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)


class EnhancedNEBValidator:
    """
    增强型 NEB 结构验证器
    
    特点：
    1. 完全基于 atoms.info 中的标签，无需手动输入
    2. 自动推断反应类型
    3. 智能修复常见问题
    4. 详细的错误报告和建议
    """
    
    # 标准共价键长（Å）
    STANDARD_BOND_LENGTHS = {
        ('C', 'H'): (1.0, 1.2),
        ('C', 'C'): (1.4, 1.6),
        ('C', 'O'): (1.2, 1.5),
        ('O', 'H'): (0.9, 1.1),
        ('C', 'N'): (1.3, 1.5),
        ('N', 'H'): (1.0, 1.1),
        ('H', 'H'): (0.7, 0.8),
    }
    
    # 原子范德华半径（用于检测重叠）
    VDW_RADIUS = {
        'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52,
        'F': 1.47, 'S': 1.8, 'P': 1.8, 'Pt': 1.75,
        'Cu': 1.4, 'Ni': 1.63, 'Fe': 1.56, 'Co': 1.52,
    }
    
    # 反应类型特征
    REACTION_SIGNATURES = {
        ReactionType.HYDROGENATION: {
            "description": "加氢反应（*X + H → *XH）",
            "atom_change": {"H": +1},
            "forming_bonds": [("C", "H"), ("O", "H"), ("N", "H")],
            "critical_distance": 2.0,  # H 应在 2Å 内
        },
        ReactionType.DEHYDROGENATION: {
            "description": "脱氢反应（*XH → *X + H）",
            "atom_change": {"H": -1},
            "breaking_bonds": [("C", "H"), ("O", "H"), ("N", "H")],
        },
        ReactionType.DISSOCIATION: {
            "description": "解离反应（*XY → *X + *Y）",
            "breaking_bonds": [("C", "O"), ("C", "C"), ("O", "O")],
        },
        ReactionType.COUPLING: {
            "description": "偶联反应（*X + *Y → *XY）",
            "atom_change": {},  # 原子数不变
            "forming_bonds": [("C", "C"), ("C", "O"), ("C", "N")],
        },
    }
    
    def __init__(self, auto_fix: bool = True, fix_attempts: int = 3):
        """
        初始化验证器
        
        Parameters
        ----------
        auto_fix : bool
            是否自动修复问题
        fix_attempts : int
            自动修复的最大尝试次数
        """
        self.auto_fix = auto_fix
        self.fix_attempts = fix_attempts
    
    def validate_endpoint_pair(
        self,
        initial: Atoms,
        final: Atoms,
        expected_reaction: Optional[str] = None
    ) -> ValidationReport:
        """
        验证一对 NEB 初末态结构
        
        Parameters
        ----------
        initial : Atoms
            初态结构（必须包含 atoms.info['adsorbate_indices']）
        final : Atoms
            末态结构（必须包含 atoms.info['adsorbate_indices']）
        expected_reaction : str, optional
            期望的反应类型（如 "hydrogenation"）
        
        Returns
        -------
        ValidationReport
            验证报告，包含问题和修复建议
        """
        issues = []
        fixes = []
        
        # 1. 检查必要的元数据
        if 'adsorbate_indices' not in initial.info:
            raise ValueError("Initial structure missing 'adsorbate_indices' in info")
        if 'adsorbate_indices' not in final.info:
            raise ValueError("Final structure missing 'adsorbate_indices' in info")
        
        # 2. 推断反应类型
        reaction_type = self._infer_reaction_type(initial, final, expected_reaction)
        logger.info(f"Detected reaction type: {reaction_type.value}")
        
        # 3. 验证初态结构
        initial_issues = self._validate_single_structure(initial, "initial")
        issues.extend(initial_issues)
        
        # 4. 验证末态结构
        final_issues = self._validate_single_structure(final, "final")
        issues.extend(final_issues)
        
        # 5. 验证初末态一致性
        consistency_issues = self._validate_consistency(initial, final, reaction_type)
        issues.extend(consistency_issues)
        
        # 6. 反应类型特定的验证
        specific_issues = self._validate_reaction_specific(initial, final, reaction_type)
        issues.extend(specific_issues)
        
        # 7. 自动修复（如果启用）
        if self.auto_fix and any(i.auto_fixable for i in issues):
            initial, final, applied_fixes = self._attempt_auto_fix(initial, final, issues)
            fixes.extend(applied_fixes)
            # 重新验证
            if applied_fixes:
                logger.info(f"Applied {len(applied_fixes)} auto-fixes, re-validating...")
                return self.validate_endpoint_pair(initial, final, expected_reaction)
        
        # 8. 确定最终状态
        critical_issues = [i for i in issues if i.severity == "critical"]
        warnings = [i for i in issues if i.severity == "warning"]
        
        if critical_issues:
            status = ValidationStatus.FAIL
        elif warnings:
            status = ValidationStatus.WARNING
        elif fixes:
            status = ValidationStatus.AUTO_FIXED
        else:
            status = ValidationStatus.PASS
        
        # 9. 生成建议
        recommendations = self._generate_recommendations(issues, reaction_type)
        
        return ValidationReport(
            status=status,
            reaction_type=reaction_type,
            initial_valid=len([i for i in initial_issues if i.severity == "critical"]) == 0,
            final_valid=len([i for i in final_issues if i.severity == "critical"]) == 0,
            issues=issues,
            fixes_applied=fixes,
            recommendations=recommendations,
            metadata={
                "initial_formula": initial.get_chemical_formula(),
                "final_formula": final.get_chemical_formula(),
                "initial_adsorbate_indices": initial.info['adsorbate_indices'],
                "final_adsorbate_indices": final.info['adsorbate_indices'],
            }
        )
    
    def _infer_reaction_type(
        self,
        initial: Atoms,
        final: Atoms,
        expected: Optional[str] = None
    ) -> ReactionType:
        """从结构变化推断反应类型"""
        if expected:
            try:
                return ReactionType(expected.lower())
            except ValueError:
                pass
        
        # 计算原子数变化
        initial_symbols = initial.get_chemical_symbols()
        final_symbols = final.get_chemical_symbols()
        
        from collections import Counter
        initial_counts = Counter(initial_symbols)
        final_counts = Counter(final_symbols)
        
        # 检测 H 原子变化
        h_change = final_counts.get('H', 0) - initial_counts.get('H', 0)
        
        if h_change == 1:
            return ReactionType.HYDROGENATION
        elif h_change == -1:
            return ReactionType.DEHYDROGENATION
        
        # 检测是否解离（吸附物原子数增加）
        initial_adsorbate = len(initial.info.get('adsorbate_indices', []))
        final_adsorbate = len(final.info.get('adsorbate_indices', []))
        
        if final_adsorbate > initial_adsorbate:
            return ReactionType.DISSOCIATION
        elif final_adsorbate < initial_adsorbate:
            return ReactionType.COUPLING
        
        return ReactionType.UNKNOWN
    
    def _validate_single_structure(
        self,
        atoms: Atoms,
        label: str
    ) -> List[GeometryIssue]:
        """验证单个结构"""
        issues = []
        
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        adsorbate_indices = atoms.info.get('adsorbate_indices', [])
        
        if not adsorbate_indices:
            issues.append(GeometryIssue(
                issue_type="no_adsorbate",
                atom_indices=(),
                description=f"{label}: No adsorbate atoms found",
                severity="critical",
                suggestion="Check adsorbate_indices in atoms.info",
                auto_fixable=False
            ))
            return issues
        
        # 1. 检查原子重叠
        for i, idx1 in enumerate(adsorbate_indices):
            for idx2 in adsorbate_indices[i+1:]:
                dist = np.linalg.norm(positions[idx1] - positions[idx2])
                elem1, elem2 = symbols[idx1], symbols[idx2]
                min_dist = self.VDW_RADIUS.get(elem1, 1.5) + self.VDW_RADIUS.get(elem2, 1.5)
                min_dist *= 0.5  # 允许一定程度的重叠
                
                if dist < min_dist:
                    issues.append(GeometryIssue(
                        issue_type="atomic_overlap",
                        atom_indices=(idx1, idx2),
                        description=f"{label}: {elem1}[{idx1}] and {elem2}[{idx2}] are too close ({dist:.2f} Å < {min_dist:.2f} Å)",
                        severity="critical",
                        suggestion="Atoms may be overlapping, need to re-relax structure",
                        auto_fixable=False
                    ))
        
        # 2. 检查吸附物-表面距离
        surface_indices = [i for i in range(len(atoms)) if i not in adsorbate_indices]
        if surface_indices:
            surface_z = max(positions[i, 2] for i in surface_indices)
            for idx in adsorbate_indices:
                height = positions[idx, 2] - surface_z
                if height > 5.0:
                    issues.append(GeometryIssue(
                        issue_type="adsorbate_too_high",
                        atom_indices=(idx,),
                        description=f"{label}: {symbols[idx]}[{idx}] is {height:.2f} Å above surface (too high)",
                        severity="warning",
                        suggestion="Adsorbate may be desorbed, check structure",
                        auto_fixable=False
                    ))
                elif height < 0.5:
                    issues.append(GeometryIssue(
                        issue_type="adsorbate_too_low",
                        atom_indices=(idx,),
                        description=f"{label}: {symbols[idx]}[{idx}] is {height:.2f} Å above surface (too low, may be embedded)",
                        severity="warning",
                        suggestion="Check if adsorbate is properly placed",
                        auto_fixable=False
                    ))
        
        # 3. 检查吸附物内部键长
        for i, idx1 in enumerate(adsorbate_indices):
            for idx2 in adsorbate_indices[i+1:]:
                dist = np.linalg.norm(positions[idx1] - positions[idx2])
                elem_pair = tuple(sorted([symbols[idx1], symbols[idx2]]))
                
                if elem_pair in self.STANDARD_BOND_LENGTHS:
                    min_len, max_len = self.STANDARD_BOND_LENGTHS[elem_pair]
                    if dist > max_len * 1.5:
                        issues.append(GeometryIssue(
                            issue_type="bond_too_long",
                            atom_indices=(idx1, idx2),
                            description=f"{label}: {elem_pair[0]}-{elem_pair[1]} distance is {dist:.2f} Å (expected {min_len}-{max_len} Å)",
                            severity="warning",
                            suggestion=f"Atoms may not be bonded, check structure",
                            auto_fixable=False
                        ))
        
        return issues
    
    def _validate_consistency(
        self,
        initial: Atoms,
        final: Atoms,
        reaction_type: ReactionType
    ) -> List[GeometryIssue]:
        """验证初末态一致性"""
        issues = []
        
        # 1. 检查表面原子数是否守恒
        initial_surface = len([i for i in range(len(initial)) 
                              if i not in initial.info.get('adsorbate_indices', [])])
        final_surface = len([i for i in range(len(final)) 
                            if i not in final.info.get('adsorbate_indices', [])])
        
        if initial_surface != final_surface:
            issues.append(GeometryIssue(
                issue_type="surface_atoms_changed",
                atom_indices=(),
                description=f"Surface atom count changed: {initial_surface} → {final_surface}",
                severity="critical",
                suggestion="Surface should not change during reaction, check indices",
                auto_fixable=False
            ))
        
        # 2. 检查是否有不合理的原子数变化
        initial_symbols = initial.get_chemical_symbols()
        final_symbols = final.get_chemical_symbols()
        
        from collections import Counter
        initial_counts = Counter(initial_symbols)
        final_counts = Counter(final_symbols)
        
        # 排除表面元素（通常是金属）
        surface_elem = max(initial_counts, key=initial_counts.get)
        
        for elem in set(initial_counts.keys()) | set(final_counts.keys()):
            if elem == surface_elem:
                continue
            change = final_counts.get(elem, 0) - initial_counts.get(elem, 0)
            if abs(change) > 2:
                issues.append(GeometryIssue(
                    issue_type="excessive_atom_change",
                    atom_indices=(),
                    description=f"{elem} atom count changed by {change:+d}",
                    severity="warning",
                    suggestion="Large atom changes may indicate incorrect pathway",
                    auto_fixable=False
                ))
        
        return issues
    
    def _validate_reaction_specific(
        self,
        initial: Atoms,
        final: Atoms,
        reaction_type: ReactionType
    ) -> List[GeometryIssue]:
        """反应类型特定的验证"""
        issues = []
        
        signature = self.REACTION_SIGNATURES.get(reaction_type)
        if not signature:
            return issues
        
        # 加氢反应特定检查
        if reaction_type == ReactionType.HYDROGENATION:
            issues.extend(self._check_hydrogenation_geometry(initial, final, signature))
        
        # 脱氢反应特定检查
        elif reaction_type == ReactionType.DEHYDROGENATION:
            issues.extend(self._check_dehydrogenation_geometry(initial, final, signature))
        
        return issues
    
    def _check_hydrogenation_geometry(
        self,
        initial: Atoms,
        final: Atoms,
        signature: Dict
    ) -> List[GeometryIssue]:
        """检查加氢反应的几何结构"""
        issues = []
        
        # 在末态中，H 应该靠近 C/O/N
        final_pos = final.get_positions()
        final_symbols = final.get_chemical_symbols()
        final_adsorbate = final.info.get('adsorbate_indices', [])
        
        # 找到新增的 H
        initial_counts = {}
        for s in initial.get_chemical_symbols():
            initial_counts[s] = initial_counts.get(s, 0) + 1
        
        final_h_atoms = [i for i in final_adsorbate if final_symbols[i] == 'H']
        
        for h_idx in final_h_atoms:
            # 找到最近的 C/O/N
            min_dist = float('inf')
            min_idx = -1
            min_elem = ""
            
            for idx in final_adsorbate:
                if idx == h_idx:
                    continue
                if final_symbols[idx] in ['C', 'O', 'N']:
                    dist = np.linalg.norm(final_pos[h_idx] - final_pos[idx])
                    if dist < min_dist:
                        min_dist = dist
                        min_idx = idx
                        min_elem = final_symbols[idx]
            
            if min_dist > signature.get("critical_distance", 2.0):
                issues.append(GeometryIssue(
                    issue_type="h_too_far_from_target",
                    atom_indices=(h_idx, min_idx),
                    description=f"Final: H[{h_idx}] is {min_dist:.2f} Å from {min_elem}[{min_idx}] (should be < 2.0 Å for hydrogenation)",
                    severity="critical",
                    suggestion="Hydrogen is too far from adsorbate. In initial state, place H within 2 Å of target atom.",
                    auto_fixable=True
                ))
        
        return issues
    
    def _check_dehydrogenation_geometry(
        self,
        initial: Atoms,
        final: Atoms,
        signature: Dict
    ) -> List[GeometryIssue]:
        """检查脱氢反应的几何结构"""
        issues = []
        
        # 在初态中，应该有一个合理的 C-H/O-H/N-H 键
        initial_pos = initial.get_positions()
        initial_symbols = initial.get_chemical_symbols()
        initial_adsorbate = initial.info.get('adsorbate_indices', [])
        
        h_atoms = [i for i in initial_adsorbate if initial_symbols[i] == 'H']
        
        for h_idx in h_atoms:
            # 找到最近的 C/O/N
            min_dist = float('inf')
            
            for idx in initial_adsorbate:
                if idx == h_idx:
                    continue
                if initial_symbols[idx] in ['C', 'O', 'N']:
                    dist = np.linalg.norm(initial_pos[h_idx] - initial_pos[idx])
                    min_dist = min(min_dist, dist)
            
            if min_dist > 1.5:  # C-H 键通常 ~1.1 Å
                issues.append(GeometryIssue(
                    issue_type="h_not_bonded",
                    atom_indices=(h_idx,),
                    description=f"Initial: H[{h_idx}] is {min_dist:.2f} Å from nearest C/O/N (should be bonded for dehydrogenation)",
                    severity="warning",
                    suggestion="Hydrogen should be bonded to the adsorbate for dehydrogenation",
                    auto_fixable=False
                ))
        
        return issues
    
    def _attempt_auto_fix(
        self,
        initial: Atoms,
        final: Atoms,
        issues: List[GeometryIssue]
    ) -> Tuple[Atoms, Atoms, List[str]]:
        """尝试自动修复问题"""
        applied_fixes = []
        
        for issue in issues:
            if not issue.auto_fixable:
                continue
            
            if issue.issue_type == "h_too_far_from_target":
                # 修复：将 H 原子移动到目标原子附近
                fixed_initial, fix_desc = self._fix_h_position(initial, issue)
                if fix_desc:
                    initial = fixed_initial
                    applied_fixes.append(fix_desc)
        
        return initial, final, applied_fixes
    
    def _fix_h_position(self, atoms: Atoms, issue: GeometryIssue) -> Tuple[Atoms, str]:
        """修复 H 原子位置"""
        h_idx, target_idx = issue.atom_indices
        
        atoms_fixed = atoms.copy()
        positions = atoms_fixed.get_positions()
        
        target_pos = positions[target_idx]
        
        # 将 H 放在目标原子附近，稍微偏上方
        new_pos = target_pos + np.array([0.5, 0.5, 1.5])  # 1.5 Å above and offset
        
        positions[h_idx] = new_pos
        atoms_fixed.set_positions(positions)
        
        return atoms_fixed, f"Moved H[{h_idx}] to {new_pos} (near target atom)"
    
    def _generate_recommendations(
        self,
        issues: List[GeometryIssue],
        reaction_type: ReactionType
    ) -> List[str]:
        """生成修复建议"""
        recommendations = []
        
        critical = [i for i in issues if i.severity == "critical"]
        warnings = [i for i in issues if i.severity == "warning"]
        
        if critical:
            recommendations.append(f"Critical issues found ({len(critical)}): Must fix before NEB")
            for issue in critical[:3]:  # 只显示前3个
                recommendations.append(f"  - {issue.suggestion}")
        
        if warnings:
            recommendations.append(f"Warnings ({len(warnings)}): Review recommended")
        
        if not critical and not warnings:
            recommendations.append("Structure validation passed. Ready for NEB calculation.")
        
        # 反应类型特定建议
        if reaction_type == ReactionType.HYDROGENATION:
            recommendations.append("Tip: For hydrogenation, ensure H is within 2 Å of the target atom in initial state")
        
        return recommendations


# =============================================================================
# 便捷函数
# =============================================================================

def validate_neb_endpoints(
    initial: Atoms,
    final: Atoms,
    expected_reaction: Optional[str] = None,
    auto_fix: bool = True
) -> ValidationReport:
    """
    便捷函数：验证 NEB 初末态
    
    Example:
        >>> from ase.io import read
        >>> initial = read("initial.vasp")
        >>> final = read("final.vasp")
        >>> report = validate_neb_endpoints(initial, final, "hydrogenation")
        >>> print(report.status)
        >>> if report.status != ValidationStatus.PASS:
        ...     print(report.recommendations)
    """
    validator = EnhancedNEBValidator(auto_fix=auto_fix)
    return validator.validate_endpoint_pair(initial, final, expected_reaction)
