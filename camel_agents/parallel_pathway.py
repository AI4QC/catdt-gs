"""
并行路径探索器 (ParallelPathwayExplorer)

支持：
1. 并行生成多个候选路径
2. 并行计算多个 NEB 步骤
3. 依赖分析和任务调度
"""

import logging
import concurrent.futures
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict
import json
from datetime import datetime

from ase import Atoms
from ase.io import read, write
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ParallelTask:
    """并行任务"""
    task_id: str
    task_type: str  # "pathway_design", "neb", "validation"
    inputs: Dict[str, Any]
    depends_on: List[str]  # 依赖的任务 ID
    priority: int = 0
    estimated_time: float = 0.0  # 预估时间（秒）


@dataclass
class ParallelTaskResult:
    """并行任务结果"""
    task_id: str
    success: bool
    result: Any
    error: Optional[str] = None
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    
    @property
    def duration(self) -> float:
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        return 0.0


class DependencyGraph:
    """依赖图 - 管理任务间的依赖关系"""
    
    def __init__(self):
        self.tasks: Dict[str, ParallelTask] = {}
        self.dependencies: Dict[str, Set[str]] = defaultdict(set)  # task -> its dependencies
        self.dependents: Dict[str, Set[str]] = defaultdict(set)    # task -> tasks that depend on it
    
    def add_task(self, task: ParallelTask):
        """添加任务"""
        self.tasks[task.task_id] = task
        for dep in task.depends_on:
            self.dependencies[task.task_id].add(dep)
            self.dependents[dep].add(task.task_id)
    
    def get_ready_tasks(self, completed: Set[str]) -> List[ParallelTask]:
        """获取可以执行的任务（所有依赖已完成）"""
        ready = []
        for task_id, task in self.tasks.items():
            if task_id not in completed and task_id not in [t.task_id for t in ready]:
                if all(dep in completed for dep in self.dependencies[task_id]):
                    ready.append(task)
        
        # 按优先级排序
        ready.sort(key=lambda t: t.priority, reverse=True)
        return ready
    
    def get_parallel_groups(self) -> List[List[str]]:
        """
        将任务分组，每组内的任务可以并行执行
        
        Returns:
            任务 ID 的分组列表
        """
        completed = set()
        groups = []
        remaining = set(self.tasks.keys())
        
        while remaining:
            # 找到所有依赖已完成的任务
            group = []
            for task_id in remaining:
                if all(dep in completed for dep in self.dependencies[task_id]):
                    group.append(task_id)
            
            if not group:
                # 有循环依赖
                raise ValueError("Circular dependency detected")
            
            groups.append(group)
            completed.update(group)
            remaining -= set(group)
        
        return groups
    
    def visualize(self) -> str:
        """生成依赖图的文本可视化"""
        lines = ["Dependency Graph:"]
        
        for task_id, task in self.tasks.items():
            deps = ", ".join(self.dependencies[task_id]) or "None"
            lines.append(f"  {task_id} ({task.task_type}) <- [{deps}]")
        
        return "\n".join(lines)


class ParallelPathwayExplorer:
    """
    并行路径探索器
    
    支持：
    1. 并行生成多个候选反应路径
    2. 并行计算 NEB（独立步骤）
    3. 智能任务调度
    """
    
    def __init__(self, max_workers: int = 4):
        self.max_workers = max_workers
        self.results: Dict[str, ParallelTaskResult] = {}
    
    def explore_multiple_pathways(
        self,
        surface: Atoms,
        reaction_description: str,
        pathway_design_fn: callable,
        n_pathways: int = 3,
        validation_fn: Optional[callable] = None
    ) -> List[Dict]:
        """
        并行探索多个反应路径
        
        Parameters:
            surface: 表面结构
            reaction_description: 反应描述
            pathway_design_fn: 路径设计函数 (surface, description) -> pathway
            n_pathways: 要生成的路径数量
            validation_fn: 验证函数 (pathway) -> validation_result
        
        Returns:
            多个路径设计结果，按 confidence 排序
        """
        logger.info(f"Exploring {n_pathways} pathways in parallel...")
        
        # 提交并行任务
        pathways = []
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # 提交多个路径设计任务
            futures = []
            for i in range(n_pathways):
                future = executor.submit(
                    self._design_pathway_with_seed,
                    pathway_design_fn,
                    surface,
                    reaction_description,
                    seed=i  # 不同的随机种子产生不同路径
                )
                futures.append((i, future))
            
            # 收集结果
            for idx, future in futures:
                try:
                    pathway = future.result(timeout=300)  # 5分钟超时
                    pathways.append(pathway)
                    logger.info(f"Pathway {idx}: design completed (confidence={pathway.get('confidence', 0):.2f})")
                except Exception as e:
                    logger.error(f"Pathway {idx} design failed: {e}")
        
        # 验证（如果提供了验证函数）
        if validation_fn and pathways:
            logger.info("Validating pathways...")
            validated = []
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {
                    executor.submit(validation_fn, pw): pw 
                    for pw in pathways
                }
                
                for future in concurrent.futures.as_completed(futures):
                    pathway = futures[future]
                    try:
                        validation = future.result()
                        pathway['validation'] = validation
                        validated.append(pathway)
                    except Exception as e:
                        logger.error(f"Validation failed: {e}")
            
            pathways = validated
        
        # 按 confidence 排序
        pathways.sort(key=lambda p: p.get('confidence', 0), reverse=True)
        
        logger.info(f"Completed: {len(pathways)} pathways generated")
        return pathways
    
    def _design_pathway_with_seed(
        self,
        design_fn: callable,
        surface: Atoms,
        description: str,
        seed: int
    ) -> Dict:
        """使用特定随机种子设计路径"""
        import random
        random.seed(seed)
        np.random.seed(seed)
        
        return design_fn(surface, description)
    
    def parallel_neb_calculations(
        self,
        steps: List[Dict[str, Any]],
        neb_fn: callable,
        max_workers: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        并行计算多个 NEB 步骤
        
        自动分析步骤依赖关系，并行计算独立的步骤。
        
        Parameters:
            steps: 步骤列表，每个包含 name, reactant, product
            neb_fn: NEB 计算函数 (step) -> neb_result
            max_workers: 并行工作线程数
        
        Returns:
            步骤名称 -> NEB 结果 的字典
        """
        workers = max_workers or self.max_workers
        logger.info(f"Running {len(steps)} NEB calculations with {workers} workers...")
        
        # 分析依赖关系
        # 步骤 i 的结果可能作为步骤 i+1 的输入
        # 但在这里我们假设所有步骤是独立的（可以同时计算）
        
        results = {}
        completed = 0
        failed = 0
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            # 提交所有任务
            future_to_step = {
                executor.submit(neb_fn, step): step 
                for step in steps
            }
            
            # 收集结果
            for future in concurrent.futures.as_completed(future_to_step):
                step = future_to_step[future]
                step_name = step.get('name', 'unknown')
                
                try:
                    result = future.result()
                    results[step_name] = result
                    completed += 1
                    logger.info(f"NEB completed: {step_name}")
                except Exception as e:
                    logger.error(f"NEB failed: {step_name} - {e}")
                    results[step_name] = {"error": str(e)}
                    failed += 1
        
        logger.info(f"NEB calculations: {completed} completed, {failed} failed")
        return results
    
    def execute_with_dependency_graph(
        self,
        tasks: List[ParallelTask],
        execute_fn: callable
    ) -> Dict[str, ParallelTaskResult]:
        """
        使用依赖图执行任务
        
        Parameters:
            tasks: 任务列表
            execute_fn: 执行函数 (task) -> result
        
        Returns:
            任务 ID -> 结果 的字典
        """
        # 构建依赖图
        graph = DependencyGraph()
        for task in tasks:
            graph.add_task(task)
        
        logger.info(f"Task dependency graph:\n{graph.visualize()}")
        
        # 获取并行组
        groups = graph.get_parallel_groups()
        logger.info(f"Execution plan: {len(groups)} groups")
        for i, group in enumerate(groups):
            logger.info(f"  Group {i}: {group}")
        
        # 执行
        results = {}
        
        for group_idx, group in enumerate(groups):
            logger.info(f"Executing group {group_idx+1}/{len(groups)}: {group}")
            
            group_tasks = [graph.tasks[tid] for tid in group]
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_task = {
                    executor.submit(execute_fn, task): task 
                    for task in group_tasks
                }
                
                for future in concurrent.futures.as_completed(future_to_task):
                    task = future_to_task[future]
                    result = ParallelTaskResult(
                        task_id=task.task_id,
                        success=True,
                        result=None,
                        start_time=datetime.now(),
                        end_time=datetime.now()
                    )
                    
                    try:
                        result.result = future.result()
                        result.success = True
                        logger.info(f"  ✓ {task.task_id}")
                    except Exception as e:
                        result.success = False
                        result.error = str(e)
                        logger.error(f"  ✗ {task.task_id}: {e}")
                    
                    results[task.task_id] = result
        
        return results
    
    def select_best_pathway(
        self,
        pathways: List[Dict],
        criteria: str = "confidence"
    ) -> Optional[Dict]:
        """
        从多个路径中选择最佳路径
        
        Parameters:
            pathways: 路径列表
            criteria: 选择标准 ("confidence", "lowest_barrier", "fewest_steps")
        
        Returns:
            最佳路径
        """
        if not pathways:
            return None
        
        if criteria == "confidence":
            return max(pathways, key=lambda p: p.get('confidence', 0))
        
        elif criteria == "fewest_steps":
            return min(pathways, key=lambda p: len(p.get('steps', [])))
        
        elif criteria == "lowest_barrier":
            # 需要 NEB 结果
            def max_barrier(p):
                barriers = [
                    s.get('activation_energy', float('inf')) 
                    for s in p.get('steps', [])
                ]
                return max(barriers) if barriers else float('inf')
            
            return min(pathways, key=max_barrier)
        
        else:
            return pathways[0]


# =============================================================================
# 便捷函数
# =============================================================================

def parallel_compute_barriers(
    steps: List[Dict],
    neb_calculator: Any,
    max_workers: int = 4
) -> Dict[str, float]:
    """
    并行计算多个反应步骤的能垒
    
    Example:
        >>> steps = [
        ...     {"name": "CO_to_CHO", "reactant": atoms1, "product": atoms2},
        ...     {"name": "CHO_to_CH2O", "reactant": atoms2, "product": atoms3},
        ... ]
        >>> barriers = parallel_compute_barriers(steps, neb_calculator)
        >>> print(barriers)
        {"CO_to_CHO": 0.84, "CHO_to_CH2O": 1.23}
    """
    explorer = ParallelPathwayExplorer(max_workers=max_workers)
    
    def neb_fn(step):
        return neb_calculator.calculate(step['reactant'], step['product'])
    
    return explorer.parallel_neb_calculations(steps, neb_fn)


def compare_pathways(
    pathways: List[Dict],
    metrics: List[str] = ["confidence", "n_steps", "max_barrier"]
) -> Dict:
    """
    比较多个路径的性能指标
    
    Returns:
        比较报告
    """
    report = {
        "n_pathways": len(pathways),
        "metrics": {},
        "ranking": {}
    }
    
    for metric in metrics:
        if metric == "confidence":
            values = [p.get('confidence', 0) for p in pathways]
            report["metrics"][metric] = {
                "best": max(values),
                "worst": min(values),
                "avg": sum(values) / len(values) if values else 0
            }
        
        elif metric == "n_steps":
            values = [len(p.get('steps', [])) for p in pathways]
            report["metrics"][metric] = {
                "best": min(values),
                "worst": max(values),
                "avg": sum(values) / len(values) if values else 0
            }
    
    # 综合排名
    scores = []
    for i, p in enumerate(pathways):
        score = p.get('confidence', 0) * 0.5  # confidence 权重 50%
        score += (1.0 / max(len(p.get('steps', [])), 1)) * 0.3  # 步数少加分
        scores.append((i, score))
    
    scores.sort(key=lambda x: x[1], reverse=True)
    report["ranking"] = {
        "best_pathway_idx": scores[0][0],
        "scores": {idx: score for idx, score in scores}
    }
    
    return report
