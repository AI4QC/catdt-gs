"""
自适应参数调整器 (AdaptiveParameterTuner)

根据中间结果自动调整计算参数，包括：
1. MC 温度调整
2. NEB 参数优化
3. 收敛判断和参数重试
"""

import logging
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass
from enum import Enum
import numpy as np

logger = logging.getLogger(__name__)


class ConvergenceStatus(Enum):
    """收敛状态"""
    CONVERGED = "converged"
    IMPROVING = "improving"
    OSCILLATING = "oscillating"
    STUCK = "stuck"
    DIVERGING = "diverging"
    UNKNOWN = "unknown"


@dataclass
class ParameterAdjustment:
    """参数调整建议"""
    parameter: str
    current_value: Any
    suggested_value: Any
    reason: str
    confidence: float  # 0-1


class MCTemperatureTuner:
    """
    MC 温度自适应调整器
    
    根据能量轨迹自动调整温度，确保：
    - 足够的探索（能量变化明显）
    - 不过度震荡
    - 合理的接受率（20-60%）
    """
    
    # 目标接受率范围
    TARGET_ACCEPTANCE_MIN = 0.20
    TARGET_ACCEPTANCE_MAX = 0.60
    
    # 温度调整因子
    TEMP_ADJUST_UP = 1.5
    TEMP_ADJUST_DOWN = 0.7
    
    def __init__(
        self,
        initial_temperature: float = 500.0,
        min_temperature: float = 100.0,
        max_temperature: float = 2000.0,
        adjustment_interval: int = 100  # 每多少步调整一次
    ):
        self.temperature = initial_temperature
        self.min_temp = min_temperature
        self.max_temp = max_temperature
        self.adjustment_interval = adjustment_interval
        
        self.energy_history: List[float] = []
        self.acceptance_history: List[bool] = []
        self.temperature_history: List[float] = [initial_temperature]
    
    def update(self, energy: float, accepted: bool, step: int):
        """更新状态"""
        self.energy_history.append(energy)
        self.acceptance_history.append(accepted)
        
        # 按间隔调整温度
        if step > 0 and step % self.adjustment_interval == 0:
            self._adjust_temperature(step)
    
    def _adjust_temperature(self, current_step: int):
        """调整温度"""
        window_size = min(self.adjustment_interval, len(self.acceptance_history))
        if window_size == 0:
            return
        
        # 计算最近窗口的接受率
        recent_acceptance = sum(self.acceptance_history[-window_size:]) / window_size
        
        # 计算能量变化
        recent_energies = self.energy_history[-window_size:]
        energy_std = np.std(recent_energies) if len(recent_energies) > 1 else 0
        
        old_temp = self.temperature
        adjustment_reason = ""
        
        # 策略 1: 接受率过低 -> 升温
        if recent_acceptance < self.TARGET_ACCEPTANCE_MIN:
            self.temperature = min(
                self.temperature * self.TEMP_ADJUST_UP,
                self.max_temp
            )
            adjustment_reason = f"acceptance rate {recent_acceptance:.2%} too low"
        
        # 策略 2: 接受率过高 -> 降温
        elif recent_acceptance > self.TARGET_ACCEPTANCE_MAX:
            self.temperature = max(
                self.temperature * self.TEMP_ADJUST_DOWN,
                self.min_temp
            )
            adjustment_reason = f"acceptance rate {recent_acceptance:.2%} too high"
        
        # 策略 3: 能量无变化（卡住）-> 大幅升温
        elif energy_std < 0.001:  # 能量几乎不变
            self.temperature = min(
                self.temperature * 2.0,
                self.max_temp
            )
            adjustment_reason = f"energy stuck (std={energy_std:.6f})"
        
        # 策略 4: 能量震荡过大 -> 小幅降温
        elif energy_std > 1.0:  # 能量变化太大
            self.temperature = max(
                self.temperature * 0.9,
                self.min_temp
            )
            adjustment_reason = f"energy oscillating (std={energy_std:.3f})"
        
        if self.temperature != old_temp:
            self.temperature_history.append(self.temperature)
            logger.info(
                f"MC temperature adjusted at step {current_step}: "
                f"{old_temp:.1f}K -> {self.temperature:.1f}K ({adjustment_reason})"
            )
    
    def get_convergence_status(self, window_size: int = 50) -> ConvergenceStatus:
        """判断收敛状态"""
        if len(self.energy_history) < window_size:
            return ConvergenceStatus.UNKNOWN
        
        recent_energies = self.energy_history[-window_size:]
        
        # 检查是否收敛（能量变化很小）
        energy_range = max(recent_energies) - min(recent_energies)
        if energy_range < 0.01:  # 能量范围 < 10 meV
            return ConvergenceStatus.CONVERGED
        
        # 检查趋势
        early_avg = np.mean(recent_energies[:window_size//2])
        late_avg = np.mean(recent_energies[window_size//2:])
        
        if late_avg < early_avg - 0.1:  # 能量明显下降
            return ConvergenceStatus.IMPROVING
        
        if abs(late_avg - early_avg) < 0.01:  # 能量几乎不变
            return ConvergenceStatus.STUCK
        
        # 检查震荡
        energy_diffs = np.diff(recent_energies)
        sign_changes = sum(1 for i in range(1, len(energy_diffs)) 
                         if energy_diffs[i] * energy_diffs[i-1] < 0)
        if sign_changes > window_size * 0.6:  # 频繁变号
            return ConvergenceStatus.OSCILLATING
        
        return ConvergenceStatus.UNKNOWN
    
    def suggest_restart_parameters(self) -> Dict[str, Any]:
        """建议重启参数（如果当前参数无效）"""
        status = self.get_convergence_status()
        
        if status == ConvergenceStatus.STUCK:
            return {
                "temperature": self.max_temp * 0.8,
                "reason": "Simulation stuck, trying higher temperature",
                "sweep_multiplier": 2.0  # 增加步数
            }
        
        elif status == ConvergenceStatus.OSCILLATING:
            return {
                "temperature": self.temperature * 0.5,
                "reason": "Energy oscillating, reducing temperature",
                "sweep_multiplier": 1.0
            }
        
        return {}


class NEBParameterOptimizer:
    """
    NEB 参数优化器
    
    根据 NEB 收敛情况自动调整参数：
    - 图片数量
    - spring constant
    - climbing image 策略
    """
    
    def __init__(self):
        self.attempt_history: List[Dict] = []
    
    def analyze_failure(
        self,
        neb_result: Any,
        initial_structure: Any,
        final_structure: Any
    ) -> List[ParameterAdjustment]:
        """
        分析 NEB 失败原因并给出调整建议
        
        Returns:
            参数调整建议列表
        """
        adjustments = []
        
        # 获取 NEB 结果信息
        fmax = getattr(neb_result, 'fmax', None)
        n_steps = getattr(neb_result, 'n_steps', 0)
        energies = getattr(neb_result, 'energies', [])
        
        # 失败类型 1: fmax 爆炸（力太大）
        if fmax and fmax > 10.0:
            adjustments.append(ParameterAdjustment(
                parameter="n_images",
                current_value=7,  # 默认值
                suggested_value=11,
                reason="fmax too large, need more images for smoother path",
                confidence=0.8
            ))
            
            adjustments.append(ParameterAdjustment(
                parameter="spring_constant",
                current_value=0.1,
                suggested_value=0.05,  # 降低弹簧常数
                reason="Reduce spring constant to allow more flexibility",
                confidence=0.7
            ))
        
        # 失败类型 2: 不收敛（步数用尽）
        elif n_steps >= 300:  # 假设最大 300 步
            adjustments.append(ParameterAdjustment(
                parameter="max_steps",
                current_value=300,
                suggested_value=500,
                reason="Not converged within max_steps, need more iterations",
                confidence=0.9
            ))
            
            adjustments.append(ParameterAdjustment(
                parameter="fmax",
                current_value=0.05,
                suggested_value=0.1,  # 放宽收敛标准
                reason="Relax convergence criterion slightly",
                confidence=0.6
            ))
        
        # 失败类型 3: 能量不单调（可能错过了 TS）
        elif energies and len(energies) > 2:
            # 检查是否有多个峰值
            peaks = sum(1 for i in range(1, len(energies)-1) 
                       if energies[i] > energies[i-1] and energies[i] > energies[i+1])
            
            if peaks > 1:
                adjustments.append(ParameterAdjustment(
                    parameter="climbing_image",
                    current_value=False,
                    suggested_value=True,
                    reason="Multiple peaks detected, use climbing image to find true TS",
                    confidence=0.8
                ))
            
            # 检查是否存在能量阱（不合理）
            valleys = sum(1 for i in range(1, len(energies)-1) 
                         if energies[i] < energies[i-1] and energies[i] < energies[i+1])
            if valleys > 0:
                adjustments.append(ParameterAdjustment(
                    parameter="interpolation_method",
                    current_value="linear",
                    suggested_value="idpp",  # Image Dependent Pair Potential
                    reason="Energy wells detected, use better interpolation",
                    confidence=0.7
                ))
        
        # 失败类型 4: 初末态结构问题（从验证器获取）
        # 这部分在 enhanced_neb_validator 中处理
        
        return adjustments
    
    def get_optimized_parameters(
        self,
        attempt_number: int,
        previous_failures: List[Dict]
    ) -> Dict[str, Any]:
        """
        根据失败历史获取优化后的参数
        
        Parameters:
            attempt_number: 当前尝试次数
            previous_failures: 之前失败的记录
        
        Returns:
            优化后的参数字典
        """
        params = {
            "n_images": 7,
            "spring_constant": 0.1,
            "max_steps": 300,
            "fmax": 0.05,
            "climbing_image": False,
            "interpolation": "linear"
        }
        
        # 根据尝试次数逐步调整
        if attempt_number == 1:
            # 第一次：标准参数
            pass
        
        elif attempt_number == 2:
            # 第二次：增加图片数，使用更好的插值
            params["n_images"] = 9
            params["interpolation"] = "idpp"
        
        elif attempt_number == 3:
            # 第三次：更多图片，CI-NEB，更宽松的收敛
            params["n_images"] = 11
            params["climbing_image"] = True
            params["max_steps"] = 500
            params["fmax"] = 0.1
        
        else:
            # 更多尝试：激进调整
            params["n_images"] = 15
            params["spring_constant"] = 0.05
            params["max_steps"] = 800
            params["climbing_image"] = True
        
        # 根据历史失败模式调整
        for failure in previous_failures:
            adjustments = self.analyze_failure(
                failure.get('result'),
                failure.get('initial'),
                failure.get('final')
            )
            
            for adj in adjustments:
                if adj.confidence > 0.7:
                    params[adj.parameter] = adj.suggested_value
        
        logger.info(f"NEB optimized parameters (attempt {attempt_number}): {params}")
        return params


class AdaptiveNEBRunner:
    """
    自适应 NEB 运行器
    
    自动尝试不同的参数组合直到成功
    """
    
    def __init__(self, max_attempts: int = 5):
        self.max_attempts = max_attempts
        self.optimizer = NEBParameterOptimizer()
        self.failure_history: List[Dict] = []
    
    def run_with_auto_retry(
        self,
        neb_fn: callable,
        initial: Any,
        final: Any,
        validate_fn: Optional[callable] = None
    ) -> Tuple[Any, Dict]:
        """
        运行 NEB，自动重试不同的参数
        
        Parameters:
            neb_fn: NEB 计算函数 (params) -> result
            initial: 初态结构
            final: 末态结构
            validate_fn: 结果验证函数
        
        Returns:
            (result, metadata)
        """
        for attempt in range(1, self.max_attempts + 1):
            logger.info(f"NEB attempt {attempt}/{self.max_attempts}")
            
            # 获取优化参数
            params = self.optimizer.get_optimized_parameters(
                attempt, self.failure_history
            )
            
            try:
                # 运行 NEB
                result = neb_fn(params)
                
                # 验证结果
                if validate_fn and not validate_fn(result):
                    raise ValueError("Validation failed")
                
                # 成功！
                metadata = {
                    "attempts": attempt,
                    "successful_params": params,
                    "failure_history": len(self.failure_history)
                }
                logger.info(f"NEB succeeded on attempt {attempt}")
                return result, metadata
            
            except Exception as e:
                logger.warning(f"NEB attempt {attempt} failed: {e}")
                self.failure_history.append({
                    "attempt": attempt,
                    "params": params,
                    "error": str(e),
                    "result": getattr(e, 'partial_result', None)
                })
        
        # 所有尝试都失败
        raise RuntimeError(f"NEB failed after {self.max_attempts} attempts")


# =============================================================================
# 集成到工作流
# =============================================================================

class AdaptiveWorkflowParameters:
    """
    自适应工作流参数管理
    
    根据前面的步骤结果自动调整后续参数
    """
    
    def __init__(self, base_config: Dict[str, Any]):
        self.base_config = base_config
        self.adjustments: List[Dict] = []
    
    def adjust_for_surface_complexity(
        self,
        surface_result: Any
    ) -> Dict[str, Any]:
        """
        根据表面复杂度调整参数
        
        复杂表面（大晶胞、多元素）需要更多 MC 步数
        """
        config = self.base_config.copy()
        
        # 获取表面信息
        n_atoms = getattr(surface_result, 'n_atoms', 0)
        n_elements = len(set(getattr(surface_result, 'symbols', [])))
        
        # 调整 MC 步数
        if n_atoms > 100:
            config['mc_total_sweeps'] = config.get('mc_total_sweeps', 100) * 2
            self.adjustments.append({
                "reason": f"Large surface ({n_atoms} atoms), doubled MC sweeps",
                "parameter": "mc_total_sweeps",
                "new_value": config['mc_total_sweeps']
            })
        
        if n_elements > 2:
            config['mc_temperature_k'] = config.get('mc_temperature_k', 500) * 1.2
            self.adjustments.append({
                "reason": f"Multi-element surface ({n_elements} elements), increased temperature",
                "parameter": "mc_temperature_k",
                "new_value": config['mc_temperature_k']
            })
        
        return config
    
    def adjust_for_reaction_complexity(
        self,
        pathway_result: Any
    ) -> Dict[str, Any]:
        """
        根据反应复杂度调整参数
        
        复杂反应（多步骤、高能量变化）需要更精细的 NEB
        """
        config = self.base_config.copy()
        
        n_steps = len(getattr(pathway_result, 'steps', []))
        
        if n_steps > 5:
            # 多步骤反应，减少每步的 NEB 图片数以节省时间
            config['neb_n_images'] = 5  # 而不是默认 7
            self.adjustments.append({
                "reason": f"Complex reaction ({n_steps} steps), reduced NEB images to save time",
                "parameter": "neb_n_images",
                "new_value": 5
            })
        
        return config
    
    def get_adjustment_summary(self) -> str:
        """获取调整摘要"""
        if not self.adjustments:
            return "No parameter adjustments made"
        
        lines = ["Parameter Adjustments:"]
        for adj in self.adjustments:
            lines.append(f"  - {adj['parameter']}: {adj['new_value']} ({adj['reason']})")
        return "\n".join(lines)


# =============================================================================
# 便捷函数
# =============================================================================

def auto_tune_mc_parameters(
    energy_history: List[float],
    acceptance_history: List[bool],
    current_temp: float
) -> Dict[str, Any]:
    """
    便捷的 MC 参数自动调优
    
    Example:
        >>> energies = [...]  # MC 能量历史
        >>> accepted = [...]  # 每一步是否接受
        >>> new_params = auto_tune_mc_parameters(energies, accepted, 500.0)
        >>> print(new_params['temperature'])
    """
    tuner = MCTemperatureTuner(initial_temperature=current_temp)
    
    for i, (e, a) in enumerate(zip(energy_history, acceptance_history)):
        tuner.update(e, a, i)
    
    status = tuner.get_convergence_status()
    restart_params = tuner.suggest_restart_parameters()
    
    return {
        "current_temperature": tuner.temperature,
        "convergence_status": status.value,
        "restart_suggested": len(restart_params) > 0,
        "restart_params": restart_params
    }
