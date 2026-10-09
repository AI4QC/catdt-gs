"""
Workflow Dashboard - 实时监控和可视化仪表板

提供：
1. 实时工作流状态监控
2. 步骤进度追踪
3. 资源使用监控
4. 结果可视化
"""

import html as html_lib
import json
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field, asdict
from datetime import datetime, timedelta
import threading
import time
from collections import deque

logger = logging.getLogger(__name__)


@dataclass
class StepMetrics:
    """步骤指标"""
    step_name: str
    status: str  # "pending", "running", "completed", "failed", "skipped"
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    progress: float = 0.0  # 0-1
    metrics: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    
    @property
    def duration(self) -> Optional[float]:
        """运行时长（秒）"""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        elif self.start_time:
            return (datetime.now() - self.start_time).total_seconds()
        return None
    
    def to_dict(self) -> Dict:
        return {
            "step_name": self.step_name,
            "status": self.status,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "end_time": self.end_time.isoformat() if self.end_time else None,
            "duration": self.duration,
            "progress": self.progress,
            "metrics": self.metrics,
            "error": self.error
        }


@dataclass
class WorkflowMetrics:
    """工作流指标"""
    run_id: str
    start_time: datetime
    end_time: Optional[datetime] = None
    status: str = "running"  # "running", "completed", "failed", "stopped"
    steps: List[StepMetrics] = field(default_factory=list)
    current_step: Optional[str] = None
    
    @property
    def total_duration(self) -> float:
        """总运行时长"""
        end = self.end_time or datetime.now()
        return (end - self.start_time).total_seconds()
    
    @property
    def overall_progress(self) -> float:
        """整体进度"""
        if not self.steps:
            return 0.0
        return sum(s.progress for s in self.steps) / len(self.steps)
    
    def to_dict(self) -> Dict:
        return {
            "run_id": self.run_id,
            "start_time": self.start_time.isoformat(),
            "end_time": self.end_time.isoformat() if self.end_time else None,
            "status": self.status,
            "total_duration": self.total_duration,
            "overall_progress": self.overall_progress,
            "current_step": self.current_step,
            "steps": [s.to_dict() for s in self.steps]
        }


class WorkflowDashboard:
    """
    工作流监控仪表板
    
    生成实时监控页面和报告
    """
    
    def __init__(
        self,
        run_id: str,
        output_dir: str = "output/dashboard",
        auto_refresh: bool = True,
        refresh_interval: int = 5  # 秒
    ):
        self.run_id = run_id
        self.output_dir = Path(output_dir) / run_id
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.metrics = WorkflowMetrics(
            run_id=run_id,
            start_time=datetime.now()
        )
        
        self.auto_refresh = auto_refresh
        self.refresh_interval = refresh_interval
        self._update_thread: Optional[threading.Thread] = None
        self._stop_event = threading.Event()
        # Guards self.metrics: update_step (caller thread) mutates it while
        # the background renderer thread iterates it.
        self._metrics_lock = threading.Lock()
        
        # 历史数据（用于趋势图）
        self.energy_history: deque = deque(maxlen=1000)
        self.progress_history: deque = deque(maxlen=100)
        
        logger.info(f"Dashboard initialized: {self.output_dir}")
    
    def start_monitoring(self):
        """开始后台监控"""
        if self.auto_refresh and not self._update_thread:
            self._stop_event.clear()
            self._update_thread = threading.Thread(target=self._background_update)
            self._update_thread.daemon = True
            self._update_thread.start()
            logger.info("Dashboard monitoring started")
    
    def stop_monitoring(self):
        """停止监控"""
        if self._update_thread:
            self._stop_event.set()
            self._update_thread.join(timeout=5)
            logger.info("Dashboard monitoring stopped")
    
    def _background_update(self):
        """后台更新线程"""
        while not self._stop_event.is_set():
            self.generate_dashboard()
            time.sleep(self.refresh_interval)
    
    def update_step(
        self,
        step_name: str,
        status: str,
        progress: Optional[float] = None,
        metrics: Optional[Dict] = None,
        error: Optional[str] = None
    ):
        """更新步骤状态"""
        with self._metrics_lock:
            # 查找或创建步骤
            step = next((s for s in self.metrics.steps if s.step_name == step_name), None)

            if step is None:
                step = StepMetrics(step_name=step_name, status=status)
                self.metrics.steps.append(step)

            # 更新状态
            old_status = step.status
            step.status = status

            # 时间戳管理
            if status == "running" and old_status != "running":
                step.start_time = datetime.now()
            elif status in ["completed", "failed"] and old_status not in ["completed", "failed"]:
                step.end_time = datetime.now()

            if progress is not None:
                step.progress = progress

            if metrics:
                step.metrics.update(metrics)

            if error:
                step.error = error

            # 更新当前步骤
            if status == "running":
                self.metrics.current_step = step_name

            # 更新工作流状态
            if all(s.status in ["completed", "failed", "skipped"] for s in self.metrics.steps):
                self.metrics.end_time = datetime.now()
                self.metrics.status = "completed" if all(s.status == "completed" for s in self.metrics.steps) else "failed"

            # 记录历史
            self.progress_history.append({
                "time": datetime.now().isoformat(),
                "progress": self.metrics.overall_progress,
                "step": step_name
            })
    
    def log_energy(self, step_name: str, energy: float, **kwargs):
        """记录能量值"""
        self.energy_history.append({
            "time": datetime.now().isoformat(),
            "step": step_name,
            "energy": energy,
            **kwargs
        })
    
    def generate_dashboard(self) -> Path:
        """生成监控页面"""
        with self._metrics_lock:
            html = self._generate_html()
            metrics_dict = self.metrics.to_dict()

        dashboard_path = self.output_dir / "index.html"
        with open(dashboard_path, 'w', encoding='utf-8') as f:
            f.write(html)

        # 同时保存 JSON 数据
        json_path = self.output_dir / "metrics.json"
        with open(json_path, 'w') as f:
            json.dump(metrics_dict, f, indent=2, default=str)

        return dashboard_path
    
    def _generate_html(self) -> str:
        """生成 HTML 页面"""
        m = self.metrics
        
        # 步骤表格行（用户来源的字符串一律 html.escape，防止注入）
        step_rows = []
        for step in m.steps:
            duration_str = f"{step.duration:.1f}s" if step.duration else "N/A"
            progress_bar = self._generate_progress_bar(step.progress)
            status_color = {
                "pending": "gray",
                "running": "blue",
                "completed": "green",
                "failed": "red",
                "skipped": "orange"
            }.get(step.status, "black")

            metrics_str = "<br>".join(
                f"{html_lib.escape(str(k))}: {html_lib.escape(str(v))}"
                for k, v in step.metrics.items()
            ) if step.metrics else "-"

            step_rows.append(f"""
                <tr>
                    <td>{html_lib.escape(str(step.step_name))}</td>
                    <td style="color: {status_color}">{html_lib.escape(str(step.status).upper())}</td>
                    <td>{progress_bar}</td>
                    <td>{duration_str}</td>
                    <td>{metrics_str}</td>
                    <td>{html_lib.escape(str(step.error)) if step.error else '-'}</td>
                </tr>
            """)
        
        step_table = "\n".join(step_rows) if step_rows else "<tr><td colspan='6'>No steps yet</td></tr>"
        
        # 整体进度
        overall_progress_bar = self._generate_progress_bar(m.overall_progress, width=400)
        
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>CatDT Workflow Dashboard - {html_lib.escape(str(m.run_id))}</title>
    <meta http-equiv="refresh" content="{self.refresh_interval}">
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
        }}
        .status-card {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 15px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .status-running {{ color: #2196F3; font-weight: bold; }}
        .status-completed {{ color: #4CAF50; font-weight: bold; }}
        .status-failed {{ color: #f44336; font-weight: bold; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th {{
            background: #333;
            color: white;
            padding: 12px;
            text-align: left;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ddd;
        }}
        tr:hover {{ background: #f9f9f9; }}
        .progress-bar {{
            background: #e0e0e0;
            border-radius: 10px;
            overflow: hidden;
            height: 20px;
        }}
        .progress-fill {{
            background: linear-gradient(90deg, #4CAF50, #8BC34A);
            height: 100%;
            transition: width 0.3s;
        }}
        .timestamp {{
            color: #666;
            font-size: 0.9em;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .metric-box {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        .metric-label {{
            color: #666;
            margin-top: 5px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧪 CatDT Workflow Dashboard</h1>
        <p>Run ID: <strong>{html_lib.escape(str(m.run_id))}</strong></p>
        <p class="timestamp">Started: {m.start_time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p class="timestamp">Last Update: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="metrics-grid">
        <div class="metric-box">
            <div class="metric-value status-{m.status}">{m.status.upper()}</div>
            <div class="metric-label">Status</div>
        </div>
        <div class="metric-box">
            <div class="metric-value">{m.total_duration/60:.1f}m</div>
            <div class="metric-label">Duration</div>
        </div>
        <div class="metric-box">
            <div class="metric-value">{len(m.steps)}</div>
            <div class="metric-label">Total Steps</div>
        </div>
        <div class="metric-box">
            <div class="metric-value">{sum(1 for s in m.steps if s.status == 'completed')}</div>
            <div class="metric-label">Completed</div>
        </div>
    </div>
    
    <div class="status-card">
        <h3>Overall Progress</h3>
        {overall_progress_bar}
        <p>{m.overall_progress*100:.1f}% complete</p>
        {f'<p>Current Step: <strong>{html_lib.escape(str(m.current_step))}</strong></p>' if m.current_step else ''}
    </div>
    
    <div class="status-card">
        <h3>Step Details</h3>
        <table>
            <thead>
                <tr>
                    <th>Step</th>
                    <th>Status</th>
                    <th>Progress</th>
                    <th>Duration</th>
                    <th>Metrics</th>
                    <th>Error</th>
                </tr>
            </thead>
            <tbody>
                {step_table}
            </tbody>
        </table>
    </div>
    
    <script>
        // Auto-scroll to current step
        const runningRow = document.querySelector('td[style*="blue"]');
        if (runningRow) {{
            runningRow.scrollIntoView({{behavior: 'smooth', block: 'center'}});
        }}
    </script>
</body>
</html>
        """
        
        return html
    
    def _generate_progress_bar(self, progress: float, width: int = 200) -> str:
        """生成进度条 HTML"""
        percentage = min(100, max(0, progress * 100))
        return f"""
            <div class="progress-bar" style="width: {width}px;">
                <div class="progress-fill" style="width: {percentage}%;"></div>
            </div>
            <span>{percentage:.0f}%</span>
        """
    
    def generate_summary_report(self) -> Path:
        """生成最终摘要报告"""
        m = self.metrics
        
        report_lines = [
            "=" * 60,
            "CatDT Workflow Summary Report",
            "=" * 60,
            f"",
            f"Run ID: {m.run_id}",
            f"Status: {m.status.upper()}",
            f"Start Time: {m.start_time.strftime('%Y-%m-%d %H:%M:%S')}",
            f"End Time: {m.end_time.strftime('%Y-%m-%d %H:%M:%S') if m.end_time else 'N/A'}",
            f"Total Duration: {m.total_duration/60:.1f} minutes",
            f"",
            f"Overall Progress: {m.overall_progress*100:.1f}%",
            f"",
            "Step Summary:",
            "-" * 60,
        ]
        
        for step in m.steps:
            duration = f"{step.duration:.1f}s" if step.duration else "N/A"
            report_lines.append(f"  {step.step_name}:")
            report_lines.append(f"    Status: {step.status.upper()}")
            report_lines.append(f"    Duration: {duration}")
            report_lines.append(f"    Progress: {step.progress*100:.0f}%")
            if step.metrics:
                for k, v in step.metrics.items():
                    report_lines.append(f"    {k}: {v}")
            if step.error:
                report_lines.append(f"    ERROR: {step.error}")
            report_lines.append("")
        
        report_lines.append("=" * 60)
        
        report_text = "\n".join(report_lines)
        
        report_path = self.output_dir / "summary_report.txt"
        with open(report_path, 'w') as f:
            f.write(report_text)
        
        return report_path


class WorkflowMonitor:
    """
    工作流监控器 - 便捷的上下文管理器
    
    Example:
        >>> with WorkflowMonitor("my_run") as monitor:
        ...     monitor.update_step("agent1", "running")
        ...     # 执行步骤
        ...     monitor.update_step("agent1", "completed", progress=1.0)
    """
    
    def __init__(
        self,
        run_id: str,
        output_dir: str = "output/dashboard",
        auto_refresh: bool = True
    ):
        self.dashboard = WorkflowDashboard(run_id, output_dir, auto_refresh)
    
    def __enter__(self):
        self.dashboard.start_monitoring()
        return self.dashboard
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:
            self.dashboard.metrics.status = "failed"
        else:
            self.dashboard.metrics.status = "completed"
        self.dashboard.metrics.end_time = datetime.now()
        
        self.dashboard.stop_monitoring()
        self.dashboard.generate_summary_report()
        
        return False


# =============================================================================
# 便捷函数
# =============================================================================

def create_dashboard_for_run(run_id: str, output_dir: str = "output/dashboard") -> WorkflowDashboard:
    """
    为运行创建仪表板
    
    Example:
        >>> dashboard = create_dashboard_for_run("co_to_ch4_001")
        >>> dashboard.update_step("agent1_surface", "running")
        >>> # ... 执行步骤
        >>> dashboard.update_step("agent1_surface", "completed", progress=1.0)
        >>> dashboard.generate_dashboard()
    """
    return WorkflowDashboard(run_id, output_dir)


def generate_workflow_summary(metrics_file: str) -> str:
    """从指标文件生成摘要"""
    with open(metrics_file, 'r') as f:
        data = json.load(f)
    
    lines = [
        "Workflow Summary",
        "=" * 50,
        f"Run ID: {data['run_id']}",
        f"Status: {data['status']}",
        f"Duration: {data['total_duration']/60:.1f} minutes",
        f"Progress: {data['overall_progress']*100:.1f}%",
        f"",
        "Steps:",
    ]
    
    for step in data['steps']:
        lines.append(f"  {step['step_name']}: {step['status']} ({step.get('duration', 'N/A')}s)")
    
    return "\n".join(lines)
