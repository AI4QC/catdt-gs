"""
最终综合测试：所有增强功能
"""

import sys
from pathlib import Path
import tempfile
import time

repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

from ase import Atoms, Atom
import numpy as np


def test_parallel_pathway():
    """测试并行路径探索"""
    print("\n" + "="*60)
    print("Test: Parallel Pathway Explorer")
    print("="*60)
    
    from camel_agents.parallel_pathway import (
        ParallelPathwayExplorer, DependencyGraph, ParallelTask,
        parallel_compute_barriers
    )
    
    explorer = ParallelPathwayExplorer(max_workers=2)
    
    # 测试依赖图
    print("\n1. Testing Dependency Graph...")
    graph = DependencyGraph()
    
    tasks = [
        ParallelTask("A", "design", {}, [], priority=1),
        ParallelTask("B", "design", {}, ["A"], priority=2),
        ParallelTask("C", "neb", {}, ["A"], priority=1),
        ParallelTask("D", "validation", {}, ["B", "C"], priority=3),
    ]
    
    for task in tasks:
        graph.add_task(task)
    
    groups = graph.get_parallel_groups()
    print(f"   Parallel groups: {groups}")
    assert groups[0] == ["A"], "First group should be A"
    assert set(groups[1]) == {"B", "C"}, "Second group should be B and C"
    assert groups[2] == ["D"], "Third group should be D"
    
    # 测试最佳路径选择
    print("\n2. Testing Best Pathway Selection...")
    pathways = [
        {"name": "path1", "confidence": 0.9, "steps": [1, 2, 3]},
        {"name": "path2", "confidence": 0.7, "steps": [1, 2]},
        {"name": "path3", "confidence": 0.8, "steps": [1, 2, 3, 4]},
    ]
    
    best = explorer.select_best_pathway(pathways, criteria="confidence")
    print(f"   Best by confidence: {best['name']}")
    assert best["name"] == "path1"
    
    best_fewest = explorer.select_best_pathway(pathways, criteria="fewest_steps")
    print(f"   Best by fewest steps: {best_fewest['name']}")
    assert best_fewest["name"] == "path2"
    
    # 测试路径比较
    print("\n3. Testing Pathway Comparison...")
    from camel_agents.parallel_pathway import compare_pathways
    
    comparison = compare_pathways(pathways)
    print(f"   Total pathways: {comparison['n_pathways']}")
    print(f"   Best pathway index: {comparison['ranking']['best_pathway_idx']}")
    
    return True


def test_adaptive_parameters():
    """测试自适应参数调整"""
    print("\n" + "="*60)
    print("Test: Adaptive Parameter Tuner")
    print("="*60)
    
    from camel_agents.adaptive_parameters import (
        MCTemperatureTuner, NEBParameterOptimizer,
        auto_tune_mc_parameters, ConvergenceStatus
    )
    
    # 测试 MC 温度调整
    print("\n1. Testing MC Temperature Tuner...")
    tuner = MCTemperatureTuner(initial_temperature=500.0)
    
    # 模拟低接受率场景
    print("   Simulating low acceptance rate...")
    for i in range(100):
        energy = -100 + i * 0.01  # 缓慢降低的能量
        accepted = i % 10 == 0  # 10% 接受率
        tuner.update(energy, accepted, i)
    
    print(f"   Initial temp: 500.0K")
    print(f"   Current temp: {tuner.temperature:.1f}K")
    print(f"   Temperature history: {[f'{t:.0f}' for t in tuner.temperature_history]}")
    
    status = tuner.get_convergence_status()
    print(f"   Convergence status: {status.value}")
    
    # 测试 NEB 参数优化
    print("\n2. Testing NEB Parameter Optimizer...")
    optimizer = NEBParameterOptimizer()
    
    # 模拟失败结果
    class MockNEBResult:
        fmax = 15.0  # 力太大
        n_steps = 50
        energies = [0, 1, 2, 1.5, 1, 0.5]  # 多个峰值
    
    adjustments = optimizer.analyze_failure(
        MockNEBResult(), None, None
    )
    
    print(f"   Adjustments for fmax=15.0:")
    for adj in adjustments:
        print(f"     - {adj.parameter}: {adj.current_value} -> {adj.suggested_value} ({adj.reason})")
    
    # 测试参数序列
    print("\n3. Testing Parameter Sequence...")
    for attempt in range(1, 4):
        params = optimizer.get_optimized_parameters(attempt, [])
        print(f"   Attempt {attempt}: n_images={params['n_images']}, CI={params['climbing_image']}")
    
    # 测试便捷的 MC 调优函数
    print("\n4. Testing Auto-tune MC Function...")
    energies = [-100 + np.sin(i/10)*0.5 for i in range(100)]
    accepted = [True if i % 3 == 0 else False for i in range(100)]
    
    result = auto_tune_mc_parameters(energies, accepted, 500.0)
    print(f"   Current temperature: {result['current_temperature']:.1f}K")
    print(f"   Convergence status: {result['convergence_status']}")
    print(f"   Restart suggested: {result['restart_suggested']}")
    
    return True


def test_workflow_dashboard():
    """测试工作流仪表板"""
    print("\n" + "="*60)
    print("Test: Workflow Dashboard")
    print("="*60)
    
    from camel_agents.tooling.dashboard import WorkflowDashboard, WorkflowMonitor
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 创建仪表板
        print("\n1. Creating Dashboard...")
        dashboard = WorkflowDashboard(
            run_id="test_dashboard_001",
            output_dir=tmpdir,
            auto_refresh=False  # 不启动后台线程
        )
        
        # 更新步骤状态
        print("\n2. Updating Step Status...")
        steps = [
            ("agent1_surface", "completed", 1.0, {"n_surfaces": 3}),
            ("agent2_adsorption", "completed", 1.0, {"best_energy": -1.5}),
            ("agent3_reconstruction", "running", 0.6, {"current_sweep": 60}),
            ("agent4_pathway", "pending", 0.0, {}),
        ]
        
        for name, status, progress, metrics in steps:
            dashboard.update_step(name, status, progress, metrics)
            print(f"   {name}: {status} ({progress*100:.0f}%)")
        
        # 生成仪表板
        print("\n3. Generating Dashboard...")
        dashboard_path = dashboard.generate_dashboard()
        print(f"   Dashboard saved: {dashboard_path}")
        
        # 检查指标
        print("\n4. Checking Metrics...")
        metrics = dashboard.metrics
        print(f"   Run ID: {metrics.run_id}")
        print(f"   Status: {metrics.status}")
        print(f"   Overall Progress: {metrics.overall_progress*100:.1f}%")
        print(f"   Total Duration: {metrics.total_duration:.1f}s")
        
        # 生成摘要报告
        print("\n5. Generating Summary Report...")
        report_path = dashboard.generate_summary_report()
        print(f"   Report saved: {report_path}")
        
        # 读取报告内容
        with open(report_path, 'r') as f:
            report_content = f.read()
        print(f"   Report length: {len(report_content)} characters")
        
        # 测试上下文管理器
        print("\n6. Testing Context Manager...")
        with WorkflowMonitor("test_context", tmpdir, auto_refresh=False) as monitor:
            monitor.update_step("step1", "running")
            time.sleep(0.1)
            monitor.update_step("step1", "completed", 1.0)
        
        print("   Context manager completed successfully")
    
    return True


def test_integration():
    """测试功能集成"""
    print("\n" + "="*60)
    print("Test: Feature Integration")
    print("="*60)
    
    # 模拟一个完整的工作流场景
    print("\nSimulating a workflow scenario...")
    
    from camel_agents.tooling.dashboard import WorkflowDashboard
    from camel_agents.adaptive_parameters import MCTemperatureTuner
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. 创建仪表板
        dashboard = WorkflowDashboard("integration_test", tmpdir, auto_refresh=False)
        
        # 2. Agent 1: 表面生成
        print("\n1. Agent 1: Surface Generation")
        dashboard.update_step("agent1", "running")
        time.sleep(0.1)
        dashboard.update_step("agent1", "completed", 1.0, {"n_surfaces": 3})
        
        # 3. Agent 2: 吸附
        print("2. Agent 2: Adsorption")
        dashboard.update_step("agent2", "running")
        time.sleep(0.1)
        dashboard.update_step("agent2", "completed", 1.0, {"sites_tested": 10})
        
        # 4. Agent 3: MC 重构（使用自适应参数）
        print("3. Agent 3: MC Reconstruction (with adaptive tuning)")
        dashboard.update_step("agent3", "running", 0.0)
        
        tuner = MCTemperatureTuner(initial_temperature=500.0)
        
        # 模拟 MC 过程
        for i in range(50):
            # 模拟能量和接受
            energy = -100 + np.random.random() * 0.1
            accepted = np.random.random() < 0.3
            
            tuner.update(energy, accepted, i)
            dashboard.update_step("agent3", "running", i/100)
            
            # 每 25 步更新一次显示
            if i % 25 == 0:
                print(f"   Step {i}: T={tuner.temperature:.0f}K, E={energy:.3f}")
        
        dashboard.update_step("agent3", "completed", 1.0, {
            "final_temp": tuner.temperature,
            "convergence": tuner.get_convergence_status().value
        })
        
        # 5. 生成最终仪表板
        dashboard.generate_dashboard()
        dashboard.generate_summary_report()
        
        print("\n✅ Integration test completed successfully!")
    
    return True


def main():
    """运行所有测试"""
    print("\n" + "="*70)
    print("Final Comprehensive Feature Tests")
    print("="*70)
    
    tests = [
        ("Parallel Pathway Explorer", test_parallel_pathway),
        ("Adaptive Parameter Tuner", test_adaptive_parameters),
        ("Workflow Dashboard", test_workflow_dashboard),
        ("Feature Integration", test_integration),
    ]
    
    results = {}
    
    for name, test_func in tests:
        try:
            print(f"\n{'='*70}")
            print(f"Running: {name}")
            print(f"{'='*70}")
            results[name] = test_func()
        except Exception as e:
            print(f"\n❌ Test failed: {e}")
            import traceback
            traceback.print_exc()
            results[name] = False
    
    # 汇总
    print("\n" + "="*70)
    print("Test Results Summary")
    print("="*70)
    
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {name}: {status}")
    
    total = len(results)
    passed = sum(results.values())
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! All features are working correctly.")
    
    return all(results.values())


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
