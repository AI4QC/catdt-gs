"""
综合测试：增强型 CatDT 功能

测试：
1. 断点续传工作流
2. 智能迭代控制
3. 跨运行缓存
4. 增强型 NEB 验证
"""

import sys
from pathlib import Path
import tempfile
import time

repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))

from ase import Atoms, Atom
from ase.io import write, read
import numpy as np


def test_smart_iteration_controller():
    """测试智能迭代控制器"""
    print("\n" + "="*60)
    print("Test 1: Smart Iteration Controller")
    print("="*60)
    
    from camel_agents.resumable import SmartIterationController
    from camel_agents.workflow import ValidationReport
    
    controller = SmartIterationController(max_iterations=5)
    
    # 模拟迭代过程
    reports = [
        ValidationReport(status="FAIL", issues=["issue1"], feedback="fix1"),
        ValidationReport(status="FAIL", issues=["issue1", "issue2"], feedback="fix2"),
        ValidationReport(status="WARNING", issues=["issue1"], feedback="fix3"),
        ValidationReport(status="PASS", issues=[], feedback=""),
    ]
    
    for i, report in enumerate(reports, 1):
        should_continue, reason = controller.should_continue(i, report)
        print(f"  Iteration {i}: {report.status} -> {'Continue' if should_continue else 'Stop'} ({reason})")
    
    # 测试循环检测
    controller2 = SmartIterationController(max_iterations=10)
    loop_report = ValidationReport(status="FAIL", issues=["same"], feedback="same_feedback")
    
    print("\n  Testing loop detection:")
    for i in range(1, 5):
        should_continue, reason = controller2.should_continue(i, loop_report)
        print(f"    Iteration {i}: {'Continue' if should_continue else 'Stop'} ({reason})")
        if not should_continue:
            break
    
    summary = controller.get_iteration_summary()
    print(f"\n  Summary: {summary['total_iterations']} iterations, "
          f"{summary['pass_count']} passes, {summary['fail_count']} fails")
    
    return True


def test_cross_run_cache():
    """测试跨运行缓存"""
    print("\n" + "="*60)
    print("Test 2: Cross-Run Cache")
    print("="*60)
    
    from camel_agents.cross_run_cache import CrossRunCache, StructureFingerprinter
    
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = CrossRunCache(cache_dir=tmpdir, max_size_gb=0.1)
        
        # 创建测试结构
        structure = Atoms('Pt4', 
                         positions=[[0,0,0], [2.8,0,0], [0,2.8,0], [2.8,2.8,0]],
                         cell=[5.6, 5.6, 20, 90, 90, 90],
                         pbc=[True, True, False])
        
        # 测试指纹生成
        print("  Testing fingerprint generation...")
        fp1 = StructureFingerprinter.generate_fingerprint(structure)
        fp2 = StructureFingerprinter.generate_fingerprint(structure)
        print(f"    Fingerprint: {fp1[:32]}...")
        print(f"    Same structure gives same fingerprint: {fp1 == fp2}")
        
        # 稍微修改结构（在容差内）
        structure2 = structure.copy()
        structure2.positions[0] += [0.005, 0, 0]  # 小变化
        fp3 = StructureFingerprinter.generate_fingerprint(structure2)
        print(f"    Similar structure (within tolerance) gives same fingerprint: {fp1 == fp3}")
        
        # 测试缓存存取
        print("\n  Testing cache operations...")
        
        def compute_mock():
            time.sleep(0.1)  # 模拟计算
            return {"energy": -123.45, "iterations": 100}
        
        # 第一次：应该计算
        start = time.time()
        result1 = cache.compute_or_cache("surff", structure, compute_mock, 
                                        params={"top_n": 5})
        time1 = time.time() - start
        print(f"    First call (compute): {time1:.3f}s")
        
        # 第二次：应该命中缓存
        start = time.time()
        result2 = cache.compute_or_cache("surff", structure, compute_mock,
                                        params={"top_n": 5})
        time2 = time.time() - start
        print(f"    Second call (cached): {time2:.3f}s")
        print(f"    Results match: {result1 == result2}")
        print(f"    Cache hit rate: {cache.get_stats()['hit_rate']}")
        
        # 测试统计
        print("\n  Cache statistics:")
        stats = cache.get_stats()
        print(f"    Total entries: {stats['total_entries']}")
        print(f"    Hits: {stats['hits']}, Misses: {stats['misses']}")
        print(f"    Saved time: {stats['saved_time_hours']*3600:.3f}s")
    
    return True


def test_resumable_basics():
    """测试可恢复工作流的基本功能"""
    print("\n" + "="*60)
    print("Test 3: Resumable Workflow Basics")
    print("="*60)
    
    from camel_agents.resumable import ResumableCatDTWorkflow
    
    # 创建工作流（不初始化 agents，只做基础测试）
    workflow = ResumableCatDTWorkflow(
        initialize_agents=False,
        enable_checkpointing=True,
        enable_smart_iteration=True
    )
    
    # 测试步骤管理
    print("  Workflow steps:")
    for i, step in enumerate(workflow.WORKFLOW_STEPS):
        print(f"    {i}: {step}")
    
    # 测试下一步计算
    print("\n  Testing _get_next_step:")
    print(f"    agent1_surface -> {workflow._get_next_step('agent1_surface')}")
    print(f"    agent3_reconstruction -> {workflow._get_next_step('agent3_reconstruction')}")
    print(f"    agent7_report -> {workflow._get_next_step('agent7_report')}")
    
    # 测试失败步骤处理
    print("\n  Testing failed step handling:")
    print(f"    agent3_reconstruction_failed -> {workflow._get_next_step('agent3_reconstruction_failed')}")
    
    # 检查是否有检查点管理器
    print(f"\n  Checkpoint manager initialized: {workflow.checkpoint_manager is not None}")
    print(f"  Smart iteration controller initialized: {workflow.iteration_controller is not None}")
    
    return True


def test_checkpoint_integration():
    """测试检查点集成"""
    print("\n" + "="*60)
    print("Test 4: Checkpoint Integration")
    print("="*60)
    
    from camel_agents.resumable import ResumableCatDTWorkflow, WorkflowState
    from camel_agents.tooling.persistence import CheckpointManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 创建检查点管理器
        checkpoint_mgr = CheckpointManager(checkpoint_dir=tmpdir)
        
        # 创建模拟状态
        state = WorkflowState(
            run_id="test_checkpoint_001",
            output_base_dir="/tmp/test_output",
            reaction_description="CO oxidation on Pt",
            surface_indices=[0, 1, 2, 3],
            adsorbate_indices=[4, 5],
            surface_path="/tmp/surface.vasp",
            intermediates=["*CO", "*O", "*CO2"]
        )
        
        # 保存检查点
        print("  Saving checkpoint...")
        checkpoint_path = checkpoint_mgr.save_checkpoint(
            run_id="test_checkpoint_001",
            step_name="agent3_reconstruction",
            step_number=3,
            state=state,
            artifacts={},
            metadata={"test": True}
        )
        print(f"    Checkpoint saved: {Path(checkpoint_path).name}")
        
        # 加载检查点
        print("  Loading checkpoint...")
        loaded = checkpoint_mgr.load_latest_checkpoint("test_checkpoint_001")
        
        if loaded:
            print(f"    Loaded: Step {loaded.step_number} - {loaded.step_name}")
            print(f"    Run ID: {loaded.run_id}")
            print(f"    Reaction: {loaded.state_data.get('reaction_description', 'N/A')}")
            
            # 恢复状态
            workflow = ResumableCatDTWorkflow(initialize_agents=False)
            restored_state = workflow._restore_state_from_checkpoint(loaded)
            
            print(f"    Restored state: {restored_state.run_id}")
            print(f"    Surface indices: {restored_state.surface_indices}")
            print(f"    Adsorbate indices: {restored_state.adsorbate_indices}")
            print(f"    Intermediates: {restored_state.intermediates}")
        
        # 列出检查点
        print("\n  Listing checkpoints:")
        checkpoints = checkpoint_mgr.list_checkpoints("test_checkpoint_001")
        for cp in checkpoints:
            print(f"    Step {cp['number']}: {cp['step']}")
    
    return True


def test_enhanced_validator_with_real_structures():
    """用真实结构测试增强验证器"""
    print("\n" + "="*60)
    print("Test 5: Enhanced Validator with Real Structures")
    print("="*60)
    
    from core.pathway.enhanced_neb_validator import (
        EnhancedNEBValidator, ValidationStatus
    )
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 创建合理的 CO 加氢结构
        surface = Atoms('Pt9',
                       positions=[[i*2.8, j*2.8, 0] for i in range(3) for j in range(3)],
                       cell=[8.4, 8.4, 20, 90, 90, 90],
                       pbc=[True, True, False])
        
        # 创建合理的初态：*CO + H（H 在表面，准备迁移到 C）
        initial = surface.copy()
        initial.append(Atom('C', position=[2.8, 2.8, 2.0]))
        initial.append(Atom('O', position=[2.8, 2.8, 3.2]))
        initial.append(Atom('H', position=[2.8, 2.0, 2.2]))  # H 靠近 C
        initial.info['adsorbate_indices'] = [9, 10, 11]
        
        # 创建合理的末态：*CHO（H 已连接到 C）
        final = surface.copy()
        final.append(Atom('C', position=[2.8, 2.8, 2.0]))
        final.append(Atom('O', position=[2.8, 2.8, 3.2]))
        final.append(Atom('H', position=[3.3, 2.8, 2.3]))  # H 连接到 C
        final.info['adsorbate_indices'] = [9, 10, 11]
        
        # 保存结构
        initial_path = Path(tmpdir) / "initial_good.vasp"
        final_path = Path(tmpdir) / "final_good.vasp"
        write(initial_path, initial)
        write(final_path, final)
        
        # 验证
        print("  Testing with reasonable structures...")
        validator = EnhancedNEBValidator(auto_fix=False)
        report = validator.validate_endpoint_pair(initial, final, "hydrogenation")
        
        print(f"    Status: {report.status.value}")
        print(f"    Issues: {len(report.issues)}")
        print(f"    Critical issues: {len([i for i in report.issues if i.severity == 'critical'])}")
        
        if report.status == ValidationStatus.PASS:
            print("    ✅ Good structures passed validation!")
        else:
            print(f"    ⚠️  Warnings found: {[i.description for i in report.issues if i.severity != 'critical']}")
        
        # 测试不合理的结构（H 离得很远）
        print("\n  Testing with problematic structures...")
        bad_initial = surface.copy()
        bad_initial.append(Atom('C', position=[2.8, 2.8, 2.0]))
        bad_initial.append(Atom('O', position=[2.8, 2.8, 3.2]))
        bad_initial.append(Atom('H', position=[5.0, 5.0, 5.0]))  # H 离得很远！
        bad_initial.info['adsorbate_indices'] = [9, 10, 11]
        
        report2 = validator.validate_endpoint_pair(bad_initial, final, "hydrogenation")
        
        print(f"    Status: {report2.status.value}")
        print(f"    Issues: {len(report2.issues)}")
        
        # 检查是否检测到 H 太远的问题
        h_far_issues = [i for i in report2.issues if 'too far' in i.description.lower()]
        if h_far_issues:
            print(f"    ✅ Correctly detected H being too far: {h_far_issues[0].description}")
        
        # 测试自动修复
        print("\n  Testing auto-fix...")
        validator_fix = EnhancedNEBValidator(auto_fix=True, fix_attempts=3)
        report3 = validator_fix.validate_endpoint_pair(bad_initial, final, "hydrogenation")
        
        print(f"    Status after fix: {report3.status.value}")
        print(f"    Fixes applied: {len(report3.fixes_applied)}")
        for fix in report3.fixes_applied:
            print(f"      - {fix}")
    
    return True


def main():
    """运行所有测试"""
    print("\n" + "="*70)
    print("Enhanced CatDT Feature Tests")
    print("="*70)
    
    tests = [
        ("Smart Iteration Controller", test_smart_iteration_controller),
        ("Cross-Run Cache", test_cross_run_cache),
        ("Resumable Workflow Basics", test_resumable_basics),
        ("Checkpoint Integration", test_checkpoint_integration),
        ("Enhanced Validator", test_enhanced_validator_with_real_structures),
    ]
    
    results = {}
    
    for name, test_func in tests:
        try:
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
        print("\n🎉 All tests passed!")
    
    return all(results.values())


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
