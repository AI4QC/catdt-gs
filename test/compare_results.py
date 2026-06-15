#!/usr/bin/env python
"""
Compare direct test results with agent system results
"""
import json
import os

print("=" * 80)
print("Comparing Direct Test vs Agent System Results")
print("=" * 80)

# Check agent system results
agent_report_path = "output/agent_auto_retry_result/agent_report.json"
if os.path.exists(agent_report_path):
    with open(agent_report_path) as f:
        agent_data = json.load(f)

    print("\n📊 Agent System Results (3rd attempt):")
    print(f"  Parameters: {agent_data['summary']['final_parameters']}")
    print(f"  Best surface: {agent_data['summary']['best_surface']}")
    print(f"  Overall ΔE: {agent_data['summary']['overall_reaction_energy']:.3f} eV")
    print(f"  Validation: {'✅ PASS' if agent_data['validation']['valid'] else '❌ FAIL'}")
    print(f"  Critical issues: {len(agent_data['validation']['critical_issues'])}")
    if agent_data['validation']['critical_issues']:
        for issue in agent_data['validation']['critical_issues'][:3]:
            print(f"    - {issue}")
else:
    print("\n❌ Agent system results not found")

# Check direct test results
direct_test_path = "output/final_test_co_oxidation"
if os.path.exists(direct_test_path):
    print(f"\n📊 Direct Test Results:")
    print(f"  Output directory exists: {direct_test_path}")

    # Try to find pathway analysis results
    import glob
    pathway_files = glob.glob(f"{direct_test_path}/**/pathway_analysis.json", recursive=True)
    if pathway_files:
        print(f"  Found {len(pathway_files)} pathway analysis files")
else:
    print(f"\n⚠️  Direct test results not found at {direct_test_path}")

print("\n" + "=" * 80)
print("Key Differences to Investigate:")
print("=" * 80)
print("1. Agent system enabled calculate_barriers=True in 3rd attempt")
print("2. Agent system used more surfaces (3 vs 1)")
print("3. Agent system used more adsorption sites (8 vs 5)")
print("4. Need to check if barrier calculation introduces errors")
print("=" * 80)
