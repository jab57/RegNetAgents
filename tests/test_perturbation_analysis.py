#!/usr/bin/env python3
"""
Test script for perturbation analysis feature
Tests the new perturbation analysis capability with TP53
"""

import asyncio
import json
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def test_perturbation_analysis():
    """Test perturbation analysis with TP53"""
    print("=" * 80)
    print("Testing Perturbation Analysis Feature")
    print("=" * 80)

    # Initialize workflow
    workflow = RegNetAgentsWorkflow()

    # Test with TP53 (heavily regulated tumor suppressor)
    print("\n1. Testing with TP53 (tumor suppressor with many regulators)...")
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    # Extract perturbation analysis
    perturbation = result.get('perturbation_analysis', {})

    if perturbation and perturbation.get('status') != 'skipped':
        print("\n[SUCCESS] Perturbation analysis completed successfully!")
        print(f"\nTarget Gene: {perturbation.get('target_gene')}")
        print(f"Cell Type: {perturbation.get('cell_type')}")
        print(f"Baseline Regulators: {perturbation.get('baseline_regulators')}")
        print(f"\n{perturbation.get('summary')}")

        # Show therapeutic insights
        insights = perturbation.get('therapeutic_insights', {})
        print(f"\n{'=' * 80}")
        print("THERAPEUTIC INSIGHTS")
        print('=' * 80)
        print(f"\nSummary: {insights.get('summary')}")

        top_target = insights.get('top_target', {})
        if top_target:
            print(f"\nTop Therapeutic Target:")
            print(f"  Regulator: {top_target.get('regulator')}")
            print(f"  Impact Score: {top_target.get('impact_score')}")
            print(f"  Potential: {top_target.get('therapeutic_potential')}")

        print(f"\nHigh Impact Targets: {insights.get('high_impact_targets', 0)}")
        print(f"\nStrategy: {insights.get('strategy')}")
        print(f"\nRecommendation: {insights.get('recommendation')}")

        # Show top 5 perturbation results
        print(f"\n{'=' * 80}")
        print("TOP 5 PERTURBATION TARGETS")
        print('=' * 80)

        results = perturbation.get('perturbation_results', [])
        for i, target in enumerate(results[:5], 1):
            print(f"\n{i}. {target.get('regulator')}")
            print(f"   Impact Score: {target.get('impact_score')}")
            print(f"   Regulatory Loss: {target.get('regulatory_loss_pct')}%")
            print(f"   Downstream Targets: {target.get('regulator_downstream_targets')}")
            print(f"   Therapeutic Potential: {target.get('therapeutic_potential')}")
            print(f"   Cascade Overlap: {target.get('cascade_overlap')} genes")
            print(f"   Note: {target.get('druggability_notes')}")

            affected = target.get('affected_cascades', [])
            if affected:
                print(f"   Affected Cascades: {', '.join(affected[:5])}")

    else:
        print("\n[FAILED] Perturbation analysis was skipped or failed")
        print(f"Reason: {perturbation.get('reason', 'Unknown')}")

    # Save full results
    output_file = "tp53_perturbation_results.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n\nFull results saved to: {output_file}")

    print("\n" + "=" * 80)
    print("Test completed!")
    print("=" * 80)

if __name__ == "__main__":
    asyncio.run(test_perturbation_analysis())
