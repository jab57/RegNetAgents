#!/usr/bin/env python3
"""
Test script for the RegNetAgents LangGraph workflow
"""

import asyncio
import sys
import os

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def test_workflow():
    """Test the LangGraph workflow with APC gene"""
    print("Initializing RegNetAgents LangGraph Workflow...")

    try:
        workflow = RegNetAgentsWorkflow()
        print("Workflow initialized successfully")

        print("\nTesting with APC gene (epithelial_cell)...")
        result = await workflow.run_analysis(
            gene="APC",
            cell_type="epithelial_cell",
            analysis_depth="comprehensive"
        )

        if result.get('status') == 'error':
            print(f"Workflow failed: {result.get('error')}")
            return

        print("Workflow completed successfully!")
        print(f"\nWorkflow Results Summary:")
        print(f"Gene: {result.get('gene_analysis_summary', {}).get('gene')}")
        print(f"Regulatory Role: {result.get('gene_analysis_summary', {}).get('regulatory_role')}")
        print(f"Analysis Steps: {len(result.get('workflow_metadata', {}).get('steps_completed', []))}")
        print(f"Total Time: {result.get('workflow_metadata', {}).get('total_analysis_time', 0):.2f}ms")

        # Show key insights
        insights = result.get('key_insights', {})
        if insights:
            print(f"\nKey Insights:")
            for key, value in insights.items():
                print(f"  {key}: {value}")

        # Show completed steps
        steps = result.get('workflow_metadata', {}).get('steps_completed', [])
        print(f"\nCompleted Analysis Steps:")
        for i, step in enumerate(steps, 1):
            print(f"  {i}. {step}")

        return True

    except Exception as e:
        print(f"Test failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = asyncio.run(test_workflow())
    if success:
        print("\nLangGraph workflow test PASSED!")
    else:
        print("\nLangGraph workflow test FAILED!")
        sys.exit(1)