#!/usr/bin/env python3
"""
Test script to run APC analysis for colon cancer screening use case
"""
import asyncio
import json
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def main():
    print("=" * 80)
    print("RegNetAgents Analysis: Colon Cancer Screening Use Case")
    print("=" * 80)
    print("\nQuestion: Analyze APC gene regulatory network in epithelial cells")
    print("\nRunning comprehensive analysis...\n")

    # Initialize workflow
    workflow = RegNetAgentsWorkflow()

    # Run analysis
    result = await workflow.run_analysis(
        gene="APC",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    # Pretty print the results
    print("\n" + "=" * 80)
    print("ANALYSIS RESULTS")
    print("=" * 80)
    print(json.dumps(result, indent=2))

    # Save to file for poster
    with open("apc_analysis_output.json", "w") as f:
        json.dumps(result, f, indent=2)

    print("\n" + "=" * 80)
    print("Analysis complete! Output saved to apc_analysis_output.json")
    print("=" * 80)

if __name__ == "__main__":
    asyncio.run(main())
