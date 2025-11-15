#!/usr/bin/env python3
"""
Test script to get actual timing for both case studies
"""
import asyncio
import json
import time
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def test_single_gene():
    """Test Case Study 2: Single Gene Deep Dive - TP53"""
    print("\n" + "=" * 80)
    print("CASE STUDY 2: Single Gene Deep Dive - TP53")
    print("=" * 80)

    workflow = RegNetAgentsWorkflow()

    start_time = time.time()

    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    execution_time = time.time() - start_time

    print(f"\nTP53 Comprehensive Analysis Complete")
    print(f"Execution Time: {execution_time:.2f} seconds")

    return execution_time

async def test_multi_gene():
    """Test Case Study 1: Multi-Gene Analysis - CRC Biomarkers"""
    print("\n" + "=" * 80)
    print("CASE STUDY 1: Multi-Gene Analysis - CRC Biomarkers")
    print("=" * 80)

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']
    print(f"Analyzing {len(genes)} genes in parallel: {', '.join(genes)}")

    workflow = RegNetAgentsWorkflow()

    start_time = time.time()

    tasks = [
        workflow.run_analysis(
            gene=gene,
            cell_type="epithelial_cell",
            analysis_depth="comprehensive"
        )
        for gene in genes
    ]

    results = await asyncio.gather(*tasks, return_exceptions=True)

    execution_time = time.time() - start_time

    success_count = sum(1 for r in results if not isinstance(r, Exception))

    print(f"\n{success_count}/{len(genes)} Genes Analyzed Successfully")
    print(f"Execution Time: {execution_time:.2f} seconds")

    return execution_time

async def main():
    print("=" * 80)
    print("RegNetAgents Conference Poster - Timing Tests")
    print("=" * 80)

    # Test both case studies
    single_time = await test_single_gene()
    multi_time = await test_multi_gene()

    # Summary
    print("\n" + "=" * 80)
    print("TIMING SUMMARY FOR POSTER")
    print("=" * 80)
    print(f"Single Gene (TP53 comprehensive):  {single_time:.2f} seconds")
    print(f"Multi-Gene (5 genes parallel):     {multi_time:.2f} seconds")
    print("=" * 80)

    # Save results
    results = {
        "single_gene_tp53_seconds": round(single_time, 2),
        "multi_gene_5_parallel_seconds": round(multi_time, 2)
    }

    with open("poster_timing_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to poster_timing_results.json")

if __name__ == "__main__":
    asyncio.run(main())
