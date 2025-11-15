#!/usr/bin/env python3
"""
Comprehensive timing test for all poster scenarios
Tests both LLM-powered and rule-based modes
"""
import asyncio
import json
import time
import os
import sys

# Add parent directory to path to import the workflow
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def test_single_gene_focused_rules():
    """Test single gene with focused analysis (rule-based)"""
    print("\n" + "=" * 80)
    print("TEST 1: Single Gene Focused Analysis (Rule-Based)")
    print("=" * 80)

    # Force rule-based mode
    os.environ['USE_LLM_AGENTS'] = 'false'

    workflow = RegNetAgentsWorkflow()
    start_time = time.time()

    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="focused"  # Skip domain agents
    )

    execution_time = time.time() - start_time
    print(f"Execution Time: {execution_time:.3f} seconds")

    return execution_time

async def test_single_gene_comprehensive_rules():
    """Test single gene comprehensive (rule-based)"""
    print("\n" + "=" * 80)
    print("TEST 2: Single Gene Comprehensive (Rule-Based)")
    print("=" * 80)

    os.environ['USE_LLM_AGENTS'] = 'false'

    workflow = RegNetAgentsWorkflow()
    start_time = time.time()

    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    execution_time = time.time() - start_time
    print(f"Execution Time: {execution_time:.3f} seconds")

    return execution_time

async def test_single_gene_comprehensive_llm():
    """Test single gene comprehensive (LLM-powered)"""
    print("\n" + "=" * 80)
    print("TEST 3: Single Gene Comprehensive (LLM-Powered)")
    print("=" * 80)

    # Check if Ollama is available
    import urllib.request
    try:
        urllib.request.urlopen('http://localhost:11434/api/tags', timeout=2)
        ollama_available = True
    except:
        ollama_available = False

    if not ollama_available:
        print("Ollama not available - skipping LLM test")
        print("To run LLM tests: Start Ollama service")
        return None

    os.environ['USE_LLM_AGENTS'] = 'true'

    workflow = RegNetAgentsWorkflow()
    start_time = time.time()

    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    execution_time = time.time() - start_time
    print(f"Execution Time: {execution_time:.3f} seconds")

    return execution_time

async def test_multi_gene_rules():
    """Test 5 genes parallel (rule-based)"""
    print("\n" + "=" * 80)
    print("TEST 4: Multi-Gene (5 genes parallel, Rule-Based)")
    print("=" * 80)

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']
    print(f"Analyzing {len(genes)} genes: {', '.join(genes)}")

    os.environ['USE_LLM_AGENTS'] = 'false'

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
    print(f"{success_count}/{len(genes)} genes completed")
    print(f"Execution Time: {execution_time:.3f} seconds")

    return execution_time

async def test_multi_gene_llm():
    """Test 5 genes parallel (LLM-powered)"""
    print("\n" + "=" * 80)
    print("TEST 5: Multi-Gene (5 genes parallel, LLM-Powered)")
    print("=" * 80)

    # Check if Ollama is available
    import urllib.request
    try:
        urllib.request.urlopen('http://localhost:11434/api/tags', timeout=2)
        ollama_available = True
    except:
        ollama_available = False

    if not ollama_available:
        print("Ollama not available - skipping LLM test")
        return None

    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']
    print(f"Analyzing {len(genes)} genes: {', '.join(genes)}")

    os.environ['USE_LLM_AGENTS'] = 'true'

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
    print(f"{success_count}/{len(genes)} genes completed")
    print(f"Execution Time: {execution_time:.3f} seconds")

    return execution_time

async def main():
    print("=" * 80)
    print("RegNetAgents - Comprehensive Timing Tests for Conference Poster")
    print("=" * 80)

    # Run all tests
    results = {}

    # Rule-based tests (always available)
    results['single_focused_rules'] = await test_single_gene_focused_rules()
    results['single_comprehensive_rules'] = await test_single_gene_comprehensive_rules()
    results['multi_5genes_rules'] = await test_multi_gene_rules()

    # LLM tests (if Ollama available)
    results['single_comprehensive_llm'] = await test_single_gene_comprehensive_llm()
    results['multi_5genes_llm'] = await test_multi_gene_llm()

    # Summary
    print("\n" + "=" * 80)
    print("TIMING SUMMARY FOR POSTER")
    print("=" * 80)
    print("\nRule-Based Mode:")
    print(f"  Single gene focused:        {results['single_focused_rules']:.3f} sec")
    print(f"  Single gene comprehensive:  {results['single_comprehensive_rules']:.3f} sec")
    print(f"  Multi-gene (5 parallel):    {results['multi_5genes_rules']:.3f} sec")

    if results['single_comprehensive_llm'] is not None:
        print("\nLLM-Powered Mode:")
        print(f"  Single gene comprehensive:  {results['single_comprehensive_llm']:.3f} sec")
        print(f"  Multi-gene (5 parallel):    {results['multi_5genes_llm']:.3f} sec")
    else:
        print("\nLLM-Powered Mode: Not tested (Ollama not running)")

    print("=" * 80)

    # Save results
    output = {
        "rule_based": {
            "single_focused_sec": round(results['single_focused_rules'], 3),
            "single_comprehensive_sec": round(results['single_comprehensive_rules'], 3),
            "multi_5genes_sec": round(results['multi_5genes_rules'], 3)
        }
    }

    if results['single_comprehensive_llm'] is not None:
        output["llm_powered"] = {
            "single_comprehensive_sec": round(results['single_comprehensive_llm'], 3),
            "multi_5genes_sec": round(results['multi_5genes_llm'], 3)
        }

    with open("comprehensive_timing_results.json", "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to comprehensive_timing_results.json")

    # Show poster-ready formatting
    print("\n" + "=" * 80)
    print("POSTER-READY FORMATTING")
    print("=" * 80)
    print("\nPerformance Comparison table values:")
    print(f"  Single gene (rule-based):        <{results['single_focused_rules']:.2f} sec")
    print(f"  Single gene (LLM-powered):       {results.get('single_comprehensive_llm', 'N/A') if results.get('single_comprehensive_llm') else 'N/A'} sec")
    print(f"  Multi-gene 5 genes (LLM):        {results.get('multi_5genes_llm', 'N/A') if results.get('multi_5genes_llm') else 'N/A'} sec")
    print(f"  Multi-gene 5 genes (rules):      {results['multi_5genes_rules']:.2f} sec")

if __name__ == "__main__":
    asyncio.run(main())
