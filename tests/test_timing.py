#!/usr/bin/env python3
"""
Comprehensive timing tests for RegNetAgents.
Tests both LLM-powered and rule-based modes, single and multi-gene analysis.
"""
import asyncio
import json
import time
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

# Standard test genes (CRC biomarker panel from manuscript)
MULTI_GENE_PANEL = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']


def check_ollama_available():
    """Check if Ollama service is running."""
    try:
        urllib.request.urlopen('http://localhost:11434/api/tags', timeout=2)
        return True
    except:
        return False


async def test_single_gene_focused_rules():
    """Test single gene with focused analysis (rule-based)."""
    print("\n" + "=" * 70)
    print("TEST 1: Single Gene Focused (Rule-Based)")
    print("=" * 70)

    os.environ['USE_LLM_AGENTS'] = 'false'
    workflow = RegNetAgentsWorkflow()

    start_time = time.time()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="focused"
    )
    execution_time = time.time() - start_time

    print(f"Execution Time: {execution_time:.3f} seconds")
    return execution_time


async def test_single_gene_comprehensive_rules():
    """Test single gene comprehensive (rule-based)."""
    print("\n" + "=" * 70)
    print("TEST 2: Single Gene Comprehensive (Rule-Based)")
    print("=" * 70)

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
    """Test single gene comprehensive (LLM-powered)."""
    print("\n" + "=" * 70)
    print("TEST 3: Single Gene Comprehensive (LLM-Powered)")
    print("=" * 70)

    if not check_ollama_available():
        print("Ollama not available - skipping")
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


async def test_multi_gene_parallel_rules():
    """Test 5 genes parallel (rule-based)."""
    print("\n" + "=" * 70)
    print("TEST 4: Multi-Gene Parallel (Rule-Based)")
    print("=" * 70)
    print(f"Genes: {', '.join(MULTI_GENE_PANEL)}")

    os.environ['USE_LLM_AGENTS'] = 'false'
    workflow = RegNetAgentsWorkflow()

    start_time = time.time()
    tasks = [
        workflow.run_analysis(
            gene=gene,
            cell_type="epithelial_cell",
            analysis_depth="comprehensive"
        )
        for gene in MULTI_GENE_PANEL
    ]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    execution_time = time.time() - start_time

    success_count = sum(1 for r in results if not isinstance(r, Exception))
    print(f"{success_count}/{len(MULTI_GENE_PANEL)} genes completed")
    print(f"Execution Time: {execution_time:.3f} seconds")
    return execution_time


async def test_multi_gene_parallel_llm():
    """Test 5 genes parallel (LLM-powered)."""
    print("\n" + "=" * 70)
    print("TEST 5: Multi-Gene Parallel (LLM-Powered)")
    print("=" * 70)
    print(f"Genes: {', '.join(MULTI_GENE_PANEL)}")

    if not check_ollama_available():
        print("Ollama not available - skipping")
        return None

    os.environ['USE_LLM_AGENTS'] = 'true'
    workflow = RegNetAgentsWorkflow()

    start_time = time.time()
    tasks = [
        workflow.run_analysis(
            gene=gene,
            cell_type="epithelial_cell",
            analysis_depth="comprehensive"
        )
        for gene in MULTI_GENE_PANEL
    ]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    execution_time = time.time() - start_time

    success_count = sum(1 for r in results if not isinstance(r, Exception))
    print(f"{success_count}/{len(MULTI_GENE_PANEL)} genes completed")
    print(f"Execution Time: {execution_time:.3f} seconds")
    return execution_time


async def test_multi_gene_sequential_llm():
    """Test 5 genes sequential (LLM-powered) - manuscript verification.

    This test runs genes sequentially to verify the ~62 second manuscript claim.
    """
    print("\n" + "=" * 70)
    print("TEST 6: Multi-Gene Sequential (LLM-Powered) - Manuscript Verification")
    print("=" * 70)
    print(f"Genes: {', '.join(MULTI_GENE_PANEL)}")

    if not check_ollama_available():
        print("Ollama not available - skipping")
        return None, []

    os.environ['USE_LLM_AGENTS'] = 'true'
    workflow = RegNetAgentsWorkflow()

    per_gene_results = []
    start_time = time.time()

    for i, gene in enumerate(MULTI_GENE_PANEL, 1):
        print(f"  [{i}/5] Analyzing {gene}...", end=" ")
        gene_start = time.time()

        try:
            result = await workflow.run_analysis(
                gene=gene,
                cell_type="epithelial_cell",
                analysis_depth="comprehensive"
            )
            gene_time = time.time() - gene_start
            per_gene_results.append({'gene': gene, 'time': gene_time, 'success': True})
            print(f"{gene_time:.2f} sec")
        except Exception as e:
            gene_time = time.time() - gene_start
            per_gene_results.append({'gene': gene, 'time': gene_time, 'success': False})
            print(f"FAILED ({e})")

    execution_time = time.time() - start_time

    print(f"\nTotal: {execution_time:.2f} seconds (manuscript claim: ~62 sec)")
    if abs(execution_time - 62) < 20:
        print("MATCH - within acceptable range")
    else:
        print(f"DISCREPANCY - {execution_time - 62:+.1f} seconds difference")

    return execution_time, per_gene_results


async def main():
    print("=" * 70)
    print("RegNetAgents - Comprehensive Timing Tests")
    print("=" * 70)

    results = {}

    # Rule-based tests (always available)
    results['single_focused_rules'] = await test_single_gene_focused_rules()
    results['single_comprehensive_rules'] = await test_single_gene_comprehensive_rules()
    results['multi_parallel_rules'] = await test_multi_gene_parallel_rules()

    # LLM tests (if Ollama available)
    results['single_comprehensive_llm'] = await test_single_gene_comprehensive_llm()
    results['multi_parallel_llm'] = await test_multi_gene_parallel_llm()
    seq_time, seq_details = await test_multi_gene_sequential_llm()
    results['multi_sequential_llm'] = seq_time

    # Summary
    print("\n" + "=" * 70)
    print("TIMING SUMMARY")
    print("=" * 70)

    print("\nRule-Based Mode:")
    print(f"  Single gene focused:        {results['single_focused_rules']:.3f} sec")
    print(f"  Single gene comprehensive:  {results['single_comprehensive_rules']:.3f} sec")
    print(f"  Multi-gene (5 parallel):    {results['multi_parallel_rules']:.3f} sec")

    if results['single_comprehensive_llm'] is not None:
        print("\nLLM-Powered Mode:")
        print(f"  Single gene comprehensive:  {results['single_comprehensive_llm']:.3f} sec")
        print(f"  Multi-gene (5 parallel):    {results['multi_parallel_llm']:.3f} sec")
        print(f"  Multi-gene (5 sequential):  {results['multi_sequential_llm']:.3f} sec")
    else:
        print("\nLLM-Powered Mode: Skipped (Ollama not running)")

    print("=" * 70)

    # Save results
    output = {
        "rule_based": {
            "single_focused_sec": round(results['single_focused_rules'], 3),
            "single_comprehensive_sec": round(results['single_comprehensive_rules'], 3),
            "multi_parallel_sec": round(results['multi_parallel_rules'], 3)
        }
    }

    if results['single_comprehensive_llm'] is not None:
        output["llm_powered"] = {
            "single_comprehensive_sec": round(results['single_comprehensive_llm'], 3),
            "multi_parallel_sec": round(results['multi_parallel_llm'], 3),
            "multi_sequential_sec": round(results['multi_sequential_llm'], 3)
        }

    with open("timing_results.json", "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to timing_results.json")


if __name__ == "__main__":
    asyncio.run(main())
