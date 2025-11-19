#!/usr/bin/env python3
"""
Test EXACT timing from manuscript:
5-gene colorectal cancer panel (MYC, CTNNB1, CCND1, TP53, KRAS)
LLM-powered mode
"""

import asyncio
import time
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

async def test_manuscript_panel():
    """Test the exact 5-gene panel from manuscript"""

    # EXACT genes from manuscript Table 2
    genes = ['MYC', 'CTNNB1', 'CCND1', 'TP53', 'KRAS']
    cell_type = 'epithelial_cell'

    print("=" * 70)
    print("MANUSCRIPT TIMING VERIFICATION")
    print("=" * 70)
    print(f"Testing 5-gene colorectal cancer panel: {', '.join(genes)}")
    print(f"Cell type: {cell_type}")
    print(f"Mode: LLM-powered (comprehensive)")
    print("=" * 70)
    print()

    # Initialize workflow
    workflow = RegNetAgentsWorkflow()

    # Check Ollama status
    print(f"[OK] Workflow initialized")
    print()

    results = []
    start_time = time.time()

    for i, gene in enumerate(genes, 1):
        print(f"[{i}/5] Analyzing {gene}...")
        gene_start = time.time()

        try:
            result = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                analysis_depth='comprehensive'  # Full analysis with LLM
            )

            gene_time = time.time() - gene_start
            results.append({
                'gene': gene,
                'time': gene_time,
                'success': True
            })
            print(f"  [OK] Complete in {gene_time:.2f} sec")

        except Exception as e:
            gene_time = time.time() - gene_start
            results.append({
                'gene': gene,
                'time': gene_time,
                'success': False,
                'error': str(e)
            })
            print(f"  [FAIL] Error: {e}")

    total_time = time.time() - start_time

    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Average per gene: {total_time/5:.2f} seconds")
    print()

    print("Per-gene breakdown:")
    for r in results:
        status = "[OK]" if r['success'] else "[FAIL]"
        print(f"  {status} {r['gene']:10s} - {r['time']:6.2f} sec")

    print()
    print("=" * 70)
    print("COMPARISON TO MANUSCRIPT")
    print("=" * 70)
    print(f"Manuscript claim: ~62 seconds")
    print(f"Test result:      {total_time:.2f} seconds")

    if abs(total_time - 62) < 20:  # Within 20 seconds
        print(f"[OK] MATCH - within acceptable range")
    else:
        diff = total_time - 62
        print(f"[WARN] DISCREPANCY - {diff:+.1f} seconds difference")

    print("=" * 70)

    return total_time, results

if __name__ == "__main__":
    total, results = asyncio.run(test_manuscript_panel())
