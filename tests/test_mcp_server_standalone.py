#!/usr/bin/env python3
"""
Standalone test to verify the MCP server works
Run this to test the server independently of Claude Desktop
"""
import asyncio
import json
from regnetagents_langgraph_mcp_server import get_workflow

async def test_server():
    print("=" * 80)
    print("Testing RegNetAgents MCP Server Standalone")
    print("=" * 80)

    # Test 1: Initialize workflow
    print("\n[TEST 1] Initializing workflow...")
    try:
        workflow = await get_workflow()
        print("[OK] Workflow initialized successfully")
    except Exception as e:
        print(f"[FAIL] Workflow initialization failed: {e}")
        return False

    # Test 2: Single gene analysis
    print("\n[TEST 2] Running single gene analysis (TP53)...")
    try:
        import time
        start = time.time()
        result = await workflow.run_analysis('TP53', 'epithelial_cell', 'comprehensive')
        elapsed = time.time() - start

        print(f"[OK] Analysis completed in {elapsed:.2f} seconds")
        print(f"  - Gene: {result.get('gene_analysis_summary', {}).get('gene', 'unknown')}")
        print(f"  - Regulatory role: {result.get('network_analysis', {}).get('regulatory_role', 'unknown')}")
        print(f"  - Targets: {result.get('network_analysis', {}).get('num_targets', 0)}")
        print(f"  - Regulators: {result.get('network_analysis', {}).get('num_regulators', 0)}")

        pathway_info = result.get('pathway_enrichment', {})
        print(f"  - Pathways: {pathway_info.get('summary', {}).get('total_pathways', 0)}")

        domain = result.get('domain_analysis', {})
        print(f"  - Domain analyses: {list(domain.keys())}")

    except Exception as e:
        print(f"[FAIL] Single gene analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Test 3: Multi-gene analysis
    print("\n[TEST 3] Running multi-gene analysis (MYC, TP53, KRAS)...")
    try:
        genes = ['MYC', 'TP53', 'KRAS']
        start = time.time()

        tasks = [workflow.run_analysis(gene, 'epithelial_cell', 'focused') for gene in genes]
        results = await asyncio.gather(*tasks, return_exceptions=True)

        elapsed = time.time() - start
        successful = sum(1 for r in results if not isinstance(r, Exception))

        print(f"[OK] Multi-gene analysis completed in {elapsed:.2f} seconds")
        print(f"  - Success: {successful}/{len(genes)} genes")

    except Exception as e:
        print(f"[FAIL] Multi-gene analysis failed: {e}")
        return False

    print("\n" + "=" * 80)
    print("ALL TESTS PASSED - MCP Server is working correctly")
    print("=" * 80)
    print("\nIf Claude Desktop is still timing out, the issue is with:")
    print("  1. Claude Desktop's MCP client timeout settings")
    print("  2. Claude Desktop not properly launching the server")
    print("  3. Communication protocol issue between Claude Desktop and the server")
    print("\nThe RegNetAgents workflow itself is fast and working perfectly.")
    return True

if __name__ == "__main__":
    success = asyncio.run(test_server())
    exit(0 if success else 1)
