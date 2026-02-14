#!/usr/bin/env python3
"""
Tests for the query_network MCP tool.

Tests network queries against the pre-computed cache, including:
- Top regulators ranked by out-degree
- Top targets ranked by in-degree
- Gene neighbors (regulators and targets)
- Network statistics
- Invalid gene handling
- top_n parameter
"""

import asyncio
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_mcp_server import get_workflow


async def run_tests():
    print("=" * 70)
    print("Testing query_network Tool")
    print("=" * 70)

    workflow = await get_workflow()
    agent = workflow.modeling_agent
    passed = 0
    failed = 0

    # Test 1: top_regulators — returns a sorted list with gene symbols
    print("\n[TEST 1] top_regulators")
    result = agent.query_network("top_regulators", "epithelial_cell")
    if (result.get("query_type") == "top_regulators"
            and isinstance(result.get("results"), list)
            and len(result["results"]) > 0
            and "gene" in result["results"][0]
            and "num_targets" in result["results"][0]):
        # Verify sorted descending by num_targets
        counts = [r["num_targets"] for r in result["results"]]
        is_sorted = all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1))
        if is_sorted:
            print(f"  [OK] {len(result['results'])} regulators, top={result['results'][0]['gene']} "
                  f"({result['results'][0]['num_targets']} targets), sorted=True")
            passed += 1
        else:
            print(f"  [FAIL] Results not sorted by num_targets: {counts}")
            failed += 1
    else:
        print(f"  [FAIL] Unexpected result: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 2: top_targets — returns a sorted list by regulator count
    print("\n[TEST 2] top_targets")
    result = agent.query_network("top_targets", "epithelial_cell")
    if (result.get("query_type") == "top_targets"
            and isinstance(result.get("results"), list)
            and len(result["results"]) > 0
            and "num_regulators" in result["results"][0]):
        counts = [r["num_regulators"] for r in result["results"]]
        is_sorted = all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1))
        if is_sorted:
            print(f"  [OK] {len(result['results'])} targets, top={result['results'][0]['gene']} "
                  f"({result['results'][0]['num_regulators']} regulators), sorted=True")
            passed += 1
        else:
            print(f"  [FAIL] Results not sorted by num_regulators: {counts}")
            failed += 1
    else:
        print(f"  [FAIL] Unexpected result: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 3: gene_neighbors for TP53 — returns regulators and targets
    print("\n[TEST 3] gene_neighbors: TP53")
    result = agent.query_network("gene_neighbors", "epithelial_cell", gene="TP53")
    if (result.get("query_type") == "gene_neighbors"
            and result.get("gene") == "TP53"
            and "regulators" in result
            and "targets" in result
            and isinstance(result["regulators"], list)
            and isinstance(result["targets"], list)):
        print(f"  [OK] TP53 has {result['num_regulators']} regulators and {result['num_targets']} targets")
        passed += 1
    else:
        print(f"  [FAIL] Unexpected result: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 4: gene_neighbors with invalid gene — graceful error
    print("\n[TEST 4] gene_neighbors: FAKEGENE (invalid)")
    result = agent.query_network("gene_neighbors", "epithelial_cell", gene="FAKEGENE")
    if result.get("error") is True and "not found" in result.get("message", "").lower():
        print(f"  [OK] error=True, message={result['message']}")
        passed += 1
    else:
        print(f"  [FAIL] Expected error, got: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 5: network_stats — returns num_genes, num_edges, density
    print("\n[TEST 5] network_stats")
    result = agent.query_network("network_stats", "epithelial_cell")
    if (result.get("query_type") == "network_stats"
            and "num_genes" in result
            and "num_edges" in result
            and "density" in result
            and result["num_genes"] > 0
            and result["num_edges"] > 0):
        print(f"  [OK] {result['num_genes']} genes, {result['num_edges']} edges, "
              f"{result['num_regulons']} regulons, density={result['density']}")
        passed += 1
    else:
        print(f"  [FAIL] Unexpected result: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 6: top_n parameter — verify list length respects it
    print("\n[TEST 6] top_n=3 parameter")
    result = agent.query_network("top_regulators", "epithelial_cell", top_n=3)
    if len(result.get("results", [])) == 3:
        print(f"  [OK] top_n=3 returned exactly 3 results")
        passed += 1
    else:
        print(f"  [FAIL] Expected 3 results, got {len(result.get('results', []))}")
        failed += 1

    # Summary
    print("\n" + "=" * 70)
    total = passed + failed
    print(f"Results: {passed}/{total} passed, {failed}/{total} failed")
    print("=" * 70)

    return failed == 0


if __name__ == "__main__":
    success = asyncio.run(run_tests())
    sys.exit(0 if success else 1)
