#!/usr/bin/env python3
"""
Tests for the validate_gene MCP tool.

Tests gene validation against the network cache, including:
- Valid gene lookup with stats
- Invalid/unknown gene handling
- Misspelled gene fuzzy matching suggestions
- Case insensitivity
"""

import asyncio
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_mcp_server import get_workflow


async def run_tests():
    print("=" * 70)
    print("Testing validate_gene Tool")
    print("=" * 70)

    workflow = await get_workflow()
    agent = workflow.modeling_agent
    passed = 0
    failed = 0

    # Test 1: Valid gene (TP53) → found=True with stats
    print("\n[TEST 1] Valid gene: TP53")
    result = agent.validate_gene("TP53", "epithelial_cell")
    if result["found"] is True and result["gene"] == "TP53" and "quick_stats" in result:
        stats = result["quick_stats"]
        print(f"  [OK] found=True, ensembl_id={result['ensembl_id']}")
        print(f"       regulators={stats['num_regulators']}, targets={stats['num_targets']}, role={stats['regulatory_role']}")
        passed += 1
    else:
        print(f"  [FAIL] Expected found=True with stats, got: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 2: Invalid gene (FAKEGENE) → found=False
    print("\n[TEST 2] Invalid gene: FAKEGENE")
    result = agent.validate_gene("FAKEGENE", "epithelial_cell")
    if result["found"] is False:
        print(f"  [OK] found=False, message={result['message']}")
        passed += 1
    else:
        print(f"  [FAIL] Expected found=False, got: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 3: Misspelled gene (TP5) → found=False with suggestions
    print("\n[TEST 3] Misspelled gene: TP5")
    result = agent.validate_gene("TP5", "epithelial_cell")
    suggestions = result.get("suggestions", [])
    if result["found"] is False and len(suggestions) > 0:
        print(f"  [OK] found=False, suggestions={suggestions}")
        if "TP53" in suggestions:
            print(f"       TP53 is among suggestions")
        passed += 1
    else:
        print(f"  [FAIL] Expected found=False with suggestions, got: {json.dumps(result, indent=2)}")
        failed += 1

    # Test 4: Case insensitivity (tp53 → same as TP53)
    print("\n[TEST 4] Case insensitivity: tp53")
    result_lower = agent.validate_gene("tp53", "epithelial_cell")
    result_upper = agent.validate_gene("TP53", "epithelial_cell")
    if result_lower["found"] == result_upper["found"] and result_lower.get("ensembl_id") == result_upper.get("ensembl_id"):
        print(f"  [OK] tp53 and TP53 resolve identically (found={result_lower['found']})")
        passed += 1
    else:
        print(f"  [FAIL] Results differ: tp53={json.dumps(result_lower)} vs TP53={json.dumps(result_upper)}")
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
