#!/usr/bin/env python3
"""
Tests for MCP Resource endpoints.

Tests the list_resources, list_resource_templates, and read_resource handlers
to verify cell-type listing, network summaries, and gene lookups.
"""

import asyncio
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pydantic import AnyUrl
from regnetagents_langgraph_mcp_server import (
    handle_list_resources,
    handle_list_resource_templates,
    handle_read_resource,
)


async def run_tests():
    print("=" * 70)
    print("Testing MCP Resources")
    print("=" * 70)

    passed = 0
    failed = 0

    # ---- Test 1: list_resources returns all 10 cell types + cell-types resource ----
    print("\n[TEST 1] list_resources returns static + per-cell-type resources")
    resources = await handle_list_resources()
    cell_type_uris = [str(r.uri) for r in resources if "network/" in str(r.uri)]
    has_cell_types_resource = any(str(r.uri) == "regnetagents://cell-types" for r in resources)
    if has_cell_types_resource and len(cell_type_uris) == 10:
        print(f"  [OK] Found cell-types resource + {len(cell_type_uris)} network resources")
        passed += 1
    else:
        print(f"  [FAIL] has_cell_types={has_cell_types_resource}, network_count={len(cell_type_uris)}")
        failed += 1

    # ---- Test 2: list_resource_templates returns gene template ----
    print("\n[TEST 2] list_resource_templates returns gene lookup template")
    templates = await handle_list_resource_templates()
    has_gene_template = any("gene" in t.uriTemplate for t in templates)
    if has_gene_template and len(templates) >= 1:
        print(f"  [OK] Found gene lookup template: {templates[0].uriTemplate}")
        passed += 1
    else:
        print(f"  [FAIL] templates={templates}")
        failed += 1

    # ---- Test 3: read cell-types resource ----
    print("\n[TEST 3] read_resource for regnetagents://cell-types")
    result = await handle_read_resource(AnyUrl("regnetagents://cell-types"))
    data = json.loads(result[0].content)
    cell_types = data.get("cell_types", [])
    if len(cell_types) == 10 and all("num_genes" in ct and "num_edges" in ct for ct in cell_types):
        print(f"  [OK] {len(cell_types)} cell types, all with gene/edge counts")
        for ct in cell_types:
            print(f"       {ct['cell_type']}: {ct['num_genes']} genes, {ct['num_edges']} edges")
        passed += 1
    else:
        print(f"  [FAIL] Got {len(cell_types)} cell types: {json.dumps(data, indent=2)}")
        failed += 1

    # ---- Test 4: read network summary for epithelial_cell ----
    print("\n[TEST 4] read_resource for regnetagents://network/epithelial_cell")
    result = await handle_read_resource(AnyUrl("regnetagents://network/epithelial_cell"))
    data = json.loads(result[0].content)
    expected_keys = {"cell_type", "num_genes", "num_edges", "density", "top_regulators_by_pagerank", "cache_version"}
    if expected_keys.issubset(data.keys()) and data["num_genes"] > 0:
        print(f"  [OK] {data['num_genes']} genes, {data['num_edges']} edges, density={data['density']}")
        print(f"       Top regulator: {data['top_regulators_by_pagerank'][0]}")
        passed += 1
    else:
        print(f"  [FAIL] Missing keys or empty: {json.dumps(data, indent=2)}")
        failed += 1

    # ---- Test 5: read gene resource for TP53 (known gene) ----
    print("\n[TEST 5] read_resource for regnetagents://gene/TP53")
    result = await handle_read_resource(AnyUrl("regnetagents://gene/TP53"))
    data = json.loads(result[0].content)
    if data.get("found") is True and data.get("present_in_cell_types", 0) > 0:
        print(f"  [OK] TP53 found in {data['present_in_cell_types']} cell types")
        for detail in data.get("cell_type_details", [])[:3]:
            print(f"       {detail['cell_type']}: {detail['num_regulators']} regulators, {detail['num_targets']} targets")
        passed += 1
    else:
        print(f"  [FAIL] {json.dumps(data, indent=2)}")
        failed += 1

    # ---- Test 6: read gene resource for invalid gene ----
    print("\n[TEST 6] read_resource for regnetagents://gene/FAKEGENE999")
    result = await handle_read_resource(AnyUrl("regnetagents://gene/FAKEGENE999"))
    data = json.loads(result[0].content)
    if data.get("found") is False:
        print(f"  [OK] found=False, message={data.get('message')}")
        passed += 1
    else:
        print(f"  [FAIL] Expected found=False, got: {json.dumps(data, indent=2)}")
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
