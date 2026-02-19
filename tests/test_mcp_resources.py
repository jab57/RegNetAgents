#!/usr/bin/env python3
"""
Tests for MCP Resource endpoints.

Tests the list_resources, list_resource_templates, and read_resource handlers
to verify cell-type listing, network summaries, and gene lookups.
"""

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


async def test_list_resources_has_all_cell_types():
    """list_resources → 1 cell-types resource + 10 per-network resources."""
    resources = await handle_list_resources()
    cell_type_uris = [str(r.uri) for r in resources if "network/" in str(r.uri)]
    has_cell_types = any(str(r.uri) == "regnetagents://cell-types" for r in resources)
    assert has_cell_types, "Missing regnetagents://cell-types resource"
    assert len(cell_type_uris) == 10, f"Expected 10 network resources, got {len(cell_type_uris)}"


async def test_list_resource_templates_has_gene_template():
    """list_resource_templates → includes a gene lookup URI template."""
    templates = await handle_list_resource_templates()
    assert len(templates) >= 1
    assert any("gene" in t.uriTemplate for t in templates)


async def test_read_cell_types_resource():
    """regnetagents://cell-types → 10 cell types each with num_genes and num_edges."""
    result = await handle_read_resource(AnyUrl("regnetagents://cell-types"))
    data = json.loads(result[0].content)
    cell_types = data.get("cell_types", [])
    assert len(cell_types) == 10, f"Expected 10 cell types, got {len(cell_types)}"
    for ct in cell_types:
        assert "num_genes" in ct and "num_edges" in ct


async def test_read_network_summary():
    """regnetagents://network/epithelial_cell → valid summary with expected fields."""
    result = await handle_read_resource(AnyUrl("regnetagents://network/epithelial_cell"))
    data = json.loads(result[0].content)
    for key in ("cell_type", "num_genes", "num_edges", "density", "top_regulators_by_pagerank"):
        assert key in data, f"Missing key: {key}"
    assert data["num_genes"] > 0
    assert data["num_edges"] > 0


async def test_read_gene_resource_valid():
    """regnetagents://gene/TP53 → found=True in at least one cell type."""
    result = await handle_read_resource(AnyUrl("regnetagents://gene/TP53"))
    data = json.loads(result[0].content)
    assert data.get("found") is True
    assert data.get("present_in_cell_types", 0) > 0


async def test_read_gene_resource_invalid():
    """regnetagents://gene/FAKEGENE999 → found=False."""
    result = await handle_read_resource(AnyUrl("regnetagents://gene/FAKEGENE999"))
    data = json.loads(result[0].content)
    assert data.get("found") is False
