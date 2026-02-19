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

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_mcp_server import get_workflow


async def test_query_network_top_regulators():
    """top_regulators → non-empty list sorted by num_targets descending."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("top_regulators", "epithelial_cell")
    assert result.get("query_type") == "top_regulators"
    assert isinstance(result.get("results"), list)
    assert len(result["results"]) > 0
    assert "gene" in result["results"][0]
    assert "num_targets" in result["results"][0]
    counts = [r["num_targets"] for r in result["results"]]
    assert all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1)), \
        "Results not sorted by num_targets descending"


async def test_query_network_top_targets():
    """top_targets → non-empty list sorted by num_regulators descending."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("top_targets", "epithelial_cell")
    assert result.get("query_type") == "top_targets"
    assert isinstance(result.get("results"), list)
    assert len(result["results"]) > 0
    assert "num_regulators" in result["results"][0]
    counts = [r["num_regulators"] for r in result["results"]]
    assert all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1)), \
        "Results not sorted by num_regulators descending"


async def test_query_network_gene_neighbors():
    """gene_neighbors for TP53 → returns regulators and targets lists."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("gene_neighbors", "epithelial_cell", gene="TP53")
    assert result.get("query_type") == "gene_neighbors"
    assert result.get("gene") == "TP53"
    assert isinstance(result.get("regulators"), list)
    assert isinstance(result.get("targets"), list)
    assert result["num_regulators"] > 0 or result["num_targets"] > 0


async def test_query_network_invalid_gene():
    """gene_neighbors for FAKEGENE → graceful error response."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("gene_neighbors", "epithelial_cell", gene="FAKEGENE")
    assert result.get("error") is True
    assert "not found" in result.get("message", "").lower()


async def test_query_network_stats():
    """network_stats → returns num_genes, num_edges, density."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("network_stats", "epithelial_cell")
    assert result.get("query_type") == "network_stats"
    assert "num_genes" in result and result["num_genes"] > 0
    assert "num_edges" in result and result["num_edges"] > 0
    assert "density" in result


async def test_query_network_top_n():
    """top_n=3 → returns exactly 3 results."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("top_regulators", "epithelial_cell", top_n=3)
    assert len(result.get("results", [])) == 3
