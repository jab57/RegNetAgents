#!/usr/bin/env python3
"""
Tests for the MCP server workflow integration.

Verifies single-gene and multi-gene analysis via the MCP server's
get_workflow() entry point, covering the fields Claude Desktop receives.
"""

import asyncio

from regnetagents_langgraph_mcp_server import get_workflow


async def test_mcp_server_initializes():
    """get_workflow() returns a valid workflow instance."""
    workflow = await get_workflow()
    assert workflow is not None
    assert hasattr(workflow, "modeling_agent")


async def test_mcp_server_single_gene_analysis():
    """TP53 comprehensive analysis → result contains gene summary and network fields."""
    workflow = await get_workflow()
    result = await workflow.run_analysis("TP53", "epithelial_cell", "comprehensive")

    assert isinstance(result, dict), "Result must be a dict"
    assert result.get("status") != "error", f"Analysis failed: {result.get('error')}"
    assert result.get("gene_analysis_summary", {}).get("gene") == "TP53"

    network = result.get("network_analysis", {})
    assert "regulatory_role" in network
    assert "num_targets" in network
    assert "num_regulators" in network


async def test_mcp_server_multi_gene_analysis():
    """Parallel analysis of MYC, TP53, KRAS → all three succeed."""
    workflow = await get_workflow()
    genes = ["MYC", "TP53", "KRAS"]
    tasks = [workflow.run_analysis(gene, "epithelial_cell", "focused") for gene in genes]
    results = await asyncio.gather(*tasks, return_exceptions=True)

    failures = [r for r in results if isinstance(r, Exception)]
    assert len(failures) == 0, f"Exceptions in multi-gene analysis: {failures}"

    successful = [r for r in results if isinstance(r, dict) and r.get("status") != "error"]
    assert len(successful) == len(genes), \
        f"Expected {len(genes)} successful results, got {len(successful)}"
