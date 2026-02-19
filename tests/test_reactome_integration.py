#!/usr/bin/env python3
"""
Tests for Reactome pathway enrichment integration.

Tests verify that the PathwayEnricherAgent returns well-structured results
and that the workflow handles Reactome API availability gracefully.
Note: Live API calls may return empty results if Reactome is unreachable.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import PathwayEnricherAgent, RegNetAgentsWorkflow


async def test_reactome_agent():
    """PathwayEnricherAgent → always returns well-structured result with required fields."""
    agent = PathwayEnricherAgent()
    result = await agent.enrich_pathways_reactome(["TP53", "APC", "BRCA1", "MYC"])

    assert isinstance(result, dict), "Result must be a dict"
    assert "status" in result, "Result missing 'status' field"
    assert "genes_analyzed" in result, "Result missing 'genes_analyzed' field"
    assert "summary" in result, "Result missing 'summary' field"

    summary = result["summary"]
    assert "total_pathways" in summary, "summary missing 'total_pathways'"
    assert "significant_pathways" in summary, "summary missing 'significant_pathways'"


async def test_workflow_integration():
    """Full workflow for TP53 (basic depth) → completes without error; pathway field present."""
    workflow = RegNetAgentsWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="basic"
    )

    assert isinstance(result, dict), "Result must be a dict"
    assert result.get("status") != "error", f"Workflow failed: {result.get('error')}"
    assert "pathway_enrichment" in result, "Result missing 'pathway_enrichment' field"

    # pathway_enrichment may be None if Reactome API is unreachable in CI
    pathway_enrichment = result.get("pathway_enrichment")
    if pathway_enrichment and isinstance(pathway_enrichment, dict):
        assert "summary" in pathway_enrichment, "pathway_enrichment missing 'summary'"
        summary = pathway_enrichment["summary"]
        assert "total_pathways" in summary
        assert "significant_pathways" in summary
