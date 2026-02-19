#!/usr/bin/env python3
"""
Tests for Reactome pathway enrichment integration.

Requires an internet connection for live Reactome API calls.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import PathwayEnricherAgent, RegNetAgentsWorkflow


async def test_reactome_agent():
    """PathwayEnricherAgent → returns enrichment result with recognized genes and pathways."""
    agent = PathwayEnricherAgent()
    result = await agent.enrich_pathways_reactome(["TP53", "APC", "BRCA1", "MYC"])

    assert result.get("status") in ("success", "partial", "error"), \
        f"Unexpected status: {result.get('status')}"
    assert "genes_analyzed" in result
    assert "summary" in result
    summary = result["summary"]
    assert "total_pathways" in summary
    assert "significant_pathways" in summary

    if result.get("status") == "success":
        assert result.get("genes_recognized", 0) > 0, \
            "Expected at least one gene recognized by Reactome"
        assert summary["total_pathways"] > 0, \
            "Expected at least one pathway returned for known cancer genes"


async def test_workflow_integration():
    """Full workflow for TP53 (basic depth) → pathway_enrichment field is populated."""
    workflow = RegNetAgentsWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="basic"
    )

    assert isinstance(result, dict), "Result must be a dict"
    assert result.get("status") != "error", f"Workflow failed: {result.get('error')}"

    pathway_enrichment = result.get("pathway_enrichment") or {}
    assert "status" in pathway_enrichment, "pathway_enrichment missing status field"
    assert "summary" in pathway_enrichment, "pathway_enrichment missing summary field"
    summary = pathway_enrichment["summary"]
    assert "total_pathways" in summary
