#!/usr/bin/env python3
"""
Tests for the RegNetAgents LangGraph workflow.

Verifies end-to-end workflow execution: gene validation, network lookup,
pathway enrichment, domain analysis, and result structure.
"""

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow


async def test_workflow():
    """Full LangGraph workflow for APC in epithelial cells → expected result structure."""
    workflow = RegNetAgentsWorkflow()
    result = await workflow.run_analysis(
        gene="APC",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    assert isinstance(result, dict), "Result must be a dict"
    assert result.get("status") != "error", f"Workflow failed: {result.get('error')}"
    assert "gene_analysis_summary" in result, "Missing gene_analysis_summary"
    assert result["gene_analysis_summary"]["gene"] == "APC"
    assert "regulatory_role" in result["gene_analysis_summary"]

    steps = result.get("workflow_metadata", {}).get("steps_completed", [])
    assert len(steps) > 0, "No workflow steps were recorded as completed"
