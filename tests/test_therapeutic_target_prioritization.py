#!/usr/bin/env python3
"""
Tests for the therapeutic target prioritization feature.

Verifies that TP53 comprehensive analysis returns ranked regulators
with the expected structure and PageRank-based scoring.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow


async def test_therapeutic_target_prioritization():
    """TP53 comprehensive analysis → ranked regulators with centrality metrics."""
    workflow = RegNetAgentsWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    assert isinstance(result, dict), "Result must be a dict"
    assert result.get("status") != "error", f"Analysis failed: {result.get('error')}"

    target_prio = result.get("therapeutic_target_prioritization", {})
    assert target_prio, "therapeutic_target_prioritization missing from result"
    assert target_prio.get("status") != "skipped", \
        f"Prioritization was skipped: {target_prio.get('reason')}"
    assert target_prio.get("target_gene") == "TP53"

    ranked = target_prio.get("ranked_regulators", [])
    assert len(ranked) > 0, "Expected at least one ranked regulator for TP53"

    top = ranked[0]
    assert "regulator" in top, "Regulator entry missing 'regulator' field"
    assert "centrality_metrics" in top, "Regulator entry missing 'centrality_metrics'"
    assert "pagerank" in top["centrality_metrics"], "Missing pagerank score"
    assert top["centrality_metrics"]["pagerank"] > 0
