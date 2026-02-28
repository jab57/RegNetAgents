#!/usr/bin/env python3
"""
Tests for the find_master_regulators MCP tool.

Tests reverse-direction ARACNe analysis: given a gene signature, identify
which TFs in the network are most significantly driving that signature.

Test strategy: TP53 is a well-characterised hub regulator with many known
targets. Supplying a list of canonical TP53 targets should rank TP53 (or a
close co-regulator) near the top of the results.
"""

from regnetagents_langgraph_mcp_server import get_workflow

# Canonical TP53 target genes used as test signature
TP53_TARGETS = ["CDKN1A", "MDM2", "BAX", "GADD45A", "PUMA", "BBC3",
                "CCNG1", "SESN2", "DDB2", "POLK"]


async def test_find_master_regulators_returns_results():
    """Known TP53 targets → master_regulators list is non-empty."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(TP53_TARGETS, cell_type="epithelial_cell", top_n=10)
    assert "master_regulators" in result, f"Unexpected result: {result}"
    assert len(result["master_regulators"]) > 0


async def test_find_master_regulators_result_structure():
    """Each master regulator entry has the expected fields."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(TP53_TARGETS, cell_type="epithelial_cell", top_n=5)
    for entry in result["master_regulators"]:
        assert "rank" in entry
        assert "gene" in entry
        assert "ensembl_id" in entry
        assert "regulon_size" in entry
        assert "overlap_count" in entry
        assert "enrichment_score" in entry
        assert "p_value" in entry
        assert "overlapping_genes" in entry
        assert entry["overlap_count"] > 0
        assert 0.0 <= entry["p_value"] <= 1.0
        assert entry["enrichment_score"] > 0.0


async def test_find_master_regulators_query_summary():
    """query_summary fields are present and internally consistent."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(TP53_TARGETS, cell_type="epithelial_cell")
    summary = result["query_summary"]
    assert summary["gene_set_size"] == len(TP53_TARGETS)
    assert summary["genes_found_in_network"] <= summary["gene_set_size"]
    assert summary["genes_found_in_network"] >= 0
    assert summary["network_size"] > 0
    assert summary["cell_type"] == "epithelial_cell"


async def test_find_master_regulators_ranked_by_pvalue():
    """Results are ordered by p-value ascending."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(TP53_TARGETS, cell_type="epithelial_cell", top_n=10)
    p_values = [r["p_value"] for r in result["master_regulators"]]
    assert p_values == sorted(p_values), "Results not sorted by p-value"


async def test_find_master_regulators_unknown_genes():
    """All-unknown gene set → error response with message."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(["FAKEGENE1", "FAKEGENE2"], cell_type="epithelial_cell")
    assert result.get("error") is True
    assert "message" in result


async def test_find_master_regulators_top_n_respected():
    """top_n parameter limits the number of returned regulators."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.find_master_regulators(TP53_TARGETS, cell_type="epithelial_cell", top_n=3)
    assert len(result["master_regulators"]) <= 3
