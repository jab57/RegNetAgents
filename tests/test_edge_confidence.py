#!/usr/bin/env python3
"""
Tests for ARACNe edge confidence filtering in query_network (issue #12).

Verifies:
- "all" returns all edges (baseline behavior unchanged)
- "medium" (MI > 0.05) returns a subset of edges
- "high" (MI > 0.1 AND bootstrap_count >= 3) returns fewer edges than "medium"
- gene_neighbors includes mi_score and bootstrap_count per edge
- high-confidence edges appear in the medium set (monotonicity)
- Invalid confidence_level returns an error
- network_stats includes has_edge_confidence flag
"""

from regnetagents_langgraph_mcp_server import get_workflow

CELL_TYPE = "epithelial_cell"
GENE = "TP53"


async def test_confidence_all_is_default_behavior():
    """confidence_level='all' produces same num_targets as original (no filter)."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result_all = agent.query_network("top_regulators", CELL_TYPE, confidence_level="all")
    result_default = agent.query_network("top_regulators", CELL_TYPE)
    assert result_all["results"][0]["num_targets"] == result_default["results"][0]["num_targets"]


async def test_high_confidence_fewer_targets_than_all():
    """'high' confidence returns <= targets compared to 'all'."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result_all = agent.query_network("top_regulators", CELL_TYPE, top_n=1, confidence_level="all")
    result_high = agent.query_network("top_regulators", CELL_TYPE, top_n=1, confidence_level="high")
    assert result_all["results"][0]["num_targets"] >= result_high["results"][0]["num_targets"]


async def test_medium_between_all_and_high():
    """'medium' confidence returns edge count between 'all' and 'high'."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    r_all = agent.query_network("top_regulators", CELL_TYPE, top_n=1, confidence_level="all")
    r_med = agent.query_network("top_regulators", CELL_TYPE, top_n=1, confidence_level="medium")
    r_high = agent.query_network("top_regulators", CELL_TYPE, top_n=1, confidence_level="high")
    all_count = r_all["results"][0]["num_targets"]
    med_count = r_med["results"][0]["num_targets"]
    high_count = r_high["results"][0]["num_targets"]
    assert high_count <= med_count <= all_count


async def test_gene_neighbors_includes_mi_score():
    """gene_neighbors targets and regulators include mi_score when cache has MI data."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("gene_neighbors", CELL_TYPE, gene=GENE)
    # Cache must be version 3 for this test to be meaningful
    assert result.get("query_type") == "gene_neighbors"
    targets = result.get("targets", [])
    if targets:
        assert "mi_score" in targets[0], "target entries should include mi_score"
        assert "bootstrap_count" in targets[0], "target entries should include bootstrap_count"
        assert isinstance(targets[0]["mi_score"], float)
        assert isinstance(targets[0]["bootstrap_count"], int)


async def test_gene_neighbors_high_confidence_subset():
    """High-confidence gene_neighbors returns fewer or equal edges than 'all'."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    r_all = agent.query_network("gene_neighbors", CELL_TYPE, gene=GENE, confidence_level="all")
    r_high = agent.query_network("gene_neighbors", CELL_TYPE, gene=GENE, confidence_level="high")
    assert r_all["num_targets"] >= r_high["num_targets"]
    assert r_all["num_regulators"] >= r_high["num_regulators"]


async def test_gene_neighbors_high_confidence_mi_threshold():
    """All edges returned by 'high' confidence have MI > 0.1."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("gene_neighbors", CELL_TYPE, gene=GENE, confidence_level="high")
    for edge in result.get("targets", []):
        assert edge.get("mi_score", 0.0) > 0.1, \
            f"High-confidence target edge has MI <= 0.1: {edge}"
    for edge in result.get("regulators", []):
        assert edge.get("mi_score", 0.0) > 0.1, \
            f"High-confidence regulator edge has MI <= 0.1: {edge}"


async def test_confidence_level_in_response():
    """Response includes confidence_level field matching what was requested."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    for level in ("all", "medium", "high"):
        result = agent.query_network("top_regulators", CELL_TYPE, confidence_level=level)
        assert result.get("confidence_level") == level, \
            f"Expected confidence_level='{level}' in response"


async def test_invalid_confidence_level_returns_error():
    """Invalid confidence_level value returns an error dict."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("top_regulators", CELL_TYPE, confidence_level="ultra")
    assert result.get("error") is True
    assert "confidence_level" in result.get("message", "").lower()


async def test_network_stats_has_edge_confidence_flag():
    """network_stats includes has_edge_confidence=True when cache is version 3."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.query_network("network_stats", CELL_TYPE)
    assert "has_edge_confidence" in result
    assert result["has_edge_confidence"] is True
