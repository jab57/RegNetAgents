#!/usr/bin/env python3
"""
Tests for the validate_gene MCP tool.

Tests gene validation against the network cache, including:
- Valid gene lookup with stats
- Invalid/unknown gene handling
- Misspelled gene fuzzy matching suggestions
- Case insensitivity
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_mcp_server import get_workflow


async def test_validate_gene_valid():
    """Valid gene (TP53) → found=True with quick_stats."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.validate_gene("TP53", "epithelial_cell")
    assert result["found"] is True
    assert result["gene"] == "TP53"
    assert "quick_stats" in result
    stats = result["quick_stats"]
    assert "num_regulators" in stats
    assert "num_targets" in stats
    assert "regulatory_role" in stats


async def test_validate_gene_invalid():
    """Invalid gene (FAKEGENE) → found=False with a message."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.validate_gene("FAKEGENE", "epithelial_cell")
    assert result["found"] is False
    assert "message" in result


async def test_validate_gene_fuzzy_suggestions():
    """Misspelled gene (TP5) → found=False with TP53 among suggestions."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result = agent.validate_gene("TP5", "epithelial_cell")
    assert result["found"] is False
    suggestions = result.get("suggestions", [])
    assert len(suggestions) > 0, "Expected fuzzy suggestions for close misspelling"
    assert "TP53" in suggestions


async def test_validate_gene_case_insensitive():
    """Case insensitivity: tp53 resolves identically to TP53."""
    workflow = await get_workflow()
    agent = workflow.modeling_agent
    result_lower = agent.validate_gene("tp53", "epithelial_cell")
    result_upper = agent.validate_gene("TP53", "epithelial_cell")
    assert result_lower["found"] == result_upper["found"]
    assert result_lower.get("ensembl_id") == result_upper.get("ensembl_id")
