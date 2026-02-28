#!/usr/bin/env python3
"""
Tests for the RegNetAgents LangGraph-powered MCP Server.
Tests the hybrid MCP-LangGraph integration via handle_call_tool.
"""

import asyncio
import json
import sys

from regnetagents_langgraph_mcp_server import handle_call_tool


async def test_comprehensive_gene_analysis():
    """comprehensive_gene_analysis for APC → gene summary with regulatory role."""
    result = await handle_call_tool("comprehensive_gene_analysis", {
        "gene": "APC",
        "cell_type": "epithelial_cell",
        "analysis_depth": "comprehensive"
    })
    assert result and len(result) > 0, "No response received"
    response = json.loads(result[0].text)
    assert not response.get("error"), f"Tool returned error: {response.get('error')}"
    summary = response.get("gene_analysis_summary", {})
    assert summary.get("gene") == "APC"
    assert "regulatory_role" in summary


async def test_multi_gene_analysis():
    """multi_gene_analysis for APC, BRCA1, TP53 → all 3 analyzed successfully."""
    result = await handle_call_tool("multi_gene_analysis", {
        "genes": ["APC", "BRCA1", "TP53"],
        "cell_type": "epithelial_cell",
        "analysis_depth": "focused"
    })
    assert result and len(result) > 0, "No response received"
    response = json.loads(result[0].text)
    assert not response.get("error"), f"Tool returned error: {response.get('error')}"
    summary = response.get("multi_gene_analysis", {})
    assert summary.get("total_genes") == 3
    assert summary.get("successful_analyses") == 3


async def test_workflow_status():
    """workflow_status → response contains status and workflow_type fields."""
    result = await handle_call_tool("workflow_status", {
        "gene": "APC",
        "show_state": True
    })
    assert result and len(result) > 0, "No response received"
    response = json.loads(result[0].text)
    assert "status" in response
    assert "workflow_type" in response


async def test_pathway_focused_analysis():
    """pathway_focused_analysis for APC/wnt_signaling → primary_gene and pathway_focus present."""
    result = await handle_call_tool("pathway_focused_analysis", {
        "gene": "APC",
        "pathway_focus": "wnt_signaling",
        "cell_type": "epithelial_cell"
    })
    assert result and len(result) > 0, "No response received"
    response = json.loads(result[0].text)
    assert not response.get("error"), f"Tool returned error: {response.get('error')}"
    pathway_info = response.get("pathway_focused_analysis", {})
    assert pathway_info.get("primary_gene") == "APC"
    assert pathway_info.get("pathway_focus") == "wnt_signaling"


async def test_workflow_insights():
    """workflow_insights (performance) → message and architecture keys present."""
    result = await handle_call_tool("workflow_insights", {
        "analysis_type": "performance"
    })
    assert result and len(result) > 0, "No response received"
    response = json.loads(result[0].text)
    assert "message" in response, f"Missing 'message' key in response: {response}"
    assert "architecture" in response, f"Missing 'architecture' key in response: {response}"


async def test_create_analysis_report():
    """create_analysis_report for APC → report contains gene, report_format, generated_by."""
    result = await handle_call_tool("create_analysis_report", {
        "gene": "APC",
        "report_format": "summary",
        "include_visualizations": True
    })
    assert result and len(result) > 0, "No response received"
    response = json.loads(result[0].text)
    report = response.get("analysis_report", {})
    assert report.get("gene") == "APC"
    assert "report_format" in report
    assert "generated_by" in report


async def demonstrate_workflow_advantages():
    """Print a summary of LangGraph MCP server advantages (for manual inspection)."""
    print("\n" + "=" * 60)
    print("LANGGRAPH MCP SERVER ADVANTAGES:")
    print("=" * 60)

    advantages = {
        "Intelligent Workflow Orchestration": [
            "Visual state management with clear transitions",
            "Smart routing based on gene characteristics",
            "Conditional execution of relevant analyses only",
            "Graceful error handling and recovery"
        ],
        "Enhanced Performance": [
            "Shared cache instances across workflow",
            "Parallel processing for multi-gene analysis",
            "Optimized resource usage",
            "Fast gene mapping integration"
        ],
        "Advanced Features": [
            "Comprehensive workflow status tracking",
            "Pathway-focused analysis capabilities",
            "Multi-gene batch processing",
            "Detailed execution insights"
        ],
        "MCP Integration Benefits": [
            "Seamless Claude Desktop compatibility",
            "Rich tool descriptions and schemas",
            "Structured JSON responses",
            "Error handling with actionable feedback"
        ]
    }

    for category, items in advantages.items():
        print(f"\n{category}:")
        for item in items:
            print(f"  • {item}")

    print(f"\n{'=' * 60}")
    print("CONFIGURATION FOR CLAUDE DESKTOP:")
    print("Add this to your claude_desktop_config.json:")
    print("=" * 60)

    config = {
        "mcpServers": {
            "regnetagents-langgraph-server": {
                "command": "python",
                "args": ["C:/path/to/RegNetAgents/regnetagents_langgraph_mcp_server.py"],
                "env": {
                    "PYTHONPATH": "C:/path/to/RegNetAgents"
                }
            }
        }
    }

    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    async def _main():
        await test_comprehensive_gene_analysis()
        await test_multi_gene_analysis()
        await test_workflow_status()
        await test_pathway_focused_analysis()
        await test_workflow_insights()
        await test_create_analysis_report()
        await demonstrate_workflow_advantages()
        print("\nLangGraph MCP Server is ready for Claude Desktop!")

    try:
        asyncio.run(_main())
    except AssertionError as e:
        print(f"\nTest FAILED: {e}")
        sys.exit(1)
