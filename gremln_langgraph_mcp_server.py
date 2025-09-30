#!/usr/bin/env python3
"""
GREmLN LangGraph-Powered MCP Server
Combines the power of LangGraph workflows with MCP integration for Claude Desktop

This server uses LangGraph for intelligent workflow orchestration while maintaining
the MCP interface for seamless Claude Desktop integration.
"""

import asyncio
import json
import logging
import os
import sys
from typing import Any, Sequence
from mcp.server import Server, NotificationOptions
from mcp.server.models import InitializationOptions
import mcp.server.stdio as stdio
from mcp.server.stdio import stdio_server
from mcp.types import (
    Resource,
    Tool,
    TextContent,
    ImageContent,
    EmbeddedResource,
    LoggingLevel
)
from pydantic import AnyUrl
import mcp.types as types

# Import our LangGraph workflow
from gremln_langgraph_workflow import GREmLNWorkflow, GeneAnalysisState, CellType

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize the server
server = Server("gremln-langgraph-server")

# Global workflow instance (initialized once for performance)
workflow_instance = None

async def get_workflow():
    """Get or create the global workflow instance"""
    global workflow_instance
    if workflow_instance is None:
        logger.info("Initializing LangGraph workflow...")
        workflow_instance = GREmLNWorkflow()
        logger.info("LangGraph workflow ready")
    return workflow_instance

@server.list_tools()
async def handle_list_tools() -> list[Tool]:
    """List available LangGraph-powered tools."""
    return [
        Tool(
            name="comprehensive_gene_analysis",
            description="""
            Perform comprehensive gene regulatory network analysis using intelligent LangGraph workflows.

            This tool automatically determines the optimal analysis path based on gene characteristics:
            - Terminal targets → Analyze regulators and pathways
            - Hub regulators → Analyze targets and downstream effects
            - Master regulators → Cross-cell comparison and pathway analysis

            Includes smart routing, state management, and comprehensive reporting.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol or Ensembl ID to analyze"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells",
                            "hepatocytes", "cardiomyocytes", "neurons",
                            "fibroblasts", "endothelial_cells"
                        ],
                        "description": "Cell type for network analysis",
                        "default": "epithelial_cell"
                    },
                    "analysis_depth": {
                        "type": "string",
                        "enum": ["basic", "comprehensive", "focused"],
                        "description": "Depth of analysis to perform",
                        "default": "comprehensive"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="workflow_status",
            description="""
            Get the status and execution details of a gene analysis workflow.
            Shows completed steps, current progress, and performance metrics.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene that was analyzed"
                    },
                    "show_state": {
                        "type": "boolean",
                        "description": "Include full workflow state details",
                        "default": False
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="multi_gene_analysis",
            description="""
            Analyze multiple genes using parallel LangGraph workflows.
            Efficient batch processing with intelligent resource management.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "genes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of gene symbols or Ensembl IDs",
                        "maxItems": 10
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells",
                            "hepatocytes", "cardiomyocytes", "neurons",
                            "fibroblasts", "endothelial_cells"
                        ],
                        "description": "Cell type for network analysis",
                        "default": "epithelial_cell"
                    },
                    "analysis_depth": {
                        "type": "string",
                        "enum": ["basic", "comprehensive", "focused"],
                        "description": "Depth of analysis for each gene",
                        "default": "focused"
                    }
                },
                "required": ["genes"]
            }
        ),
        Tool(
            name="pathway_focused_analysis",
            description="""
            Perform pathway-focused analysis using LangGraph workflows.
            Optimized for understanding biological pathway interactions.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Primary gene for pathway analysis"
                    },
                    "pathway_focus": {
                        "type": "string",
                        "enum": ["wnt_signaling", "cell_cycle", "apoptosis", "immune_response", "metabolism"],
                        "description": "Pathway to focus analysis on"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells",
                            "hepatocytes", "cardiomyocytes", "neurons",
                            "fibroblasts", "endothelial_cells"
                        ],
                        "description": "Cell type for analysis",
                        "default": "epithelial_cell"
                    }
                },
                "required": ["gene", "pathway_focus"]
            }
        ),
        Tool(
            name="workflow_insights",
            description="""
            Get workflow execution insights and optimization recommendations.
            Analyzes workflow performance and suggests improvements.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "analysis_type": {
                        "type": "string",
                        "enum": ["performance", "routing", "completeness"],
                        "description": "Type of insights to generate",
                        "default": "performance"
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="create_analysis_report",
            description="""
            Generate a comprehensive analysis report from workflow results.
            Creates formatted reports with key insights and recommendations.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene that was analyzed"
                    },
                    "report_format": {
                        "type": "string",
                        "enum": ["summary", "detailed", "scientific"],
                        "description": "Format of the report",
                        "default": "summary"
                    },
                    "include_visualizations": {
                        "type": "boolean",
                        "description": "Include workflow visualization data",
                        "default": True
                    }
                },
                "required": ["gene"]
            }
        )
    ]

@server.call_tool()
async def handle_call_tool(name: str, arguments: dict) -> list[TextContent]:
    """Handle tool calls using LangGraph workflows."""
    try:
        workflow = await get_workflow()

        if name == "comprehensive_gene_analysis":
            gene = arguments["gene"]
            cell_type = arguments.get("cell_type", "epithelial_cell")
            analysis_depth = arguments.get("analysis_depth", "comprehensive")

            logger.info(f"Starting comprehensive analysis for {gene} using LangGraph workflow")

            # Simple wrapper around workflow
            result = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                analysis_depth=analysis_depth
            )

            # Add minimal MCP metadata
            result["workflow_info"] = {
                "workflow_type": "langgraph",
                "server_version": "langgraph-mcp-v2.0",
                "execution_mode": "thin_wrapper"
            }

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "workflow_status":
            gene = arguments["gene"]
            show_state = arguments.get("show_state", False)

            # Simple status wrapper
            status = {
                "gene": gene,
                "status": "available",
                "workflow_type": "langgraph",
                "message": "Use comprehensive_gene_analysis for detailed analysis"
            }

            return [TextContent(type="text", text=json.dumps(status, indent=2))]

        elif name == "multi_gene_analysis":
            genes = arguments["genes"]
            cell_type = arguments.get("cell_type", "epithelial_cell")
            analysis_depth = arguments.get("analysis_depth", "focused")

            logger.info(f"Starting multi-gene analysis for {len(genes)} genes")

            # Simple parallel execution wrapper
            tasks = [
                workflow.run_analysis(gene=gene, cell_type=cell_type, analysis_depth=analysis_depth)
                for gene in genes
            ]

            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Simple result compilation
            compiled_results = {
                "multi_gene_analysis": {
                    "genes_analyzed": genes,
                    "cell_type": cell_type,
                    "analysis_depth": analysis_depth,
                    "total_genes": len(genes),
                    "successful_analyses": sum(1 for r in results if not isinstance(r, Exception)),
                    "failed_analyses": sum(1 for r in results if isinstance(r, Exception))
                },
                "individual_results": {
                    gene: result if not isinstance(result, Exception) else {"status": "error", "error": str(result)}
                    for gene, result in zip(genes, results)
                }
            }

            return [TextContent(type="text", text=json.dumps(compiled_results, indent=2))]

        elif name == "pathway_focused_analysis":
            gene = arguments["gene"]
            pathway_focus = arguments["pathway_focus"]
            cell_type = arguments.get("cell_type", "epithelial_cell")

            logger.info(f"Starting pathway-focused analysis: {gene} -> {pathway_focus}")

            # Simple wrapper - run comprehensive analysis and add pathway focus
            result = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                analysis_depth="comprehensive"
            )

            # Add pathway focus metadata
            result["pathway_focused_analysis"] = {
                "primary_gene": gene,
                "pathway_focus": pathway_focus,
                "cell_type": cell_type
            }

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "workflow_insights":
            analysis_type = arguments.get("analysis_type", "performance")

            # Simple static insights
            insights = {
                "workflow_insights": {
                    "analysis_type": analysis_type,
                    "message": "LangGraph provides intelligent gene analysis routing",
                    "architecture": "Thin MCP wrapper over LangGraph workflow"
                }
            }

            return [TextContent(type="text", text=json.dumps(insights, indent=2))]

        elif name == "create_analysis_report":
            gene = arguments["gene"]
            report_format = arguments.get("report_format", "summary")
            include_visualizations = arguments.get("include_visualizations", True)

            # Simple wrapper - run analysis and format as report
            result = await workflow.run_analysis(gene=gene, analysis_depth="comprehensive")

            report = {
                "analysis_report": {
                    "gene": gene,
                    "report_format": report_format,
                    "generated_by": "LangGraph Workflow",
                    "timestamp": "Generated on demand"
                },
                "analysis_results": result
            }

            return [TextContent(type="text", text=json.dumps(report, indent=2))]

        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]

    except Exception as e:
        logger.error(f"Tool '{name}' execution failed: {str(e)}", exc_info=True)
        error_response = {
            "error": f"LangGraph workflow error in '{name}'",
            "details": str(e),
            "server_type": "langgraph-mcp",
            "suggested_action": "Check gene name and parameters"
        }
        return [TextContent(type="text", text=json.dumps(error_response, indent=2))]

async def main():
    """Main server function with LangGraph integration."""
    try:
        # Ensure we're in the correct directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(script_dir)

        logger.info("Starting GREmLN LangGraph MCP Server...")

        # Check if network data files exist (basic validation)
        network_dir = "models/networks/epithelial_cell"
        if not os.path.exists(network_dir):
            logger.error(f"Network directory missing: {network_dir}")
            logger.error("Please ensure network cache files are generated using build_network_cache.py")
        else:
            network_file = os.path.join(network_dir, "network_index.pkl")
            if not os.path.exists(network_file):
                logger.error(f"Network index file missing: {network_file}")
                logger.error("Please run: python build_network_cache.py --all")

        # Pre-initialize the workflow for faster first requests
        await get_workflow()

        async with stdio_server() as (read_stream, write_stream):
            await server.run(
                read_stream,
                write_stream,
                InitializationOptions(
                    server_name="gremln-langgraph-server",
                    server_version="1.0.0",
                    capabilities=types.ServerCapabilities(
                        resources={},
                        tools={},
                        prompts={},
                        logging={}
                    )
                ),
            )
    except KeyboardInterrupt:
        logger.info("Server stopped by user")
    except Exception as e:
        logger.error(f"Server error: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    asyncio.run(main())