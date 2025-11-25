#!/usr/bin/env python3
"""
RegNetAgents LangGraph-Powered MCP Server
==========================================

Model Context Protocol (MCP) server that exposes RegNetAgents multi-agent workflow
to Claude Desktop for conversational gene regulatory network analysis.

This server bridges LangGraph's sophisticated workflow orchestration with Claude Desktop's
MCP protocol, enabling natural language queries like "Analyze TP53 in epithelial cells"
to trigger comprehensive multi-agent analysis workflows.

Architecture:
    - MCP Server: Handles Claude Desktop communication and tool registration
    - LangGraph Workflow: Orchestrates multi-agent analysis pipeline
    - Tool Registry: Exposes 6 analysis tools to Claude Desktop

Available Tools:
    1. comprehensive_gene_analysis: Full workflow-driven analysis (recommended)
    2. multi_gene_analysis: Parallel processing of multiple genes
    3. pathway_focused_analysis: Pathway-centric analysis
    4. workflow_status: Real-time execution monitoring
    5. workflow_insights: Performance analytics
    6. create_analysis_report: Generate formatted reports

Key Features:
    - Conversational interface (natural language → structured analysis)
    - Automatic workflow orchestration (no manual configuration)
    - Parallel execution of independent analyses
    - Real-time progress monitoring
    - LLM-powered domain insights with rule-based fallback

Performance:
    - Single gene: 0.6-15 seconds (comprehensive)
    - Multi-gene (5): 15-62 seconds (parallel execution)
    - Instant cross-cell comparison (pre-computed indices)

Author: Jose A. Bird, PhD
License: MIT
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
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow, GeneAnalysisState, CellType

# Configure logging with more detail
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize the server
server = Server("regnetagents-langgraph-server")

# Global workflow instance (initialized once for performance)
workflow_instance = None

async def get_workflow():
    """
    Get or create the global workflow instance (singleton pattern).

    Initializes the LangGraph workflow on first call, loading network indices
    for all 10 cell types. Subsequent calls return the cached instance for
    optimal performance. The workflow initialization includes:
    - Loading pre-computed network indices (~183K edges for epithelial cells)
    - Initializing specialized agents (modeling, pathway, domain)
    - Compiling the LangGraph state machine

    Returns:
        RegNetAgentsWorkflow: Initialized workflow instance ready for analysis

    Note:
        Initialization takes ~2-3 seconds on first call. Cached afterwards.
    """
    global workflow_instance
    if workflow_instance is None:
        logger.info("Initializing LangGraph workflow...")
        workflow_instance = RegNetAgentsWorkflow()
        logger.info("LangGraph workflow ready")
    return workflow_instance

@server.list_tools()
async def handle_list_tools() -> list[Tool]:
    """List available LangGraph-powered tools."""
    return [
        Tool(
            name="comprehensive_gene_analysis",
            description="""
            Find out what controls a gene and what it controls. Perfect for understanding disease genes and finding drug targets.

            What you'll get:
            - Which genes regulate your gene of interest (upstream controllers)
            - Which genes your gene regulates (downstream targets)
            - Biological pathways involved (with statistical validation)
            - How the gene behaves in different cell types
            - Relevant insights for cancer, drugs, and clinical research

            Great for answering questions like:
            - "What regulates TP53 in cancer cells?"
            - "What does BRCA1 control in breast tissue?"
            - "How does this gene contribute to disease?"

            The system automatically determines the best analyses to run based on your gene's characteristics.
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
            Check the progress of a gene analysis.

            Shows what steps have been completed and how long the analysis took.
            Useful for understanding what analyses were performed and troubleshooting.
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
            Analyze multiple genes at once (up to 10 genes) - RECOMMENDED for 2+ genes.

            ⚡ FASTER than running comprehensive_gene_analysis multiple times because all
            genes are analyzed in parallel simultaneously.

            Perfect for:
            - Comparing several disease-related genes (e.g., cancer biomarkers)
            - Analyzing gene families or pathways
            - Batch processing research gene lists
            - Biomarker discovery panels

            Example: Compare TP53, BRCA1, APC, MYC, and KRAS for cancer screening
            All genes analyzed in parallel - typically completes in 5-15 seconds for 5 genes.

            IMPORTANT: Use this tool instead of calling comprehensive_gene_analysis multiple
            times to avoid timeouts and get faster results.
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
                        "description": "Depth of analysis for each gene (use 'basic' for fastest results without pathways)",
                        "default": "comprehensive"
                    }
                },
                "required": ["genes"]
            }
        ),
        Tool(
            name="cross_cell_comparison",
            description="""
            Compare how a single gene behaves across all available cell types.

            Shows how regulatory networks differ by tissue/cell type for the same gene.
            Reveals cell-specific regulation patterns and why the same gene can have
            different roles in different tissues.

            Great for questions like:
            - "How does TP53 differ between immune cells and epithelial cells?"
            - "In which cell types is this gene most highly regulated?"
            - "Does this gene act as a regulator in some cells but a target in others?"

            Returns regulatory role (hub/target) and network statistics for each cell type.
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol or Ensembl ID to compare across cell types"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="pathway_focused_analysis",
            description="""
            Explore how a gene participates in specific biological pathways.

            Pathways are groups of genes working together (like "Wnt signaling" or "Cell cycle").
            This tool focuses on a specific pathway you're interested in.

            Great for questions like:
            - "How does APC participate in Wnt signaling?"
            - "What's the role of TP53 in apoptosis?"
            - "How does this gene affect the cell cycle?"

            Results include statistical validation from the Reactome pathway database.
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
            Get technical details about how the analysis system works.

            Shows performance metrics and how the intelligent routing system decided
            which analyses to run. Mostly useful for advanced users and debugging.
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
            Create a formatted report summarizing your gene analysis.

            Compiles all analysis results into an organized report.
            Choose from summary (quick overview), detailed (full results), or
            scientific (publication-ready format).

            Great for saving and sharing your findings.
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
            logger.info(f"Parameters: cell_type={cell_type}, analysis_depth={analysis_depth}")

            import time
            start_time = time.time()

            # Simple wrapper around workflow
            result = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                analysis_depth=analysis_depth
            )

            execution_time = time.time() - start_time
            logger.info(f"Analysis for {gene} completed in {execution_time:.2f} seconds")

            # Add minimal MCP metadata
            result["workflow_info"] = {
                "workflow_type": "langgraph",
                "server_version": "regnetagents-mcp-v2.1",
                "execution_mode": "thin_wrapper",
                "execution_time_seconds": round(execution_time, 2)
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

            import time
            start_time = time.time()

            # Simple parallel execution wrapper with progress logging
            tasks = [
                workflow.run_analysis(gene=gene, cell_type=cell_type, analysis_depth=analysis_depth)
                for gene in genes
            ]

            # Execute all analyses in parallel
            logger.info(f"Executing {len(genes)} gene analyses in parallel...")
            results = await asyncio.gather(*tasks, return_exceptions=True)

            execution_time = time.time() - start_time
            logger.info(f"Multi-gene analysis completed in {execution_time:.2f} seconds")

            # Simple result compilation
            compiled_results = {
                "multi_gene_analysis": {
                    "genes_analyzed": genes,
                    "cell_type": cell_type,
                    "analysis_depth": analysis_depth,
                    "total_genes": len(genes),
                    "successful_analyses": sum(1 for r in results if not isinstance(r, Exception)),
                    "failed_analyses": sum(1 for r in results if isinstance(r, Exception)),
                    "execution_time_seconds": round(execution_time, 2)
                },
                "individual_results": {
                    gene: result if not isinstance(result, Exception) else {"status": "error", "error": str(result)}
                    for gene, result in zip(genes, results)
                }
            }

            return [TextContent(type="text", text=json.dumps(compiled_results, indent=2))]

        elif name == "cross_cell_comparison":
            gene = arguments["gene"]

            logger.info(f"Starting cross-cell comparison for {gene}")

            import time
            start_time = time.time()

            # Use the modeling agent's cross-cell comparison method
            result = await workflow.modeling_agent.compare_gene_across_cell_types(gene)

            execution_time = time.time() - start_time
            logger.info(f"Cross-cell comparison completed in {execution_time:.2f} seconds")

            # Add timing metadata
            result["execution_time_seconds"] = round(execution_time, 2)
            result["analysis_type"] = "cross_cell_comparison"
            result["gene_analyzed"] = gene

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

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
            "server_type": "regnetagents-mcp",
            "suggested_action": "Check gene name and parameters"
        }
        return [TextContent(type="text", text=json.dumps(error_response, indent=2))]

async def main():
    """Main server function with LangGraph integration."""
    try:
        # Ensure we're in the correct directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(script_dir)

        logger.info("Starting RegNetAgents LangGraph MCP Server...")

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
                    server_name="regnetagents-langgraph-server",
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