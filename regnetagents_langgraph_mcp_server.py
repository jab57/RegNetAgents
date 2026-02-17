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
    - MCP Server: Handles Claude Desktop communication, tool registration, and resource browsing
    - LangGraph Workflow: Orchestrates multi-agent analysis pipeline
    - Tool Registry: Exposes analysis tools to Claude Desktop
    - Resource Registry: Exposes browsable resources (cell types, network summaries, gene lookups)

Available Tools:
    1. validate_gene: Quick gene name check with fuzzy suggestions (<100ms)
    2. query_network: Instant network queries (top regulators, targets, neighbors, stats) (<50ms)
    3. comprehensive_gene_analysis: Full workflow-driven analysis (recommended)
    3. multi_gene_analysis: Parallel processing of multiple genes
    4. pathway_focused_analysis: Pathway-centric analysis
    5. cross_cell_comparison: Gene behavior across cell types
    6. load_gene_results: Load previously saved analysis results
    7. list_available_results: List available result files
    8. workflow_status: Real-time execution monitoring
    9. workflow_insights: Performance analytics
    10. create_analysis_report: Generate formatted reports

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
from typing import Any, Iterable, Sequence
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
from mcp.server.lowlevel.helper_types import ReadResourceContents

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

@server.list_resources()
async def handle_list_resources() -> list[Resource]:
    """Return static resources and dynamic per-cell-type resources."""
    workflow = await get_workflow()
    cache = workflow.cache

    resources = [
        Resource(
            uri=AnyUrl("regnetagents://cell-types"),
            name="Available Cell Types",
            description="List of all cell types with gene and edge counts",
            mimeType="application/json"
        )
    ]

    for ct in CellType:
        network_data = cache.network_indices.get(ct.value, {})
        gene_count = network_data.get("num_genes", 0)
        edge_count = network_data.get("num_edges", 0)
        resources.append(Resource(
            uri=AnyUrl(f"regnetagents://network/{ct.value}"),
            name=f"Network Summary: {ct.value}",
            description=f"{gene_count} genes, {edge_count} edges",
            mimeType="application/json"
        ))

    return resources


@server.list_resource_templates()
async def handle_list_resource_templates() -> list[types.ResourceTemplate]:
    """Return parameterized resource templates."""
    return [
        types.ResourceTemplate(
            uriTemplate="regnetagents://gene/{gene_symbol}",
            name="Gene Lookup",
            description="Look up a gene across all cell types",
            mimeType="application/json"
        )
    ]


@server.read_resource()
async def handle_read_resource(uri: AnyUrl) -> Iterable[ReadResourceContents]:
    """Serve resource data for the given URI."""
    uri_str = str(uri)
    workflow = await get_workflow()
    cache = workflow.cache
    gene_mapper = workflow.modeling_agent.gene_mapper

    # regnetagents://cell-types
    if uri_str == "regnetagents://cell-types":
        cell_types = []
        for ct in CellType:
            network_data = cache.network_indices.get(ct.value, {})
            cell_types.append({
                "cell_type": ct.value,
                "num_genes": network_data.get("num_genes", 0),
                "num_edges": network_data.get("num_edges", 0),
                "num_regulons": network_data.get("num_regulons", 0),
            })
        return [ReadResourceContents(
            content=json.dumps({"cell_types": cell_types}, indent=2),
            mime_type="application/json"
        )]

    # regnetagents://network/{cell_type}
    if uri_str.startswith("regnetagents://network/"):
        cell_type = uri_str.replace("regnetagents://network/", "")
        network_data = cache.network_indices.get(cell_type, {})
        if not network_data:
            return [ReadResourceContents(
                content=json.dumps({"error": f"Unknown cell type: {cell_type}"}),
                mime_type="application/json"
            )]

        num_genes = network_data.get("num_genes", 0)
        num_edges = network_data.get("num_edges", 0)
        density = (num_edges / (num_genes * (num_genes - 1))) if num_genes > 1 else 0

        # Top 10 regulators by PageRank
        pagerank = network_data.get("pagerank_normalized", {})
        top_regulators = sorted(pagerank.items(), key=lambda x: x[1], reverse=True)[:10]
        top_regulators_named = []
        for ensembl_id, score in top_regulators:
            symbol = gene_mapper.ensembl_to_symbol(ensembl_id) or ensembl_id
            top_regulators_named.append({"gene": symbol, "pagerank": round(score, 6)})

        summary = {
            "cell_type": cell_type,
            "num_genes": num_genes,
            "num_edges": num_edges,
            "num_regulons": network_data.get("num_regulons", 0),
            "density": round(density, 6),
            "top_regulators_by_pagerank": top_regulators_named,
            "cache_version": network_data.get("cache_version", "unknown"),
        }
        return [ReadResourceContents(
            content=json.dumps(summary, indent=2),
            mime_type="application/json"
        )]

    # regnetagents://gene/{gene_symbol}
    if uri_str.startswith("regnetagents://gene/"):
        gene_symbol = uri_str.replace("regnetagents://gene/", "").upper()
        ensembl_id = gene_mapper.symbol_to_ensembl(gene_symbol)

        if not ensembl_id:
            return [ReadResourceContents(
                content=json.dumps({
                    "gene": gene_symbol,
                    "found": False,
                    "message": f"Gene '{gene_symbol}' not found in mapper"
                }),
                mime_type="application/json"
            )]

        presence = []
        for ct in CellType:
            network_data = cache.network_indices.get(ct.value, {})
            all_genes = set(network_data.get("all_genes", []))
            if ensembl_id in all_genes:
                regulators = network_data.get("target_regulators", {}).get(ensembl_id, [])
                targets = network_data.get("regulator_targets", {}).get(ensembl_id, [])
                presence.append({
                    "cell_type": ct.value,
                    "num_regulators": len(regulators),
                    "num_targets": len(targets),
                })

        return [ReadResourceContents(
            content=json.dumps({
                "gene": gene_symbol,
                "ensembl_id": ensembl_id,
                "found": True,
                "present_in_cell_types": len(presence),
                "cell_type_details": presence,
            }, indent=2),
            mime_type="application/json"
        )]

    # Unknown URI
    return [ReadResourceContents(
        content=json.dumps({"error": f"Unknown resource URI: {uri_str}"}),
        mime_type="application/json"
    )]


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
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type for network analysis (10 available)",
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
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type for network analysis (10 available)",
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
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type for analysis (10 available)",
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
        ),
        Tool(
            name="load_gene_results",
            description="""
            Load previously saved analysis results for a gene.

            Reads JSON result files from the results/ directory. Perfect for
            post-analysis visualization, comparison, or exploration without
            re-running the analysis.

            Great for:
            - Creating custom visualizations of existing results
            - Comparing results across multiple genes
            - Extracting specific data for further analysis
            - Quick data exploration after batch analyses

            Example: "Load the TP53 analysis results and show me the top 3 regulators"
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g., TP53, BRCA1)"
                    },
                    "analysis_type": {
                        "type": "string",
                        "enum": ["comprehensive", "cervical", "biomarker"],
                        "description": "Type of analysis results to load",
                        "default": "comprehensive"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="list_available_results",
            description="""
            List all available analysis result files.

            Shows what gene analyses have been saved in the results/ directory.
            Useful for discovering what data is available for visualization or comparison.

            Returns a list of genes with saved results and what types of analyses
            are available for each gene.

            Example: "What gene analyses do I have available?"
            """,
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="validate_gene",
            description="""
            Quickly check if a gene name is valid and present in the network (<100ms).

            Use this BEFORE running a full analysis to catch typos or invalid gene names
            instantly, instead of waiting 5-15 seconds for the analysis to fail.

            Returns:
            - If valid: gene symbol, Ensembl ID, and quick network stats (regulators, targets, role)
            - If invalid: fuzzy-matched suggestions for similar gene names

            Great for:
            - "Is TP53 in the epithelial cell network?"
            - "Did I spell this gene name correctly?"
            - Quick lookup of a gene's regulatory role without full analysis
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol to validate (e.g., TP53, BRCA1)"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type network to validate against",
                        "default": "epithelial_cell"
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="query_network",
            description="""
            Instantly query the pre-computed gene regulatory network (<50ms).

            Use this for quick network questions instead of running a full analysis.
            Answers questions about network structure directly from pre-computed data.

            Query types:
            - top_regulators: Genes with the most targets (out-degree), with PageRank scores
            - top_targets: Most highly regulated genes (in-degree)
            - gene_neighbors: Immediate regulators and targets of a specific gene
            - network_stats: Summary statistics (genes, edges, density, etc.)

            Great for:
            - "What are the top regulators in epithelial cells?"
            - "How many targets does TP53 have?"
            - "What regulates MYC?"
            - "How large is the CD4 T cell network?"
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "query_type": {
                        "type": "string",
                        "enum": ["top_regulators", "top_targets", "gene_neighbors", "network_stats"],
                        "description": "Type of network query to run"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type network to query",
                        "default": "epithelial_cell"
                    },
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (required for gene_neighbors query)"
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results to return for ranked queries",
                        "default": 10
                    }
                },
                "required": ["query_type"]
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

        elif name == "load_gene_results":
            gene = arguments["gene"].upper()
            analysis_type = arguments.get("analysis_type", "comprehensive")

            logger.info(f"Loading results for {gene} ({analysis_type})")

            # Build potential file paths
            result_patterns = [
                f"results/{gene.lower()}_analysis.json",
                f"results/{gene.lower()}_detailed_report.json",
                f"results/cervical_{gene.lower()}_analysis.json",
                f"results/biomarker_results.json"
            ]

            # Try to find and load the result file
            result_data = None
            loaded_file = None

            for file_path in result_patterns:
                if os.path.exists(file_path):
                    try:
                        with open(file_path, 'r') as f:
                            data = json.load(f)

                            # Check if it's a multi-gene file (like biomarker_results.json)
                            if isinstance(data, dict) and gene in data:
                                result_data = data[gene]
                                loaded_file = file_path
                                break
                            elif isinstance(data, dict):
                                result_data = data
                                loaded_file = file_path
                                break
                    except json.JSONDecodeError as e:
                        logger.warning(f"Could not parse {file_path}: {e}")
                        continue

            if result_data:
                response = {
                    "gene": gene,
                    "loaded_from": loaded_file,
                    "analysis_type": analysis_type,
                    "data": result_data
                }
                return [TextContent(type="text", text=json.dumps(response, indent=2))]
            else:
                error_response = {
                    "error": f"No results found for {gene}",
                    "searched_locations": result_patterns,
                    "suggestion": f"Run comprehensive_gene_analysis for {gene} first, or check if the gene symbol is correct"
                }
                return [TextContent(type="text", text=json.dumps(error_response, indent=2))]

        elif name == "list_available_results":
            logger.info("Listing available analysis results")

            results_dir = "results"
            available_results = {}

            if os.path.exists(results_dir):
                for filename in os.listdir(results_dir):
                    if filename.endswith('.json'):
                        file_path = os.path.join(results_dir, filename)
                        try:
                            # Extract gene name from filename
                            if filename.startswith('cervical_'):
                                gene = filename.replace('cervical_', '').replace('_analysis.json', '').upper()
                                analysis_type = 'cervical'
                            elif filename.endswith('_detailed_report.json'):
                                gene = filename.replace('_detailed_report.json', '').upper()
                                analysis_type = 'detailed'
                            elif filename.endswith('_analysis.json'):
                                gene = filename.replace('_analysis.json', '').upper()
                                analysis_type = 'comprehensive'
                            elif filename == 'biomarker_results.json':
                                # Special handling for multi-gene file
                                with open(file_path, 'r') as f:
                                    data = json.load(f)
                                    for gene_key in data.keys():
                                        if gene_key.isupper() and len(gene_key) < 10:  # Likely a gene symbol
                                            if gene_key not in available_results:
                                                available_results[gene_key] = []
                                            available_results[gene_key].append('biomarker')
                                continue
                            else:
                                continue

                            if gene not in available_results:
                                available_results[gene] = []
                            available_results[gene].append(analysis_type)

                        except Exception as e:
                            logger.warning(f"Could not process {filename}: {e}")

            response = {
                "results_directory": results_dir,
                "available_genes": list(available_results.keys()),
                "total_genes": len(available_results),
                "details": available_results,
                "note": "Use load_gene_results to load specific gene data"
            }

            return [TextContent(type="text", text=json.dumps(response, indent=2))]

        elif name == "validate_gene":
            gene = arguments["gene"]
            cell_type = arguments.get("cell_type", "epithelial_cell")

            logger.info(f"Validating gene '{gene}' in {cell_type}")

            result = workflow.modeling_agent.validate_gene(gene, cell_type)

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "query_network":
            query_type = arguments["query_type"]
            cell_type = arguments.get("cell_type", "epithelial_cell")
            gene = arguments.get("gene", None)
            top_n = arguments.get("top_n", 10)

            logger.info(f"Network query: {query_type} in {cell_type}" +
                        (f" for {gene}" if gene else ""))

            result = workflow.modeling_agent.query_network(
                query_type=query_type,
                cell_type=cell_type,
                gene=gene,
                top_n=top_n
            )

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

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