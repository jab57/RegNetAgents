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
    - Prompt Registry: Exposes guided prompt templates for common analysis workflows

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
import time
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
# IMPORTANT: MCP uses stdout for JSON-RPC transport, so all logging MUST go to stderr
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)

# Initialize the server
server = Server("regnetagents-langgraph-server")

# Global workflow instance (initialized once for performance)
workflow_instance = None

# In-memory result cache: {(gene, cell_type, analysis_depth): (timestamp, result)}
_result_cache = {}
_CACHE_TTL_SECONDS = 600  # 10 minutes

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


# ---------------------------------------------------------------------------
# Cancer gene panels used by the cancer_biomarker_panel prompt
# ---------------------------------------------------------------------------
_CANCER_PANELS = {
    "colorectal":   ["APC", "TP53", "KRAS", "MYC", "CTNNB1"],
    "breast":       ["BRCA1", "BRCA2", "TP53", "ERBB2", "ESR1"],
    "lung":         ["KRAS", "EGFR", "TP53", "ALK", "STK11"],
    "prostate":     ["AR", "TP53", "PTEN", "MYC", "ERG"],
    "melanoma":     ["BRAF", "CDKN2A", "TP53", "PTEN", "MITF"],
    "pancreatic":   ["KRAS", "TP53", "SMAD4", "CDKN2A", "BRCA2"],
    "glioblastoma": ["EGFR", "TP53", "PTEN", "IDH1", "CDKN2A"],
    "general":      ["TP53", "MYC", "KRAS", "BRCA1", "EGFR"],
}

_CELL_TYPES = [
    "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
    "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
    "erythrocytes", "nk_cells", "nkt_cells",
    "monocyte-derived_dendritic_cells"
]


@server.list_prompts()
async def handle_list_prompts() -> list[types.Prompt]:
    """Return available prompt templates."""
    return [
        types.Prompt(
            name="gene_deep_dive",
            description=(
                "Comprehensive guided analysis of a single gene. "
                "Validates the gene, runs full network analysis, queries immediate neighbors, "
                "and summarizes regulatory role, pathways, and clinical relevance."
            ),
            arguments=[
                types.PromptArgument(
                    name="gene",
                    description="Gene symbol to analyze (e.g. TP53, MYC, BRCA1)",
                    required=True
                ),
                types.PromptArgument(
                    name="cell_type",
                    description=f"Cell type network to use. Options: {', '.join(_CELL_TYPES)}",
                    required=False
                ),
            ]
        ),
        types.Prompt(
            name="cancer_biomarker_panel",
            description=(
                "Analyze a pre-defined panel of cancer-related genes in parallel. "
                "Available panels: colorectal, breast, lung, prostate, general."
            ),
            arguments=[
                types.PromptArgument(
                    name="cancer_type",
                    description="Cancer type panel to use: colorectal, breast, lung, prostate, or general",
                    required=True
                ),
                types.PromptArgument(
                    name="cell_type",
                    description=f"Cell type network to use. Options: {', '.join(_CELL_TYPES)}",
                    required=False
                ),
            ]
        ),
        types.Prompt(
            name="cross_cell_comparison",
            description=(
                "Compare how a single gene behaves across all 10 cell-type networks. "
                "Highlights differences in regulatory role between immune and epithelial contexts."
            ),
            arguments=[
                types.PromptArgument(
                    name="gene",
                    description="Gene symbol to compare across cell types (e.g. TP53, MYC)",
                    required=True
                ),
            ]
        ),
        types.Prompt(
            name="tumor_context_analysis",
            description=(
                "Full domain analysis of a gene against a TCGA tumor-state ARACNe network. "
                "Mirrors gene_deep_dive but uses tumor-state topology instead of population-averaged GREmLN networks. "
                "Includes MoA breakdown (activating vs. repressive regulators), master regulator identification, "
                "druggability, and clinical actionability."
            ),
            arguments=[
                types.PromptArgument(
                    name="gene",
                    description="Gene symbol to analyze (e.g. YAP1, ESR1, MYC)",
                    required=True
                ),
                types.PromptArgument(
                    name="cancer_type",
                    description="TCGA cancer type: brca, coad, hnsc, luad, lusc, ov, prad, ucec",
                    required=True
                ),
            ]
        ),
        types.Prompt(
            name="network_context_comparison",
            description=(
                "Compare a gene's regulatory context between population-averaged (GREmLN epithelial) "
                "and tumor-state (TCGA) networks. Returns conserved regulators, population-averaged-only "
                "regulators, and tumor-state-only regulators with biological interpretation."
            ),
            arguments=[
                types.PromptArgument(
                    name="gene",
                    description="Gene symbol to compare across network contexts (e.g. MYC, TP53)",
                    required=True
                ),
                types.PromptArgument(
                    name="cancer_type",
                    description="TCGA cancer type: brca, coad, hnsc, luad, lusc, ov, prad, ucec",
                    required=True
                ),
            ]
        ),
        types.Prompt(
            name="candidate_prioritization",
            description=(
                "Two-step regulatory candidate prioritization workflow from the RegNetAgents NAR paper. "
                "Step 1: identifies source-labeled regulatory candidates (TCGA-only / GREmLN-only / Both) "
                "filtered against OncoKB. Step 2: runs comprehensive domain analysis per candidate with "
                "source-driven network routing. Returns a structured summary table."
            ),
            arguments=[
                types.PromptArgument(
                    name="gene",
                    description="Focal gene to find regulatory candidates for (e.g. CTNNB1, MYC)",
                    required=True
                ),
                types.PromptArgument(
                    name="cancer_type",
                    description="TCGA cancer type: brca, coad, hnsc, luad, lusc, ov, prad, ucec",
                    required=True
                ),
            ]
        ),
    ]


@server.get_prompt()
async def handle_get_prompt(name: str, arguments: dict | None) -> types.GetPromptResult:
    """Generate a prompt message for the requested template."""
    args = arguments or {}

    if name == "gene_deep_dive":
        gene = args.get("gene", "").upper() or "TP53"
        cell_type = args.get("cell_type", "epithelial_cell")

        text = (
            f"Please perform a comprehensive deep-dive analysis of **{gene}** "
            f"in **{cell_type}** using the following steps:\n\n"
            f"1. Use `validate_gene` to confirm **{gene}** exists in the {cell_type} network "
            f"and get basic stats.\n"
            f"2. Run `comprehensive_gene_analysis` for the full multi-agent report "
            f"(network position, pathways, therapeutic targets, domain insights).\n"
            f"3. Use `query_network` with `query_type=\"gene_neighbors\"` and `gene=\"{gene}\"` "
            f"to show the immediate regulatory context.\n"
            f"4. Summarize:\n"
            f"   - Regulatory role (hub, intermediate, heavily regulated?)\n"
            f"   - Top upstream regulators and downstream targets\n"
            f"   - Key enriched pathways\n"
            f"   - Cancer relevance and druggability\n"
            f"   - Any clinical significance\n"
        )
        return types.GetPromptResult(
            description=f"Deep-dive analysis of {gene} in {cell_type}",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(type="text", text=text)
                )
            ]
        )

    elif name == "cancer_biomarker_panel":
        cancer_type = args.get("cancer_type", "general").lower()
        cell_type = args.get("cell_type", "epithelial_cell")
        panel = _CANCER_PANELS.get(cancer_type, _CANCER_PANELS["general"])
        gene_list_str = ", ".join(panel)

        text = (
            f"Please analyze the **{cancer_type} cancer biomarker panel** "
            f"in **{cell_type}** using the following steps:\n\n"
            f"Gene panel: {gene_list_str}\n\n"
            f"1. Use `multi_gene_analysis` with genes={panel} and cell_type=\"{cell_type}\" "
            f"to run all analyses in parallel.\n"
            f"2. For each gene, briefly summarize:\n"
            f"   - Regulatory role in {cell_type}\n"
            f"   - Number of regulators and targets\n"
            f"   - Cancer relevance\n"
            f"3. Identify which gene in the panel has the **highest network centrality** "
            f"(most targets or highest PageRank) — this is the likely master regulator.\n"
            f"4. Note any genes that are **absent from the network** "
            f"(not expressed in {cell_type}).\n"
            f"5. Provide a one-paragraph summary of the panel's collective regulatory landscape.\n"
        )
        return types.GetPromptResult(
            description=f"{cancer_type.title()} cancer biomarker panel in {cell_type}",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(type="text", text=text)
                )
            ]
        )

    elif name == "cross_cell_comparison":
        gene = args.get("gene", "").upper() or "TP53"

        text = (
            f"Please compare **{gene}** across all available cell-type networks "
            f"using the following steps:\n\n"
            f"1. Use `cross_cell_comparison` with `gene=\"{gene}\"` to get network statistics "
            f"across all 10 cell types.\n"
            f"2. Summarize the results in a table with columns: "
            f"Cell Type | Regulatory Role | # Regulators | # Targets | In Network.\n"
            f"3. Highlight:\n"
            f"   - Which cell type shows the **highest activity** for {gene} "
            f"(most regulators + targets)\n"
            f"   - Differences between **immune cells** (T cells, B cells, monocytes, NK cells) "
            f"and **epithelial cells**\n"
            f"   - Any cell types where {gene} is **absent from the network**\n"
            f"4. Provide a brief biological interpretation of why {gene} might differ "
            f"across these contexts.\n"
        )
        return types.GetPromptResult(
            description=f"Cross-cell-type comparison of {gene}",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(type="text", text=text)
                )
            ]
        )

    elif name == "tumor_context_analysis":
        gene = args.get("gene", "").upper() or "YAP1"
        cancer_type = args.get("cancer_type", "brca").lower()

        text = (
            f"Please perform a comprehensive tumor-context analysis of **{gene}** "
            f"in the **TCGA {cancer_type.upper()} tumor-state network** using the following steps:\n\n"
            f"1. Use `query_network` with `network_source=\"tcga\"`, `tcga_network=\"{cancer_type}\"`, "
            f"`gene=\"{gene}\"`, and `query_type=\"gene_neighbors\"` to retrieve tumor-state regulatory neighbors, "
            f"including MoA (mode of action: activating vs. repressive).\n"
            f"2. Use `find_master_regulators` with `network_source=\"tcga\"` and `tcga_network=\"{cancer_type}\"` "
            f"to identify transcription factors driving {gene}'s regulon in the tumor context.\n"
            f"3. Run `comprehensive_gene_analysis` with `gene=\"{gene}\"` and `tcga_network=\"{cancer_type}\"` "
            f"for full domain analysis (oncogenic potential, druggability, clinical actionability, network vulnerability).\n"
            f"4. Summarize:\n"
            f"   - Regulatory role in the {cancer_type.upper()} tumor network\n"
            f"   - MoA breakdown: how many upstream regulators are activating vs. repressive\n"
            f"   - Top upstream regulators by PageRank\n"
            f"   - Druggability and clinical actionability assessment\n"
            f"   - Network vulnerability (is {gene} a critical hub in the tumor network?)\n"
        )
        return types.GetPromptResult(
            description=f"Tumor-context analysis of {gene} in TCGA {cancer_type.upper()}",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(type="text", text=text)
                )
            ]
        )

    elif name == "network_context_comparison":
        gene = args.get("gene", "").upper() or "MYC"
        cancer_type = args.get("cancer_type", "brca").lower()

        text = (
            f"Please compare the regulatory context of **{gene}** between the "
            f"population-averaged GREmLN epithelial network and the TCGA {cancer_type.upper()} "
            f"tumor-state network using the following steps:\n\n"
            f"1. Use `compare_network_contexts` with `gene=\"{gene}\"` and `cancer_type=\"{cancer_type}\"` "
            f"to retrieve conserved, population-averaged-only, and tumor-state-only regulator sets.\n"
            f"2. Summarize the results:\n"
            f"   - Conserved regulators (present in both networks)\n"
            f"   - Population-averaged-only regulators (present in GREmLN epithelial, absent in TCGA {cancer_type.upper()})\n"
            f"   - Tumor-state-only regulators (present in TCGA {cancer_type.upper()}, absent in GREmLN epithelial)\n"
            f"   - Conserved fraction (what proportion of regulators are shared?)\n"
            f"3. Provide a biological interpretation:\n"
            f"   - What does the conserved regulator set suggest about {gene}'s core regulatory program?\n"
            f"   - Are the tumor-state-only regulators known cancer drivers or oncogenes?\n"
            f"   - Note: GREmLN networks represent population-averaged mixed cell states "
            f"(not purely normal tissue) — interpret differences as reflecting both technical "
            f"(single-cell vs. bulk) and biological (cell state) factors.\n"
        )
        return types.GetPromptResult(
            description=f"Network context comparison of {gene}: GREmLN epithelial vs. TCGA {cancer_type.upper()}",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(type="text", text=text)
                )
            ]
        )

    elif name == "candidate_prioritization":
        gene = args.get("gene", "").upper() or "CTNNB1"
        cancer_type = args.get("cancer_type", "brca").lower()

        text = (
            f"Please run the regulatory candidate prioritization workflow for **{gene}** "
            f"in **TCGA {cancer_type.upper()}** using the following two-step process:\n\n"
            f"**Step 1 — Identify source-labeled candidates:**\n"
            f"Use `compare_network_contexts` with `gene=\"{gene}\"` and `cancer_type=\"{cancer_type}\"` "
            f"to retrieve the OncoKB-filtered candidate shortlist. Note the source label for each candidate "
            f"(TCGA-only, GREmLN-only, or Both) and MoA where available.\n\n"
            f"**Step 2 — Comprehensive domain analysis per candidate:**\n"
            f"For each candidate returned in Step 1, run `comprehensive_gene_analysis` with "
            f"source-driven network routing:\n"
            f"   - TCGA-only candidates → use `tcga_network=\"{cancer_type}\"`\n"
            f"   - GREmLN-only candidates → use `cell_type=\"epithelial_cell\"`\n"
            f"   - Both-source candidates → use `tcga_network=\"{cancer_type}\"` (tumor-state default)\n\n"
            f"**Step 3 — Summarize as a table:**\n"
            f"Present results with columns: "
            f"Candidate | Source | MoA | Oncogenic Potential | Druggability | "
            f"Clinical Actionability | Network Vulnerability | PageRank\n\n"
            f"**Step 4 — Interpret:**\n"
            f"   - Which candidates are highest priority for therapeutic follow-up "
            f"(activating regulators upstream of oncogenes, or repressive regulators upstream of tumor suppressors)?\n"
            f"   - Do TCGA-only candidates differ meaningfully from Both-source candidates in their domain profiles?\n"
        )
        return types.GetPromptResult(
            description=f"Candidate prioritization for {gene} in TCGA {cancer_type.upper()}",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(type="text", text=text)
                )
            ]
        )

    else:
        return types.GetPromptResult(
            description="Unknown prompt",
            messages=[
                types.PromptMessage(
                    role="user",
                    content=types.TextContent(
                        type="text",
                        text=f"Unknown prompt '{name}'. Available prompts: gene_deep_dive, cancer_biomarker_panel, cross_cell_comparison, tumor_context_analysis, network_context_comparison, candidate_prioritization."
                    )
                )
            ]
        )


def _format_markdown(gene: str, cell_type: str, result: dict, sections: list) -> str:
    """Format a run_analysis result as a markdown report."""
    include_all = "all" in sections

    lines = [f"# {gene} Regulatory Network Analysis — {cell_type.replace('_', ' ').title()}", ""]

    # --- Summary ---
    if include_all or "summary" in sections:
        network = result.get("network_analysis") or {}
        gene_sum = result.get("gene_analysis_summary") or {}
        key = result.get("key_insights") or {}

        lines += [
            "## Summary",
            "| Property | Value |",
            "|----------|-------|",
            f"| Gene | {gene} |",
            f"| Cell Type | {cell_type} |",
            f"| Regulatory Role | {network.get('regulatory_role') or gene_sum.get('regulatory_role', 'N/A')} |",
            f"| Regulators | {network.get('num_regulators', 'N/A')} |",
            f"| Targets | {network.get('num_targets', 'N/A')} |",
            f"| In Network | {network.get('in_network', 'N/A')} |",
            f"| Cancer Relevance | {key.get('cancer_relevance', 'N/A')} |",
            f"| Druggability | {key.get('druggability', 'N/A')} |",
            f"| Clinical Significance | {key.get('clinical_significance', 'N/A')} |",
            f"| Biomarker Potential | {key.get('biomarker_potential', 'N/A')} |",
            "",
        ]

    # --- Regulators ---
    if include_all or "regulators" in sections:
        reg_data = result.get("regulatory_analysis") or {}
        hub_regs = reg_data.get("hub_regulators") or []
        lines += [
            "## Top Regulators",
            "| Gene | Ensembl ID | Regulatory Strength |",
            "|------|-----------|---------------------|",
        ]
        if hub_regs:
            for r in hub_regs:
                lines.append(
                    f"| {r.get('gene_symbol', 'N/A')} | {r.get('ensembl_id', 'N/A')} | {r.get('regulatory_strength', 'N/A')} |"
                )
        else:
            lines.append("| N/A | — | — |")
        lines.append("")

    # --- Targets ---
    if include_all or "targets" in sections:
        tgt_data = result.get("target_analysis") or {}
        cascade = tgt_data.get("cascade_targets") or []
        lines += [
            "## Top Targets",
            "| Gene | Ensembl ID | Cascade Level |",
            "|------|-----------|--------------|",
        ]
        if cascade:
            for t in cascade:
                lines.append(
                    f"| {t.get('gene_symbol', 'N/A')} | {t.get('ensembl_id', 'N/A')} | {t.get('cascade_level', 'N/A')} |"
                )
        else:
            lines.append("| N/A | — | — |")
        lines.append("")

    # --- Pathways ---
    if include_all or "pathways" in sections:
        pw_data = result.get("pathway_enrichment") or {}
        pathways = pw_data.get("enriched_pathways") or []
        lines += [
            "## Enriched Pathways",
            "| Pathway | p-value | FDR | Genes Found | Genes Total |",
            "|---------|---------|-----|-------------|-------------|",
        ]
        if pathways and isinstance(pathways, list):
            for p in pathways:
                pval = p.get("p_value")
                fdr = p.get("fdr")
                lines.append(
                    f"| {p.get('pathway_name', 'N/A')} | "
                    f"{f'{pval:.2e}' if pval is not None else 'N/A'} | "
                    f"{f'{fdr:.2e}' if fdr is not None else 'N/A'} | "
                    f"{p.get('genes_found', 'N/A')} | "
                    f"{p.get('genes_total', 'N/A')} |"
                )
        else:
            lines.append("| N/A | — | — | — | — |")
        lines.append("")

    return "\n".join(lines)


def _format_csv(gene: str, cell_type: str, result: dict, sections: list) -> str:
    """Format a run_analysis result as CSV (three sections separated by blank lines)."""
    include_all = "all" in sections
    parts = []

    # --- Summary ---
    if include_all or "summary" in sections:
        network = result.get("network_analysis") or {}
        gene_sum = result.get("gene_analysis_summary") or {}
        key = result.get("key_insights") or {}
        rows = [
            "# Summary",
            "property,value",
            f"gene,{gene}",
            f"cell_type,{cell_type}",
            f"regulatory_role,{network.get('regulatory_role') or gene_sum.get('regulatory_role', '')}",
            f"num_regulators,{network.get('num_regulators', '')}",
            f"num_targets,{network.get('num_targets', '')}",
            f"in_network,{network.get('in_network', '')}",
            f"cancer_relevance,{key.get('cancer_relevance', '')}",
            f"druggability,{key.get('druggability', '')}",
            f"clinical_significance,{key.get('clinical_significance', '')}",
            f"biomarker_potential,{key.get('biomarker_potential', '')}",
        ]
        parts.append("\n".join(rows))

    # --- Regulators ---
    if include_all or "regulators" in sections:
        reg_data = result.get("regulatory_analysis") or {}
        hub_regs = reg_data.get("hub_regulators") or []
        rows = ["# Regulators", "gene,regulator,ensembl_id,regulatory_strength"]
        for r in hub_regs:
            rows.append(
                f"{gene},{r.get('gene_symbol', '')},{r.get('ensembl_id', '')},{r.get('regulatory_strength', '')}"
            )
        parts.append("\n".join(rows))

    # --- Targets ---
    if include_all or "targets" in sections:
        tgt_data = result.get("target_analysis") or {}
        cascade = tgt_data.get("cascade_targets") or []
        rows = ["# Targets", "gene,target,ensembl_id,cascade_level"]
        for t in cascade:
            rows.append(
                f"{gene},{t.get('gene_symbol', '')},{t.get('ensembl_id', '')},{t.get('cascade_level', '')}"
            )
        parts.append("\n".join(rows))

    # --- Pathways ---
    if include_all or "pathways" in sections:
        pw_data = result.get("pathway_enrichment") or {}
        pathways = pw_data.get("enriched_pathways") or []
        rows = ["# Pathways", "pathway_id,pathway_name,p_value,fdr,genes_found,genes_total"]
        if isinstance(pathways, list):
            for p in pathways:
                pval = p.get("p_value", "")
                fdr = p.get("fdr", "")
                name = str(p.get("pathway_name", "")).replace(",", ";")
                rows.append(
                    f"{p.get('pathway_id', '')},{name},{pval},{fdr},"
                    f"{p.get('genes_found', '')},{p.get('genes_total', '')}"
                )
        parts.append("\n".join(rows))

    return "\n\n".join(parts)


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
                    },
                    "tcga_network": {
                        "type": "string",
                        "enum": ["brca", "coad", "hnsc", "luad", "lusc", "ov", "prad", "ucec"],
                        "description": "TCGA cancer-type network to use for topology analysis (e.g. 'brca'). "
                                       "When provided, TCGA network topology is used instead of the GREmLN cell-type network. "
                                       "Useful for analyzing candidates identified from TCGA ARACNe networks."
                    },
                    "use_cache": {
                        "type": "boolean",
                        "description": "Return cached results if available (faster). Set to false to force fresh analysis.",
                        "default": True
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
                    },
                    "use_cache": {
                        "type": "boolean",
                        "description": "Return cached results for genes already analyzed recently. Set to false to force fresh analysis.",
                        "default": True
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
            - gene_neighbors: Immediate regulators and targets of a specific gene (includes MI score per edge)
            - network_stats: Summary statistics (genes, edges, density, etc.)

            Edge confidence filtering (confidence_level):
            - "all" (default): all edges, no filter
            - "medium": edges with MI score > 0.05
            - "high": edges with MI score > 0.1 AND bootstrap count >= 3

            Great for:
            - "What are the top regulators in epithelial cells?"
            - "How many high-confidence targets does TP53 have?"
            - "What regulates MYC with high confidence?"
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
                    "network_source": {
                        "type": "string",
                        "enum": ["cell_type", "tcga"],
                        "description": "'cell_type' (default) queries population-averaged GREmLN networks. 'tcga' queries tumor-state ARACNe networks — also set tcga_network.",
                        "default": "cell_type"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type network to query (used when network_source='cell_type')",
                        "default": "epithelial_cell"
                    },
                    "tcga_network": {
                        "type": "string",
                        "enum": ["brca", "coad", "hnsc", "luad", "lusc", "ov", "prad", "ucec"],
                        "description": "TCGA cancer type to query (required when network_source='tcga'). brca=breast, coad=colon, hnsc=head/neck, luad=lung adeno, lusc=lung squamous, ov=ovarian, prad=prostate, ucec=uterine."
                    },
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (required for gene_neighbors query)"
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of results to return for ranked queries",
                        "default": 10
                    },
                    "confidence_level": {
                        "type": "string",
                        "enum": ["all", "medium", "high"],
                        "description": "Edge confidence filter. 'all' = no filter (default), 'medium' = likelihood>0.05, 'high' = likelihood>0.1. For TCGA networks bootstrap counts are unavailable; 'high' uses likelihood only.",
                        "default": "all"
                    }
                },
                "required": ["query_type"]
            }
        ),
        Tool(
            name="list_prompts",
            description="List the available MCP prompt templates and how to use them. Call this when a user asks what prompts are available or how to get started.",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="export_results",
            description="Export gene analysis results as markdown or CSV for sharing and manuscripts.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol to analyze and export"
                    },
                    "format": {
                        "type": "string",
                        "enum": ["markdown", "csv"],
                        "description": "Output format: markdown (renders in Claude Desktop) or csv (for spreadsheets)",
                        "default": "markdown"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells"
                        ],
                        "default": "epithelial_cell"
                    },
                    "sections": {
                        "type": "array",
                        "items": {
                            "type": "string",
                            "enum": ["summary", "regulators", "targets", "pathways", "all"]
                        },
                        "description": "Which sections to include",
                        "default": ["all"]
                    }
                },
                "required": ["gene"]
            }
        ),
        Tool(
            name="find_master_regulators",
            description="""
            Identify which transcription factors drive a gene signature (reverse ARACNe analysis).

            Given a list of differentially expressed genes (e.g., from an RNA-seq experiment),
            finds which TFs in the regulatory network have regulons most enriched in that gene
            set. This is the core reverse-direction use case ARACNe networks were designed for.

            Uses Fisher's exact test to rank TFs by statistical enrichment significance.

            Returns the top N master regulators with:
            - Overlap count and regulon size
            - Enrichment score (fold enrichment over random)
            - p-value (Fisher's exact test, one-sided)
            - Overlapping gene symbols

            IMPORTANT - follow-up guidance: After calling this tool, if the user wants more
            detail on a specific returned TF, use query_network with query_type="gene_neighbors"
            (fast, <50ms). Do NOT automatically call comprehensive_gene_analysis on the results
            — that tool takes 15-60 seconds per gene and is not needed to interpret master
            regulator output.

            Example: "Which TFs drive this set of upregulated genes from my experiment?"
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_set": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of gene symbols (e.g., differentially expressed genes)"
                    },
                    "network_source": {
                        "type": "string",
                        "enum": ["cell_type", "tcga"],
                        "description": "'cell_type' (default) uses population-averaged GREmLN networks. 'tcga' uses tumor-state ARACNe networks — also set tcga_network.",
                        "default": "cell_type"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "Cell type network to use (when network_source='cell_type')",
                        "default": "epithelial_cell"
                    },
                    "tcga_network": {
                        "type": "string",
                        "enum": ["brca", "coad", "hnsc", "luad", "lusc", "ov", "prad", "ucec"],
                        "description": "TCGA cancer type (required when network_source='tcga')"
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of top master regulators to return",
                        "default": 10
                    }
                },
                "required": ["gene_set"]
            }
        ),
        Tool(
            name="compare_network_contexts",
            description="""
            Compare regulatory wiring for a gene across population-averaged (GREmLN) and
            tumor-state (TCGA) network contexts in a single call.

            Instead of making two separate query_network calls and comparing manually, this
            tool runs both queries internally, computes regulator and target overlap, and
            returns a structured comparison showing:
            - Which regulators are conserved across both contexts
            - Which regulators appear only in population-averaged wiring (lost in tumor)
            - Which regulators appear only in tumor-state wiring (emerge in tumor)
            - A regulatory rewiring classification (low / moderate / high)

            Rewiring is classified by Jaccard overlap of regulator sets:
            - >= 0.6 conserved fraction → "low" rewiring (stable regulatory program)
            - 0.3–0.6 → "moderate" rewiring
            - < 0.3 → "high" rewiring (substantially different tumor-state program)

            Example questions this answers:
            - "Which MYC regulators are conserved in colorectal vs. epithelial context?"
            - "How different is ESR1's regulatory wiring in BRCA vs. normal epithelium?"
            - "What regulators emerge specifically in ovarian cancer for TP53?"
            """,
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g. 'MYC', 'TP53', 'ESR1')"
                    },
                    "cancer_type": {
                        "type": "string",
                        "enum": ["brca", "coad", "hnsc", "luad", "lusc", "ov", "prad", "ucec"],
                        "description": "TCGA cancer type for the tumor-state context"
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": [
                            "epithelial_cell", "cd14_monocytes", "cd16_monocytes",
                            "cd20_b_cells", "cd4_t_cells", "cd8_t_cells",
                            "erythrocytes", "nk_cells", "nkt_cells",
                            "monocyte-derived_dendritic_cells"
                        ],
                        "description": "GREmLN reference context (default: epithelial_cell)",
                        "default": "epithelial_cell"
                    }
                },
                "required": ["gene", "cancer_type"]
            }
        ),
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
            tcga_network = arguments.get("tcga_network", None)
            use_cache = arguments.get("use_cache", True)

            cache_key = (gene, cell_type, analysis_depth, tcga_network)
            if use_cache and cache_key in _result_cache:
                cached_ts, cached_result = _result_cache[cache_key]
                if time.time() - cached_ts < _CACHE_TTL_SECONDS:
                    logger.info(f"Cache hit for {gene} ({cell_type}, {analysis_depth}, tcga={tcga_network})")
                    cached_result["workflow_info"]["source"] = "cache"
                    return [TextContent(type="text", text=json.dumps(cached_result, indent=2))]

            logger.info(f"Starting comprehensive analysis for {gene} using LangGraph workflow")
            logger.info(f"Parameters: cell_type={cell_type}, analysis_depth={analysis_depth}, tcga_network={tcga_network}")

            start_time = time.time()

            # Simple wrapper around workflow
            result = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                tcga_network=tcga_network,
                analysis_depth=analysis_depth
            )

            execution_time = time.time() - start_time
            logger.info(f"Analysis for {gene} completed in {execution_time:.2f} seconds")

            # Add minimal MCP metadata
            result["workflow_info"] = {
                "workflow_type": "langgraph",
                "server_version": "regnetagents-mcp-v2.1",
                "execution_mode": "thin_wrapper",
                "execution_time_seconds": round(execution_time, 2),
                "source": "live"
            }

            _result_cache[cache_key] = (time.time(), result)

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
            use_cache = arguments.get("use_cache", True)

            logger.info(f"Starting multi-gene analysis for {len(genes)} genes")

            start_time = time.time()

            # Check cache per gene; only run analysis for misses
            gene_results = {}
            genes_to_run = []
            for gene in genes:
                cache_key = (gene, cell_type, analysis_depth)
                if use_cache and cache_key in _result_cache:
                    cached_ts, cached_result = _result_cache[cache_key]
                    if time.time() - cached_ts < _CACHE_TTL_SECONDS:
                        logger.info(f"Cache hit for {gene}")
                        cached_result.setdefault("workflow_info", {})["source"] = "cache"
                        gene_results[gene] = cached_result
                        continue
                genes_to_run.append(gene)

            if genes_to_run:
                tasks = [
                    workflow.run_analysis(gene=gene, cell_type=cell_type, analysis_depth=analysis_depth)
                    for gene in genes_to_run
                ]
                logger.info(f"Executing {len(genes_to_run)} gene analyses in parallel...")
                live_results = await asyncio.gather(*tasks, return_exceptions=True)
                for gene, result in zip(genes_to_run, live_results):
                    if not isinstance(result, Exception):
                        result.setdefault("workflow_info", {})["source"] = "live"
                        _result_cache[(gene, cell_type, analysis_depth)] = (time.time(), result)
                    gene_results[gene] = result

            results = [gene_results[gene] for gene in genes]

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
            network_source = arguments.get("network_source", "cell_type")
            cell_type = arguments.get("cell_type", "epithelial_cell")
            tcga_network = arguments.get("tcga_network", None)
            gene = arguments.get("gene", None)
            top_n = arguments.get("top_n", 10)
            confidence_level = arguments.get("confidence_level", "all")

            context = f"tcga/{tcga_network}" if network_source == "tcga" else cell_type
            logger.info(f"Network query: {query_type} in {context}" +
                        (f" for {gene}" if gene else "") +
                        (f" [{confidence_level} confidence]" if confidence_level != "all" else ""))

            result = workflow.modeling_agent.query_network(
                query_type=query_type,
                cell_type=cell_type,
                gene=gene,
                top_n=top_n,
                confidence_level=confidence_level,
                network_source=network_source,
                tcga_network=tcga_network,
            )

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "list_prompts":
            result = {
                "available_prompts": [
                    {
                        "name": "gene_deep_dive",
                        "description": "Comprehensive guided analysis of a single gene.",
                        "arguments": {
                            "gene": "required — gene symbol, e.g. TP53, MYC, BRCA1",
                            "cell_type": "optional — default: epithelial_cell"
                        },
                        "what_it_does": (
                            "Validates the gene, runs full comprehensive analysis, "
                            "queries immediate network neighbors, then summarizes "
                            "regulatory role, top regulators, top targets, key pathways, "
                            "and clinical/cancer relevance."
                        ),
                        "example": "Use the gene_deep_dive prompt for TP53 in epithelial cells"
                    },
                    {
                        "name": "cancer_biomarker_panel",
                        "description": "Parallel analysis of a pre-defined cancer gene panel.",
                        "arguments": {
                            "cancer_type": "required — one of: colorectal, breast, lung, prostate, general",
                            "cell_type": "optional — default: epithelial_cell"
                        },
                        "panels": {
                            "colorectal": ["APC", "TP53", "KRAS", "MYC", "CTNNB1"],
                            "breast":     ["BRCA1", "BRCA2", "TP53", "ERBB2", "ESR1"],
                            "lung":       ["KRAS", "EGFR", "TP53", "ALK", "STK11"],
                            "prostate":   ["AR", "TP53", "PTEN", "MYC", "ERG"],
                            "general":    ["TP53", "MYC", "KRAS", "BRCA1", "EGFR"]
                        },
                        "what_it_does": (
                            "Runs multi_gene_analysis on the selected panel in parallel, "
                            "identifies the most central regulator, flags genes absent "
                            "from the network, and summarizes the panel's regulatory landscape."
                        ),
                        "example": "Use the cancer_biomarker_panel prompt for colorectal cancer"
                    },
                    {
                        "name": "cross_cell_comparison",
                        "description": "Compare a gene across all 10 cell-type networks.",
                        "arguments": {
                            "gene": "required — gene symbol, e.g. TP53, MYC"
                        },
                        "what_it_does": (
                            "Runs cross_cell_comparison, summarizes results in a table "
                            "(cell type / regulatory role / regulators / targets / in network), "
                            "highlights the most active cell type, and compares immune cells "
                            "vs. epithelial cells."
                        ),
                        "example": "Use the cross_cell_comparison prompt for MYC"
                    },
                    {
                        "name": "tumor_context_analysis",
                        "description": "Full domain analysis of a gene against a TCGA tumor-state network.",
                        "arguments": {
                            "gene": "required — gene symbol, e.g. YAP1, ESR1",
                            "cancer_type": "required — one of: brca, coad, hnsc, luad, lusc, ov, prad, ucec"
                        },
                        "what_it_does": (
                            "Queries tumor-state regulatory neighbors (with MoA), finds master regulators "
                            "in tumor context, runs comprehensive_gene_analysis with tcga_network, and "
                            "summarizes regulatory role, MoA breakdown, druggability, and clinical actionability."
                        ),
                        "example": "Use the tumor_context_analysis prompt for YAP1 in brca"
                    },
                    {
                        "name": "network_context_comparison",
                        "description": "Compare a gene's regulatory context between population-averaged GREmLN and TCGA tumor-state networks.",
                        "arguments": {
                            "gene": "required — gene symbol, e.g. MYC, TP53",
                            "cancer_type": "required — one of: brca, coad, hnsc, luad, lusc, ov, prad, ucec"
                        },
                        "what_it_does": (
                            "Runs compare_network_contexts, returns conserved regulators, "
                            "population-averaged-only regulators, and tumor-state-only regulators "
                            "with biological interpretation. Frames output as population-averaged "
                            "vs. tumor-state context (not 'rewiring')."
                        ),
                        "example": "Use the network_context_comparison prompt for MYC in coad"
                    },
                    {
                        "name": "candidate_prioritization",
                        "description": "Two-step regulatory candidate prioritization workflow.",
                        "arguments": {
                            "gene": "required — focal gene, e.g. CTNNB1, MYC",
                            "cancer_type": "required — one of: brca, coad, hnsc, luad, lusc, ov, prad, ucec"
                        },
                        "what_it_does": (
                            "Step 1: compare_network_contexts → source-labeled OncoKB-filtered candidate shortlist. "
                            "Step 2: comprehensive_gene_analysis per candidate with source-driven network routing "
                            "(TCGA-only→tcga_network, GREmLN-only→cell_type, Both→tcga_network). "
                            "Returns structured summary table: candidate | source | MoA | oncogenic potential | "
                            "druggability | clinical actionability | network vulnerability | PageRank."
                        ),
                        "example": "Use the candidate_prioritization prompt for CTNNB1 in brca"
                    }
                ],
                "how_to_use": (
                    "These prompts scaffold common workflows. "
                    "You can trigger them by describing what you want in plain language — "
                    "for example: 'Do a deep dive on TP53', "
                    "'Analyze the colorectal cancer panel', "
                    "'Compare BRCA1 across all cell types', "
                    "'Analyze YAP1 in the BRCA tumor network', or "
                    "'Prioritize CTNNB1 regulatory candidates in BRCA'."
                )
            }
            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "find_master_regulators":
            gene_set = arguments["gene_set"]
            network_source = arguments.get("network_source", "cell_type")
            cell_type = arguments.get("cell_type", "epithelial_cell")
            tcga_network = arguments.get("tcga_network", None)
            top_n = arguments.get("top_n", 10)

            context = f"tcga/{tcga_network}" if network_source == "tcga" else cell_type
            logger.info(f"Finding master regulators for {len(gene_set)}-gene set in {context}")

            result = await asyncio.to_thread(
                workflow.modeling_agent.find_master_regulators,
                gene_set,
                cell_type,
                top_n,
                network_source,
                tcga_network,
            )

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "export_results":
            gene = arguments["gene"]
            fmt = arguments.get("format", "markdown")
            cell_type = arguments.get("cell_type", "epithelial_cell")
            sections = arguments.get("sections", ["all"])

            logger.info(f"Exporting results for {gene} as {fmt} (cell_type={cell_type})")

            result = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                analysis_depth="comprehensive"
            )

            if fmt == "markdown":
                text = _format_markdown(gene.upper(), cell_type, result, sections)
            else:
                text = _format_csv(gene.upper(), cell_type, result, sections)

            return [TextContent(type="text", text=text)]

        elif name == "compare_network_contexts":
            from regnetagents.context_comparison import compare_network_contexts
            gene        = arguments["gene"]
            cancer_type = arguments["cancer_type"]
            cell_type   = arguments.get("cell_type", "epithelial_cell")

            logger.info(f"Comparing network contexts for {gene}: {cell_type} vs tcga/{cancer_type}")

            result = compare_network_contexts(
                agent=workflow.modeling_agent,
                gene=gene,
                cancer_type=cancer_type,
                cell_type=cell_type,
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