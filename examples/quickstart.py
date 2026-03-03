#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RegNetAgents Quick Start
========================

Minimal single-gene example — the simplest way to use RegNetAgents
programmatically without Claude Desktop, Ollama, or any API keys.

What this shows:
----------------
- Import and initialize RegNetAgentsWorkflow
- Analyze one gene in one cell type
- Read the key output fields

Expected runtime: 2-5 seconds (rule-based mode, no LLM required)

Usage:
------
    python examples/quickstart.py

For multi-gene panel analysis, see:
    python demo_biomarker_analysis.py
"""

import asyncio
import os
import sys

# Allow running from the examples/ subdirectory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

# ─── Configure your analysis here ──────────────────────────────────────────

GENE = "TP53"           # Gene symbol to analyze
CELL_TYPE = "epithelial_cell"  # Available cell types: epithelial_cell,
                               # cd4_t_cells, cd8_t_cells, cd14_monocytes,
                               # cd16_monocytes, cd20_b_cells, nk_cells,
                               # nkt_cells, erythrocytes,
                               # monocyte-derived_dendritic_cells

# ───────────────────────────────────────────────────────────────────────────


async def main():
    print(f"RegNetAgents Quick Start — analyzing {GENE} in {CELL_TYPE}")
    print("-" * 60)

    # Initialize the workflow (loads network cache from models/networks/)
    workflow = RegNetAgentsWorkflow()

    # Run comprehensive analysis for the gene
    # analysis_depth options:
    #   "focused"       — network positions and pathway enrichment only
    #   "basic"         — adds rule-based domain insights
    #   "comprehensive" — adds therapeutic target prioritization (default)
    result = await workflow.run_analysis(
        gene=GENE,
        cell_type=CELL_TYPE,
        analysis_depth="comprehensive",
    )

    # ── Network topology ──────────────────────────────────────────────────
    network = result.get("network_analysis", {})
    print(f"\nNetwork topology ({CELL_TYPE}):")
    print(f"  Regulatory role : {network.get('regulatory_role', 'unknown')}")
    print(f"  Upstream regulators : {network.get('num_regulators', 0)}")
    print(f"  Downstream targets  : {network.get('num_targets', 0)}")

    # ── Top upstream regulator by PageRank ────────────────────────────────
    targets = result.get("therapeutic_target_prioritization", {})
    ranked = targets.get("ranked_regulators", [])
    if ranked:
        top = ranked[0]
        pagerank = top.get("centrality_metrics", {}).get("pagerank", 0)
        print(f"  Top upstream regulator : {top['regulator']} "
              f"(PageRank: {pagerank:.3f})")

    # ── Pathway enrichment ────────────────────────────────────────────────
    pathway = result.get("pathway_enrichment", {})
    n_pathways = pathway.get("summary", {}).get("total_pathways", 0)
    n_sig = pathway.get("summary", {}).get("significant_pathways", 0)
    print(f"\nPathway enrichment:")
    print(f"  Total pathways     : {n_pathways}")
    print(f"  Significant (FDR<0.05) : {n_sig}")

    # ── Domain insights (rule-based) ──────────────────────────────────────
    domain = result.get("domain_analysis", {})
    cancer = domain.get("cancer_analysis", {}).get("insights", {})
    drug = domain.get("drug_analysis", {}).get("insights", {})
    if cancer or drug:
        print(f"\nDomain insights (rule-based):")
        if cancer:
            print(f"  Oncogenic potential    : "
                  f"{cancer.get('oncogenic_potential', 'unknown')}")
        if drug:
            print(f"  Druggability           : "
                  f"{drug.get('druggability_assessment', 'unknown')}")

    print("\n" + "-" * 60)
    print("Done.")
    print("\nNext steps:")
    print("  Multi-gene analysis  : python demo_biomarker_analysis.py")
    print("  MCP / Claude Desktop : see docs/REGNETAGENTS_MCP_SETUP.md")


if __name__ == "__main__":
    asyncio.run(main())
