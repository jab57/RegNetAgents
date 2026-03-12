#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RegNetAgents TCGA Query Example
================================

Query tumor-state ARACNe networks from TCGA alongside the population-averaged
GREmLN networks already bundled with RegNetAgents.

What this shows:
----------------
- How to use network_source='tcga' with query_network()
- How to use network_source='tcga' with find_master_regulators()
- The unique MoA (mode-of-action) field present in TCGA edges
- Side-by-side comparison: population-averaged vs. tumor-state context

Prerequisites:
--------------
1. Download the 8 TCGA ARACNe CSVs from Figshare:
       https://figshare.com/s/5d1ffd9f8b2e86e37ed6

2. Place each file (renamed to network.csv) at:
       models/networks/tcga/{cancer_type}/network.csv
   e.g. models/networks/tcga/brca/network.csv

3. Build the PKL caches:
       python scripts/build_tcga_cache.py --all

Usage:
------
    python examples/tcga_query.py
"""

import asyncio
import os
import sys

# Allow running from the examples/ subdirectory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

# ─── Configure your analysis here ──────────────────────────────────────────

GENE = "TP53"                  # Gene to query
TCGA_CANCER = "brca"           # TCGA cancer type: brca, coad, hnsc, luad,
                               # lusc, ov, prad, ucec
CELL_TYPE = "epithelial_cell"  # GREmLN population-averaged reference

GENE_SIGNATURE = [             # Example DEG signature for master regulator test
    "CDKN1A", "MDM2", "BAX", "GADD45A", "PUMA",
    "BBC3",   "CCNG1", "SESN2", "DDB2",  "POLK",
]

# ───────────────────────────────────────────────────────────────────────────


def _check_tcga_available(workflow, cancer_type: str) -> bool:
    return bool(workflow.tcga_cache.tcga_indices.get(cancer_type))


async def main():
    print("RegNetAgents TCGA Query Example")
    print("=" * 60)

    workflow = RegNetAgentsWorkflow()
    agent = workflow.modeling_agent

    if not _check_tcga_available(workflow, TCGA_CANCER):
        print(
            f"\nNo TCGA cache found for '{TCGA_CANCER}'.\n"
            "To use this example:\n"
            "  1. Download CSVs from https://figshare.com/s/5d1ffd9f8b2e86e37ed6\n"
            "  2. Place each as  models/networks/tcga/{cancer_type}/network.csv\n"
            "  3. Run:  python scripts/build_tcga_cache.py --all\n"
        )
        return

    # ── 1. Network stats comparison ───────────────────────────────────────
    print(f"\n1. Network stats — {CELL_TYPE} (population-averaged) vs TCGA {TCGA_CANCER}")
    print("-" * 60)

    gremln_stats = agent.query_network("network_stats", cell_type=CELL_TYPE)
    tcga_stats   = agent.query_network("network_stats",
                                       network_source="tcga",
                                       tcga_network=TCGA_CANCER)

    print(f"  {'Metric':<22} {'GREmLN (pop-avg)':>18}  {'TCGA ' + TCGA_CANCER.upper():>18}")
    print(f"  {'-'*22} {'-'*18}  {'-'*18}")
    for key in ("num_genes", "num_edges", "num_regulons"):
        g = gremln_stats.get(key, "n/a")
        t = tcga_stats.get(key, "n/a")
        print(f"  {key:<22} {g:>18,}  {t:>18,}")

    # ── 2. Gene neighbors with MoA (TCGA only) ───────────────────────────
    print(f"\n2. {GENE} targets in TCGA {TCGA_CANCER} (with Mode of Action)")
    print("-" * 60)

    tcga_neighbors = agent.query_network(
        "gene_neighbors",
        gene=GENE,
        network_source="tcga",
        tcga_network=TCGA_CANCER,
    )

    if tcga_neighbors.get("error"):
        print(f"  {tcga_neighbors['message']}")
    else:
        targets = tcga_neighbors.get("targets", [])
        print(f"  {GENE} has {tcga_neighbors['num_targets']} targets in TCGA {TCGA_CANCER}")
        if targets:
            print(f"  Top 5 targets (likelihood-sorted):")
            sorted_targets = sorted(targets,
                                    key=lambda x: x.get("likelihood", 0),
                                    reverse=True)
            for t in sorted_targets[:5]:
                moa = t.get("moa")
                moa_label = (
                    "activation" if moa == 1.0 else
                    "repression" if moa == -1.0 else
                    "unknown"
                )
                print(
                    f"    {t['gene']:<12}  "
                    f"likelihood={t.get('likelihood', 0):.4f}  "
                    f"MoA={moa_label}"
                )

    # ── 3. Side-by-side gene neighbors (GREmLN vs TCGA) ──────────────────
    print(f"\n3. {GENE} in {CELL_TYPE} vs TCGA {TCGA_CANCER} — regulator counts")
    print("-" * 60)

    gremln_nb = agent.query_network(
        "gene_neighbors", gene=GENE, cell_type=CELL_TYPE
    )
    print(
        f"  GREmLN ({CELL_TYPE}) : "
        f"{gremln_nb.get('num_targets', 0)} targets, "
        f"{gremln_nb.get('num_regulators', 0)} regulators"
    )
    if not tcga_neighbors.get("error"):
        print(
            f"  TCGA {TCGA_CANCER.upper():<17}: "
            f"{tcga_neighbors.get('num_targets', 0)} targets, "
            f"{tcga_neighbors.get('num_regulators', 0)} regulators"
        )

    # ── 4. Master regulators in TCGA network ─────────────────────────────
    print(f"\n4. Master regulators for {len(GENE_SIGNATURE)}-gene signature in TCGA {TCGA_CANCER}")
    print("-" * 60)

    mr_result = agent.find_master_regulators(
        gene_set=GENE_SIGNATURE,
        network_source="tcga",
        tcga_network=TCGA_CANCER,
        top_n=5,
    )

    if mr_result.get("error"):
        print(f"  {mr_result['message']}")
    else:
        summary = mr_result.get("query_summary", {})
        print(
            f"  Signature: {summary.get('genes_found_in_network')}/{summary.get('gene_set_size')} "
            f"genes found in TCGA {TCGA_CANCER} network"
        )
        print(f"  Top master regulators:")
        for r in mr_result.get("master_regulators", []):
            print(
                f"    #{r['rank']}  {r['gene']:<12}  "
                f"overlap={r['overlap_count']}/{r['regulon_size']}  "
                f"p={r['p_value']:.2e}"
            )

    print("\n" + "=" * 60)
    print("Done.")
    print("\nNext steps:")
    print("  Compare two contexts : see issue #17 (compare_network_contexts)")
    print("  MCP / Claude Desktop : see docs/REGNETAGENTS_MCP_SETUP.md")


if __name__ == "__main__":
    asyncio.run(main())
