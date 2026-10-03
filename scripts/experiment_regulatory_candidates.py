#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
experiment_regulatory_candidates.py
====================================
Cross-context regulatory candidate identification for the RegNetAgents research
paper (Briefings in Bioinformatics).

For each focal gene in BRCA and COAD, queries both the GREmLN epithelial
network and the TCGA tumor network, then generates a source-labeled candidate
therapeutic target list filtered by IntOGen cancer driver annotation.

Reference dataset (committed to the repo, no download):
  - regnetagents/reference_data/intogen_drivers.tsv : IntOGen Compendium of
    Mutational Cancer Driver Genes, release 2024.09.20 (CC0 1.0), pan-cancer
    (all genes, any role). v1/v2 of the paper used IntOGen; v3 uses IntOGen.

Outputs:
  results/
  - experiment_rewiring_results.json           : full statistics
  - target_list_brca.png                       : candidate counts by source, BRCA (NAR Fig 2A)
  - target_list_coad.png                       : candidate counts by source, COAD (NAR Fig 2B)
  - experiment_rewiring_barchart_brca.png      : regulator count bar chart, BRCA
  - experiment_rewiring_barchart_coad.png      : regulator count bar chart, COAD

  supplementary/  (named directly in the paper's Data Availability section)
  - table_s1_brca_candidates.csv                : source-labeled candidate list, BRCA
  - table_s2_coad_candidates.csv                : source-labeled candidate list, COAD
  - table_s3_tier_specificity.csv               : focal genes vs random genes, per tier
                                                  (the "tier_specificity" key of the JSON
                                                  above holds the per-gene detail)

  manuscript/  (NAR paper figures — overwrite in place)
  - figure_heatmap_brca.png                    : OR enrichment heatmap, BRCA (NAR Fig 3A)
  - figure_heatmap_coad.png                    : OR enrichment heatmap, COAD (NAR Fig 3B)
  - figure_negcontrol_brca.png                 : negative controls, BRCA (NAR Fig 4A)
  - figure_negcontrol_coad.png                 : negative controls, COAD (NAR Fig 4B)

Background: all genes in each cancer type's TCGA network (symbol-native PKL).

Usage:
  python scripts/experiment_regulatory_candidates.py

Dependencies (all in requirements.txt): scipy, numpy, matplotlib, seaborn
"""

import csv
import json
import math
import os
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

# Allow running from scripts/ or repo root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents.context_comparison import compare_network_contexts
from regnetagents.driver_gene_client import load_driver_roles
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

# ── Configuration ──────────────────────────────────────────────────────────────

BRCA_GENES  = [
    "TP53", "MYC", "CTNNB1", "CCND1",          # original panel
    "BRCA1", "BRCA2", "PIK3CA", "PTEN",         # BRCA hallmarks
    "RB1", "ERBB2", "ESR1", "GATA3",            # additional BRCA drivers
    "CDH1",                                     # lobular BRCA driver
]
COAD_GENES  = [
    "TP53", "MYC", "CTNNB1", "CCND1",           # original panel
    "KRAS", "APC", "SMAD4",                      # original COAD panel
    "BRAF", "PIK3CA", "PTEN",                    # additional CRC drivers
    "FBXW7", "TCF7L2", "RNF43",                 # additional CRC drivers
]
HOUSEKEEPING_GENES = ["ACTB", "GAPDH", "HPRT1", "LDHA", "TUBB"]
# PABPC1 was dropped in v3: it is an IntOGen driver (LIHC, LUSC), so it no longer
# meets the "non-driver" definition of this panel.
NEUTRAL_GENES      = ["FASN", "PCNA", "PKM", "VIM"]
CELL_TYPE   = "epithelial_cell"
N_PERMUTATIONS = 1000
RANDOM_SEED    = 42

RESULTS_DIR    = "results"
MANUSCRIPT_DIR = "manuscript"
SUPPLEMENTARY_DIR = "supplementary"
# Maps cancer_type -> supplementary table number, per the paper's Data Availability
# section (only BRCA/COAD have a named supplementary CSV; other cancer types fall
# back to writing under RESULTS_DIR).
SUPPLEMENTARY_TABLE_NUMBER = {"brca": 1, "coad": 2}

# Tier specificity (Table S3): focal panel vs random non-IntOGen genes
TIER_NULL_N_RANDOM     = 300
TIER_NULL_SEEDS        = [RANDOM_SEED, 1, 2, 3, 4]
TIER_NULL_RESAMPLES    = 2000   # random panels drawn per seed for the Stouffer Z null
TIER_LOGOR_PERMUTATIONS = 5000  # label permutations per seed for the log-OR test
# Floor for Stouffer Z: a Fisher p that underflows to 0 would give Z=+inf.
TIER_P_FLOOR           = 1e-15
TIER_SPECIFICITY_TABLE = "table_s3_tier_specificity.csv"

# ── Reference set loader ───────────────────────────────────────────────────────

def load_intogen_roles() -> dict:
    """Return {SYMBOL: role} for every gene in the committed IntOGen snapshot.

    Pan-cancer, any role (oncogene / tumor_suppressor / mixed / ambiguous); the
    reference set is the full key set. Uses the same loader as the tool's driver
    annotation (regnetagents.driver_gene_client).
    """
    return load_driver_roles()


# ── Statistics ─────────────────────────────────────────────────────────────────

def fisher_enrichment(query: set, reference: set, background: set) -> dict:
    """
    One-tailed Fisher's exact test: is query enriched in reference?

    Contingency table:
                    in ref    not in ref
    in query          a           b
    not in query      c           d
    """
    bg = background | query          # ensure query is a subset of background
    not_query = bg - query
    a = len(query & reference)
    b = len(query - reference)
    c = len(not_query & reference)
    d = len(not_query - reference)
    or_val, p_val = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    or_val = float(or_val)  # type: ignore[arg-type]
    p_val  = float(p_val)   # type: ignore[arg-type]
    return {
        "a": a, "b": b, "c": c, "d": d,
        "odds_ratio":  round(or_val, 4),
        # Full precision: rounding to 6 decimals turned p < 5e-7 into 0, which made
        # stouffer_z return +inf and distorted BH-FDR.
        "p_value":     p_val,
        "query_size":  len(query),
        "ref_overlap": a,
    }


def permutation_test(
    query: set,
    reference: set,
    background: set,
    n: int = 1000,
    seed: int = 42,
) -> dict:
    """
    Empirical null: draw n random gene sets of the same size as query from
    background and compute their odds ratios. Returns empirical p-value.
    """
    rng = random.Random(seed)
    bg_list = sorted(background)
    size    = min(len(query), len(bg_list))
    obs_or  = fisher_enrichment(query, reference, background)["odds_ratio"]
    null_ors = [
        fisher_enrichment(
            set(rng.sample(bg_list, size)), reference, background
        )["odds_ratio"]
        for _ in range(n)
    ]
    emp_p = sum(1 for x in null_ors if x >= obs_or) / n
    return {
        "observed_or":   obs_or,
        "empirical_p":   round(emp_p, 4),
        "null_or_mean":  round(float(np.mean(null_ors)), 4),
        "null_or_std":   round(float(np.std(null_ors)),  4),
    }


def bh_fdr(p_values: list) -> list:
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    n = len(p_values)
    if n == 0:
        return []
    order  = sorted(range(n), key=lambda i: p_values[i])
    adj    = [0.0] * n
    prev   = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        adj[i] = min(p_values[i] * n / (n - rank + 1), prev)
        prev   = adj[i]
    return [min(v, 1.0) for v in adj]


def stouffer_z(p_values: list, weights: list = None) -> dict:
    """Stouffer's weighted combined Z-score from a list of one-tailed p-values."""
    if not p_values:
        return {"combined_z": float("nan"), "combined_p": float("nan")}
    # isf(p) == ppf(1 - p) but keeps precision for tiny p; floor keeps p = 0 finite.
    zs = np.array([stats.norm.isf(max(p, 1e-300)) if p < 1.0 else 0.0 for p in p_values])
    w  = np.ones(len(zs)) if weights is None else np.array(weights, dtype=float)
    z  = float(np.dot(w, zs) / np.sqrt(np.dot(w, w)))
    p  = float(1 - stats.norm.cdf(z))
    return {"combined_z": round(z, 4), "combined_p": round(p, 6)}


# ── Figures ────────────────────────────────────────────────────────────────────

def plot_workflow_figure() -> None:
    """Figure 1: NAR paper analysis pipeline schematic."""
    from matplotlib.patches import FancyBboxPatch

    C_INPUT    = "#E8F4F8"   # light blue
    C_AGENT    = "#B8E6B8"   # light green
    C_CLASS    = "#FFE6B8"   # light orange
    C_CAND     = "#E8D4F8"   # light purple
    C_ENRICH   = "#FFF0B8"   # light yellow
    C_OUTPUT   = "#F4E8E8"   # light pink

    fig, ax = plt.subplots(figsize=(7, 13))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 18)
    ax.axis("off")

    def box(y, h, color, title, lines=(), title_size=11):
        patch = FancyBboxPatch((1, y), 8, h, boxstyle="round,pad=0.15",
                               facecolor=color, edgecolor="#555555", linewidth=1.2)
        ax.add_patch(patch)
        ax.text(5, y + h - 0.38, title, ha="center", va="top",
                fontsize=title_size, fontweight="bold")
        for i, line in enumerate(lines):
            ax.text(5, y + h - 0.78 - i * 0.42, line, ha="center", va="top",
                    fontsize=9, color="#333333")

    def arrow(y_top, y_bot, label=""):
        ax.annotate("", xy=(5, y_bot), xytext=(5, y_top),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1.5))
        if label:
            ax.text(5.25, (y_top + y_bot) / 2, label, ha="left", va="center",
                    fontsize=8, color="#555555", style="italic")

    def sub_boxes(y, h, labels, colors):
        n = len(labels)
        w = 7.0 / n
        for i, (lbl, c) in enumerate(zip(labels, colors)):
            x0 = 1.5 + i * w
            patch = FancyBboxPatch((x0, y), w - 0.15, h,
                                   boxstyle="round,pad=0.1",
                                   facecolor=c, edgecolor="#888888", linewidth=0.8)
            ax.add_patch(patch)
            ax.text(x0 + (w - 0.15) / 2, y + h / 2, lbl, ha="center", va="center",
                    fontsize=9, fontweight="bold")

    # ── boxes top→bottom ───────────────────────────────────────────────────────
    box(15.5, 1.9, C_INPUT, "INPUT",
        ("Focal gene  ·  GREmLN cell type  ·  TCGA cancer type",))

    arrow(15.5, 14.7)

    box(12.9, 1.9, C_AGENT, "compare_network_contexts  (RegNetAgents)",
        ("Query focal gene in GREmLN epithelial_cell",
         "and TCGA ARACNe tumor network"))

    arrow(12.9, 12.1)

    box(9.8, 2.8, C_CLASS, "REGULATOR CLASSIFICATION", ())
    sub_boxes(10.15, 0.9,
              ["Both", "GREmLN-only", "TCGA-only"],
              ["#c8e6c8", "#ffe0b2", "#ffcccc"])
    ax.text(5, 10.05,
            "Context-specificity = 1 − J,  J = |GREmLN ∩ TCGA| / |GREmLN ∪ TCGA|",
            ha="center", va="top", fontsize=8.5, color="#444444")

    arrow(9.8, 9.0, "TCGA-only + GREmLN-only candidates")

    box(7.0, 2.3, C_CAND, "SOURCE-LABELED CANDIDATE LIST",
        ("Filter all candidates against IntOGen",
         "Source (TCGA-only / GREmLN-only / Both)  ·  IntOGen role",
         "MoA direction (activating / repressive)  for TCGA-only"))

    arrow(7.0, 6.2, "TCGA-only candidates")

    box(4.2, 2.1, C_ENRICH, "ENRICHMENT VALIDATION",
        ("Fisher's exact test  ·  IntOGen reference",
         "Permutation control (n=1,000)  ·  BH-FDR correction"))

    arrow(4.2, 3.4)

    box(1.5, 2.0, C_OUTPUT, "OUTPUT",
        ("Source-labeled candidate list  ·  OR · BH-FDR per gene",
         "Stouffer Z across panel"))

    plt.tight_layout()
    out = os.path.join(MANUSCRIPT_DIR, "figure_workflow.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"            {out}  (NAR Fig 1)")



def plot_or_heatmap(
    results: dict, genes: list, cancer_type: str
) -> None:
    """Heatmap: rows = focal genes, columns = reference sets, values = OR."""
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    ref_keys   = ["intogen"]
    ref_labels = ["IntOGen\n(pan-cancer)"]

    or_matrix    = []
    annot_matrix = []
    for gene in genes:
        row_or, row_annot = [], []
        for ref in ref_keys:
            enr   = results[gene]["enrichment"].get(ref, {})
            or_v  = enr.get("odds_ratio", float("nan"))
            adj_p = enr.get("fdr_adjusted_p", enr.get("p_value", 1.0))
            star  = "**" if adj_p < 0.01 else ("*" if adj_p < 0.05 else "")
            safe_or = or_v if (or_v != float("inf") and not math.isnan(or_v)) else 5.0
            row_or.append(safe_or)
            row_annot.append(
                f"{or_v:.1f}{star}" if not math.isnan(or_v) else "n/a"
            )
        or_matrix.append(row_or)
        annot_matrix.append(row_annot)

    fig, ax = plt.subplots(figsize=(9, len(genes) * 0.9 + 1.5))
    sns.heatmap(
        or_matrix,
        annot=annot_matrix,
        fmt="",
        xticklabels=ref_labels,
        yticklabels=genes,
        cmap="RdBu",
        center=1.0,
        vmin=0,
        vmax=16,
        ax=ax,
        linewidths=0.5,
        cbar_kws={"label": "Odds Ratio"},
    )
    ax.set_title(
        f"Enrichment of {CT}-specific regulators in cancer driver gene set\n"
        "(* BH-FDR < 0.05, ** BH-FDR < 0.01; primary test = IntOGen)",
        fontsize=11,
    )
    ax.set_xlabel("Reference gene set", fontsize=10)
    ax.set_ylabel("Focal gene", fontsize=10)
    plt.tight_layout()
    plt.savefig(
        os.path.join(MANUSCRIPT_DIR, f"figure_heatmap_{ct}.png"), dpi=150
    )
    plt.close()


def plot_gremln_heatmap(
    gremln_comparison: dict, cancer_type: str, genes: list
) -> None:
    """Heatmap of GREmLN-only enrichment ORs (Figure 3C/3D)."""
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    ct_data = gremln_comparison.get(ct, {})
    gene_results = ct_data.get("gene_results", {})

    or_matrix, annot_matrix = [], []
    for gene in genes:
        gr = gene_results.get(gene, {})
        if gr.get("skipped") or gene not in gene_results:
            or_matrix.append([float("nan")])
            annot_matrix.append(["n/a"])
        else:
            or_v  = gr.get("gremln_only_or", float("nan"))
            adj_p = gr.get("gremln_only_fdr", gr.get("gremln_only_p", 1.0))
            star  = "**" if adj_p < 0.01 else ("*" if adj_p < 0.05 else "")
            safe_or = or_v if (or_v != float("inf") and not math.isnan(or_v)) else 5.0
            or_matrix.append([safe_or])
            annot_matrix.append([f"{or_v:.1f}{star}" if not math.isnan(or_v) else "n/a"])

    fig, ax = plt.subplots(figsize=(9, len(genes) * 0.9 + 1.5))
    sns.heatmap(
        or_matrix,
        annot=annot_matrix,
        fmt="",
        xticklabels=["IntOGen\n(pan-cancer)"],
        yticklabels=genes,
        cmap="RdBu",
        center=1.0,
        vmin=0,
        vmax=16,
        ax=ax,
        linewidths=0.5,
        cbar_kws={"label": "Odds Ratio"},
        mask=np.array([[math.isnan(row[0])] for row in or_matrix]),
    )
    ax.set_title(
        f"Enrichment of {CT} GREmLN-only candidates in cancer driver gene set\n"
        "(* BH-FDR < 0.05, ** BH-FDR < 0.01; GREmLN epithelial_cell background)",
        fontsize=11,
    )
    ax.set_xlabel("Reference gene set", fontsize=10)
    ax.set_ylabel("Focal gene", fontsize=10)
    plt.tight_layout()
    plt.savefig(
        os.path.join(MANUSCRIPT_DIR, f"figure_heatmap_gremln_{ct}.png"), dpi=150
    )
    plt.close()


def plot_regulator_counts(
    comparisons: dict,
    genes: list,
    results: dict,
    cancer_type: str,
    out_dir: str,
) -> None:
    """Bar chart: conserved vs. cancer-specific regulators, with IntOGen overlap."""
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    conserved_n = [comparisons[g]["regulators"]["conserved_count"] for g in genes]
    specific_n  = [len(comparisons[g]["regulators"]["tumor_state_only"]) for g in genes]
    intogen_n    = [
        results[g]["enrichment"].get("intogen", {}).get("ref_overlap", 0)
        for g in genes
    ]

    x     = np.arange(len(genes))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x - width / 2, conserved_n, width,
           label="Conserved regulators", color="#4e79a7")
    ax.bar(x + width / 2, specific_n, width,
           label=f"{CT}-specific regulators", color="#f28e2b")
    ax.bar(x + width / 2, intogen_n, width,
           label=f"{CT}-specific + IntOGen", color="#e15759", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=11)
    ax.set_ylabel("Number of regulators", fontsize=10)
    ax.set_title(
        f"Conserved vs. {CT}-specific regulators per focal gene\n"
        "(red overlay = overlap with IntOGen cancer drivers)",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(
        os.path.join(out_dir, f"experiment_rewiring_barchart_{ct}.png"), dpi=150
    )
    plt.close()


def plot_neg_controls(
    focal_results: dict,
    neg_results: dict,
    focal_genes: list,
    neg_genes: list,
    cancer_type: str,
) -> None:
    """Bar chart: IntOGen OR for focal cancer genes vs. housekeeping negative controls."""
    ct = cancer_type.upper()

    focal_ors  = []
    focal_lbls = []
    for g in focal_genes:
        if g in focal_results and not focal_results[g].get("skipped"):
            v = focal_results[g]["enrichment"].get("intogen", {}).get("odds_ratio", 0)
            focal_ors.append(min(float(v), 20.0) if not math.isnan(float(v)) else 0)
            focal_lbls.append(g)

    neg_ors  = []
    neg_lbls = []
    for g in neg_genes:
        if g in neg_results and not neg_results[g].get("skipped", True):
            v = neg_results[g]["enrichment"].get("intogen", {}).get("odds_ratio", 0)
            neg_ors.append(min(float(v), 20.0) if not math.isnan(float(v)) else 0)
            neg_lbls.append(g)

    all_ors  = focal_ors  + [None] + neg_ors
    all_lbls = focal_lbls + [""]   + neg_lbls
    colors   = (["#e15759"] * len(focal_ors)) + ["white"] + (["#76b7b2"] * len(neg_ors))

    x   = np.arange(len(all_ors))
    fig, ax = plt.subplots(figsize=(max(9, len(all_ors) * 0.9), 4))
    for i, (v, c) in enumerate(zip(all_ors, colors)):
        if v is not None:
            ax.bar(i, v, color=c)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(all_lbls, fontsize=10)
    ax.set_ylabel("Odds Ratio vs. IntOGen", fontsize=10)
    ax.set_title(
        f"Negative control validation ({ct}): cancer driver genes vs. housekeeping genes\n"
        "Red = cancer focal genes; teal = housekeeping negative controls (expected OR ~1)",
        fontsize=11,
    )
    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(color="#e15759", label="Cancer driver genes"),
                 Patch(color="#76b7b2", label="Housekeeping (negative control)")],
        fontsize=9,
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(MANUSCRIPT_DIR, f"figure_negcontrol_{cancer_type.lower()}.png"), dpi=150
    )
    plt.close()


def plot_neutral_controls(
    focal_results: dict,
    neutral_results: dict,
    focal_genes: list,
    neutral_genes: list,
    cancer_type: str,
) -> None:
    """Bar chart: IntOGen OR for focal cancer genes vs. tumor-expressed neutral controls."""
    ct = cancer_type.upper()

    focal_ors  = []
    focal_lbls = []
    for g in focal_genes:
        if g in focal_results and not focal_results[g].get("skipped"):
            v = focal_results[g]["enrichment"].get("intogen", {}).get("odds_ratio", 0)
            focal_ors.append(min(float(v), 20.0) if not math.isnan(float(v)) else 0)
            focal_lbls.append(g)

    neutral_ors  = []
    neutral_lbls = []
    for g in neutral_genes:
        if g in neutral_results and not neutral_results[g].get("skipped", True):
            v = neutral_results[g]["enrichment"].get("intogen", {}).get("odds_ratio", 0)
            neutral_ors.append(min(float(v), 20.0) if not math.isnan(float(v)) else 0)
            neutral_lbls.append(g)

    all_ors  = focal_ors  + [None] + neutral_ors
    all_lbls = focal_lbls + [""]   + neutral_lbls
    colors   = (["#e15759"] * len(focal_ors)) + ["white"] + (["#f28e2b"] * len(neutral_ors))

    x   = np.arange(len(all_ors))
    fig, ax = plt.subplots(figsize=(max(9, len(all_ors) * 0.9), 4))
    for i, (v, c) in enumerate(zip(all_ors, colors)):
        if v is not None:
            ax.bar(i, v, color=c)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(all_lbls, fontsize=10)
    ax.set_ylabel("Odds Ratio vs. IntOGen", fontsize=10)
    ax.set_title(
        f"Neutral control validation ({ct}): cancer driver genes vs. tumor-expressed non-driver genes\n"
        "Red = cancer focal genes; orange = neutral controls (tumor-expressed, non-IntOGen; expected OR ~1)",
        fontsize=11,
    )
    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(color="#e15759", label="Cancer driver genes"),
                 Patch(color="#f28e2b", label="Neutral (tumor-expressed, non-IntOGen)")],
        fontsize=9,
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(MANUSCRIPT_DIR, f"figure_neutralcontrol_{cancer_type.lower()}.png"), dpi=150
    )
    plt.close()


# ── Target list (source-labeled) ───────────────────────────────────────────────

def generate_target_list(comparison: dict, intogen: set, moa_map: dict, intogen_roles: dict) -> list:
    """
    Return all IntOGen-overlapping regulators from either network, labeled by source.

    Source values:
      "Both"         — present in both GREmLN and TCGA networks (highest confidence)
      "TCGA-only"    — tumor-selective; MoA available
      "GREmLN-only"  — present in normal epithelium only
    """
    conserved   = set(comparison["regulators"]["conserved"])
    tumor_only  = set(comparison["regulators"]["tumor_state_only"])
    normal_only = set(comparison["regulators"]["population_averaged_only"])

    rows = []
    for g in sorted(conserved & intogen):
        rows.append({"regulator": g, "source": "Both",
                     "moa": moa_map.get(g), "direction": _direction(moa_map.get(g)),
                     "intogen_role": intogen_roles.get(g, "")})
    for g in sorted(tumor_only & intogen):
        rows.append({"regulator": g, "source": "TCGA-only",
                     "moa": moa_map.get(g), "direction": _direction(moa_map.get(g)),
                     "intogen_role": intogen_roles.get(g, "")})
    for g in sorted(normal_only & intogen):
        rows.append({"regulator": g, "source": "GREmLN-only",
                     "moa": None, "direction": "",
                     "intogen_role": intogen_roles.get(g, "")})

    # sort: Both first, then TCGA-only (activating before repressive), then GREmLN-only
    order = {"Both": 0, "TCGA-only": 1, "GREmLN-only": 2}
    rows.sort(key=lambda r: (order[r["source"]], -(r["moa"] or 0)))
    return rows


def _direction(moa) -> str:
    if moa is None:
        return ""
    return "activating" if moa > 0 else ("repressive" if moa < 0 else "unknown")


def save_target_table(all_targets: dict, cancer_type: str, out_dir: str) -> None:
    """Write source-labeled target list to CSV.

    BRCA/COAD write directly to their named supplementary table
    (supplementary/table_sN_<cancer_type>_candidates.csv), matching the paper's
    Data Availability section. Other cancer types fall back to out_dir.
    """
    ct = cancer_type.upper()
    table_n = SUPPLEMENTARY_TABLE_NUMBER.get(cancer_type.lower())
    if table_n is not None:
        path = os.path.join(SUPPLEMENTARY_DIR, f"table_s{table_n}_{cancer_type.lower()}_candidates.csv")
    else:
        path = os.path.join(out_dir, f"target_list_{cancer_type.lower()}.csv")
    fields = ["focal_gene", "regulator", "source", "intogen_role", "moa", "direction"]
    rows = []
    for focal_gene, targets in all_targets.items():
        for entry in targets:
            rows.append({
                "focal_gene":  focal_gene,
                "regulator":   entry["regulator"],
                "source":      entry["source"],
                "intogen_role": entry.get("intogen_role", ""),
                "moa":         round(entry["moa"], 3) if entry["moa"] is not None else "",
                "direction":   entry["direction"],
            })
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[{ct}] Target list ({len(rows)} entries) -> {path}")



def plot_target_list(all_targets: dict, cancer_type: str) -> None:
    """Stacked bar: IntOGen targets per source per focal gene."""
    ct    = cancer_type.upper()
    genes = [g for g, t in all_targets.items() if t]
    if not genes:
        return

    def _count(g, src):
        return sum(1 for r in all_targets[g] if r["source"] == src)

    n_both  = [_count(g, "Both")        for g in genes]
    n_tcga  = [_count(g, "TCGA-only")   for g in genes]
    n_grem  = [_count(g, "GREmLN-only") for g in genes]
    x = range(len(genes))

    fig, ax = plt.subplots(figsize=(max(8, len(genes) * 0.9), 4))
    ax.bar(x, n_both, label="Both networks",   color="#2166ac")
    ax.bar(x, n_tcga, bottom=n_both,           label="TCGA-only",   color="#d6604d")
    ax.bar(x, n_grem, bottom=[a+b for a,b in zip(n_both, n_tcga)],
           label="GREmLN-only", color="#4dac26")
    ax.set_xticks(list(x))
    ax.set_xticklabels(genes, fontsize=11)
    ax.set_ylabel("IntOGen-overlapping regulators", fontsize=10)
    ax.set_title(
        f"{ct}: Candidate therapeutic regulators by network source\n"
        "(IntOGen-filtered; source indicates which network context identified each regulator)",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(MANUSCRIPT_DIR, f"target_list_{cancer_type.lower()}.png"), dpi=150)
    plt.close()


# ── Per-cancer analysis helper ─────────────────────────────────────────────────

def run_cancer_analysis(
    agent,
    workflow,
    cancer_type: str,
    focal_genes: list,
    intogen_raw: set,
    intogen_roles: dict,
    out_dir: str,
) -> tuple:
    """
    Run cross-context comparison + enrichment tests for one cancer type.

    Returns: (comparisons, results, background, combined_stouffer, testable_genes)
    """
    ct = cancer_type.lower()
    CT = cancer_type.upper()

    background = set(workflow.tcga_cache.tcga_indices[ct].get("all_genes", []))
    print(f"\n[{CT}] Background (TCGA gene universe): {len(background):,} genes")

    intogen = intogen_raw & background
    print(f"[{CT}] IntOGen in background: {len(intogen):,} genes")

    ref_sets = {"intogen": intogen}

    # Cross-context comparisons
    print(f"[{CT}] Running compare_network_contexts for {focal_genes} ...")
    comparisons: dict = {}
    for gene in focal_genes:
        print(f"  {gene} ...", end=" ", flush=True)
        result = compare_network_contexts(agent, gene, ct, CELL_TYPE)
        if result.get("error"):
            print(f"SKIP ({result.get('message', 'unknown error')})")
            continue
        comparisons[gene] = result
        r = result["regulators"]
        print(
            f"both={r['conserved_count']}  "
            f"gremln_only={len(r['population_averaged_only'])}  "
            f"{ct}_only={len(r['tumor_state_only'])}  "
            f"context_specificity={result['interpretation']['regulatory_rewiring']}"
        )

    # Enrichment tests
    print(f"[{CT}] Fisher's exact tests + permutation controls ...")
    results: dict = {}
    for gene, comp in comparisons.items():
        specific = set(comp["regulators"]["tumor_state_only"])
        gene_res: dict = {
            "comparison_summary": {
                "conserved_count":       comp["regulators"]["conserved_count"],
                f"{ct}_specific_count":  len(specific),
                "epithelial_only_count": len(comp["regulators"]["population_averaged_only"]),
                "rewiring":              comp["interpretation"]["regulatory_rewiring"],
                "conserved_fraction":    comp["regulators"]["conserved_fraction"],
            },
            "enrichment":    {},
            "permutation":   {},
            "moa_extension": {},
            "skipped": False,
        }

        # Populate MoA for tumor-specific regulators from TCGA network
        tcga_neighbors = agent.query_network(
            "gene_neighbors", gene=gene, network_source="tcga", tcga_network=ct
        )
        moa_map = {r["gene"]: r.get("moa") for r in tcga_neighbors.get("regulators", [])}
        gene_res["moa_extension"] = {g: moa_map[g] for g in specific if g in moa_map}
        gene_res["target_list"]   = generate_target_list(comp, intogen, moa_map, intogen_roles)

        if len(specific) < 3:
            print(f"  {gene}: skipping -- only {len(specific)} {ct}-specific regulators")
            gene_res["skipped"] = True
            results[gene] = gene_res
            continue

        for ref_name, ref_set in ref_sets.items():
            fisher = fisher_enrichment(specific, ref_set, background)
            perm   = permutation_test(
                specific, ref_set, background,
                n=N_PERMUTATIONS, seed=RANDOM_SEED,
            )
            gene_res["enrichment"][ref_name]  = fisher
            gene_res["permutation"][ref_name] = perm
            print(
                f"  {gene:7s} vs {ref_name:<32}: "
                f"OR={fisher['odds_ratio']:5.2f}  p={fisher['p_value']:.4f}  "
                f"emp_p={perm['empirical_p']:.4f}  "
                f"overlap={fisher['ref_overlap']}/{fisher['query_size']}"
            )

        results[gene] = gene_res

    # BH-FDR across primary tests
    testable = [g for g in focal_genes if g in results and not results[g]["skipped"]]
    primary_ps = [results[g]["enrichment"]["intogen"]["p_value"] for g in testable]
    fdr_vals   = bh_fdr(primary_ps)
    for gene, adj_p in zip(testable, fdr_vals):
        results[gene]["enrichment"]["intogen"]["fdr_adjusted_p"] = round(adj_p, 6)

    # Stouffer combined Z
    combined_stats: dict = {}
    for ref_name in ref_sets:
        ps = [results[g]["enrichment"][ref_name]["p_value"]
              for g in testable if ref_name in results[g].get("enrichment", {})]
        ws = [results[g]["enrichment"][ref_name]["query_size"]
              for g in testable if ref_name in results[g].get("enrichment", {})]
        combined_stats[ref_name] = stouffer_z(ps, ws)
    print(f"\n[{CT}] Stouffer combined Z (all testable focal genes):")
    for ref_name, s in combined_stats.items():
        print(f"  {ref_name:<32}: Z={s['combined_z']:6.3f}  p={s['combined_p']:.4f}")

    # Source-labeled target table
    all_targets = {g: results[g]["target_list"] for g in testable if "target_list" in results[g]}
    save_target_table(all_targets, ct, out_dir)

    # Figures
    plot_or_heatmap(results, testable, ct)
    plot_regulator_counts(comparisons, testable, results, ct, out_dir)
    plot_target_list(all_targets, ct)

    return comparisons, results, background, combined_stats, testable


# ── Negative controls ──────────────────────────────────────────────────────────

def run_negative_controls(
    agent,
    cancer_type: str,
    intogen_raw: set,
    background: set,
) -> dict:
    """
    Run the same enrichment test on housekeeping genes.
    Expected result: OR ≈ 1 (no enrichment) — validates that the enrichment
    seen for cancer driver genes is specific, not a general network property.
    """
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    intogen = intogen_raw & background
    ref_sets = {"intogen": intogen}

    print(f"\n[NEG CTRL / {CT}] Running housekeeping gene controls: {HOUSEKEEPING_GENES}")
    neg_results: dict = {}
    for gene in HOUSEKEEPING_GENES:
        print(f"  {gene} ...", end=" ", flush=True)
        result = compare_network_contexts(agent, gene, ct, CELL_TYPE)
        if result.get("error"):
            print(f"SKIP ({result.get('message', 'unknown error')})")
            continue
        specific = set(result["regulators"]["tumor_state_only"])
        r = result["regulators"]
        print(
            f"both={r['conserved_count']}  "
            f"{ct}_only={len(specific)}  "
            f"context_specificity={result['interpretation']['regulatory_rewiring']}"
        )
        if len(specific) < 3:
            neg_results[gene] = {"skipped": True, "specific_count": len(specific)}
            continue
        gene_res: dict = {"skipped": False, "specific_count": len(specific), "enrichment": {}}
        for ref_name, ref_set in ref_sets.items():
            fisher = fisher_enrichment(specific, ref_set, background)
            gene_res["enrichment"][ref_name] = fisher
            print(
                f"    vs {ref_name:<32}: "
                f"OR={fisher['odds_ratio']:5.2f}  p={fisher['p_value']:.4f}  "
                f"overlap={fisher['ref_overlap']}/{fisher['query_size']}"
            )
        neg_results[gene] = gene_res
    return neg_results


def run_neutral_controls(
    agent,
    cancer_type: str,
    intogen_raw: set,
    background: set,
) -> dict:
    """
    Run the same enrichment test on tumor-expressed, non-IntOGen, non-housekeeping
    genes (FASN, PCNA, PKM, PABPC1, VIM).

    These genes are present in the TCGA network (tumor-expressed) and have
    substantial network connectivity, but have no cancer-driver annotation in
    IntOGen. Expected result: OR ≈ 1 — validates that the enrichment seen for
    cancer driver genes requires cancer-specific biology, not merely tumor-network
    membership or high network degree.
    """
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    intogen = intogen_raw & background
    ref_sets = {"intogen": intogen}

    print(f"\n[NEUTRAL CTRL / {CT}] Running tumor-expressed neutral gene controls: {NEUTRAL_GENES}")
    neutral_results: dict = {}
    for gene in NEUTRAL_GENES:
        print(f"  {gene} ...", end=" ", flush=True)
        result = compare_network_contexts(agent, gene, ct, CELL_TYPE)
        if result.get("error"):
            print(f"SKIP ({result.get('message', 'unknown error')})")
            continue
        specific = set(result["regulators"]["tumor_state_only"])
        r = result["regulators"]
        print(
            f"both={r['conserved_count']}  "
            f"{ct}_only={len(specific)}  "
            f"context_specificity={result['interpretation']['regulatory_rewiring']}"
        )
        if len(specific) < 3:
            neutral_results[gene] = {"skipped": True, "specific_count": len(specific)}
            continue
        gene_res: dict = {"skipped": False, "specific_count": len(specific), "enrichment": {}}
        for ref_name, ref_set in ref_sets.items():
            fisher = fisher_enrichment(specific, ref_set, background)
            gene_res["enrichment"][ref_name] = fisher
            print(
                f"    vs {ref_name:<32}: "
                f"OR={fisher['odds_ratio']:5.2f}  p={fisher['p_value']:.4f}  "
                f"overlap={fisher['ref_overlap']}/{fisher['query_size']}"
            )
        neutral_results[gene] = gene_res
    return neutral_results


# ── GREmLN formal enrichment analysis ─────────────────────────────────────────

def run_gremln_comparison(
    agent,
    brca_comparisons: dict,
    coad_comparisons: dict,
    intogen_raw: set,
    brca_bg: set,
    coad_bg: set,
) -> dict:
    """
    Formal enrichment analysis of GREmLN-only candidates vs IntOGen, mirroring
    the TCGA analysis in run_cancer_analysis().

    Background: GREmLN epithelial_cell gene universe (ENSG->symbol via pre-built
    gene ID cache; run build_network_cache.py --enrich-gene-cache for full coverage).
    Reference: IntOGen cancer genes intersected with GREmLN background.

    Runs per-gene Fisher's exact test + permutation control (n=1,000), BH-FDR
    correction across testable genes, and Stouffer Z per cancer type.

    Note: the GREmLN epithelial_cell network is not cancer-type-specific (it is a
    pan-tissue healthy epithelial network), so the background is the same for BRCA
    and COAD. Results should be interpreted as measuring enrichment of normal
    epithelial regulatory candidates in IntOGen, not tumor-specific enrichment.
    """
    print("\n" + "=" * 70)
    print("GREmLN-only IntOGen enrichment analysis (formal statistics)")
    print("Background: GREmLN epithelial_cell gene universe")
    print("=" * 70)

    import pickle as _pickle
    cache_path = "cache/gene_id_cache.pkl"
    try:
        with open(cache_path, "rb") as f:
            id_cache = _pickle.load(f)
        e2s = id_cache.get("ensembl_to_symbol", {})
    except Exception as exc:
        print(f"WARNING: Could not load gene ID cache ({exc}). Skipping GREmLN analysis.")
        return {}

    gremln_idx = agent.cache.network_indices.get(CELL_TYPE, {})
    ensg_ids   = gremln_idx.get("all_genes", set())
    gremln_bg  = {e2s[e].upper() for e in ensg_ids if e in e2s}
    intogen_g   = intogen_raw & gremln_bg
    coverage   = 100 * len(gremln_bg) / max(len(ensg_ids), 1)

    print(f"\nGREmLN background: {len(gremln_bg):,} symbols from {len(ensg_ids):,} "
          f"ENSG IDs ({coverage:.0f}% coverage)")
    print(f"IntOGen in GREmLN background: {len(intogen_g):,}")

    # Gene universe overlap between GREmLN and each TCGA network
    print("\nGene universe overlap (GREmLN epithelial_cell vs TCGA):")
    universe_overlap: dict = {}
    for cancer, tcga_bg in [("brca", brca_bg), ("coad", coad_bg)]:
        overlap = gremln_bg & tcga_bg
        pct_tcga = 100 * len(overlap) / max(len(tcga_bg), 1)
        pct_gremln = 100 * len(overlap) / max(len(gremln_bg), 1)
        universe_overlap[cancer] = {
            "gremln_size":    len(gremln_bg),
            "tcga_size":      len(tcga_bg),
            "overlap":        len(overlap),
            "pct_tcga_in_gremln":  round(pct_tcga, 1),
            "pct_gremln_in_tcga":  round(pct_gremln, 1),
        }
        print(f"  {cancer.upper()}: overlap={len(overlap):,}  "
              f"{pct_tcga:.1f}% of TCGA in GREmLN  "
              f"{pct_gremln:.1f}% of GREmLN in TCGA")

    out: dict = {}
    for cancer, comparisons, tcga_bg, focal_genes in [
        ("brca", brca_comparisons, brca_bg, BRCA_GENES),
        ("coad", coad_comparisons, coad_bg, COAD_GENES),
    ]:
        CT = cancer.upper()
        intogen_t = intogen_raw & tcga_bg
        out[cancer] = {"universe_overlap": universe_overlap[cancer]}

        print(f"\n[GREmLN / {CT}] Fisher's exact tests + permutation controls ...")
        gene_results: dict = {}

        for gene in focal_genes:
            if gene not in comparisons:
                continue
            comp        = comparisons[gene]
            tcga_only   = set(comp["regulators"]["tumor_state_only"])
            gremln_only = set(comp["regulators"]["population_averaged_only"])

            # TCGA-only enrichment (against TCGA background, for reference)
            if len(tcga_only) >= 3:
                t = fisher_enrichment(tcga_only, intogen_t, tcga_bg)
                t_or, t_n, t_ov = t["odds_ratio"], t["query_size"], t["ref_overlap"]
            else:
                t_or, t_n, t_ov = 0.0, len(tcga_only), 0

            # GREmLN-only enrichment (against GREmLN background)
            gene_res: dict = {
                "tcga_only_or":               t_or,
                "tcga_only_n":                t_n,
                "tcga_only_intogen_overlap":   t_ov,
                "gremln_only_n":              len(gremln_only),
                "gremln_only_intogen_overlap": 0,
                "gremln_only_or":             0.0,
                "gremln_only_p":              1.0,
                "gremln_only_emp_p":          1.0,
                "gremln_only_fdr":            1.0,
                "skipped": False,
            }

            if len(gremln_only) < 3:
                print(f"  {gene}: skipping -- only {len(gremln_only)} GREmLN-only regulators")
                gene_res["skipped"] = True
                gene_results[gene] = gene_res
                continue

            g_fisher = fisher_enrichment(gremln_only, intogen_g, gremln_bg)
            g_perm   = permutation_test(
                gremln_only, intogen_g, gremln_bg,
                n=N_PERMUTATIONS, seed=RANDOM_SEED,
            )
            gene_res.update({
                "gremln_only_or":             g_fisher["odds_ratio"],
                "gremln_only_p":              g_fisher["p_value"],
                "gremln_only_emp_p":          g_perm["empirical_p"],
                "gremln_only_intogen_overlap": g_fisher["ref_overlap"],
            })
            gene_results[gene] = gene_res

            print(
                f"  {gene:7s} vs intogen (GREmLN)              : "
                f"OR={g_fisher['odds_ratio']:5.2f}  p={g_fisher['p_value']:.4f}  "
                f"emp_p={g_perm['empirical_p']:.4f}  "
                f"overlap={g_fisher['ref_overlap']}/{g_fisher['query_size']}"
            )

        # BH-FDR across testable GREmLN genes
        testable = [g for g in focal_genes
                    if g in gene_results and not gene_results[g]["skipped"]]
        if testable:
            gremln_ps  = [gene_results[g]["gremln_only_p"] for g in testable]
            fdr_vals   = bh_fdr(gremln_ps)
            for gene, adj_p in zip(testable, fdr_vals):
                gene_results[gene]["gremln_only_fdr"] = round(adj_p, 6)

            # Stouffer Z for GREmLN-only panel
            ws = [gene_results[g]["gremln_only_n"] for g in testable]
            gremln_stouffer = stouffer_z(gremln_ps, ws)
            print(f"\n[GREmLN / {CT}] Stouffer combined Z (GREmLN-only, all testable):")
            print(f"  intogen (GREmLN background)      : "
                  f"Z={gremln_stouffer['combined_z']:6.3f}  "
                  f"p={gremln_stouffer['combined_p']:.4f}")
        else:
            gremln_stouffer = {"combined_z": float("nan"), "combined_p": float("nan")}

        out[cancer] = {
            "gene_results":      gene_results,
            "combined_stouffer": gremln_stouffer,
            "background_size":   len(gremln_bg),
            "intogen_in_bg":      len(intogen_g),
            "coverage_pct":      round(coverage, 1),
        }

    print(f"\nGREmLN background note: {len(gremln_bg):,}/{len(ensg_ids):,} ENSG IDs "
          f"resolved via pre-built cache ({coverage:.0f}% coverage). "
          f"epithelial_cell network is pan-tissue (not cancer-type-specific).")
    out["universe_overlap"] = universe_overlap
    return out


# ── Tier specificity: focal genes vs random genes (Table S3) ──────────────────

def _tier_sets(symbol: str, gremln_idx: dict, e2s: dict, sym2id: dict,
               tcga_idx: dict, tcga_bg: set):
    """(tcga_only, gremln_only) regulator sets read straight from the network
    indices, or None if the gene is not in both networks. Gives the same sets as
    compare_network_contexts (checked against it for every focal gene)."""
    ensg = sym2id.get(symbol)
    if ensg is None or symbol not in tcga_bg:
        return None
    g = {(e2s.get(r) or r).upper() for r in gremln_idx["target_regulators"].get(ensg, [])}
    t = {r.upper() for r in tcga_idx["target_regulators"].get(symbol, [])}
    return t - g, g - t


def _haldane_log_or(query: set, reference: set, background: set) -> float:
    """Log odds ratio with +0.5 in every cell (finite when a cell is zero)."""
    bg = background | query
    a = len(query & reference)
    b = len(query) - a
    c = len((bg - query) & reference)
    d = len(bg) - len(query) - c
    return float(np.log(((a + .5) * (d + .5)) / ((b + .5) * (c + .5))))


def _tier_test(query: set, reference: set, background: set) -> dict:
    if len(query) < 3:
        return {"skipped": True, "n": len(query)}
    f = fisher_enrichment(query, reference, background)
    return {"skipped": False, "n": len(query), "overlap": f["ref_overlap"],
            "or": f["odds_ratio"], "p": f["p_value"],
            "log_or": _haldane_log_or(query, reference, background)}


def _tier_rows_testable(rows: dict) -> list:
    return [r for r in rows.values() if not r.get("skipped", True)]


def _tier_add_fdr(rows: dict) -> None:
    names = [g for g, r in rows.items() if not r.get("skipped", True)]
    for g, q in zip(names, bh_fdr([rows[g]["p"] for g in names])):
        rows[g]["fdr"] = q


def _tier_stouffer(rows: list) -> float:
    ps = [max(r["p"], TIER_P_FLOOR) for r in rows]
    return stouffer_z(ps, [r["n"] for r in rows])["combined_z"]


def _tier_summary(rows: dict) -> dict:
    t = _tier_rows_testable(rows)
    if not t:
        return {"testable": 0}
    return {
        "testable":       len(t),
        "median_or":      float(np.median([r["or"] for r in t])),
        "frac_p_lt_0.05": sum(r["p"] < 0.05 for r in t) / len(t),
        "n_fdr_lt_0.05":  sum(r.get("fdr", 1.0) < 0.05 for r in t),
        "stouffer_z":     _tier_stouffer(t),
    }


def _focal_vs_random(focal_rows: dict, random_rows: dict, seed: int) -> dict:
    """
    Compare the focal panel with random genes two ways:
      - Stouffer Z (the paper's statistic) vs TIER_NULL_RESAMPLES random panels of
        the same size; grows with candidate-set size, so read it with care.
      - Mean log OR, focal minus random, permutation test; size-independent.
    """
    f, n = _tier_rows_testable(focal_rows), _tier_rows_testable(random_rows)
    focal_z = _tier_stouffer(f)
    rng = np.random.default_rng(seed)
    null_zs = np.array([
        _tier_stouffer([n[i] for i in rng.choice(len(n), len(f), replace=False)])
        for _ in range(TIER_NULL_RESAMPLES)
    ])
    f_lor = np.array([r["log_or"] for r in f])
    n_lor = np.array([r["log_or"] for r in n])
    observed = f_lor.mean() - n_lor.mean()
    pooled, k = np.concatenate([f_lor, n_lor]), len(f_lor)
    perm = np.empty(TIER_LOGOR_PERMUTATIONS)
    for i in range(TIER_LOGOR_PERMUTATIONS):
        x = rng.permutation(pooled)
        perm[i] = x[:k].mean() - x[k:].mean()
    return {
        "seed":                 seed,
        "n_focal":              len(f),
        "n_random":             len(n),
        "focal_z":              focal_z,
        "random_z_median":      float(np.median(null_zs)),
        "random_z_95th":        float(np.percentile(null_zs, 95)),
        "stouffer_empirical_p": float((np.sum(null_zs >= focal_z) + 1) / (len(null_zs) + 1)),
        "focal_geomean_or":     float(np.exp(f_lor.mean())),
        "random_geomean_or":    float(np.exp(n_lor.mean())),
        "logor_permutation_p":  float((np.sum(perm >= observed) + 1) / (len(perm) + 1)),
    }


def run_tier_specificity(agent, workflow, intogen_raw: set,
                         brca_comparisons: dict, coad_comparisons: dict) -> dict:
    """
    Do focal cancer genes' candidates beat random genes' candidates?

    Regulators are themselves enriched in cancer genes relative to all genes, so
    any gene's candidates look enriched against an all-gene background. This runs
    the TCGA-only and GREmLN-only tests on TIER_NULL_N_RANDOM random genes (in both
    networks, not IntOGen, not focal/control) per seed, and compares the focal panel
    with them — under the paper's all-gene background and a regulator-only one.
    Control genes are also run through both tiers for reference.
    """
    print("\n" + "=" * 70)
    print("Tier specificity: focal genes vs random genes (Table S3)")
    print("=" * 70)

    import pickle as _pickle
    with open("cache/gene_id_cache.pkl", "rb") as f:
        e2s = _pickle.load(f).get("ensembl_to_symbol", {})
    gremln_idx = agent.cache.network_indices[CELL_TYPE]
    in_net = set(gremln_idx["all_genes"])
    sym2id: dict = {}  # symbol -> in-network Ensembl ID
    for ensg, sym in e2s.items():
        if ensg in in_net:
            sym2id.setdefault(sym.upper(), ensg)
    g_bg   = {e2s[e].upper() for e in gremln_idx["all_genes"] if e in e2s}
    g_regs = {e2s[e].upper() for e in gremln_idx["regulator_targets"] if e in e2s}

    def regulator_bias(regs: set, bg: set) -> dict:
        a, b = len(intogen_raw & regs) / len(regs), len(intogen_raw & bg) / len(bg)
        return {"n_regulators": len(regs), "n_background": len(bg),
                "intogen_frac_regulators": a, "intogen_frac_background": b, "ratio": a / b}

    out: dict = {
        "config": {"n_random": TIER_NULL_N_RANDOM, "seeds": TIER_NULL_SEEDS,
                   "n_null_resamples": TIER_NULL_RESAMPLES,
                   "n_logor_permutations": TIER_LOGOR_PERMUTATIONS},
        "regulator_bias": {"gremln": regulator_bias(g_regs, g_bg)},
    }
    focal = {"brca": BRCA_GENES, "coad": COAD_GENES}
    pipeline = {"brca": brca_comparisons, "coad": coad_comparisons}
    excluded = set(BRCA_GENES) | set(COAD_GENES) | set(HOUSEKEEPING_GENES) | set(NEUTRAL_GENES)

    for ct in ["brca", "coad"]:
        CT = ct.upper()
        tcga_idx = workflow.tcga_cache.tcga_indices[ct]
        t_bg   = {g.upper() for g in tcga_idx["all_genes"]}
        t_regs = {g.upper() for g in tcga_idx["regulator_targets"]}
        out["regulator_bias"][ct] = regulator_bias(t_regs, t_bg)

        backgrounds = {
            "all_gene_bg":  {"tcga": (intogen_raw & t_bg, t_bg),
                             "gremln": (intogen_raw & g_bg, g_bg)},
            "regulator_bg": {"tcga": (intogen_raw & t_regs, t_regs),
                             "gremln": (intogen_raw & g_regs, g_regs)},
        }
        groups = {"focal": focal[ct], "housekeeping": HOUSEKEEPING_GENES,
                  "neutral": NEUTRAL_GENES}

        sets = {g: _tier_sets(g, gremln_idx, e2s, sym2id, tcga_idx, t_bg)
                for genes in groups.values() for g in genes}
        for g, comp in pipeline[ct].items():  # must match compare_network_contexts
            r = comp["regulators"]
            via_pipeline = ({x.upper() for x in r["tumor_state_only"]},
                            {x.upper() for x in r["population_averaged_only"]})
            if sets.get(g) != via_pipeline:
                raise RuntimeError(f"[{CT}] {g}: tier sets differ from compare_network_contexts")

        pool = sorted((g_bg & t_bg) - intogen_raw - excluded)
        random_sets = {}
        for seed in TIER_NULL_SEEDS:
            sample = random.Random(seed).sample(pool, min(TIER_NULL_N_RANDOM, len(pool)))
            random_sets[seed] = {
                g: s for g in sample
                if (s := _tier_sets(g, gremln_idx, e2s, sym2id, tcga_idx, t_bg)) is not None}

        out[ct] = {}
        for bg_name, tiers in backgrounds.items():
            out[ct][bg_name] = {}
            for tier, (ref, bg) in tiers.items():
                ix = 0 if tier == "tcga" else 1
                res: dict = {}
                for grp, genes in groups.items():
                    rows = {g: ({"skipped": True, "error": "not in both networks"}
                                if sets[g] is None else _tier_test(sets[g][ix], ref, bg))
                            for g in genes}
                    _tier_add_fdr(rows)
                    res[grp] = {"genes": rows, "summary": _tier_summary(rows)}
                per_seed = []
                for seed, rs in random_sets.items():
                    random_rows = {g: _tier_test(s[ix], ref, bg) for g, s in rs.items()}
                    _tier_add_fdr(random_rows)
                    if seed == RANDOM_SEED:
                        res["random_genes"] = {"seed": seed, "summary": _tier_summary(random_rows)}
                    per_seed.append(_focal_vs_random(res["focal"]["genes"], random_rows, seed))
                res["focal_vs_random_by_seed"] = per_seed
                out[ct][bg_name][tier] = res
                lo = min(r["logor_permutation_p"] for r in per_seed)
                hi = max(r["logor_permutation_p"] for r in per_seed)
                print(f"[{CT} {bg_name} {tier}-only] focal Z={per_seed[0]['focal_z']:.2f}  "
                      f"focal vs random log-OR p={lo:.3f}-{hi:.3f} ({len(per_seed)} seeds)")
    return out


def save_tier_specificity_table(tier_results: dict) -> None:
    """One row per cancer type x tier x background x seed -> supplementary Table S3."""
    path = os.path.join(SUPPLEMENTARY_DIR, TIER_SPECIFICITY_TABLE)
    fields = ["cancer_type", "tier", "background", "seed", "n_focal", "n_random",
              "focal_geomean_or", "random_geomean_or", "logor_permutation_p",
              "focal_stouffer_z", "random_stouffer_z_median", "stouffer_empirical_p"]
    tier_label = {"tcga": "TCGA-only", "gremln": "GREmLN-only"}
    bg_label = {"all_gene_bg": "all genes", "regulator_bg": "regulators only"}
    rows = []
    for ct in ["brca", "coad"]:
        for bg_name in ["all_gene_bg", "regulator_bg"]:
            for tier in ["tcga", "gremln"]:
                for r in tier_results[ct][bg_name][tier]["focal_vs_random_by_seed"]:
                    rows.append({
                        "cancer_type":              ct.upper(),
                        "tier":                     tier_label[tier],
                        "background":               bg_label[bg_name],
                        "seed":                     r["seed"],
                        "n_focal":                  r["n_focal"],
                        "n_random":                 r["n_random"],
                        "focal_geomean_or":         round(r["focal_geomean_or"], 3),
                        "random_geomean_or":        round(r["random_geomean_or"], 3),
                        "logor_permutation_p":      round(r["logor_permutation_p"], 4),
                        "focal_stouffer_z":         round(r["focal_z"], 3),
                        "random_stouffer_z_median": round(r["random_z_median"], 3),
                        "stouffer_empirical_p":     round(r["stouffer_empirical_p"], 4),
                    })
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Tier specificity table ({len(rows)} rows) -> {path}")


# ── Main ───────────────────────────────────────────────────────────────────────

def run_experiment() -> None:
    for d in (RESULTS_DIR, MANUSCRIPT_DIR, SUPPLEMENTARY_DIR):
        os.makedirs(d, exist_ok=True)
    random.seed(RANDOM_SEED)

    print("Loading RegNetAgents workflow...")
    workflow = RegNetAgentsWorkflow()
    agent    = workflow.modeling_agent

    for ct in ["brca", "coad"]:
        if not workflow.tcga_cache.tcga_indices.get(ct):
            print(f"ERROR: No TCGA cache for '{ct}'. Run build_tcga_cache.py first.")
            sys.exit(1)

    intogen_roles = load_intogen_roles()
    intogen_raw   = set(intogen_roles)
    print(f"IntOGen drivers: {len(intogen_raw):,} total")

    # ── BRCA analysis ──────────────────────────────────────────────────────────
    (brca_comp, brca_res, brca_bg,
     brca_stouffer, brca_testable) = run_cancer_analysis(
        agent, workflow, "brca", BRCA_GENES, intogen_raw, intogen_roles, RESULTS_DIR,
    )
    brca_neg     = run_negative_controls(agent, "brca", intogen_raw, brca_bg)
    brca_neutral = run_neutral_controls(agent, "brca", intogen_raw, brca_bg)
    plot_neg_controls(brca_res, brca_neg, brca_testable, HOUSEKEEPING_GENES, "brca")
    plot_neutral_controls(brca_res, brca_neutral, brca_testable, NEUTRAL_GENES, "brca")

    # ── COAD analysis ──────────────────────────────────────────────────────────
    (coad_comp, coad_res, coad_bg,
     coad_stouffer, coad_testable) = run_cancer_analysis(
        agent, workflow, "coad", COAD_GENES, intogen_raw, intogen_roles, RESULTS_DIR,
    )
    coad_neg     = run_negative_controls(agent, "coad", intogen_raw, coad_bg)
    coad_neutral = run_neutral_controls(agent, "coad", intogen_raw, coad_bg)
    plot_neg_controls(coad_res, coad_neg, coad_testable, HOUSEKEEPING_GENES, "coad")
    plot_neutral_controls(coad_res, coad_neutral, coad_testable, NEUTRAL_GENES, "coad")

    # ── Exploratory: GREmLN-only vs TCGA-only enrichment comparison ────────────
    gremln_comparison = run_gremln_comparison(
        agent,
        brca_comparisons=brca_comp,
        coad_comparisons=coad_comp,
        intogen_raw=intogen_raw,
        brca_bg=brca_bg,
        coad_bg=coad_bg,
    )

    # ── Tier specificity: focal genes vs random genes (Table S3) ──────────────
    tier_specificity = run_tier_specificity(agent, workflow, intogen_raw, brca_comp, coad_comp)
    save_tier_specificity_table(tier_specificity)

    # ── Save combined JSON ─────────────────────────────────────────────────────
    output = {
        "config": {
            "brca_focal_genes": BRCA_GENES,
            "coad_focal_genes": COAD_GENES,
            "cell_type":        CELL_TYPE,
            "n_permutations":   N_PERMUTATIONS,
        },
        "brca": {
            "background_size":   len(brca_bg),
            "gene_results":      brca_res,
            "combined_stouffer": brca_stouffer,
        },
        "coad": {
            "background_size":   len(coad_bg),
            "gene_results":      coad_res,
            "combined_stouffer": coad_stouffer,
        },
        "negative_controls": {
            "brca": brca_neg,
            "coad": coad_neg,
        },
        "neutral_controls": {
            "brca": brca_neutral,
            "coad": coad_neutral,
        },
        "gremln_comparison": gremln_comparison,
        "tier_specificity":  tier_specificity,
    }
    plot_gremln_heatmap(gremln_comparison, "brca", BRCA_GENES)
    plot_gremln_heatmap(gremln_comparison, "coad", COAD_GENES)
    plot_workflow_figure()

    out_json = os.path.join(RESULTS_DIR, "experiment_rewiring_results.json")
    with open(out_json, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults -> {out_json}")
    print(f"Figures  -> {MANUSCRIPT_DIR}/figure_workflow.png  (NAR Fig 1)")
    print(f"            {MANUSCRIPT_DIR}/figure_heatmap_brca.png  (NAR Fig 3A)")
    print(f"            {MANUSCRIPT_DIR}/figure_heatmap_coad.png  (NAR Fig 3B)")
    print(f"            {MANUSCRIPT_DIR}/figure_heatmap_gremln_brca.png  (NAR Fig 3C)")
    print(f"            {MANUSCRIPT_DIR}/figure_heatmap_gremln_coad.png  (NAR Fig 3D)")
    print(f"            {MANUSCRIPT_DIR}/figure_negcontrol_brca.png  (NAR Fig 4A)")
    print(f"            {MANUSCRIPT_DIR}/figure_negcontrol_coad.png  (NAR Fig 4B)")
    print(f"            {MANUSCRIPT_DIR}/figure_neutralcontrol_brca.png  (neutral ctrl BRCA)")
    print(f"            {MANUSCRIPT_DIR}/figure_neutralcontrol_coad.png  (neutral ctrl COAD)")
    print(f"            {MANUSCRIPT_DIR}/target_list_brca.png  (NAR Fig 2A)")
    print(f"            {MANUSCRIPT_DIR}/target_list_coad.png  (NAR Fig 2B)")
    print(f"            {RESULTS_DIR}/experiment_rewiring_barchart_brca.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_barchart_coad.png")
    print(f"            {SUPPLEMENTARY_DIR}/{TIER_SPECIFICITY_TABLE}  (Table S3)")
    print("\nDone.")


if __name__ == "__main__":
    run_experiment()
