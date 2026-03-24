#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
experiment_regulatory_candidates.py
====================================
Cross-context regulatory candidate identification for the RegNetAgents research
paper (Briefings in Bioinformatics).

For each focal gene in BRCA and COAD, queries both the GREmLN epithelial
network and the TCGA tumor network, then generates a source-labeled candidate
therapeutic target list filtered by OncoKB cancer driver annotation.

Reference datasets (pre-downloaded to data/):
  - data/oncokb_cancer_genes.tsv   : OncoKB pan-cancer driver list

Outputs:
  results/
  - experiment_rewiring_results.json           : full statistics
  - target_list_brca.csv                       : source-labeled candidate list, BRCA
  - target_list_coad.csv                       : source-labeled candidate list, COAD
  - target_list_brca.png                       : candidate counts by source, BRCA (NAR Fig 2A)
  - target_list_coad.png                       : candidate counts by source, COAD (NAR Fig 2B)
  - experiment_rewiring_barchart_brca.png      : regulator count bar chart, BRCA
  - experiment_rewiring_barchart_coad.png      : regulator count bar chart, COAD

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
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow

# ── Configuration ──────────────────────────────────────────────────────────────

BRCA_GENES  = [
    "TP53", "MYC", "CTNNB1", "CCND1",          # original panel
    "BRCA1", "BRCA2", "PIK3CA", "PTEN",         # BRCA hallmarks
    "RB1", "ERBB2", "ESR1", "GATA3",            # additional BRCA drivers
]  # CDH1 not in GREmLN epithelial network
COAD_GENES  = [
    "TP53", "MYC", "CTNNB1", "CCND1",           # original panel
    "KRAS", "APC", "SMAD4",                      # original COAD panel
    "BRAF", "PIK3CA", "PTEN",                    # additional CRC drivers
    "FBXW7", "TCF7L2", "RNF43",                 # additional CRC drivers
]
HOUSEKEEPING_GENES = ["ACTB", "GAPDH", "HPRT1", "LDHA", "TUBB"]
CELL_TYPE   = "epithelial_cell"
N_PERMUTATIONS = 1000
RANDOM_SEED    = 42

ONCOKB_PATH    = "data/oncokb_cancer_genes.tsv"
RESULTS_DIR    = "results"
MANUSCRIPT_DIR = "manuscript"

# ── Reference set loaders ──────────────────────────────────────────────────────

def download_oncokb_if_missing(path: str) -> None:
    """Download OncoKB cancer gene list from public API if not present locally."""
    if os.path.exists(path):
        return
    import urllib.request
    url = "https://www.oncokb.org/api/v1/utils/cancerGeneList"
    print(f"OncoKB file not found. Downloading from {url} ...")
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["Hugo Symbol", "Gene Type"])
            for entry in data:
                writer.writerow([
                    entry.get("hugoSymbol", ""),
                    entry.get("geneType", "OTHER"),
                ])
        print(f"OncoKB cancer gene list saved to {path} ({len(data)} genes)")
    except Exception as exc:
        print(f"ERROR: Could not download OncoKB data: {exc}")
        print("Please download manually from https://www.oncokb.org/cancer-genes")
        sys.exit(1)


def load_oncokb(path: str) -> set:
    """Return symbol set for ONCOGENE / TSG / ONCOGENE_AND_TSG genes."""
    genes = set()
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["Gene Type"] in {"ONCOGENE", "TSG", "ONCOGENE_AND_TSG"}:
                genes.add(row["Hugo Symbol"].strip().upper())
    return genes


def load_oncokb_roles(path: str) -> dict:
    """Return {symbol: role} for ONCOGENE / TSG / ONCOGENE_AND_TSG genes."""
    roles = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gt = row["Gene Type"]
            if gt in {"ONCOGENE", "TSG", "ONCOGENE_AND_TSG"}:
                roles[row["Hugo Symbol"].strip().upper()] = gt
    return roles


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
        "p_value":     round(p_val, 6),
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
    zs = np.array([stats.norm.ppf(1 - p) if p < 1.0 else 0.0 for p in p_values])
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
        ("Filter all candidates against OncoKB",
         "Source (TCGA-only / GREmLN-only / Both)  ·  OncoKB role",
         "MoA direction (activating / repressive)  for TCGA-only"))

    arrow(7.0, 6.2, "TCGA-only candidates")

    box(4.2, 2.1, C_ENRICH, "ENRICHMENT VALIDATION",
        ("Fisher's exact test  ·  OncoKB reference",
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
    ref_keys   = ["oncokb"]
    ref_labels = ["OncoKB\n(pan-cancer)"]

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
        "(* BH-FDR < 0.05, ** BH-FDR < 0.01; primary test = OncoKB)",
        fontsize=11,
    )
    ax.set_xlabel("Reference gene set", fontsize=10)
    ax.set_ylabel("Focal gene", fontsize=10)
    plt.tight_layout()
    plt.savefig(
        os.path.join(MANUSCRIPT_DIR, f"figure_heatmap_{ct}.png"), dpi=150
    )
    plt.close()


def plot_regulator_counts(
    comparisons: dict,
    genes: list,
    results: dict,
    cancer_type: str,
    out_dir: str,
) -> None:
    """Bar chart: conserved vs. cancer-specific regulators, with OncoKB overlap."""
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    conserved_n = [comparisons[g]["regulators"]["conserved_count"] for g in genes]
    specific_n  = [len(comparisons[g]["regulators"]["tumor_state_only"]) for g in genes]
    oncokb_n    = [
        results[g]["enrichment"].get("oncokb", {}).get("ref_overlap", 0)
        for g in genes
    ]

    x     = np.arange(len(genes))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x - width / 2, conserved_n, width,
           label="Conserved regulators", color="#4e79a7")
    ax.bar(x + width / 2, specific_n, width,
           label=f"{CT}-specific regulators", color="#f28e2b")
    ax.bar(x + width / 2, oncokb_n, width,
           label=f"{CT}-specific + OncoKB", color="#e15759", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=11)
    ax.set_ylabel("Number of regulators", fontsize=10)
    ax.set_title(
        f"Conserved vs. {CT}-specific regulators per focal gene\n"
        "(red overlay = overlap with OncoKB cancer drivers)",
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
    """Bar chart: OncoKB OR for focal cancer genes vs. housekeeping negative controls."""
    ct = cancer_type.upper()

    focal_ors  = []
    focal_lbls = []
    for g in focal_genes:
        if g in focal_results and not focal_results[g].get("skipped"):
            v = focal_results[g]["enrichment"].get("oncokb", {}).get("odds_ratio", 0)
            focal_ors.append(min(float(v), 20.0) if not math.isnan(float(v)) else 0)
            focal_lbls.append(g)

    neg_ors  = []
    neg_lbls = []
    for g in neg_genes:
        if g in neg_results and not neg_results[g].get("skipped", True):
            v = neg_results[g]["enrichment"].get("oncokb", {}).get("odds_ratio", 0)
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
    ax.set_ylabel("Odds Ratio vs. OncoKB", fontsize=10)
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


# ── Target list (source-labeled) ───────────────────────────────────────────────

def generate_target_list(comparison: dict, oncokb: set, moa_map: dict, oncokb_roles: dict) -> list:
    """
    Return all OncoKB-overlapping regulators from either network, labeled by source.

    Source values:
      "Both"         — present in both GREmLN and TCGA networks (highest confidence)
      "TCGA-only"    — tumor-selective; MoA available
      "GREmLN-only"  — present in normal epithelium only
    """
    conserved   = set(comparison["regulators"]["conserved"])
    tumor_only  = set(comparison["regulators"]["tumor_state_only"])
    normal_only = set(comparison["regulators"]["population_averaged_only"])

    rows = []
    for g in sorted(conserved & oncokb):
        rows.append({"regulator": g, "source": "Both",
                     "moa": moa_map.get(g), "direction": _direction(moa_map.get(g)),
                     "oncokb_role": oncokb_roles.get(g, "")})
    for g in sorted(tumor_only & oncokb):
        rows.append({"regulator": g, "source": "TCGA-only",
                     "moa": moa_map.get(g), "direction": _direction(moa_map.get(g)),
                     "oncokb_role": oncokb_roles.get(g, "")})
    for g in sorted(normal_only & oncokb):
        rows.append({"regulator": g, "source": "GREmLN-only",
                     "moa": None, "direction": "",
                     "oncokb_role": oncokb_roles.get(g, "")})

    # sort: Both first, then TCGA-only (activating before repressive), then GREmLN-only
    order = {"Both": 0, "TCGA-only": 1, "GREmLN-only": 2}
    rows.sort(key=lambda r: (order[r["source"]], -(r["moa"] or 0)))
    return rows


def _direction(moa) -> str:
    if moa is None:
        return ""
    return "activating" if moa > 0 else ("repressive" if moa < 0 else "unknown")


def save_target_table(all_targets: dict, cancer_type: str, out_dir: str) -> None:
    """Write source-labeled target list to CSV."""
    ct   = cancer_type.upper()
    path = os.path.join(out_dir, f"target_list_{cancer_type.lower()}.csv")
    fields = ["focal_gene", "regulator", "source", "oncokb_role", "moa", "direction"]
    rows = []
    for focal_gene, targets in all_targets.items():
        for entry in targets:
            rows.append({
                "focal_gene":  focal_gene,
                "regulator":   entry["regulator"],
                "source":      entry["source"],
                "oncokb_role": entry.get("oncokb_role", ""),
                "moa":         round(entry["moa"], 3) if entry["moa"] is not None else "",
                "direction":   entry["direction"],
            })
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[{ct}] Target list ({len(rows)} entries) -> {path}")


def plot_target_list(all_targets: dict, cancer_type: str, out_dir: str) -> None:
    """Stacked bar: OncoKB targets per source per focal gene."""
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
    ax.set_ylabel("OncoKB-overlapping regulators", fontsize=10)
    ax.set_title(
        f"{ct}: Candidate therapeutic regulators by network source\n"
        "(OncoKB-filtered; source indicates which network context identified each regulator)",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"target_list_{cancer_type.lower()}.png"), dpi=150)
    plt.close()


# ── Per-cancer analysis helper ─────────────────────────────────────────────────

def run_cancer_analysis(
    agent,
    workflow,
    cancer_type: str,
    focal_genes: list,
    oncokb_raw: set,
    oncokb_roles: dict,
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

    oncokb = oncokb_raw & background
    print(f"[{CT}] OncoKB in background: {len(oncokb):,} genes")

    ref_sets = {"oncokb": oncokb}

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
        gene_res["target_list"]   = generate_target_list(comp, oncokb, moa_map, oncokb_roles)

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
    primary_ps = [results[g]["enrichment"]["oncokb"]["p_value"] for g in testable]
    fdr_vals   = bh_fdr(primary_ps)
    for gene, adj_p in zip(testable, fdr_vals):
        results[gene]["enrichment"]["oncokb"]["fdr_adjusted_p"] = round(adj_p, 6)

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
    plot_target_list(all_targets, ct, out_dir)

    return comparisons, results, background, combined_stats, testable


# ── Negative controls ──────────────────────────────────────────────────────────

def run_negative_controls(
    agent,
    cancer_type: str,
    oncokb_raw: set,
    background: set,
) -> dict:
    """
    Run the same enrichment test on housekeeping genes.
    Expected result: OR ≈ 1 (no enrichment) — validates that the enrichment
    seen for cancer driver genes is specific, not a general network property.
    """
    ct = cancer_type.lower()
    CT = cancer_type.upper()
    oncokb = oncokb_raw & background
    ref_sets = {"oncokb": oncokb}

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


# ── Exploratory: GREmLN-only vs TCGA-only enrichment comparison ────────────────

def run_gremln_comparison(
    agent,
    brca_comparisons: dict,
    coad_comparisons: dict,
    oncokb_raw: set,
    brca_bg: set,
    coad_bg: set,
) -> dict:
    """
    Exploratory comparison: GREmLN-only vs TCGA-only OncoKB enrichment.

    Uses the GREmLN epithelial_cell gene universe as background (translated from
    ENSG via the pre-built gene ID cache). Only ~35% of ENSG IDs resolve to gene
    symbols via the cache, so this is presented as qualitative support rather than
    a formal statistical comparison (see Discussion). Results are included in the
    output JSON for reproducibility.
    """
    print("\n" + "=" * 70)
    print("EXPLORATORY: GREmLN-only vs TCGA-only OncoKB enrichment comparison")
    print("NOTE: GREmLN background derived from ~35% ENSG-to-symbol resolution.")
    print("      Results are qualitative support, not formal statistics.")
    print("=" * 70)

    import pickle as _pickle
    cache_path = "cache/gene_id_cache.pkl"
    try:
        with open(cache_path, "rb") as f:
            id_cache = _pickle.load(f)
        e2s = id_cache.get("ensembl_to_symbol", {})
    except Exception as exc:
        print(f"WARNING: Could not load gene ID cache ({exc}). Skipping GREmLN comparison.")
        return {}

    gremln_idx = agent.cache.network_indices.get(CELL_TYPE, {})
    ensg_ids   = gremln_idx.get("all_genes", set())
    gremln_bg  = {e2s[e].upper() for e in ensg_ids if e in e2s}
    oncokb_g   = oncokb_raw & gremln_bg
    coverage   = 100 * len(gremln_bg) / max(len(ensg_ids), 1)

    print(f"\nGREmLN background: {len(gremln_bg):,} symbols from {len(ensg_ids):,} "
          f"ENSG IDs ({coverage:.0f}% coverage)")
    print(f"OncoKB in GREmLN background: {len(oncokb_g):,}")

    print(f"\n{'Gene':<10} {'Cancer':<6} {'TCGA-OR':>8} {'TCGA-n':>7} {'TCGA-ov':>8} "
          f"{'GREmLN-OR':>10} {'GREmLN-n':>9} {'GREmLN-ov':>10}")
    print("-" * 80)

    out: dict = {}
    for cancer, comparisons, bg, focal_genes in [
        ("brca", brca_comparisons, brca_bg, BRCA_GENES),
        ("coad", coad_comparisons, coad_bg, COAD_GENES),
    ]:
        oncokb_t = oncokb_raw & bg
        out[cancer] = {}
        for gene in focal_genes:
            if gene not in comparisons:
                continue
            comp        = comparisons[gene]
            tcga_only   = set(comp["regulators"]["tumor_state_only"])
            gremln_only = set(comp["regulators"]["population_averaged_only"])

            if len(tcga_only) >= 3:
                t    = fisher_enrichment(tcga_only, oncokb_t, bg)
                t_or, t_n, t_ov = t["odds_ratio"], t["query_size"], t["ref_overlap"]
            else:
                t_or, t_n, t_ov = 0.0, len(tcga_only), 0

            if len(gremln_only) >= 3:
                g    = fisher_enrichment(gremln_only, oncokb_g, gremln_bg)
                g_or, g_n, g_ov = g["odds_ratio"], g["query_size"], g["ref_overlap"]
            else:
                g_or, g_n, g_ov = 0.0, len(gremln_only), 0

            print(f"{gene:<10} {cancer:<6} {t_or:>8.2f} {t_n:>7} {t_ov:>8} "
                  f"{g_or:>10.2f} {g_n:>9} {g_ov:>10}")
            out[cancer][gene] = {
                "tcga_only_or":             t_or,
                "tcga_only_n":              t_n,
                "tcga_only_oncokb_overlap": t_ov,
                "gremln_only_or":           g_or,
                "gremln_only_n":            g_n,
                "gremln_only_oncokb_overlap": g_ov,
            }

    print(f"\nCoverage note: {len(gremln_bg):,}/{len(ensg_ids):,} ENSG IDs resolved to gene "
          f"symbols via pre-built cache. Incomplete coverage may affect GREmLN background "
          f"size and OR estimates; treat as exploratory.")
    return out


# ── Main ───────────────────────────────────────────────────────────────────────

def run_experiment() -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    random.seed(RANDOM_SEED)

    print("Loading RegNetAgents workflow...")
    workflow = RegNetAgentsWorkflow()
    agent    = workflow.modeling_agent

    for ct in ["brca", "coad"]:
        if not workflow.tcga_cache.tcga_indices.get(ct):
            print(f"ERROR: No TCGA cache for '{ct}'. Run build_tcga_cache.py first.")
            sys.exit(1)

    download_oncokb_if_missing(ONCOKB_PATH)
    oncokb_raw   = load_oncokb(ONCOKB_PATH)
    oncokb_roles = load_oncokb_roles(ONCOKB_PATH)
    print(f"OncoKB drivers: {len(oncokb_raw):,} total")

    # ── BRCA analysis ──────────────────────────────────────────────────────────
    (brca_comp, brca_res, brca_bg,
     brca_stouffer, brca_testable) = run_cancer_analysis(
        agent, workflow, "brca", BRCA_GENES, oncokb_raw, oncokb_roles, RESULTS_DIR,
    )
    brca_neg = run_negative_controls(agent, "brca", oncokb_raw, brca_bg)
    plot_neg_controls(brca_res, brca_neg, brca_testable, HOUSEKEEPING_GENES, "brca")

    # ── COAD analysis ──────────────────────────────────────────────────────────
    (coad_comp, coad_res, coad_bg,
     coad_stouffer, coad_testable) = run_cancer_analysis(
        agent, workflow, "coad", COAD_GENES, oncokb_raw, oncokb_roles, RESULTS_DIR,
    )
    coad_neg = run_negative_controls(agent, "coad", oncokb_raw, coad_bg)
    plot_neg_controls(coad_res, coad_neg, coad_testable, HOUSEKEEPING_GENES, "coad")

    # ── Exploratory: GREmLN-only vs TCGA-only enrichment comparison ────────────
    gremln_comparison = run_gremln_comparison(
        agent,
        brca_comparisons=brca_comp,
        coad_comparisons=coad_comp,
        oncokb_raw=oncokb_raw,
        brca_bg=brca_bg,
        coad_bg=coad_bg,
    )

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
        "gremln_comparison_exploratory": gremln_comparison,
    }
    plot_workflow_figure()

    out_json = os.path.join(RESULTS_DIR, "experiment_rewiring_results.json")
    with open(out_json, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults -> {out_json}")
    print(f"Figures  -> {MANUSCRIPT_DIR}/figure_workflow.png  (NAR Fig 1)")
    print(f"            {MANUSCRIPT_DIR}/figure_heatmap_brca.png  (NAR Fig 3A)")
    print(f"            {MANUSCRIPT_DIR}/figure_heatmap_coad.png  (NAR Fig 3B)")
    print(f"            {MANUSCRIPT_DIR}/figure_negcontrol_brca.png  (NAR Fig 4A)")
    print(f"            {MANUSCRIPT_DIR}/figure_negcontrol_coad.png  (NAR Fig 4B)")
    print(f"            {RESULTS_DIR}/target_list_brca.png  (NAR Fig 2A)")
    print(f"            {RESULTS_DIR}/target_list_coad.png  (NAR Fig 2B)")
    print(f"            {RESULTS_DIR}/experiment_rewiring_barchart_brca.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_barchart_coad.png")
    print("\nDone.")


if __name__ == "__main__":
    run_experiment()
