#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
experiment_tcga_rewiring.py
===========================
Validation experiment for the RegNetAgents research paper (Briefings in
Bioinformatics).

Two-cancer validation:
  1. BRCA: tumor-specific regulators enriched in OncoKB
  2. COAD: same test for colorectal cancer context (replication)

Reference datasets (pre-downloaded to data/):
  - data/oncokb_cancer_genes.tsv   : OncoKB pan-cancer driver list

Outputs (written to results/):
  - results/experiment_rewiring_results.json       : full statistics
  - results/experiment_rewiring_heatmap_brca.png   : OR heatmap, BRCA
  - results/experiment_rewiring_heatmap_coad.png   : OR heatmap, COAD
  - results/experiment_rewiring_barchart_brca.png  : regulator count bar chart, BRCA
  - results/experiment_rewiring_barchart_coad.png  : regulator count bar chart, COAD
  - results/experiment_rewiring_negcontrol_brca.png: negative controls, BRCA
  - results/experiment_rewiring_negcontrol_coad.png: negative controls, COAD

Background: all genes in each cancer type's TCGA network (symbol-native PKL).

Usage:
  python scripts/experiment_tcga_rewiring.py

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

ONCOKB_PATH = "data/oncokb_cancer_genes.tsv"
RESULTS_DIR = "results"

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

def plot_or_heatmap(
    results: dict, genes: list, cancer_type: str, out_dir: str
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
        cmap="RdYlGn",
        center=1.0,
        vmin=0,
        vmax=5,
        ax=ax,
        linewidths=0.5,
        cbar_kws={"label": "Odds Ratio"},
    )
    ax.set_title(
        f"Enrichment of {CT}-specific regulators in cancer driver gene sets\n"
        "(* BH-FDR < 0.05, ** BH-FDR < 0.01; primary test = OncoKB)",
        fontsize=11,
    )
    ax.set_xlabel("Reference gene set", fontsize=10)
    ax.set_ylabel("Focal gene", fontsize=10)
    plt.tight_layout()
    plt.savefig(
        os.path.join(out_dir, f"experiment_rewiring_heatmap_{ct}.png"), dpi=150
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
    out_dir: str,
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
        os.path.join(out_dir, f"experiment_rewiring_negcontrol_{cancer_type.lower()}.png"), dpi=150
    )
    plt.close()


# ── Per-cancer analysis helper ─────────────────────────────────────────────────

def run_cancer_analysis(
    agent,
    workflow,
    cancer_type: str,
    focal_genes: list,
    oncokb_raw: set,
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
            f"conserved={r['conserved_count']}  "
            f"epithelial_only={len(r['population_averaged_only'])}  "
            f"{ct}_only={len(r['tumor_state_only'])}  "
            f"rewiring={result['interpretation']['regulatory_rewiring']}"
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

    # Figures
    plot_or_heatmap(results, testable, ct, out_dir)
    plot_regulator_counts(comparisons, testable, results, ct, out_dir)

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
            f"conserved={r['conserved_count']}  "
            f"{ct}_only={len(specific)}  "
            f"rewiring={result['interpretation']['regulatory_rewiring']}"
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
    oncokb_raw = load_oncokb(ONCOKB_PATH)
    print(f"OncoKB drivers: {len(oncokb_raw):,} total")

    # ── BRCA analysis ──────────────────────────────────────────────────────────
    (_, brca_res, brca_bg,
     brca_stouffer, brca_testable) = run_cancer_analysis(
        agent, workflow, "brca", BRCA_GENES, oncokb_raw, RESULTS_DIR,
    )
    brca_neg = run_negative_controls(agent, "brca", oncokb_raw, brca_bg)
    plot_neg_controls(brca_res, brca_neg, brca_testable, HOUSEKEEPING_GENES, "brca", RESULTS_DIR)

    # ── COAD analysis ──────────────────────────────────────────────────────────
    (_, coad_res, coad_bg,
     coad_stouffer, coad_testable) = run_cancer_analysis(
        agent, workflow, "coad", COAD_GENES, oncokb_raw, RESULTS_DIR,
    )
    coad_neg = run_negative_controls(agent, "coad", oncokb_raw, coad_bg)
    plot_neg_controls(coad_res, coad_neg, coad_testable, HOUSEKEEPING_GENES, "coad", RESULTS_DIR)

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
    }
    out_json = os.path.join(RESULTS_DIR, "experiment_rewiring_results.json")
    with open(out_json, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults -> {out_json}")
    print(f"Figures  -> {RESULTS_DIR}/experiment_rewiring_heatmap_brca.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_heatmap_coad.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_barchart_brca.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_barchart_coad.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_negcontrol_brca.png")
    print(f"            {RESULTS_DIR}/experiment_rewiring_negcontrol_coad.png")
    print("\nDone.")


if __name__ == "__main__":
    run_experiment()
