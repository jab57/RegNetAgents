#!/usr/bin/env python3
"""
Convert TCGA ARACNe network CSVs to pickle cache format for RegNetAgents.

Source data: Figshare pre-exported CSVs from the Bioconductor aracne.networks
package (Lim & Califano).  Download from:

    https://figshare.com/s/5d1ffd9f8b2e86e37ed6

Rename each file (e.g. brca_regul.csv → network.csv) and place at:

    models/networks/tcga/{cancer_type}/network.csv

Then run this script:

    python scripts/build_tcga_cache.py --all
    python scripts/build_tcga_cache.py --cancer-type brca

CSV column format (Figshare aracne.networks):
    Regulator, Target, MoA, Likelihood
    - MoA      : mode of action  — +1 activation, -1 repression, 0 unknown
    - Likelihood: edge weight, 0–1

Gene ID strategy:
    - Networks are symbol-native (nodes keyed by uppercase gene symbol).
    - MyGene.info batch POST is used at build time to validate that symbols
      are real human genes (≥90% must resolve to canonical Ensembl IDs).
    - PKL is stored as symbol-keyed (id_type: "symbol").
    - Do NOT use GeneIDMapper — it creates synthetic IDs that mask true resolution.
"""

import argparse
import csv
import json
import os
import pickle
import sys
import time
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Set, Tuple

import networkx as nx
import urllib.request
import urllib.parse
import urllib.error

# Allow importing from the same scripts/ directory
sys.path.insert(0, os.path.dirname(__file__))
from build_network_cache import compute_thresholds, update_threshold_config

# Ensure project root is on path for regnetagents imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from regnetagents.tcga_registry import TCGA_NETWORK_REGISTRY, TCGA_CANCER_TYPES

CACHE_VERSION = 4
MYGENE_URL = "https://mygene.info/v3/gene"


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

def load_tcga_csv(
    csv_path: str,
) -> Tuple[
    Dict[str, List[str]],
    Dict[str, List[str]],
    Set[str],
    Dict[str, Dict[str, float]],
    Dict[str, Dict[str, float]],
]:
    """
    Parse a TCGA ARACNe CSV (columns: Regulator, Target, MoA, Likelihood).

    Returns:
        (regulator_targets, target_regulators, all_genes,
         regulator_target_mi, regulator_target_moa)
    where:
        regulator_target_mi[reg][tgt]  = Likelihood score (0–1)
        regulator_target_moa[reg][tgt] = MoA value (-1, 0, or +1)
    """
    regulator_targets: Dict[str, List[str]] = defaultdict(list)
    target_regulators: Dict[str, List[str]] = defaultdict(list)
    all_genes: Set[str] = set()
    regulator_target_mi: Dict[str, Dict[str, float]] = defaultdict(dict)
    regulator_target_moa: Dict[str, Dict[str, float]] = defaultdict(dict)

    print(f"Loading TCGA network from {csv_path} ...")

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"Regulator", "Target", "MoA", "Likelihood"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"CSV must have columns {required}. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            reg = row["Regulator"].strip().upper()
            tgt = row["Target"].strip().upper()
            if not reg or not tgt:
                continue

            likelihood = float(row.get("Likelihood") or 0)
            moa = float(row.get("MoA") or 0)

            regulator_targets[reg].append(tgt)
            target_regulators[tgt].append(reg)
            all_genes.add(reg)
            all_genes.add(tgt)
            regulator_target_mi[reg][tgt] = likelihood
            regulator_target_moa[reg][tgt] = moa

    print(
        f"  Loaded {len(all_genes):,} genes, "
        f"{sum(len(v) for v in regulator_targets.values()):,} edges, "
        f"{len(regulator_targets):,} regulons"
    )
    return (
        dict(regulator_targets),
        dict(target_regulators),
        all_genes,
        dict(regulator_target_mi),
        dict(regulator_target_moa),
    )


# ---------------------------------------------------------------------------
# MyGene.info symbol validation
# ---------------------------------------------------------------------------

def validate_symbols_mygene(
    symbols: List[str],
    min_resolution_rate: float = 0.90,
) -> Dict[str, str]:
    """
    POST to MyGene.info to validate that gene symbols are real human genes.

    Classifies each symbol into:
        - resolved   : canonical Ensembl ENSG... returned
        - alias      : multiple hits, at least one canonical match
        - unresolved : no confident hit (logged but kept in network as-is)

    Aborts if (resolved + alias) / total < min_resolution_rate.

    Returns:
        Dict mapping symbol → canonical Ensembl ID (for resolved/alias only).
        Unresolved symbols are not present in the returned dict.
    """
    print(f"Validating {len(symbols):,} unique symbols via MyGene.info ...")

    # Batch into chunks of 1000 (API limit)
    batch_size = 1000
    resolved: Dict[str, str] = {}
    unresolved: List[str] = []

    for i in range(0, len(symbols), batch_size):
        batch = symbols[i : i + batch_size]
        payload = json.dumps({
            "q": batch,
            "scopes": "symbol",
            "fields": "ensembl.gene",
            "species": "human",
        }).encode()

        req = urllib.request.Request(
            MYGENE_URL,
            data=payload,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                hits = json.loads(resp.read())
        except urllib.error.URLError as exc:
            print(f"  WARNING: MyGene.info request failed ({exc}). Skipping validation.")
            return {}

        for hit in hits:
            symbol = hit.get("query", "").upper()
            if hit.get("notfound"):
                unresolved.append(symbol)
                continue

            ensg = hit.get("ensembl", {})
            if isinstance(ensg, list):
                # Multiple Ensembl entries — pick first canonical ENSG
                ensg_ids = [
                    e["gene"] for e in ensg
                    if isinstance(e, dict) and e.get("gene", "").startswith("ENSG")
                ]
                if ensg_ids:
                    resolved[symbol] = ensg_ids[0]
                else:
                    unresolved.append(symbol)
            elif isinstance(ensg, dict):
                gene_id = ensg.get("gene", "")
                if gene_id.startswith("ENSG"):
                    resolved[symbol] = gene_id
                else:
                    unresolved.append(symbol)
            else:
                unresolved.append(symbol)

        # Polite pacing between batches
        if i + batch_size < len(symbols):
            time.sleep(0.5)

    total = len(symbols)
    n_resolved = len(resolved)
    n_unresolved = len(unresolved)
    rate = n_resolved / total if total else 0

    print(
        f"  Symbol resolution: {n_resolved:,} resolved, "
        f"{n_unresolved:,} unresolved  ({rate:.1%})"
    )

    if unresolved:
        print(f"  Unresolved symbols (first 20): {unresolved[:20]}")

    if rate < min_resolution_rate:
        raise RuntimeError(
            f"Symbol resolution rate {rate:.1%} is below the required "
            f"{min_resolution_rate:.0%} threshold. "
            "Check that the CSV is a valid human TCGA ARACNe network."
        )

    return resolved


# ---------------------------------------------------------------------------
# PageRank
# ---------------------------------------------------------------------------

def compute_pagerank(
    regulator_targets: Dict[str, List[str]],
    all_genes: Set[str],
) -> Dict[str, float]:
    """Compute max-normalised PageRank for all genes."""
    print("  Computing PageRank ...")
    try:
        G = nx.DiGraph()
        G.add_nodes_from(all_genes)
        for reg, targets in regulator_targets.items():
            for tgt in targets:
                G.add_edge(reg, tgt)

        raw = nx.pagerank(G, alpha=0.85, max_iter=100, tol=1e-6)
        max_val = max(raw.values()) if raw else 1.0
        if max_val == 0:
            max_val = 1.0
        return {g: round(v / max_val, 6) for g, v in raw.items()}
    except Exception as exc:
        print(f"  WARNING: PageRank failed ({exc}). Returning empty dict.")
        return {}


# ---------------------------------------------------------------------------
# Main build function
# ---------------------------------------------------------------------------

def build_tcga_cache(
    cancer_type: str,
    output_dir: str = "models/networks/tcga",
    skip_validation: bool = False,
) -> None:
    """
    Build the pickle cache for a single TCGA cancer type.

    Args:
        cancer_type   : One of the keys in TCGA_NETWORK_REGISTRY
        output_dir    : Root directory for TCGA caches
        skip_validation: Skip MyGene.info symbol validation (for offline use)
    """
    if cancer_type not in TCGA_NETWORK_REGISTRY:
        raise ValueError(
            f"Unknown cancer type '{cancer_type}'. "
            f"Valid types: {TCGA_CANCER_TYPES}"
        )

    entry = TCGA_NETWORK_REGISTRY[cancer_type]
    csv_path = entry["csv"]
    pkl_dir = os.path.join(output_dir, cancer_type)
    pkl_path = os.path.join(pkl_dir, "network_index.pkl")
    threshold_config_path = os.path.join("models", "threshold_config.json")

    print(f"\n{'='*60}")
    print(f"Building TCGA cache: {cancer_type} ({entry['label']})")
    print(f"{'='*60}")

    if not os.path.exists(csv_path):
        print(
            f"  SKIP — CSV not found: {csv_path}\n"
            "  Download from https://figshare.com/s/5d1ffd9f8b2e86e37ed6 "
            "and place at the path above."
        )
        return

    # 1. Load CSV
    (
        regulator_targets,
        target_regulators,
        all_genes,
        regulator_target_mi,
        regulator_target_moa,
    ) = load_tcga_csv(csv_path)

    # 2. Validate symbols with MyGene.info
    if skip_validation:
        print("  Skipping symbol validation (--skip-validation flag set)")
    else:
        validate_symbols_mygene(sorted(all_genes))

    # 3. PageRank
    pagerank_normalized = compute_pagerank(regulator_targets, all_genes)

    # 4. Thresholds
    print("  Computing empirical thresholds ...")
    thresholds = compute_thresholds(
        regulator_targets,
        target_regulators,
        all_genes,
        pagerank_normalized,
    )
    threshold_key = f"tcga_{cancer_type}"
    update_threshold_config(threshold_key, thresholds, threshold_config_path)

    # 5. Build PKL
    # bootstrap counts are not present in TCGA CSVs — store zeros
    regulator_target_count = {
        reg: {tgt: 0 for tgt in targets}
        for reg, targets in regulator_target_mi.items()
    }

    cache_data = {
        "regulator_targets":      dict(regulator_targets),
        "target_regulators":      dict(target_regulators),
        "all_genes":              sorted(all_genes),
        "regulator_target_mi":    dict(regulator_target_mi),
        "regulator_target_count": regulator_target_count,
        "regulator_target_moa":   dict(regulator_target_moa),
        "pagerank_normalized":    pagerank_normalized,
        "num_genes":              len(all_genes),
        "num_edges":              sum(len(v) for v in regulator_targets.values()),
        "num_regulons":           len(regulator_targets),
        "id_type":                "symbol",
        "source":                 "tcga_aracne",
        "cancer_type":            cancer_type,
        "cache_version":          CACHE_VERSION,
        "created":                datetime.now().isoformat(),
    }

    os.makedirs(pkl_dir, exist_ok=True)
    with open(pkl_path, "wb") as fh:
        pickle.dump(cache_data, fh, protocol=4)

    print(
        f"  Cache written: {pkl_path}\n"
        f"  Genes: {cache_data['num_genes']:,}  "
        f"Edges: {cache_data['num_edges']:,}  "
        f"Regulons: {cache_data['num_regulons']:,}"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build TCGA ARACNe network caches for RegNetAgents."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--all",
        action="store_true",
        help="Build caches for all 8 TCGA cancer types",
    )
    group.add_argument(
        "--cancer-type",
        choices=TCGA_CANCER_TYPES,
        metavar="CANCER_TYPE",
        help=f"Build cache for a single cancer type. "
             f"Choices: {', '.join(TCGA_CANCER_TYPES)}",
    )
    parser.add_argument(
        "--output-dir",
        default="models/networks/tcga",
        help="Root directory for output caches (default: models/networks/tcga)",
    )
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip MyGene.info symbol validation (useful for offline/testing)",
    )
    args = parser.parse_args()

    cancer_types = TCGA_CANCER_TYPES if args.all else [args.cancer_type]

    for ct in cancer_types:
        build_tcga_cache(
            cancer_type=ct,
            output_dir=args.output_dir,
            skip_validation=args.skip_validation,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
