#!/usr/bin/env python3
"""
Script to convert TSV network files to pickle cache format for RegNetAgents.

This script processes the original .tsv network files mentioned in the README
and converts them into the optimized .pkl cache format used by the system.

Expected TSV format:
- ARACNe network output format with header
- Columns: regulator.values \t target.values \t mi.values \t scc.values \t count.values \t log.p.values
- Uses columns 0-1 (regulator, target) and columns 2, 4 (mi.values, count.values) for edge confidence
- Gene IDs should be Ensembl format (ENSG...)
- Header line is automatically detected and skipped

Usage:
    python build_network_cache.py [cell_type] [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]

    # Process all cell types from models/networks/
    python build_network_cache.py --all

    # Process specific cell type
    python build_network_cache.py epithelial_cell

    # Custom input/output directories
    python build_network_cache.py epithelial_cell --input-dir /path/to/tsv --output-dir /path/to/output

    # Enrich gene_id_cache.pkl only (no PKL rebuild, requires processed PKLs already present)
    python build_network_cache.py --enrich-gene-cache

Note: --all automatically enriches cache/gene_id_cache.pkl with ENSG→symbol mappings
for all GREmLN genes via MyGene.info (requires internet access).
"""

import os
import pickle
import argparse
import json
import time
import urllib.request
import urllib.error
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Set, Tuple
import sys
import numpy as np
import networkx as nx

MYGENE_URL = "https://mygene.info/v3/query"

# Cell types supported by RegNetAgents (with available network data)
SUPPORTED_CELL_TYPES = [
    'cd14_monocytes',
    'cd16_monocytes',
    'cd20_b_cells',
    'cd4_t_cells',
    'cd8_t_cells',
    'erythrocytes',
    'nk_cells',
    'nkt_cells',
    'epithelial_cell',
    'monocyte-derived_dendritic_cells'
]

def load_tsv_network(tsv_file: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Set[str], Dict[str, Dict[str, float]], Dict[str, Dict[str, int]]]:
    """
    Load network from TSV file and build regulator-target mappings.

    Args:
        tsv_file: Path to TSV file with regulator-target pairs

    Returns:
        Tuple of (regulator_targets, target_regulators, all_genes,
                  regulator_target_mi, regulator_target_count)
        where regulator_target_mi[reg][tgt] = mutual information score
        and   regulator_target_count[reg][tgt] = bootstrap reproducibility count
    """
    regulator_targets = defaultdict(list)
    target_regulators = defaultdict(list)
    all_genes = set()
    regulator_target_mi = defaultdict(dict)
    regulator_target_count = defaultdict(dict)

    print(f"Loading network from {tsv_file}...")

    if not os.path.exists(tsv_file):
        raise FileNotFoundError(f"TSV file not found: {tsv_file}")

    with open(tsv_file, 'r') as f:
        header_skipped = False
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                print(f"Warning: Invalid line {line_num} in {tsv_file}: {line}")
                continue

            # Skip header line in ARACNe format
            if not header_skipped and (parts[0] == 'regulator.values' or 'regulator' in parts[0].lower()):
                print(f"Skipping header line: {line}")
                header_skipped = True
                continue

            regulator = parts[0].strip()
            target = parts[1].strip()

            # Basic validation for Ensembl IDs
            if not regulator.startswith('ENSG') or not target.startswith('ENSG'):
                if line_num <= 10:  # Only show warnings for first 10 lines to avoid spam
                    print(f"Warning: Non-Ensembl IDs on line {line_num}: {regulator} -> {target}")
                continue

            regulator_targets[regulator].append(target)
            target_regulators[target].append(regulator)
            all_genes.add(regulator)
            all_genes.add(target)

            # Parse MI (col 2) and bootstrap count (col 4)
            try:
                mi_score = float(parts[2]) if len(parts) > 2 else 0.0
            except ValueError:
                mi_score = 0.0
            try:
                boot_count = int(parts[4]) if len(parts) > 4 else 0
            except ValueError:
                boot_count = 0
            regulator_target_mi[regulator][target] = mi_score
            regulator_target_count[regulator][target] = boot_count

    # Convert defaultdicts to regular dicts
    regulator_targets = dict(regulator_targets)
    target_regulators = dict(target_regulators)
    regulator_target_mi = dict(regulator_target_mi)
    regulator_target_count = dict(regulator_target_count)

    print(f"Loaded {len(regulator_targets)} regulators, {len(target_regulators)} targets, {len(all_genes)} total genes")

    return regulator_targets, target_regulators, all_genes, regulator_target_mi, regulator_target_count

def calculate_stats(regulator_targets: Dict[str, List[str]],
                   target_regulators: Dict[str, List[str]],
                   all_genes: Set[str]) -> Dict[str, int]:
    """Calculate network statistics."""
    num_edges = sum(len(targets) for targets in regulator_targets.values())
    num_genes = len(all_genes)
    num_regulons = len(regulator_targets)

    return {
        'num_edges': num_edges,
        'num_genes': num_genes,
        'num_regulons': num_regulons
    }

def calculate_pagerank(regulator_targets: Dict[str, List[str]],
                      num_genes: int) -> Dict[str, float]:
    """
    Calculate PageRank for all genes in the network.

    Args:
        regulator_targets: Dictionary mapping regulators to their targets
        num_genes: Total number of genes in network

    Returns:
        Dictionary of normalized PageRank scores for each gene
    """
    print(f"  Calculating PageRank for {num_genes} genes...")

    # Build NetworkX directed graph
    G = nx.DiGraph()
    edge_count = 0
    for regulator, targets in regulator_targets.items():
        for target in targets:
            G.add_edge(regulator, target)
            edge_count += 1

    print(f"  Built graph with {G.number_of_nodes()} nodes and {edge_count} edges")

    try:
        # Calculate PageRank with standard parameters
        # alpha=0.85: damping factor (standard value from PageRank algorithm)
        # max_iter=100: maximum iterations for convergence
        # tol=1e-06: convergence tolerance
        pagerank_scores = nx.pagerank(G, alpha=0.85, max_iter=100, tol=1e-06)

        # Normalize by maximum value for interpretability
        max_pagerank = max(pagerank_scores.values()) if pagerank_scores else 1.0
        pagerank_normalized = {gene_id: score / max_pagerank
                              for gene_id, score in pagerank_scores.items()}

        print(f"  PageRank calculation complete (max score: {max_pagerank:.6f})")
        return pagerank_normalized

    except Exception as e:
        print(f"  WARNING: PageRank calculation failed: {e}")
        print(f"  Returning empty PageRank dict (will fall back to on-demand calculation)")
        return {}

def compute_thresholds(
    regulator_targets: Dict[str, List[str]],
    target_regulators: Dict[str, List[str]],
    all_genes: Set[str],
    pagerank_normalized: Dict[str, float]
) -> Dict[str, float]:
    """
    Compute empirical threshold defaults from the degree distributions of a
    cell-type network.

    Thresholds are set at the 90th percentile (high) and 75th percentile
    (moderate) of target count, regulator count, and normalized PageRank
    across all genes in the network. This ensures 'high' evidence reflects
    genuine outliers within each tissue context rather than fixed absolute
    values. A minimum of 1 is applied to target and regulator thresholds so
    that the moderate tier always requires at least one connection.

    These values serve as biologically informed defaults and can be overridden
    by the user at runtime via ThresholdConfig. They are empirically motivated,
    not optimized, and should be adjusted when applying RegNetAgents to
    networks with substantially different topology (e.g., bulk RNA-seq derived
    networks).

    Args:
        regulator_targets: Mapping of regulator gene IDs to their target lists
        target_regulators: Mapping of target gene IDs to their regulator lists
        all_genes: Set of all gene IDs in the network
        pagerank_normalized: Normalized PageRank scores (max-normalized to 1.0)

    Returns:
        Dict with keys: target_high, target_moderate, regulator_high,
        regulator_moderate, pagerank_high
    """
    target_counts = [len(regulator_targets.get(g, [])) for g in all_genes]
    regulator_counts = [len(target_regulators.get(g, [])) for g in all_genes]
    pagerank_values = [pagerank_normalized.get(g, 0.0) for g in all_genes]

    def pct(values, p):
        return float(np.percentile(values, p))

    return {
        "target_high":      max(1, round(pct(target_counts, 90))),
        "target_moderate":  max(1, round(pct(target_counts, 75))),
        "regulator_high":   max(1, round(pct(regulator_counts, 90))),
        "regulator_moderate": max(1, round(pct(regulator_counts, 75))),
        "pagerank_high":    round(pct(pagerank_values, 90), 4),
        "percentiles_used": {"high": 90, "moderate": 75},
        "n_genes": len(all_genes),
    }


def update_threshold_config(cell_type: str, thresholds: Dict, config_path: str) -> None:
    """
    Update models/threshold_config.json with thresholds for a cell type.

    Creates the file if it does not exist. Existing entries for other cell
    types are preserved. Called automatically by build_network_cache so that
    adding a new cell type always produces a matching threshold entry.

    Args:
        cell_type: Cell type name (e.g. 'epithelial_cell')
        thresholds: Dict returned by compute_thresholds()
        config_path: Path to threshold_config.json
    """
    config = {}
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            config = json.load(f)

    config[cell_type] = thresholds

    os.makedirs(os.path.dirname(config_path), exist_ok=True)
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    print(f"  Thresholds saved to {config_path}")


def build_network_cache(tsv_file: str, output_file: str) -> dict:
    """
    Convert TSV network file to pickle cache format.

    Args:
        tsv_file: Path to input TSV file
        output_file: Path to output pickle file

    Returns:
        Thresholds dict computed from this network's topology distributions
    """
    # Load network data
    regulator_targets, target_regulators, all_genes, regulator_target_mi, regulator_target_count = load_tsv_network(tsv_file)

    # Calculate statistics
    stats = calculate_stats(regulator_targets, target_regulators, all_genes)

    # Calculate PageRank (pre-compute for performance)
    print("Calculating PageRank centrality...")
    pagerank_normalized = calculate_pagerank(regulator_targets, stats['num_genes'])

    # Compute empirical thresholds from this network's topology distributions
    print("Computing empirical thresholds from network topology distributions...")
    thresholds = compute_thresholds(regulator_targets, target_regulators, all_genes, pagerank_normalized)
    print(f"  target_high={thresholds['target_high']}, target_moderate={thresholds['target_moderate']}, "
          f"regulator_high={thresholds['regulator_high']}, regulator_moderate={thresholds['regulator_moderate']}, "
          f"pagerank_high={thresholds['pagerank_high']}")

    # Build cache data structure
    cache_data = {
        'regulator_targets': regulator_targets,
        'target_regulators': target_regulators,
        'all_genes': sorted(list(all_genes)),  # Sort for consistency
        'num_edges': stats['num_edges'],
        'num_genes': stats['num_genes'],
        'num_regulons': stats['num_regulons'],
        'pagerank_normalized': pagerank_normalized,  # Pre-computed PageRank
        'pagerank_params': {  # Store parameters for reproducibility
            'alpha': 0.85,
            'max_iter': 100,
            'tol': 1e-06
        },
        'regulator_target_mi': regulator_target_mi,        # {reg: {tgt: mi_score}}
        'regulator_target_count': regulator_target_count,  # {reg: {tgt: bootstrap_count}}
        'cache_version': 3,  # Version 3: adds ARACNe edge confidence (MI score, bootstrap count)
        'created': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Save to pickle file
    print(f"Saving cache to {output_file}...")
    with open(output_file, 'wb') as f:
        pickle.dump(cache_data, f)

    print(f"Cache created successfully:")
    print(f"  - {stats['num_regulons']} regulators")
    print(f"  - {len(target_regulators)} targets")
    print(f"  - {stats['num_genes']} total genes")
    print(f"  - {stats['num_edges']} total edges")
    print(f"  - Output: {output_file}")

    return thresholds

def validate_cache(cache_file: str) -> bool:
    """
    Validate a generated cache file.

    Args:
        cache_file: Path to the cache file

    Returns:
        True if valid, False otherwise
    """
    try:
        print(f"Validating cache: {cache_file}")
        with open(cache_file, 'rb') as f:
            data = pickle.load(f)

        # Required keys for all cache versions
        required_keys = ['regulator_targets', 'target_regulators', 'all_genes', 'num_edges', 'num_genes', 'num_regulons']
        for key in required_keys:
            if key not in data:
                print(f"ERROR: Missing key '{key}' in cache")
                return False

        # Basic sanity checks
        if len(data['all_genes']) != data['num_genes']:
            print(f"ERROR: Gene count mismatch: {len(data['all_genes'])} vs {data['num_genes']}")
            return False

        if len(data['regulator_targets']) != data['num_regulons']:
            print(f"ERROR: Regulator count mismatch: {len(data['regulator_targets'])} vs {data['num_regulons']}")
            return False

        # Check cache version and optional feature data
        cache_version = data.get('cache_version', 1)  # Default to version 1 if not present
        if cache_version >= 2:
            if 'pagerank_normalized' in data:
                pagerank_count = len(data['pagerank_normalized'])
                print(f"  Cache version {cache_version} with {pagerank_count} PageRank scores")
            else:
                print(f"  WARNING: Cache version {cache_version} missing PageRank data")
        if cache_version >= 3:
            if 'regulator_target_mi' in data and 'regulator_target_count' in data:
                edge_count = sum(len(v) for v in data['regulator_target_mi'].values())
                print(f"  Edge confidence data present ({edge_count} MI scores)")
            else:
                print(f"  WARNING: Cache version {cache_version} missing edge confidence data")
        if cache_version < 2:
            print(f"  Cache version 1 (legacy, no PageRank - will calculate on-demand)")

        print(f"Cache validation passed")
        return True

    except Exception as e:
        print(f"ERROR validating cache: {e}")
        return False

def process_cell_type(cell_type: str, input_dir: str, output_dir: str) -> bool:
    """
    Process a single cell type.

    Builds the network index pickle and computes empirical thresholds from
    that network's topology distributions, then updates
    models/threshold_config.json automatically.

    Args:
        cell_type: Name of cell type
        input_dir: Directory containing TSV files
        output_dir: Directory for output pickle files

    Returns:
        True if successful, False otherwise
    """
    print(f"\n=== Processing {cell_type} ===")

    tsv_file = os.path.join(input_dir, cell_type, "network.tsv")
    output_file = os.path.join(output_dir, cell_type, "network_index.pkl")
    config_path = os.path.join(output_dir, "..", "threshold_config.json")

    try:
        thresholds = build_network_cache(tsv_file, output_file)

        # Update threshold config with this cell type's empirical thresholds
        update_threshold_config(cell_type, thresholds, os.path.normpath(config_path))

        # Validate the generated cache
        if validate_cache(output_file):
            print(f"SUCCESS: {cell_type} processed successfully")
            return True
        else:
            print(f"FAILED: {cell_type} cache validation failed")
            return False

    except Exception as e:
        print(f"ERROR processing {cell_type}: {e}")
        return False

def ensg_to_symbol_batch(ensg_ids: list, batch_size: int = 500) -> dict:
    """
    Convert Ensembl gene IDs to gene symbols via MyGene.info batch POST.
    Returns dict: ensg_id -> symbol (upper-cased).
    Requires internet access.
    """
    print(f"  Converting {len(ensg_ids):,} ENSG IDs to symbols via MyGene.info ...")
    id_to_symbol = {}
    unresolved = []

    for i in range(0, len(ensg_ids), batch_size):
        batch = ensg_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(ensg_ids) + batch_size - 1) // batch_size
        print(f"    Batch {batch_num}/{total_batches} ({len(batch)} IDs) ...", end=" ", flush=True)

        payload = json.dumps({
            "q": batch,
            "scopes": "ensembl.gene",
            "fields": "symbol",
            "species": "human",
        }).encode("utf-8")

        req = urllib.request.Request(
            MYGENE_URL,
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            method="POST",
        )

        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                hits = json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, Exception) as e:
            print(f"ERROR: {e}")
            unresolved.extend(batch)
            continue

        resolved_this_batch = 0
        for hit in hits:
            ensg = hit.get("query", "")
            if hit.get("notfound") or not ensg:
                unresolved.append(ensg)
                continue
            symbol = hit.get("symbol", "")
            if symbol:
                id_to_symbol[ensg] = symbol.upper()
                resolved_this_batch += 1
            else:
                unresolved.append(ensg)

        print(f"resolved {resolved_this_batch}/{len(batch)}")

        if i + batch_size < len(ensg_ids):
            time.sleep(0.3)

    total = len(ensg_ids)
    n_resolved = len(id_to_symbol)
    print(f"  Total resolved: {n_resolved:,}/{total:,} ({n_resolved/total:.1%})")
    return id_to_symbol


def update_gene_id_cache(output_dir: str, cache_path: str = "cache/gene_id_cache.pkl") -> None:
    """
    Collect all ENSG IDs from processed GREmLN network PKLs, resolve any that
    are missing from gene_id_cache.pkl via MyGene.info, and update the cache.

    Run automatically after --all, or manually with --enrich-gene-cache.
    """
    print("\n=== Enriching gene_id_cache.pkl with GREmLN ENSG IDs ===")

    # Load existing cache
    id_cache = {"symbol_to_ensembl": {}, "ensembl_to_symbol": {}}
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            id_cache = pickle.load(f)
    e2s = id_cache.get("ensembl_to_symbol", {})
    s2e = id_cache.get("symbol_to_ensembl", {})
    print(f"  Existing cache: {len(e2s):,} ENSG->symbol entries")

    # Collect all ENSG IDs from all processed PKLs
    all_ensg: set = set()
    for ct in SUPPORTED_CELL_TYPES:
        pkl_path = os.path.join(output_dir, ct, "network_index.pkl")
        if not os.path.exists(pkl_path):
            continue
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
        genes = data.get("all_genes", [])
        ensg_in_ct = {g for g in genes if str(g).startswith("ENSG")}
        all_ensg.update(ensg_in_ct)

    print(f"  Total unique ENSG IDs across all cell types: {len(all_ensg):,}")

    missing = sorted(all_ensg - set(e2s.keys()))
    print(f"  Missing from cache: {len(missing):,}")

    if not missing:
        print("  Cache already complete — nothing to do.")
        return

    new_mappings = ensg_to_symbol_batch(missing)

    # Update both directions
    for ensg, symbol in new_mappings.items():
        e2s[ensg] = symbol
        if symbol not in s2e:
            s2e[symbol] = ensg

    id_cache["ensembl_to_symbol"] = e2s
    id_cache["symbol_to_ensembl"] = s2e

    os.makedirs(os.path.dirname(cache_path) if os.path.dirname(cache_path) else ".", exist_ok=True)
    with open(cache_path, "wb") as f:
        pickle.dump(id_cache, f)

    print(f"  Updated cache saved to {cache_path}")
    print(f"  New totals: {len(e2s):,} ENSG->symbol, {len(s2e):,} symbol->ENSG")


def main():
    parser = argparse.ArgumentParser(description="Convert TSV network files to pickle cache format")
    parser.add_argument('cell_type', nargs='?', help='Cell type to process, or --all for all types')
    parser.add_argument('--all', action='store_true', help='Process all supported cell types')
    parser.add_argument('--input-dir', default='models/networks', help='Input directory containing TSV files')
    parser.add_argument('--output-dir', default='models/networks', help='Output directory for pickle files')
    parser.add_argument('--list-types', action='store_true', help='List supported cell types and exit')
    parser.add_argument('--enrich-gene-cache', action='store_true',
                        help='Bulk-resolve all GREmLN ENSG IDs to symbols and update cache/gene_id_cache.pkl')
    parser.add_argument('--gene-cache-path', default='cache/gene_id_cache.pkl',
                        help='Path to gene_id_cache.pkl (default: cache/gene_id_cache.pkl)')

    args = parser.parse_args()

    if args.list_types:
        print("Supported cell types:")
        for cell_type in SUPPORTED_CELL_TYPES:
            print(f"  - {cell_type}")
        return

    if args.enrich_gene_cache and not (args.all or args.cell_type):
        # Standalone: just update the gene ID cache without rebuilding PKLs
        update_gene_id_cache(args.output_dir, args.gene_cache_path)
        return

    if args.all or args.cell_type == '--all':
        # Process all cell types
        print(f"Processing all cell types from {args.input_dir} to {args.output_dir}")
        success_count = 0
        for cell_type in SUPPORTED_CELL_TYPES:
            if process_cell_type(cell_type, args.input_dir, args.output_dir):
                success_count += 1

        print(f"\n=== Summary ===")
        print(f"Successfully processed: {success_count}/{len(SUPPORTED_CELL_TYPES)} cell types")

        # Always enrich gene ID cache after rebuilding all cell types
        update_gene_id_cache(args.output_dir, args.gene_cache_path)

    elif args.cell_type:
        # Process single cell type
        if args.cell_type not in SUPPORTED_CELL_TYPES:
            print(f"ERROR: Unsupported cell type '{args.cell_type}'")
            print("Supported types:", ", ".join(SUPPORTED_CELL_TYPES))
            sys.exit(1)

        success = process_cell_type(args.cell_type, args.input_dir, args.output_dir)
        if not success:
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == '__main__':
    main()