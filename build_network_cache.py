#!/usr/bin/env python3
"""
Script to convert TSV network files to pickle cache format for GREmLN.

This script processes the original .tsv network files mentioned in the README
and converts them into the optimized .pkl cache format used by the system.

Expected TSV format:
- ARACNe network output format with header
- Columns: regulator.values \t target.values \t mi.values \t scc.values \t count.values \t log.p.values
- Uses only first two columns (regulator.values, target.values)
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
"""

import os
import pickle
import argparse
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Set, Tuple
import sys

# Cell types supported by GREmLN
SUPPORTED_CELL_TYPES = [
    # Existing cell types (10)
    'cd14_monocytes',
    'cd16_monocytes',
    'cd20_b_cells',
    'cd4_t_cells',
    'cd8_t_cells',
    'erythrocytes',
    'nk_cells',
    'nkt_cells',
    'epithelial_cell',
    'monocyte-derived_dendritic_cells',

    # New priority cell types (5)
    'hepatocytes',
    'cardiomyocytes',
    'neurons',
    'fibroblasts',
    'endothelial_cells'
]

def load_tsv_network(tsv_file: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Set[str]]:
    """
    Load network from TSV file and build regulator-target mappings.

    Args:
        tsv_file: Path to TSV file with regulator-target pairs

    Returns:
        Tuple of (regulator_targets, target_regulators, all_genes)
    """
    regulator_targets = defaultdict(list)
    target_regulators = defaultdict(list)
    all_genes = set()

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

    # Convert defaultdicts to regular dicts
    regulator_targets = dict(regulator_targets)
    target_regulators = dict(target_regulators)

    print(f"Loaded {len(regulator_targets)} regulators, {len(target_regulators)} targets, {len(all_genes)} total genes")

    return regulator_targets, target_regulators, all_genes

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

def build_network_cache(tsv_file: str, output_file: str) -> None:
    """
    Convert TSV network file to pickle cache format.

    Args:
        tsv_file: Path to input TSV file
        output_file: Path to output pickle file
    """
    # Load network data
    regulator_targets, target_regulators, all_genes = load_tsv_network(tsv_file)

    # Calculate statistics
    stats = calculate_stats(regulator_targets, target_regulators, all_genes)

    # Build cache data structure
    cache_data = {
        'regulator_targets': regulator_targets,
        'target_regulators': target_regulators,
        'all_genes': sorted(list(all_genes)),  # Sort for consistency
        'num_edges': stats['num_edges'],
        'num_genes': stats['num_genes'],
        'num_regulons': stats['num_regulons'],
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

        print(f"Cache validation passed")
        return True

    except Exception as e:
        print(f"ERROR validating cache: {e}")
        return False

def process_cell_type(cell_type: str, input_dir: str, output_dir: str) -> bool:
    """
    Process a single cell type.

    Args:
        cell_type: Name of cell type
        input_dir: Directory containing TSV files
        output_dir: Directory for output pickle files

    Returns:
        True if successful, False otherwise
    """
    print(f"\n=== Processing {cell_type} ===")

    # Input TSV file
    tsv_file = os.path.join(input_dir, cell_type, "network.tsv")

    # Output pickle file
    output_file = os.path.join(output_dir, cell_type, "network_index.pkl")

    try:
        build_network_cache(tsv_file, output_file)

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

def main():
    parser = argparse.ArgumentParser(description="Convert TSV network files to pickle cache format")
    parser.add_argument('cell_type', nargs='?', help='Cell type to process, or --all for all types')
    parser.add_argument('--all', action='store_true', help='Process all supported cell types')
    parser.add_argument('--input-dir', default='models/networks', help='Input directory containing TSV files')
    parser.add_argument('--output-dir', default='models/networks', help='Output directory for pickle files')
    parser.add_argument('--list-types', action='store_true', help='List supported cell types and exit')

    args = parser.parse_args()

    if args.list_types:
        print("Supported cell types:")
        for cell_type in SUPPORTED_CELL_TYPES:
            print(f"  - {cell_type}")
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