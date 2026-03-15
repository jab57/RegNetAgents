"""
Network loader for RegNetAgents.

Provides load_network() to load a TCGA ARACNe CSV network as a NetworkX
DiGraph with gene symbol nodes and edge attributes (moa, likelihood).
"""

import csv
import os

import networkx as nx

from regnetagents.tcga_registry import TCGA_NETWORK_REGISTRY


def load_network(name: str) -> nx.DiGraph:
    """
    Load a TCGA ARACNe network by cancer type name.

    Nodes are gene symbols (e.g., "TP53").  Edge attributes:
        - moa (float)       : Mode of action — +1 activation, -1 repression,
                              0 unknown/not determined
        - likelihood (float): ARACNe edge weight (0–1)

    Args:
        name: TCGA cancer type key — one of brca, coad, hnsc, luad, lusc,
              ov, prad, ucec

    Returns:
        nx.DiGraph with gene symbol nodes and moa / likelihood edge attributes

    Raises:
        KeyError: if name is not a recognised TCGA cancer type
        FileNotFoundError: if the source CSV has not been downloaded yet
    """
    if name not in TCGA_NETWORK_REGISTRY:
        raise KeyError(
            f"Unknown TCGA cancer type '{name}'. "
            f"Valid types: {sorted(TCGA_NETWORK_REGISTRY)}"
        )

    csv_path = TCGA_NETWORK_REGISTRY[name]["csv"]

    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"TCGA network CSV not found: {csv_path}\n"
            "Extract CSVs from the Bioconductor aracne.networks tarball:\n"
            "  python scripts/extract_tcga_networks.py \\\n"
            "      --tarball /tmp/aracne.networks.tar.gz \\\n"
            "      --output-dir models/networks/tcga\n"
            f"  python scripts/build_tcga_cache.py --cancer-type {name}\n"
            "See docs/DATA_SOURCES.md for full instructions."
        )

    G = nx.DiGraph()
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            regulator = row["Regulator"].strip().upper()
            target = row["Target"].strip().upper()
            moa = float(row.get("MoA") or 0)
            likelihood = float(row.get("Likelihood") or 0)
            G.add_edge(regulator, target, moa=moa, likelihood=likelihood)

    return G
