#!/usr/bin/env python3
"""
Visualization Templates for RegNetAgents Results
================================================

Pre-made plotting functions for common visualizations of gene analysis results.
Designed to be easy for Claude Desktop's code execution tool to use.

These templates work with JSON results from RegNetAgents analyses and generate
publication-quality matplotlib figures.

Author: Jose A. Bird, PhD
License: MIT
"""

import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def load_gene_results(gene: str, results_dir: str = "results") -> Dict:
    """
    Load analysis results for a gene.

    Args:
        gene: Gene symbol (e.g., "TP53", "BRCA1")
        results_dir: Directory containing result files

    Returns:
        Dictionary with analysis results

    Example:
        >>> data = load_gene_results("TP53")
        >>> print(data['gene_info']['regulatory_role'])
    """
    gene_lower = gene.lower()

    # Try different file patterns
    patterns = [
        f"{results_dir}/{gene_lower}_analysis.json",
        f"{results_dir}/{gene_lower}_detailed_report.json",
        f"{results_dir}/cervical_{gene_lower}_analysis.json"
    ]

    for file_path in patterns:
        if Path(file_path).exists():
            with open(file_path, 'r') as f:
                return json.load(f)

    raise FileNotFoundError(f"No results found for {gene}. Tried: {patterns}")


def plot_regulator_rankings(
    gene: str,
    results_dir: str = "results",
    top_n: int = 10,
    save_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6)
) -> plt.Figure:
    """
    Create a bar chart of top regulators ranked by PageRank score.

    Args:
        gene: Gene symbol
        results_dir: Directory with results
        top_n: Number of top regulators to show
        save_path: Optional path to save figure
        figsize: Figure size (width, height)

    Returns:
        Matplotlib figure object

    Example:
        >>> fig = plot_regulator_rankings("TP53", top_n=5)
        >>> plt.show()
    """
    data = load_gene_results(gene, results_dir)

    # Extract regulator data
    if 'perturbation_analysis' in data:
        regulators = data['perturbation_analysis'].get('regulators', [])
    elif 'network_analysis' in data and 'regulators' in data['network_analysis']:
        regulators = data['network_analysis']['regulators']
    else:
        raise ValueError(f"No regulator data found for {gene}")

    # Sort by PageRank and get top N
    sorted_regs = sorted(regulators, key=lambda x: x.get('pagerank', 0), reverse=True)[:top_n]

    # Extract data for plotting
    symbols = [r['symbol'] for r in sorted_regs]
    pageranks = [r['pagerank'] for r in sorted_regs]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create horizontal bar chart
    y_pos = np.arange(len(symbols))
    bars = ax.barh(y_pos, pageranks, color='steelblue', alpha=0.8)

    # Customize appearance
    ax.set_yticks(y_pos)
    ax.set_yticklabels(symbols)
    ax.invert_yaxis()  # Top = highest PageRank
    ax.set_xlabel('PageRank Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Regulator Gene', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {len(symbols)} Regulators of {gene.upper()}',
                 fontsize=14, fontweight='bold')

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, pageranks)):
        ax.text(val + 0.01, i, f'{val:.3f}', va='center', fontsize=10)

    # Grid and styling
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig


def plot_pathway_enrichment(
    gene: str,
    results_dir: str = "results",
    top_n: int = 10,
    fdr_threshold: float = 0.05,
    save_path: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 6)
) -> plt.Figure:
    """
    Create a bar chart of enriched pathways with FDR significance.

    Args:
        gene: Gene symbol
        results_dir: Directory with results
        top_n: Number of top pathways to show
        fdr_threshold: FDR cutoff for significance
        save_path: Optional path to save figure
        figsize: Figure size

    Returns:
        Matplotlib figure object

    Example:
        >>> fig = plot_pathway_enrichment("TP53", top_n=5)
        >>> plt.show()
    """
    data = load_gene_results(gene, results_dir)

    # Extract pathway data
    if 'pathway_enrichment' not in data:
        raise ValueError(f"No pathway enrichment data found for {gene}")

    pathways = data['pathway_enrichment'].get('enriched_pathways', [])

    # Filter by FDR and get top N
    significant = [p for p in pathways if p.get('fdr', 1.0) < fdr_threshold]
    sorted_pathways = sorted(significant, key=lambda x: x.get('p_value', 1.0))[:top_n]

    if not sorted_pathways:
        print(f"No pathways with FDR < {fdr_threshold} found for {gene}")
        return None

    # Extract data
    names = [p['pathway_name'] for p in sorted_pathways]
    p_values = [p['p_value'] for p in sorted_pathways]
    fdrs = [p['fdr'] for p in sorted_pathways]

    # Transform p-values to -log10 scale
    neg_log_p = [-np.log10(p) if p > 0 else 10 for p in p_values]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create horizontal bar chart
    y_pos = np.arange(len(names))
    bars = ax.barh(y_pos, neg_log_p, color='coral', alpha=0.8)

    # Customize appearance
    ax.set_yticks(y_pos)
    # Truncate long pathway names
    truncated_names = [n[:50] + '...' if len(n) > 50 else n for n in names]
    ax.set_yticklabels(truncated_names, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('-log10(p-value)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Pathway', fontsize=12, fontweight='bold')
    ax.set_title(f'Top Enriched Pathways for {gene.upper()} (FDR < {fdr_threshold})',
                 fontsize=14, fontweight='bold')

    # Add FDR threshold line
    threshold_line = -np.log10(fdr_threshold)
    ax.axvline(threshold_line, color='red', linestyle='--', linewidth=2,
               label=f'FDR = {fdr_threshold}')
    ax.legend()

    # Add FDR labels
    for i, (bar, fdr) in enumerate(zip(bars, fdrs)):
        ax.text(bar.get_width() + 0.1, i, f'FDR={fdr:.3f}',
                va='center', fontsize=9)

    # Grid and styling
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig


def compare_genes_regulators(
    genes: List[str],
    results_dir: str = "results",
    save_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6)
) -> plt.Figure:
    """
    Compare the number of regulators and targets across multiple genes.

    Args:
        genes: List of gene symbols
        results_dir: Directory with results
        save_path: Optional path to save figure
        figsize: Figure size

    Returns:
        Matplotlib figure object

    Example:
        >>> fig = compare_genes_regulators(["TP53", "BRCA1", "MYC"])
        >>> plt.show()
    """
    # Load data for all genes
    genes_data = []
    for gene in genes:
        try:
            data = load_gene_results(gene, results_dir)
            gene_info = data.get('gene_info', {})
            genes_data.append({
                'gene': gene.upper(),
                'regulators': gene_info.get('num_regulators', 0),
                'targets': gene_info.get('num_targets', 0),
                'role': gene_info.get('regulatory_role', 'unknown')
            })
        except FileNotFoundError:
            print(f"Warning: No results found for {gene}, skipping")

    if not genes_data:
        raise ValueError("No gene data found")

    # Extract data for plotting
    gene_names = [d['gene'] for d in genes_data]
    regulators = [d['regulators'] for d in genes_data]
    targets = [d['targets'] for d in genes_data]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create grouped bar chart
    x = np.arange(len(gene_names))
    width = 0.35

    bars1 = ax.bar(x - width/2, regulators, width, label='Regulators (Upstream)',
                   color='steelblue', alpha=0.8)
    bars2 = ax.bar(x + width/2, targets, width, label='Targets (Downstream)',
                   color='coral', alpha=0.8)

    # Customize appearance
    ax.set_xlabel('Gene', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Connections', fontsize=12, fontweight='bold')
    ax.set_title('Regulatory Network Comparison', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(gene_names, fontsize=11)
    ax.legend()

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}', ha='center', va='bottom', fontsize=10)

    # Grid and styling
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig


def plot_network_overview(
    gene: str,
    results_dir: str = "results",
    save_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """
    Create a comprehensive overview figure for a gene.

    Shows: regulatory role, network statistics, and top regulators.

    Args:
        gene: Gene symbol
        results_dir: Directory with results
        save_path: Optional path to save figure
        figsize: Figure size

    Returns:
        Matplotlib figure object

    Example:
        >>> fig = plot_network_overview("TP53")
        >>> plt.show()
    """
    data = load_gene_results(gene, results_dir)
    gene_info = data.get('gene_info', {})

    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # 1. Network Statistics (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    stats = {
        'Regulators': gene_info.get('num_regulators', 0),
        'Targets': gene_info.get('num_targets', 0)
    }
    ax1.bar(stats.keys(), stats.values(), color=['steelblue', 'coral'], alpha=0.8)
    ax1.set_ylabel('Count', fontweight='bold')
    ax1.set_title('Network Statistics', fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    for i, (k, v) in enumerate(stats.items()):
        ax1.text(i, v, str(v), ha='center', va='bottom', fontweight='bold')

    # 2. Regulatory Role (top right)
    ax2 = fig.add_subplot(gs[0, 1])
    role = gene_info.get('regulatory_role', 'unknown').replace('_', ' ').title()
    ax2.text(0.5, 0.5, role, fontsize=20, fontweight='bold',
             ha='center', va='center', transform=ax2.transAxes,
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    ax2.set_title('Regulatory Role', fontweight='bold')
    ax2.axis('off')

    # 3. Top Regulators (bottom, spans both columns)
    ax3 = fig.add_subplot(gs[1, :])

    if 'perturbation_analysis' in data or 'network_analysis' in data:
        if 'perturbation_analysis' in data:
            regulators = data['perturbation_analysis'].get('regulators', [])
        else:
            regulators = data['network_analysis'].get('regulators', [])

        sorted_regs = sorted(regulators, key=lambda x: x.get('pagerank', 0), reverse=True)[:5]

        if sorted_regs:
            symbols = [r['symbol'] for r in sorted_regs]
            pageranks = [r['pagerank'] for r in sorted_regs]

            y_pos = np.arange(len(symbols))
            ax3.barh(y_pos, pageranks, color='green', alpha=0.6)
            ax3.set_yticks(y_pos)
            ax3.set_yticklabels(symbols)
            ax3.invert_yaxis()
            ax3.set_xlabel('PageRank Score', fontweight='bold')
            ax3.set_title('Top 5 Regulators', fontweight='bold')
            ax3.grid(axis='x', alpha=0.3)

            for i, val in enumerate(pageranks):
                ax3.text(val + 0.01, i, f'{val:.3f}', va='center')

    # Main title
    fig.suptitle(f'{gene.upper()} Network Overview', fontsize=16, fontweight='bold', y=0.98)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig


# Example usage functions
def example_single_gene_visualization():
    """Example: Visualize a single gene analysis."""
    gene = "TP53"

    # Create all visualizations
    fig1 = plot_regulator_rankings(gene, top_n=10)
    fig2 = plot_pathway_enrichment(gene, top_n=8)
    fig3 = plot_network_overview(gene)

    plt.show()


def example_multi_gene_comparison():
    """Example: Compare multiple genes."""
    genes = ["TP53", "BRCA1", "MYC", "KRAS"]

    fig = compare_genes_regulators(genes)
    plt.show()


if __name__ == "__main__":
    # Quick test
    print("RegNetAgents Visualization Templates")
    print("=" * 50)
    print("\nAvailable functions:")
    print("  - plot_regulator_rankings(gene)")
    print("  - plot_pathway_enrichment(gene)")
    print("  - compare_genes_regulators([genes])")
    print("  - plot_network_overview(gene)")
    print("\nUse these in Claude Desktop for quick visualizations!")
