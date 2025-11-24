"""
Generate figures for bioRxiv preprint from analysis results
"""

import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

# Set style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'

# Load data from results folder
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
results_dir = os.path.join(project_root, 'results')
output_dir = os.path.join(project_root, 'biorxiv')  # Output directory for all figures and tables

with open(os.path.join(results_dir, 'biomarker_results.json'), 'r') as f:
    biomarker_data = json.load(f)

# Use the standard centrality version with PageRank, degree, etc.
with open(os.path.join(results_dir, 'tp53_perturbation_standard_centrality.json'), 'r') as f:
    tp53_data = json.load(f)

#############################################################################
# FIGURE 1: Multi-Agent Architecture
#############################################################################

def create_figure1():
    """Create Figure 1: RegNetAgents Multi-Agent Architecture"""

    fig, ax = plt.subplots(1, 1, figsize=(12, 14))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 16)
    ax.axis('off')

    # Define colors
    color_input = '#E8F4F8'      # Light blue for input
    color_core = '#B8E6B8'        # Light green for core agents
    color_domain = '#FFE6B8'      # Light orange for domain agents
    color_output = '#E8D4F8'      # Light purple for output
    color_data = '#F0F0F0'        # Light gray for data sources

    # Helper function to create boxes
    def create_box(ax, x, y, width, height, text, color, fontsize=10, fontweight='normal'):
        box = FancyBboxPatch((x-width/2, y-height/2), width, height,
                             boxstyle="round,pad=0.1",
                             edgecolor='black', facecolor=color, linewidth=2)
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
                fontweight=fontweight, wrap=True)
        return box

    # Helper function to create arrows
    def create_arrow(ax, x1, y1, x2, y2, style='solid', color='black', linewidth=2):
        arrow = FancyArrowPatch((x1, y1), (x2, y2),
                               arrowstyle='->', mutation_scale=20,
                               linestyle=style, color=color, linewidth=linewidth,
                               connectionstyle="arc3,rad=0")
        ax.add_patch(arrow)
        return arrow

    # Title
    ax.text(5, 15.5, 'RegNetAgents Multi-Agent Architecture',
            ha='center', va='center', fontsize=16, fontweight='bold')

    # 1. User Input (Top)
    create_box(ax, 5, 14.5, 3, 0.6, 'User Query via MCP Server\n(Natural Language)',
               color_input, fontsize=10, fontweight='bold')

    # Arrow down
    create_arrow(ax, 5, 14.2, 5, 13.6)

    # 2. Initialization/Validation Agent
    create_box(ax, 5, 13, 3.5, 0.8, 'Initialization & Validation Agent\nGene Symbol → Ensembl ID',
               color_core, fontsize=10, fontweight='bold')

    # Arrow down
    create_arrow(ax, 5, 12.6, 5, 11.9)

    # 3. Network Modeling Agent
    create_box(ax, 5, 11.3, 4, 1, 'Network Modeling Agent\n• Identify Regulators (upstream)\n• Identify Targets (downstream)\n• Classify Regulatory Role',
               color_core, fontsize=9, fontweight='bold')

    # Data source annotation (left side)
    create_box(ax, 1.5, 11.3, 2, 0.8, 'ARACNe Networks\n(Pre-computed\nNetworkX cache)',
               color_data, fontsize=8)
    create_arrow(ax, 2.5, 11.3, 3, 11.3, style='solid', color='gray')

    # Arrow down to decision point
    create_arrow(ax, 5, 10.8, 5, 10.2)

    # 4. Intelligent Routing (Decision Diamond)
    # Create diamond shape
    diamond_x = [5, 5.7, 5, 4.3, 5]
    diamond_y = [10, 9.5, 9, 9.5, 10]
    diamond = plt.Polygon(list(zip(diamond_x, diamond_y)),
                         edgecolor='black', facecolor='#FFFACD', linewidth=2)
    ax.add_patch(diamond)
    ax.text(5, 9.5, '>5\nregulators?', ha='center', va='center',
            fontsize=9, fontweight='bold')

    # Conditional arrow to perturbation (dashed, right)
    ax.text(6.2, 9.5, 'Yes', fontsize=8, style='italic')
    create_arrow(ax, 5.7, 9.5, 7.2, 9.5, style='dashed', color='red', linewidth=2)

    # 5. Perturbation Analysis (right side, conditional)
    create_box(ax, 8.2, 9.5, 1.6, 1.2, 'Perturbation\nAnalysis\n• PageRank\n• Centrality\n• Target\nRanking',
               color_core, fontsize=8, fontweight='bold')

    # Arrow from perturbation back to main flow
    create_arrow(ax, 8.2, 8.9, 8.2, 7.5, style='dashed', color='red')
    create_arrow(ax, 8.2, 7.5, 5.5, 7.5, style='dashed', color='red')

    # Main flow continues down (No path)
    ax.text(4.2, 8.6, 'No', fontsize=8, style='italic')
    create_arrow(ax, 5, 9, 5, 8.3)

    # 6. Pathway Enrichment Agent
    create_box(ax, 5, 7.7, 3.5, 0.8, 'Pathway Enrichment Agent\nReactome API Query',
               color_core, fontsize=10, fontweight='bold')

    # Data source annotation (right side for Reactome)
    create_box(ax, 8.5, 7.7, 1.6, 0.6, 'Reactome\nDatabase',
               color_data, fontsize=8)
    create_arrow(ax, 7.75, 7.7, 8.2, 7.7, style='solid', color='gray')

    # Arrow down to parallel agents
    create_arrow(ax, 5, 7.3, 5, 6.7)

    # 7. Parallel Domain Analysis (4 agents)
    # Header box
    create_box(ax, 5, 6.3, 4.5, 0.5, 'Parallel Domain Analysis (LLM-Powered)',
               '#FFD700', fontsize=10, fontweight='bold')

    # Gene annotation data source (feeds into all domain agents via header)
    create_box(ax, 1.15, 5.5, 1.0, 0.6, 'NCBI+UniProt\n(15,347)',
               color_data, fontsize=6.5)
    # Arrow to domain analysis header (shows all agents get annotations)
    create_arrow(ax, 1.65, 5.5, 2.5, 6.15, style='solid', color='gray')

    # Four parallel agent boxes (smaller and better spaced)
    agent_positions = [2.8, 4.5, 6.2, 7.9]
    agent_labels = ['Cancer\nBiology', 'Drug\nDiscovery',
                    'Clinical\nRelevance', 'Systems\nBiology']

    # Arrows from header to each agent
    for pos in agent_positions:
        create_arrow(ax, 5, 6.05, pos, 5.3, style='solid')

    # Create agent boxes (smaller width and height)
    for i, (pos, label) in enumerate(zip(agent_positions, agent_labels)):
        create_box(ax, pos, 4.8, 1.15, 0.8, label, color_domain, fontsize=8, fontweight='bold')
        # Arrows from agents to integration
        create_arrow(ax, pos, 4.4, 5, 3.6, style='solid')

    # 8. Integration Agent
    create_box(ax, 5, 3, 4, 1, 'Integration & Reporting Agent\n• Synthesize multi-domain insights\n• Generate comprehensive JSON report\n• Rank therapeutic targets',
               color_core, fontsize=9, fontweight='bold')

    # Arrow down
    create_arrow(ax, 5, 2.5, 5, 1.9)

    # 9. Output
    create_box(ax, 5, 1.4, 3.5, 0.7, 'Comprehensive Analysis Report\n(Structured JSON)',
               color_output, fontsize=10, fontweight='bold')

    # Legend
    legend_y = 0.4
    ax.text(1, legend_y, 'Legend:', fontsize=9, fontweight='bold')
    # Solid arrow
    create_arrow(ax, 2, legend_y, 2.7, legend_y, style='solid', linewidth=1.5)
    ax.text(3.2, legend_y, 'Required step', fontsize=8, va='center')
    # Dashed arrow
    create_arrow(ax, 5, legend_y, 5.7, legend_y, style='dashed', color='red', linewidth=1.5)
    ax.text(6.3, legend_y, 'Conditional execution', fontsize=8, va='center')

    # Add workflow annotations
    ax.text(0.3, 13, 'Input', fontsize=9, fontweight='bold', rotation=90, va='center')
    ax.text(0.3, 9.5, 'Analysis', fontsize=9, fontweight='bold', rotation=90, va='center')
    ax.text(0.3, 4.8, 'Domain\nInsights', fontsize=9, fontweight='bold', rotation=90, va='center')
    ax.text(0.3, 2, 'Output', fontsize=9, fontweight='bold', rotation=90, va='center')

    # Add execution time annotation
    ax.text(9.5, 1.4, 'Execution Time:\n~0.6s (rule-based)\n~15s (LLM-powered)',
            fontsize=7, ha='right', va='center', style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()

    # Save figure
    plt.savefig(os.path.join(output_dir, 'figure1_architecture.png'),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, 'figure1_architecture.pdf'),
                bbox_inches='tight', facecolor='white')
    print(f"[OK] Figure 1 saved: {output_dir}/figure1_architecture.png and .pdf")
    plt.close()

#############################################################################
# FIGURE 2: Biomarker Panel Analysis
#############################################################################

def create_figure2():
    """Create Figure 2: Colorectal Cancer Biomarker Panel Analysis"""

    fig = plt.figure(figsize=(12, 8))

    # Extract data
    genes = [r['gene'] for r in biomarker_data['biomarker_results']]
    num_targets = [r['num_targets'] for r in biomarker_data['biomarker_results']]
    num_regulators = [r['num_regulators'] for r in biomarker_data['biomarker_results']]
    biomarker_types = [r['biomarker_type'] for r in biomarker_data['biomarker_results']]
    regulatory_roles = [r['regulatory_role'] for r in biomarker_data['biomarker_results']]

    # A) Bar chart: Regulators and Targets
    ax1 = plt.subplot(2, 2, 1)
    x = np.arange(len(genes))
    width = 0.35

    bars1 = ax1.bar(x - width/2, num_regulators, width, label='Regulators (upstream)',
                    color='#3498db', alpha=0.8)
    bars2 = ax1.bar(x + width/2, num_targets, width, label='Targets (downstream)',
                    color='#e74c3c', alpha=0.8)

    ax1.set_xlabel('Gene', fontweight='bold')
    ax1.set_ylabel('Number of Connections', fontweight='bold')
    ax1.set_title('A) Regulatory Network Architecture', fontweight='bold', fontsize=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(genes, fontweight='bold')
    ax1.legend(frameon=True, loc='upper left')
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_ylim(0, max(max(num_targets), max(num_regulators)) * 1.15)

    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        if height > 0:
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom', fontsize=8)

    # B) Regulatory role classification
    ax2 = plt.subplot(2, 2, 2)
    role_map = {'hub_regulator': 'Hub Regulator', 'terminal_target': 'Terminal Target',
                'target': 'Target', 'master_regulator': 'Master Regulator'}
    role_colors = {'hub_regulator': '#e74c3c', 'terminal_target': '#f39c12',
                   'target': '#3498db', 'master_regulator': '#9b59b6'}

    for i, (gene, role) in enumerate(zip(genes, regulatory_roles)):
        color = role_colors.get(role, '#95a5a6')
        ax2.barh(i, 1, color=color, alpha=0.8)
        ax2.text(0.5, i, role_map.get(role, role.replace('_', ' ').title()),
                ha='center', va='center', fontweight='bold', fontsize=10, color='white')

    ax2.set_yticks(range(len(genes)))
    ax2.set_yticklabels(genes, fontweight='bold')
    ax2.set_xlim(0, 1)
    ax2.set_xticks([])
    ax2.set_title('B) Network Regulatory Role', fontweight='bold', fontsize=12)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)

    # C) Biomarker type classification
    ax3 = plt.subplot(2, 2, 3)
    biomarker_type_map = {'diagnostic': 'Diagnostic', 'prognostic': 'Prognostic',
                          'predictive': 'Predictive'}
    type_colors = {'diagnostic': '#3498db', 'prognostic': '#e74c3c',
                   'predictive': '#f39c12'}

    for i, (gene, btype) in enumerate(zip(genes, biomarker_types)):
        color = type_colors[btype]
        ax3.barh(i, 1, color=color, alpha=0.8)
        ax3.text(0.5, i, biomarker_type_map[btype], ha='center', va='center',
                fontweight='bold', fontsize=10, color='white')

    ax3.set_yticks(range(len(genes)))
    ax3.set_yticklabels(genes, fontweight='bold')
    ax3.set_xlim(0, 1)
    ax3.set_xticks([])
    ax3.set_title('C) Biomarker Classification', fontweight='bold', fontsize=12)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)

    # D) Pathway enrichment summary
    ax4 = plt.subplot(2, 2, 4)

    # Extract real pathway counts from results
    pathway_counts = [r['num_pathways'] for r in biomarker_data['biomarker_results']]

    bars = ax4.barh(genes, pathway_counts, color='#9b59b6', alpha=0.8)
    ax4.set_xlabel('Number of Significant Pathways', fontweight='bold')
    ax4.set_title('D) Pathway Enrichment (FDR<0.05)', fontweight='bold', fontsize=12)
    ax4.set_xlim(0, max(pathway_counts) * 1.15)
    ax4.grid(axis='x', alpha=0.3)

    # Add value labels
    for i, (bar, count) in enumerate(zip(bars, pathway_counts)):
        ax4.text(count + 0.5, i, f'{count}', va='center', fontsize=9)

    plt.tight_layout()
    # Save to biorxiv folder
    plt.savefig(os.path.join(output_dir, 'figure2_biomarker_panel.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'figure2_biomarker_panel.pdf'), bbox_inches='tight')
    print(f"[OK] Figure 2 saved: {output_dir}/figure2_biomarker_panel.png and .pdf")
    plt.close()

#############################################################################
# FIGURE 3: TP53 Perturbation Analysis
#############################################################################

def create_figure3():
    """Create Figure 3: TP53 Perturbation Analysis and Therapeutic Target Ranking"""

    fig = plt.figure(figsize=(14, 6))

    # Extract perturbation data (updated format - removed redundant metrics)
    perturbations = tp53_data['perturbation_results']
    regulators = [p['regulator'] for p in perturbations]
    # Use PageRank as primary ranking metric (standard centrality)
    pagerank_scores = [p['centrality_metrics']['pagerank'] for p in perturbations]
    downstream_targets = [p['regulator_downstream_targets'] for p in perturbations]

    # Sort by PageRank (descending) - primary ranking metric
    sorted_data = sorted(zip(regulators, pagerank_scores, downstream_targets),
                        key=lambda x: x[1], reverse=True)
    regulators_sorted = [x[0] for x in sorted_data]
    pagerank_sorted = [x[1] for x in sorted_data]
    targets_sorted = [x[2] for x in sorted_data]

    # A) Network diagram (simplified)
    ax1 = plt.subplot(1, 3, 1)

    # TP53 in center
    tp53_x, tp53_y = 0.5, 0.5
    circle_tp53 = plt.Circle((tp53_x, tp53_y), 0.08, color='#e74c3c', alpha=0.8, zorder=3)
    ax1.add_patch(circle_tp53)
    ax1.text(tp53_x, tp53_y, 'TP53', ha='center', va='center',
            fontweight='bold', fontsize=11, color='white', zorder=4)

    # Regulators around TP53
    n_regs = len(regulators_sorted)
    angles = np.linspace(0, 2*np.pi, n_regs, endpoint=False)
    radius = 0.35

    for i, (angle, reg, targets) in enumerate(zip(angles, regulators_sorted, targets_sorted)):
        x = tp53_x + radius * np.cos(angle)
        y = tp53_y + radius * np.sin(angle)

        # Scale circle size by downstream targets (normalized)
        size = 0.04 + 0.04 * (targets / max(targets_sorted))

        # Highlight top 3
        color = '#2ecc71' if i < 3 else '#3498db'
        alpha = 0.9 if i < 3 else 0.6

        circle = plt.Circle((x, y), size, color=color, alpha=alpha, zorder=2)
        ax1.add_patch(circle)

        # Arrow from regulator to TP53
        ax1.annotate('', xy=(tp53_x, tp53_y), xytext=(x, y),
                    arrowprops=dict(arrowstyle='->', lw=2 if i < 3 else 1,
                                  color='#2c3e50', alpha=0.6))

        # Label
        ax1.text(x, y, reg, ha='center', va='center',
                fontweight='bold' if i < 3 else 'normal',
                fontsize=9 if i < 3 else 7, color='white', zorder=3)

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('A) TP53 Regulatory Network\n(7 upstream regulators)',
                 fontweight='bold', fontsize=11)

    # Add legend
    legend_elements = [
        mpatches.Circle((0, 0), 1, fc='#2ecc71', alpha=0.9, label='Top 3 targets'),
        mpatches.Circle((0, 0), 1, fc='#3498db', alpha=0.6, label='Other regulators')
    ]
    ax1.legend(handles=legend_elements, loc='upper left', fontsize=8, frameon=True)

    # B) PageRank scores ranked (standard centrality metric)
    ax2 = plt.subplot(1, 3, 2)

    colors = ['#2ecc71' if i < 3 else '#3498db' for i in range(len(regulators_sorted))]
    bars = ax2.barh(range(len(regulators_sorted)), pagerank_sorted, color=colors, alpha=0.8)

    ax2.set_yticks(range(len(regulators_sorted)))
    ax2.set_yticklabels(regulators_sorted, fontweight='bold')
    ax2.set_xlabel('PageRank Score', fontweight='bold')
    ax2.set_title('B) Therapeutic Target Ranking\n(by PageRank - standard centrality)',
                 fontweight='bold', fontsize=11)
    ax2.set_xlim(0, max(pagerank_sorted) * 1.15)
    ax2.grid(axis='x', alpha=0.3)
    ax2.invert_yaxis()

    # Add value labels
    for i, (bar, score) in enumerate(zip(bars, pagerank_sorted)):
        ax2.text(score + 0.01, i, f'{score:.3f}', va='center', fontsize=8)

    # C) Alternative ranking by degree centrality (downstream targets)
    ax3 = plt.subplot(1, 3, 3)

    # Sort by downstream targets for degree centrality ranking
    sorted_by_degree = sorted(zip(regulators, downstream_targets, pagerank_scores),
                             key=lambda x: x[1], reverse=True)
    regulators_by_degree = [x[0] for x in sorted_by_degree]
    targets_by_degree = [x[1] for x in sorted_by_degree]
    pagerank_by_degree = [x[2] for x in sorted_by_degree]

    # Color by whether in top 3 by PageRank
    colors_degree = []
    for reg in regulators_by_degree:
        if reg in regulators_sorted[:3]:
            colors_degree.append('#2ecc71')  # Top 3 by PageRank
        else:
            colors_degree.append('#3498db')

    bars = ax3.barh(range(len(regulators_by_degree)), targets_by_degree,
                    color=colors_degree, alpha=0.8)

    ax3.set_yticks(range(len(regulators_by_degree)))
    ax3.set_yticklabels(regulators_by_degree, fontweight='bold')
    ax3.set_xlabel('Downstream Targets (degree centrality)', fontweight='bold')
    ax3.set_title('C) Alternative Ranking by Degree\n(off-target effect estimate)',
                 fontweight='bold', fontsize=11)
    ax3.set_xlim(0, max(targets_by_degree) * 1.15)
    ax3.grid(axis='x', alpha=0.3)
    ax3.invert_yaxis()

    # Add value labels
    for i, (bar, targets) in enumerate(zip(bars, targets_by_degree)):
        ax3.text(targets + 10, i, f'{targets}', va='center', fontsize=8)

    # Add validation notes (bottom right corner for clarity)
    validation_text = "Green = Top 3 by PageRank\n[V] WWTR1/YAP1: Validated\n[?] RBPMS: Novel"
    ax3.text(0.98, 0.05, validation_text, transform=ax3.transAxes,
            ha='right', va='bottom', fontsize=8, style='italic',
            bbox=dict(boxstyle='round', facecolor='#f39c12', alpha=0.3))

    plt.tight_layout()
    # Save to biorxiv folder
    plt.savefig(os.path.join(output_dir, 'figure3_tp53_perturbation.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'figure3_tp53_perturbation.pdf'), bbox_inches='tight')
    print(f"[OK] Figure 3 saved: {output_dir}/figure3_tp53_perturbation.png and .pdf")
    plt.close()

#############################################################################
# TABLE 2: Biomarker Results (CSV for inclusion)
#############################################################################

def create_table2():
    """Create Table 2 as CSV"""

    df = pd.DataFrame({
        'Gene': [r['gene'] for r in biomarker_data['biomarker_results']],
        'Regulatory Role': [r['regulatory_role'] for r in biomarker_data['biomarker_results']],
        'Targets': [r['num_targets'] for r in biomarker_data['biomarker_results']],
        'Regulators': [r['num_regulators'] for r in biomarker_data['biomarker_results']],
        'Biomarker Type': [r['biomarker_type'] for r in biomarker_data['biomarker_results']]
    })

    # Save to biorxiv folder
    df.to_csv(os.path.join(output_dir, 'table2_biomarker_results.csv'), index=False)
    print(f"[OK] Table 2 saved: {output_dir}/table2_biomarker_results.csv")

    # Also create formatted text version
    with open(os.path.join(output_dir, 'table2_biomarker_results.txt'), 'w') as f:
        f.write(df.to_string(index=False))
    print(f"[OK] Table 2 saved: {output_dir}/table2_biomarker_results.txt")

#############################################################################
# TABLE 3: TP53 Perturbation Results (CSV for inclusion)
#############################################################################

def create_table3():
    """Create Table 3 as CSV - Updated to remove redundant metrics"""

    perturbations = tp53_data['perturbation_results']

    # Validation status based on literature review
    # WWTR1, CHD4, YAP1 are experimentally validated TP53 regulators (Hippo pathway)
    # Others are novel hypotheses
    validation_map = {
        'WWTR1': '✓ Validated',
        'CHD4': '✓ Validated',
        'YAP1': '✓ Validated',
        'RBPMS': 'Novel',
        'PRRX2': 'Novel',
        'THRA': 'Novel',
        'IKZF2': 'Novel'
    }

    df = pd.DataFrame({
        'Rank': range(1, len(perturbations) + 1),
        'Regulator': [p['regulator'] for p in perturbations],
        'PageRank': [p['centrality_metrics']['pagerank'] for p in perturbations],
        'Out-Degree Centrality': [p['centrality_metrics']['out_degree_centrality'] for p in perturbations],
        'Downstream Targets': [p['regulator_downstream_targets'] for p in perturbations],
        'Validation Status': [validation_map.get(p['regulator'], 'Novel') for p in perturbations]
    })

    # Sort by PageRank (primary ranking metric)
    df = df.sort_values('PageRank', ascending=False).reset_index(drop=True)
    df['Rank'] = range(1, len(df) + 1)

    # Save to biorxiv folder
    df.to_csv(os.path.join(output_dir, 'table3_tp53_perturbation.csv'), index=False)
    print(f"[OK] Table 3 saved: {output_dir}/table3_tp53_perturbation.csv")

    # Also create formatted text version (UTF-8 for checkmark symbols)
    with open(os.path.join(output_dir, 'table3_tp53_perturbation.txt'), 'w', encoding='utf-8') as f:
        f.write(df.to_string(index=False))
    print(f"[OK] Table 3 saved: {output_dir}/table3_tp53_perturbation.txt")

#############################################################################
# Main execution
#############################################################################

if __name__ == "__main__":
    print("\n" + "="*60)
    print("Generating figures and tables for bioRxiv preprint")
    print("="*60 + "\n")

    print("Creating Figure 1: Multi-Agent Architecture...")
    create_figure1()

    print("\nCreating Figure 2: Biomarker Panel Analysis...")
    create_figure2()

    print("\nCreating Figure 3: TP53 Perturbation Analysis...")
    create_figure3()

    print("\nCreating Table 2: Biomarker Results...")
    create_table2()

    print("\nCreating Table 3: TP53 Perturbation Results...")
    create_table3()

    print("\n" + "="*60)
    print("[SUCCESS] All figures and tables generated successfully!")
    print("="*60)
    print("\nGenerated files in biorxiv/ folder:")
    print("  - figure1_architecture.png (and .pdf)")
    print("  - figure2_biomarker_panel.png (and .pdf)")
    print("  - figure3_tp53_perturbation.png (and .pdf)")
    print("  - table2_biomarker_results.csv (and .txt)")
    print("  - table3_tp53_perturbation.csv (and .txt)")
    print("\nReady for your bioRxiv submission!")
