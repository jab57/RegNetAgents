"""
Create Figure 1: RegNetAgents Multi-Agent Architecture
Workflow schematic showing directed acyclic graph of agent execution
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Create figure
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

# Four parallel agent boxes
agent_positions = [2, 4, 6, 8]
agent_labels = ['Cancer\nBiology\nAgent', 'Drug\nDiscovery\nAgent',
                'Clinical\nRelevance\nAgent', 'Systems\nBiology\nAgent']

# Arrows from header to each agent
for pos in agent_positions:
    create_arrow(ax, 5, 6.05, pos, 5.4, style='solid')

# Create agent boxes
for i, (pos, label) in enumerate(zip(agent_positions, agent_labels)):
    create_box(ax, pos, 4.8, 1.4, 1, label, color_domain, fontsize=9, fontweight='bold')
    # Arrows from agents to integration
    create_arrow(ax, pos, 4.3, 5, 3.6, style='solid')

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
plt.savefig('C:/Dev/RegNetAgents/biorxiv/figure1_architecture.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('C:/Dev/RegNetAgents/biorxiv/figure1_architecture.pdf',
            bbox_inches='tight', facecolor='white')

print("Figure 1 created successfully!")
print("Files saved:")
print("  - figure1_architecture.png (300 DPI)")
print("  - figure1_architecture.pdf")
