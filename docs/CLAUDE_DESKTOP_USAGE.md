# Using Claude Desktop with RegNetAgents

## Overview

This guide shows how to use Claude Desktop's code execution capabilities with your RegNetAgents analysis results. This lets you explore results interactively, create custom visualizations, and compare genes using natural language.

**Key Benefit**: After running your Ollama-powered analyses (103s for 7 genes), use Claude Desktop to explore results without writing Python code manually.

---

## Prerequisites

1. **Claude Desktop** installed (with Pro subscription for code execution)
2. **MCP Server** configured (see Setup section below)
3. **Analysis Results** available in `results/` directory

**Cost**: $0 additional (uses existing Ollama + Claude Desktop Pro subscription)

---

## Setup

### 1. Configure MCP Server in Claude Desktop

Edit your Claude Desktop MCP configuration file:

**Location**:
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`
- Mac: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

**Add RegNetAgents server**:
```json
{
  "mcpServers": {
    "regnetagents": {
      "command": "python",
      "args": ["C:/Dev/RegNetAgents/regnetagents_langgraph_mcp_server.py"],
      "env": {}
    }
  }
}
```

**Note**: Adjust the path to match your installation directory.

### 2. Restart Claude Desktop

After editing the config, restart Claude Desktop to load the MCP server.

### 3. Verify Connection

In Claude Desktop, type:
```
List available analysis results
```

Claude should respond with genes that have saved results.

---

## Available MCP Tools

Claude Desktop can now use these tools (11 total):

### 1. `validate_gene`
Quick gene name check (<100ms). Returns basic stats if found, or fuzzy-matched suggestions for misspelled names. Use before full analysis to catch typos.

**Example**:
```
Is TP53 in the epithelial cell network?
```

### 2. `query_network`
Instant network queries from pre-computed data (<50ms). Supports four query types: top_regulators, top_targets, gene_neighbors, and network_stats.

**Example**:
```
What are the top regulators in epithelial cells?
How many targets does TP53 have?
What regulates MYC?
How large is the epithelial cell network?
```

### 3. `comprehensive_gene_analysis`
Run full analysis for a gene (triggers Ollama workflow).

### 4. `load_gene_results`
Load previously saved analysis results without re-running analysis.

**Example**:
```
Load the TP53 analysis results
```

### 5. `list_available_results`
See what gene analyses are available.

**Example**:
```
What gene analyses do I have?
```

### 6. `multi_gene_analysis`
Analyze multiple genes in parallel.

### 7-11. Other tools
- `pathway_focused_analysis`
- `cross_cell_comparison`
- `workflow_status`
- `workflow_insights`
- `create_analysis_report`

---

## Example Workflows

### Workflow 1: Quick Data Exploration

**Step 1**: Run analysis (using Ollama - free)
```bash
python scripts/analyze_cervical_genes.py
```
→ Takes 103 seconds, saves JSON results

**Step 2**: Explore in Claude Desktop
```
You: Load the TP53 results and show me the top 5 regulators

Claude: [Loads results/cervical_tp53_analysis.json]
        [Shows: WWTR1 (0.473), GPSM3 (0.184), ...]
```

### Workflow 2: Custom Visualization

```
You: Create a bar chart of PageRank scores for TP53's regulators

Claude: [Executes Python code]
        import matplotlib.pyplot as plt
        # ... generates plot
        [Shows visualization]
```

### Workflow 3: Compare Multiple Genes

```
You: Compare the number of regulators across all 7 cervical cancer genes

Claude: [Loads all 7 result files]
        [Creates grouped bar chart]
        - TP53: 7 regulators
        - BRCA1: 23 regulators
        - EGFR: 37 regulators
        ...
```

### Workflow 4: Pathway Analysis

```
You: Which pathways are enriched in both TP53 and MYC?

Claude: [Loads both gene results]
        [Compares pathway lists]

        Overlapping pathways:
        1. Signaling by Hippo (both genes)
        2. NOTCH signaling (both genes)
        ...
```

### Workflow 5: Statistical Analysis

```
You: What's the correlation between regulator count and PageRank
     for cervical cancer genes?

Claude: [Loads 7 gene results]
        [Runs scipy.stats.pearsonr]

        Correlation: r = 0.45, p = 0.03
        [Shows scatter plot with regression line]
```

---

## Using Visualization Templates

The `scripts/visualization_templates.py` file provides pre-made plotting functions.

**Example in Claude Desktop**:
```
You: Use the visualization templates to plot regulator rankings for BRCA1

Claude: [Executes code]
        from scripts.visualization_templates import plot_regulator_rankings
        fig = plot_regulator_rankings("BRCA1", top_n=10)
        plt.show()

        [Shows publication-quality bar chart]
```

**Available template functions**:
- `plot_regulator_rankings(gene, top_n=10)`
- `plot_pathway_enrichment(gene, top_n=10)`
- `compare_genes_regulators([genes])`
- `plot_network_overview(gene)`

---

## Common Questions

### Q: What's the difference between Claude Desktop and the Python scripts?

**Python scripts** (using Ollama):
- Run the full multi-agent analysis
- 103 seconds for 7 genes
- Generates JSON results
- Free (uses local Ollama)

**Claude Desktop**:
- Reads existing JSON results
- Interactive exploration (~5 seconds)
- Natural language queries
- Uses Claude Pro subscription (already have it)

**They work together**: Python does the heavy lifting, Claude helps you explore.

---

### Q: Do I need to change my bioRxiv manuscript?

**NO**. Your manuscript documents:
- "Analysis uses Ollama llama3.1:8b" ✓ (still true)
- "103 seconds for 7 genes" ✓ (still true)
- All timing and methodology ✓ (unchanged)

Claude Desktop is just a post-analysis exploration tool.

---

### Q: Will this add API costs?

**NO**. This uses:
- Ollama for analysis ($0)
- Claude Desktop Pro subscription (you likely already have this)
- Claude's code execution is included in Pro subscription

Total additional cost: $0

---

### Q: Can I export figures for publication?

**YES**. Tell Claude to save the figure:

```
You: Save that plot as a high-resolution PDF

Claude: [Executes code]
        plt.savefig('tp53_regulators.pdf', dpi=300, bbox_inches='tight')

        Figure saved to tp53_regulators.pdf
```

---

## Example Natural Language Queries

### Data Exploration
- "What genes have the most regulators?"
- "Show me all hub regulators in my results"
- "Which gene has the highest PageRank regulator?"

### Visualization
- "Create a heatmap of PageRank scores across genes"
- "Plot pathway enrichment for MYC"
- "Show me a network diagram for TP53"

### Comparison
- "Compare BRCA1 and BRCA2 regulatory profiles"
- "Which pathways are unique to TP53?"
- "Show regulator overlap between all genes"

### Analysis
- "Calculate average PageRank by regulatory role"
- "Find genes with similar network patterns"
- "Identify potential hub regulators for drug targeting"

### Export
- "Save the last plot as publication-quality PDF"
- "Create a table of top regulators for all genes"
- "Export pathway data to CSV"

---

## Troubleshooting

### Issue: Claude Desktop can't connect to MCP server

**Solution**:
1. Check `claude_desktop_config.json` path is correct
2. Restart Claude Desktop
3. Check server logs in Claude Desktop settings

### Issue: "No results found for gene X"

**Solution**:
1. Run analysis first: `python scripts/analyze_cervical_genes.py`
2. Check `results/` directory has JSON files
3. Verify gene symbol spelling (case-insensitive)

### Issue: Visualization shows empty plot

**Solution**:
1. Load results first: "Load TP53 results"
2. Check if gene has regulator/pathway data
3. Try a different visualization type

---

## Best Practices

### 1. Run Analysis First
Always run Python analysis before using Claude Desktop:
```bash
python scripts/analyze_cervical_genes.py
```

### 2. Start with Simple Queries
Begin with: "List available results"
Then: "Load TP53 results"
Finally: "Show me the top regulators"

### 3. Iterate on Visualizations
```
You: Plot regulator rankings for TP53
     [Claude shows plot]
You: Make that horizontal and show only top 5
     [Claude adjusts and re-plots]
You: Perfect, save as PDF
     [Claude exports figure]
```

### 4. Combine Multiple Operations
```
You: For each cervical cancer gene:
     1. Load the results
     2. Extract the top regulator
     3. Create a summary table
     4. Save as CSV
```

Claude can chain these operations automatically.

---

## Advanced Usage

### Custom Analysis Scripts

Create your own analysis in Claude Desktop:

```
You: Write a Python script that:
     1. Loads all cervical gene results
     2. Calculates average PageRank by regulatory role
     3. Performs ANOVA test
     4. Creates publication figure

Claude: [Writes and executes multi-step script]
        [Shows statistical results and figure]
```

### Integration with Other Tools

```
You: Load TP53 results and export regulator list
     in a format compatible with Cytoscape

Claude: [Converts data to Cytoscape format]
        [Saves as .txt file]
```

### Batch Processing

```
You: For all genes in results:
     - Create regulator ranking plot
     - Save as [gene]_regulators.pdf
     - Generate summary statistics CSV

Claude: [Processes all genes automatically]
        Generated 7 plots and summary.csv
```

---

## Summary

**What You Get**:
- Interactive data exploration (natural language)
- On-demand visualization (no manual coding)
- Cross-gene comparisons (automatic)
- Publication-quality figures (300 DPI)

**What It Costs**:
- $0 (uses existing Ollama + Claude Desktop Pro)

**What Stays the Same**:
- Analysis workflow (still 103 seconds)
- bioRxiv documentation (no changes needed)
- GitHub repos (all working)

**When to Use**:
- After running Python analysis
- For exploring results interactively
- When you need custom visualizations
- For answering ad-hoc questions about data

---

## Support

For questions or issues:
1. Check this documentation
2. Review visualization templates: `scripts/visualization_templates.py`
3. See MCP server code: `regnetagents_langgraph_mcp_server.py`

**GitHub**: https://github.com/jab57/RegNetAgents

---

*Last updated: December 2025*
*Author: Jose A. Bird, PhD*
