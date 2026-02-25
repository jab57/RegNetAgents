# RegNetAgents MCP Server Setup Guide

This guide helps you set up RegNetAgents with Claude Desktop using the Model Context Protocol (MCP).

---

## What is the MCP Server?

The RegNetAgents Agents MCP (Model Context Protocol) server enables conversational access to gene regulatory network analysis through Claude Desktop. It provides a natural language interface to the powerful multi-agent gene analysis workflow.

---

## Installation

### 1. Add to Claude Desktop Configuration

Add the following to your Claude Desktop MCP configuration file:

**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "regnetagents": {
      "command": "python",
      "args": ["C:\\Dev\\RegNetAgents\\regnetagents_langgraph_mcp_server.py"]
    }
  }
}
```

**Note**: Update the path to match your RegNetAgents installation directory.

### 2. Verify Installation

Restart Claude Desktop and check that the RegNetAgents tools are available. You should see tools like:
- `validate_gene`
- `query_network`
- `comprehensive_gene_analysis`
- `multi_gene_analysis`
- `pathway_focused_analysis`
- `cross_cell_comparison`

---

## Requirements

### System Requirements
- Python 3.10 or higher
- 8GB+ RAM recommended
- Network cache files (pre-computed ARACNe networks)
- **Optional**: LLM for narrative rationales — Ollama (local/cloud), OpenAI, Anthropic, or any OpenAI-compatible API

### Python Dependencies

Install required packages:

```bash
pip install -r requirements.txt
```

This installs all dependencies including:
- langgraph (workflow orchestration)
- mcp (Model Context Protocol)
- ollama (local LLM inference)
- networkx, pandas, numpy (data processing)
- requests (API calls)

### Optional: LLM Setup (for Narrative Rationales)

**Default behaviour**: `USE_LLM_AGENTS=false` — the system runs entirely on rule-based categorical assessments (high/moderate/low with evidence factors) plus automatic cross-domain contradiction detection. No LLM is needed and no model is loaded. This is the correct setting for MCP clients (Claude Desktop, Cursor, Zed).

**What enabling LLM adds**: Natural-language scientific rationales explaining each domain assessment (cancer, drug, clinical, systems biology).

**Option A — Local Ollama** (free, no API key):
1. Install from https://ollama.com/download
2. Run: `ollama pull llama3.1:8b`
3. In `.env`: set `USE_LLM_AGENTS=true` (keep `LLM_PROVIDER=ollama`)

**Option B — Cloud provider** (OpenAI, Anthropic, or compatible):
```bash
USE_LLM_AGENTS=true
LLM_PROVIDER=openai          # or: anthropic | openai_compatible
LLM_API_KEY=your-key
LLM_MODEL=gpt-4o-mini        # or claude-haiku-4-5-20251001, etc.
# LLM_API_BASE=https://api.groq.com/openai/v1  # for openai_compatible
```

**Performance**:
- Rule-based (default): <1 second per gene
- With local Ollama: ~4 seconds per gene
- With cloud provider: ~5–25 seconds per gene (network latency)

#### Alternative: Ollama Cloud (No Local Installation Required)

**For users without local GPU resources**, Ollama Cloud provides an alternative deployment option using the same methodology as the published results.

**What this is:**
- Cloud-hosted inference using Ollama's free tier
- Access to models like `gemini-3-flash-preview`
- Same analytical workflow as local Ollama
- No local GPU or model downloads required

**Setup:**

1. Get API key from https://ollama.com/settings (free tier available)
2. Edit your `.env` file:

```bash
# Use Ollama Cloud instead of local
OLLAMA_API_KEY=your-api-key-here
OLLAMA_MODEL=gemini-3-flash-preview
USE_LLM_AGENTS=true
```

3. Restart Claude Desktop

**Performance:**
- Ollama Cloud: ~20-25 seconds per gene (network latency included)
- Local Ollama: ~4 seconds per gene (published results)
- Rule-based: <1 second (instant fallback)

**Free Tier Limits:**
- 5 premium requests/month for Gemini models
- Hourly/weekly rate limits

**Note:** This is a **deployment option**, not a methodological change. The analytical workflow, domain agents, and network analysis remain identical to the published approach. Ollama Cloud simply provides an alternative infrastructure for users without local compute resources.

### Network Cache Files

Generate network cache files for 10 cell types:

```bash
# Build all networks
python scripts/build_network_cache.py --all

# Or build specific cell types
python scripts/build_network_cache.py --cell-type epithelial_cell
```

**Cell types available:**
- epithelial_cell
- cd4_t_cells
- cd8_t_cells
- cd14_monocytes
- cd16_monocytes
- cd20_b_cells
- nk_cells
- nkt_cells
- erythrocytes
- monocyte-derived_dendritic_cells

---

## Getting Started

### Example Queries

Once installed, you can ask Claude Desktop natural language questions:

#### **Example 1: Single Gene Analysis**
> *"Analyze TP53 in epithelial cells - show me what it regulates, what regulates it, and its role in cancer pathways"*

**Returns**: Network position, 163 targets, 7 regulators, therapeutic target prioritization (candidate regulators), 16 pathways, **LLM-generated cancer/drug/clinical insights with rationales** (~4 seconds with Ollama, <1 second without)

#### **Example 2: Multi-Gene Characterization**
> *"Characterize these candidate genes for CRC biomarker potential: MYC, CTNNB1, CCND1, TP53, KRAS"*

**Returns**: Regulatory networks, biomarker types, network connectivity metrics, pathway enrichment (FDR < 0.05), standard centrality rankings (1.34 seconds for 5 genes)

#### **Example 3: Cross-Cell Comparison**
> *"Compare TP53 regulatory role across CD4 T cells, B cells, and epithelial cells"*

**Returns**: Cell-type-specific network positions, differential regulation patterns

#### **Example 4: Pathway-Focused**
> *"What pathways is BRCA1 involved in and what genes does it regulate?"*

**Returns**: Reactome enrichment, DNA repair pathways, downstream cascade targets

#### **Example 5: Therapeutic Target Prioritization**
> *"What happens if I inhibit TP53 upstream regulators? Which would be the best candidate for experimental validation?"*

**Returns**: Standard centrality rankings (PageRank, degree centrality) for each regulator, regulatory input loss percentages, network connectivity analysis, ranked candidates for experimental validation using dual ranking approaches

**Note**: Results are topology-based using standard network metrics and prioritize regulators for investigation; does not predict gene expression changes

### Query Tips

- **Be specific about cell type**: "in epithelial cells", "across T cells and B cells"
- **Specify analysis depth**: "comprehensive analysis", "quick overview", "focus on cancer pathways"
- **Ask for specific domains**: "cancer insights", "drug development perspective", "clinical relevance"

---

## MCP Tools Available

### Validation Tools

1. **validate_gene**
   - Quick gene name check against the network (<100ms)
   - Returns basic stats (regulators, targets, role) if found
   - Returns fuzzy-matched suggestions for misspelled gene names
   - Use before full analysis to catch typos instantly
   - Parameters: gene, cell_type

2. **query_network**
   - Instant network queries from pre-computed data (<50ms)
   - Query types: top_regulators, top_targets, gene_neighbors, network_stats
   - Returns gene symbols, counts, PageRank scores, and network statistics
   - Use for quick questions like "top regulators?" or "what regulates MYC?"
   - Parameters: query_type, cell_type, gene (for gene_neighbors), top_n

### Core Analysis Tools

3. **comprehensive_gene_analysis**
   - Full multi-agent analysis (network + regulators + targets + therapeutic target prioritization + pathways + 4 LLM-powered domain agents)
   - **LLM Mode**: AI-generated insights with scientific rationales (Ollama/llama3.1:8b)
   - **Fallback Mode**: Rule-based heuristics if LLM unavailable
   - Parameters: gene, cell_type, analysis_depth
   - Execution: ~4 seconds per gene (LLM mode), <1 second (rule-based mode)
   - Includes therapeutic target prioritization for genes with >5 regulators

4. **multi_gene_analysis**
   - Parallel processing of multiple genes (up to 10)
   - Faster than sequential calls

5. **pathway_focused_analysis**
   - Reactome pathway enrichment analysis
   - Statistical validation (FDR < 0.05)

6. **cross_cell_comparison**
   - Compare gene across all 10 cell types
   - Identifies cell-type-specific regulation

7. **load_gene_results**
   - Load previously saved analysis results from disk

8. **list_available_results**
   - List available result files in the results directory

### Utility Tools

9. **create_analysis_report**
   - Generate formatted JSON/markdown reports

10. **workflow_status**
   - Check analysis progress (for long-running queries)

11. **workflow_insights**
    - Technical performance metrics

12. **export_results**
    - Export analysis results as markdown (renders as tables) or CSV (for spreadsheets)
    - Parameters: gene, format (markdown/csv), cell_type, sections (summary/regulators/targets/pathways/all)
    - Returns plain text, not JSON — markdown renders natively in Claude Desktop

---

## MCP Resources Available

In addition to tools, the server exposes **browsable MCP Resources** for lightweight data discovery without running an analysis.

### 1. `regnetagents://cell-types`
**List all available cell types** with gene and edge counts.

Example query in Claude Desktop:
```
What cell types are available and how large are the networks?
```

### 2. `regnetagents://network/{cell_type}`
**Network summary** for a specific cell type — gene/edge counts, density, and top 10 regulators by PageRank.

Example query:
```
Show me the epithelial_cell network summary
```

### 3. `regnetagents://gene/{gene_symbol}`
**Gene lookup** across all cell types — shows which networks contain the gene and its regulator/target counts in each.

Example query:
```
Where does TP53 appear across cell types?
```

---

## Architecture Overview

### Thin-Wrapper Design

```
┌─────────────────────────────────────────┐
│   MCP Server (regnetagents_langgraph_mcp_server.py)
│   • Protocol translation                │
│   • Tool registration (10 tools)        │
│   • Claude Desktop integration          │
│   • ~466 lines                           │
└─────────────┬───────────────────────────┘
              ↓
┌─────────────────────────────────────────┐
│   LangGraph Workflow (regnetagents_langgraph_workflow.py)
│   • Multi-agent orchestration           │
│   • Intelligent routing                 │
│   • State management                    │
│   • Domain analysis logic               │
│   • ~1370 lines                          │
└─────────────────────────────────────────┘
```

**Benefits**:
- ✅ Testable: Workflow independent of MCP protocol
- ✅ Portable: Same workflow can power CLI, API, notebooks
- ✅ Maintainable: Update workflow without changing MCP layer

---

## Troubleshooting

### "Network cache not found"

**Solution**: Run `python scripts/build_network_cache.py --all` to generate network files.

### "Reactome API timeout"

**Solution**: Check internet connection. Reactome API requires network access.

### "Gene not found in network"

**Solution**: Gene may not be present in the selected cell type network. Try a different cell type or verify gene symbol is correct.

### MCP Server Not Appearing in Claude Desktop

**Solutions**:
1. Verify path in `claude_desktop_config.json` is correct
2. Restart Claude Desktop completely
3. Check Python is in PATH: `python --version`
4. Check server logs in Claude Desktop settings

---

## Performance Notes

- **First query**: May take 2-3 seconds (network cache loading)
- **Subsequent queries**: 2-5 seconds typical
- **Multi-gene (5 genes)**: 3-8 seconds with parallel processing
- **Cross-cell comparison**: +1-2 seconds per additional cell type

---

## Credits & Attribution

### Network Data Sources

This MCP server uses **pre-computed gene regulatory networks** derived from single-cell RNA-seq data:

- **Data Source**: CellxGene Data Portal (Chan Zuckerberg Initiative)
- **Processing Method**: ARACNe algorithm (statistical mutual information method for network inference)
- **Processing Lab**: Califano Lab, Columbia University
- **Network Format**: Pre-computed graph structures (not neural network model weights)

**Note**: The RegNetAgents framework performs network analysis on pre-computed regulatory networks using NetworkX graph algorithms. It does NOT use neural networks or machine learning for network inference.

### Data Sources

The regulatory networks were derived from publicly available single-cell RNA-seq data:

- **Data Source**: CellxGene Data Portal (Chan Zuckerberg Initiative)
- **Network Inference**: ARACNe algorithm (Califano Lab, Columbia University)
- **Processing**: Standard bioinformatics pipeline (mutual information-based network reconstruction)
- **License**: All data sources are publicly available

### Additional Data Sources

- **Pathway Analysis**: Reactome Pathway Database (manually curated, peer-reviewed)
- **Gene Annotations**: Local NCBI and UniProt databases; Ensembl REST API for gene ID mapping

### RegNetAgents Agents Framework

Multi-agent orchestration, intelligent workflow routing, parallel gene analysis, LLM-powered domain-specific analysis agents, and Claude Desktop MCP integration built on top of the pre-computed regulatory networks.

---

## Support

For issues, questions, or support:
- **Contact**: jbird@birdaisolutions.com
- **Documentation**: See `README.md` for general overview

---

**Last Updated**: 2025-10-22
