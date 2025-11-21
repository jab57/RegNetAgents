# RegNetAgents Agents MCP Server Setup Guide

> **Note:** This guide is for users who have received access to the RegNetAgents software. The code is available from the corresponding author upon reasonable request for academic research purposes. Contact jbird@birdaisolutions.com for access inquiries.

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
- `comprehensive_gene_analysis`
- `analyze_regulators`
- `analyze_targets`
- `pathway_enrichment`
- `cross_cell_comparison`

---

## Requirements

### System Requirements
- Python 3.8 or higher
- 8GB+ RAM recommended
- Network cache files (pre-computed ARACNe networks)
- **Optional**: Ollama for LLM-powered domain insights (recommended)

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

### Optional: Ollama Setup (for LLM-Powered Insights)

**What this adds**: AI-generated scientific rationales for all domain analyses (cancer, drug, clinical, systems biology)

**Setup (5 minutes)**:

1. Install Ollama: https://ollama.com/download
2. Pull model: `ollama pull llama3.1:8b`
3. Configure: Create `.env` file in RegNetAgents directory:

```bash
# Copy from example
cp .env.example .env
```

**Performance**:
- With Ollama: ~4 seconds per gene (LLM-generated insights with rationales)
- Without Ollama: Instant (rule-based heuristics, automatic fallback)
- Still 100-200× faster than manual literature review

**Note**: System works perfectly without Ollama using fast rule-based heuristics. LLM mode adds scientific rationales and interpretations.

### Network Cache Files

Generate network cache files for 10 cell types:

```bash
# Build all networks
python build_network_cache.py --all

# Or build specific cell types
python build_network_cache.py --cell-type epithelial_cell
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

**Returns**: Network position, 163 targets, 7 regulators, perturbation analysis (candidate regulators), 16 pathways, **LLM-generated cancer/drug/clinical insights with rationales** (~4 seconds with Ollama, <1 second without)

#### **Example 2: Multi-Gene Characterization**
> *"Characterize these candidate genes for CRC biomarker potential: MYC, CTNNB1, CCND1, TP53, KRAS"*

**Returns**: Regulatory networks, biomarker types, network connectivity metrics, pathway enrichment (FDR < 0.05), standard centrality rankings (1.34 seconds for 5 genes)

#### **Example 3: Cross-Cell Comparison**
> *"Compare TP53 regulatory role across CD4 T cells, B cells, and epithelial cells"*

**Returns**: Cell-type-specific network positions, differential regulation patterns

#### **Example 4: Pathway-Focused**
> *"What pathways is BRCA1 involved in and what genes does it regulate?"*

**Returns**: Reactome enrichment, DNA repair pathways, downstream cascade targets

#### **Example 5: Perturbation Analysis**
> *"What happens if I inhibit TP53 upstream regulators? Which would be the best candidate for experimental validation?"*

**Returns**: Standard centrality rankings (PageRank, degree centrality) for each regulator, regulatory input loss percentages, network connectivity analysis, ranked candidates for experimental validation using dual ranking approaches

**Note**: Results are topology-based using standard network metrics and prioritize regulators for investigation; does not predict gene expression changes

### Query Tips

- **Be specific about cell type**: "in epithelial cells", "across T cells and B cells"
- **Specify analysis depth**: "comprehensive analysis", "quick overview", "focus on cancer pathways"
- **Ask for specific domains**: "cancer insights", "drug development perspective", "clinical relevance"

---

## MCP Tools Available

### Core Analysis Tools

1. **comprehensive_gene_analysis**
   - Full multi-agent analysis (network + regulators + targets + perturbation + pathways + 4 LLM-powered domain agents)
   - **LLM Mode**: AI-generated insights with scientific rationales (Ollama/llama3.1:8b)
   - **Fallback Mode**: Rule-based heuristics if LLM unavailable
   - Parameters: gene, cell_type, analysis_depth
   - Execution: ~4 seconds per gene (LLM mode), <1 second (rule-based mode)
   - Includes perturbation analysis for genes with >5 regulators

2. **analyze_regulators**
   - Detailed upstream regulator analysis
   - Shows genes that control the target gene
   - Includes perturbation analysis (simulates inhibiting each regulator)

3. **analyze_targets**
   - Detailed downstream target analysis
   - Shows genes controlled by the regulator

4. **pathway_enrichment**
   - Reactome pathway enrichment analysis
   - Statistical validation (FDR < 0.05)

5. **cross_cell_comparison**
   - Compare gene across all 10 cell types
   - Identifies cell-type-specific regulation

### Utility Tools

6. **create_analysis_report**
   - Generate formatted JSON/markdown reports

7. **workflow_status**
   - Check analysis progress (for long-running queries)

8. **workflow_insights**
   - Technical performance metrics

---

## Architecture Overview

### Thin-Wrapper Design

```
┌─────────────────────────────────────────┐
│   MCP Server (regnetagents_langgraph_mcp_server.py)
│   • Protocol translation                │
│   • Tool registration (8 tools)         │
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

**Solution**: Run `python build_network_cache.py --all` to generate network files.

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
- **Conference Poster**: See `REGNETAGENTS_CONFERENCE_POSTER.md` for technical details

---

**Last Updated**: 2025-10-22
