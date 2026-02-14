# RegNetAgents

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18500028.svg)](https://doi.org/10.5281/zenodo.18500028)

**Multi-Agent LLM Framework for Gene Regulatory Network Analysis**

RegNetAgents automates gene regulatory network analysis through Claude Desktop. Built on pre-computed ARACNe networks from single-cell RNA-seq data and powered by LangGraph workflow orchestration, it deploys four specialized domain agents (cancer biology, drug discovery, clinical relevance, systems biology) that generate scientific insights using Ollama inference (local or cloud), with support for a range of open-weight models.

Analyze **10 cell types** including immune cells, blood cells, and epithelial tissue—with AI-generated rationales and interpretations.

---

## Data Sources

- **Primary Source**: [GREmLN Foundation Model](https://github.com/czi-ai/GREmLN) (Zhang et al. 2025)
- **Networks**: 10 cell-type-specific ARACNe networks from 500K+ single cells
- **Underlying Data**: CELLxGENE Census 2024-07-01
- **Development**: CZ Biohub NY / Columbia University (Califano Lab)

See [DATA_SOURCES.md](docs/DATA_SOURCES.md) for complete details.

---

## Quick Start

**Time required:** 5-10 minutes

### Prerequisites

- Python 3.10+
- Claude Desktop
- Internet connection (for Reactome API)

### Installation

```bash
# Clone and setup
git clone https://github.com/jab57/RegNetAgents.git
cd RegNetAgents
python -m venv env

# Activate (Windows)
env\Scripts\activate
# Activate (Linux/Mac)
source env/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Install Ollama (Optional but Recommended)

```bash
# Download from https://ollama.com/download
ollama pull llama3.1:8b
cp .env.example .env
```

If Ollama unavailable, system uses fast rule-based heuristics automatically.

### Configure Claude Desktop

Add to `claude_desktop_config.json`:
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **Mac**: `~/Library/Application Support/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "regnetagents": {
      "command": "python",
      "args": ["C:\\path\\to\\RegNetAgents\\regnetagents_langgraph_mcp_server.py"]
    }
  }
}
```

Replace path with your actual installation path. Restart Claude Desktop.

### Verify Installation

```bash
python verify_installation.py
```

Expected: `7/7 checks passed` (or `6/7` without Ollama)

### Test It

In Claude Desktop:
```
Analyze the TP53 gene in epithelial cells
```

### Programmatic Usage

```bash
python demo_biomarker_analysis.py
```

---

## Features

- **Network Analysis**: Gene positions in regulatory networks
- **Therapeutic Target Prioritization**: PageRank-based ranking of upstream regulators
- **10 Cell Types**: Immune cells, blood cells, epithelial tissue
- **LLM-Powered Insights**: 4 domain agents with scientific rationales
- **Pathway Enrichment**: Reactome API with FDR correction
- **Dual Mode**: LLM-powered or deterministic rule-based fallback

### MCP Tools

| Tool | Description |
|------|-------------|
| `validate_gene` | Quick gene name check with fuzzy suggestions (<100ms) |
| `query_network` | Instant network queries — top regulators, targets, neighbors, stats (<50ms) |
| `comprehensive_gene_analysis` | Full analysis with domain insights |
| `multi_gene_analysis` | Parallel processing of multiple genes |
| `pathway_focused_analysis` | Reactome pathway enrichment |
| `cross_cell_comparison` | Gene behavior across cell types |
| `load_gene_results` | Load previously saved analysis results |
| `list_available_results` | List available result files |
| `workflow_status` | Execution monitoring |
| `workflow_insights` | Performance analytics |
| `create_analysis_report` | Generate reports |

---

## Cell Types (10 Available)

**Immune & Blood**: CD4/CD8 T cells, CD14/CD16 Monocytes, CD20 B cells, NK/NKT cells, Erythrocytes, Dendritic cells

**Epithelial**: Epithelial cells (183,247 edges - largest network)

See [ADDING_NEW_CELL_TYPES.md](docs/ADDING_NEW_CELL_TYPES.md) to expand.

---

## Example Queries

```
What regulates BRCA1 in epithelial cells?
Compare MYC, TP53, and KRAS across different cell types
What pathways is IL6 involved in for immune cells?
```

**Popular genes**: TP53, BRCA1, APC, KRAS, MYC, IL6, TNF

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Network index missing | Run `python scripts/build_network_cache.py --all` |
| Tools not in Claude Desktop | Check path in config, use absolute path, restart |
| Module not found | Activate venv, run `pip install -r requirements.txt` |
| Empty results | Use `validate_gene` to check gene name first |
| Pathway analysis fails | Check internet connection |

---

## Performance

| Mode | Single Gene | 5 Genes |
|------|-------------|---------|
| Rule-based | 0.6 sec | 15 sec |
| LLM-powered | 15 sec | 62 sec |

---

## Project Structure

```
RegNetAgents/
├── regnetagents_langgraph_mcp_server.py  # MCP server
├── regnetagents_langgraph_workflow.py    # Core workflow
├── regnetagents/                          # Package
├── models/networks/                       # Network data (10 cell types)
├── tests/                                 # Test suite
├── docs/                                  # Documentation
└── scripts/                               # Utilities
```

---

## Documentation

- [REGNETAGENTS_MCP_SETUP.md](docs/REGNETAGENTS_MCP_SETUP.md) - Setup guide
- [DATA_SOURCES.md](docs/DATA_SOURCES.md) - Cell types and sources
- [REGNETAGENTS_Analysis_Pipeline.md](docs/REGNETAGENTS_Analysis_Pipeline.md) - Architecture
- [docs/README.md](docs/README.md) - Full index

---

## Development Context

This project was developed with assistance from Claude Code (Anthropic). It is a research prototype for hypothesis generation—users should validate findings experimentally. See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines.

---

## License

MIT License - see [LICENSE](LICENSE)
