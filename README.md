# GREmLN MCP Server

Gene Regulatory network-informed Large Language Model (GREmLN) MCP Server for analyzing regulatory networks and gene relationships.

## Setup Instructions

### 1. Environment Setup
```bash
# Create and activate virtual environment
python -m venv gremln-env
source gremln-env/Scripts/activate  # Windows
# source gremln-env/bin/activate    # Linux/Mac

# Install dependencies
pip install -r requirements.txt
```

### 2. Install GREmLN Dependencies

**⚠️ scGraphLLM Installation Required**

The scGraphLLM package is not available on PyPI. You need to:

1. **Find the source repository** - Check the GREmLN paper or contact authors
2. **Install from source**:
   ```bash
   # If available on GitHub:
   pip install git+https://github.com/[author]/scGraphLLM.git
   
   # Or clone and install:
   git clone https://github.com/[author]/scGraphLLM.git
   cd scGraphLLM
   pip install -e .
   ```

### 3. Download Model and Data

You'll need:
- **GDTransformer model checkpoint**: `models/model.ckpt`
- **Regulatory networks**: Cell-type specific `.tsv` files in `models/networks/`

**Data Sources:** The regulatory networks are derived from CellxGene single-cell RNA-seq data processed through a comprehensive ARACNe-based pipeline. See `gremln_source/scripts/README.md` for complete preprocessing details and `GREmLN_Analysis_Pipeline.md` for data lineage documentation.

Expected structure:
```
GREmLN/
├── models/
│   ├── model.ckpt
│   └── networks/
│       ├── cd14_monocytes/
│       │   ├── network.tsv          # Original ARACNe output
│       │   └── network_index.pkl    # Optimized cache (generated)
│       ├── cd16_monocytes/
│       ├── cd4_t_cells/
│       ├── cd8_t_cells/
│       ├── epithelial_cell/
│       └── ... (15 total cell types)
```

**Data Processing Pipeline:**
```
CellxGene Data → QC/Metacells → ARACNe Networks → TSV Files → Pickle Caches
```

### 4. Test Installation

```bash
python gremln_mcp_server.py
```

## Features

### Current Capabilities
- **Network Analysis**: Analyze gene positions in regulatory networks
- **Cell-Type Specific**: Support for 15 different cell types
- **Gene Similarity**: Find genes with similar regulatory patterns
- **Network Statistics**: Get detailed network information

### Cell Types Supported
**15 Major Human Cell Types** (derived from CellxGene datasets):

**Immune & Blood Cell Types:**
- **CD14 Monocytes** - Classical circulating monocytes
- **CD16 Monocytes** - Non-classical patrolling monocytes
- **CD20 B Cells** - B lymphocytes
- **CD4 T Cells** - Helper T cells
- **CD8 T Cells** - Cytotoxic T cells
- **Erythrocytes** - Red blood cells
- **NK Cells** - Natural killer cells
- **NKT Cells** - Natural killer T cells
- **Monocyte-derived Dendritic Cells** - Antigen-presenting cells

**Tissue & Organ Cell Types:**
- **Epithelial Cells** - Barrier tissue cells
- **Hepatocytes** - Liver cells (drug metabolism)
- **Cardiomyocytes** - Heart muscle cells
- **Neurons** - Brain and nervous system cells
- **Fibroblasts** - Connective tissue cells
- **Endothelial Cells** - Blood vessel lining cells

Each cell type contains ARACNe-generated gene regulatory networks with 403-183,247 regulatory edges based on the top 1024 highly variable genes.

### MCP Tools Available

1. **analyze_gene_network**
   - Analyze gene's regulatory context
   - Parameters: gene, cell_type

2. **find_similar_genes**
   - Find genes with similar regulatory patterns
   - Parameters: gene, cell_type, top_n

3. **get_network_info**
   - Get network statistics
   - Parameters: cell_type

## Architecture

- **gremln_mcp_server.py**: Main MCP server
- **agents/**: Multi-agent framework (planned)
- **utils/**: Utility functions
- **models/**: Model checkpoints and networks
- **complete_gene_service.py**: Gene annotation service

## Status

🟡 **Partial Implementation** 
- ✅ Basic MCP server structure
- ✅ Core dependencies installed
- ⚠️ Requires scGraphLLM installation
- ⚠️ Requires model/network data
- 🔄 Multi-agent framework pending

## Next Steps

1. Locate and install scGraphLLM package
2. Download GREmLN model checkpoint
3. Download regulatory network data
4. Implement advanced multi-agent features
5. Add visualization capabilities