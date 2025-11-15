# Data Sources for RegNetAgents Cell Types

This document provides comprehensive information about the data sources used for the 10 currently available cell types in RegNetAgents, as well as detailed instructions for adding new cell types.

---

## Current Cell Types (10 Available)

RegNetAgents currently supports **10 human cell types** with pre-computed gene regulatory networks, derived from publicly available single-cell RNA-seq data processed through the ARACNe algorithm.

### Data Source Overview

**Primary Source**: [GREmLN Foundation Model](https://github.com/czi-ai/GREmLN) (Chan Zuckerberg Initiative / CZ Biohub NY)
**Download Location**: [GREmLN Quickstart Tutorial](https://virtualcellmodels.cziscience.com/quickstart/gremln-quickstart) (networks available via Google Drive)
**Underlying Data**: [CellxGene Data Portal](https://cellxgene.cziscience.com/) (11M scRNA-seq profiles, 162 cell types from Census release 2024-07-01)
**Networks Used**: 10 cell types (500,000+ cells subset)
**Processing Method**: ARACNe algorithm via RegNetAgents framework
**Development Team**: Zhang et al. (2025), Califano Lab (Columbia University / CZ Biohub NY)
**Format**: Pre-computed networks as TSV files (network.tsv), converted to optimized pickle caches (network_index.pkl)
**Publication**: bioRxiv 2025.07.03.663009

### How Networks Were Obtained

The 10 cell types were obtained from the **RegNetAgents Quickstart Tutorial**, which provides:

1. **`human_immune_cells.h5ad`** - 9 immune and blood cell types (2,473 cells from bone marrow and peripheral blood)
2. **`epithelial_cells.h5ad`** - 1 epithelial cell type (1,000 cells from CELLxGENE)
3. **Pre-computed ARACNe networks** - TSV files at `GREmLN_Tutorial/networks/{cell_type}/network.tsv`

**Download Instructions**: Available via gdown from Google Drive at the [GREmLN Quickstart Tutorial](https://virtualcellmodels.cziscience.com/quickstart/gremln-quickstart)

---

## Detailed Cell Type Information

### Immune & Blood Cell Types (9 types)

These cell types were extracted from the `human_immune_cells.h5ad` dataset, which aggregated multiple human immune cell populations from CellxGene collections.

#### 1. **CD14 Monocytes** (`cd14_monocytes`)
- **Description**: Classical circulating monocytes from peripheral blood
- **Network Size**: 2,009 regulatory edges
- **Cell Type Marker**: CD14 (high expression)
- **Research Applications**: Inflammation, innate immunity, monocyte differentiation
- **Typical Source Tissue**: Peripheral blood, bone marrow

#### 2. **CD16 Monocytes** (`cd16_monocytes`)
- **Description**: Non-classical patrolling monocytes
- **Network Size**: 1,236 regulatory edges
- **Cell Type Marker**: CD16 (FCGR3A)
- **Research Applications**: Vascular patrolling, tissue repair, inflammation resolution
- **Typical Source Tissue**: Peripheral blood

#### 3. **CD20 B Cells** (`cd20_b_cells`)
- **Description**: B lymphocytes (antibody-producing cells)
- **Network Size**: 1,128 regulatory edges
- **Cell Type Markers**: CD20 (MS4A1), CD19
- **Research Applications**: Adaptive immunity, antibody production, B cell malignancies
- **Typical Source Tissue**: Lymph nodes, spleen, peripheral blood

#### 4. **CD4 T Cells** (`cd4_t_cells`)
- **Description**: Helper T lymphocytes
- **Network Size**: 1,371 regulatory edges
- **Cell Type Marker**: CD4
- **Research Applications**: Adaptive immunity, helper T cell responses, autoimmunity
- **Typical Source Tissue**: Lymph nodes, spleen, thymus, peripheral blood

#### 5. **CD8 T Cells** (`cd8_t_cells`)
- **Description**: Cytotoxic T lymphocytes
- **Network Size**: 3,154 regulatory edges
- **Cell Type Marker**: CD8A, CD8B
- **Research Applications**: Cell-mediated immunity, cancer immunotherapy, viral immunity
- **Typical Source Tissue**: Lymph nodes, spleen, thymus, peripheral blood

#### 6. **Erythrocytes** (`erythrocytes`)
- **Description**: Red blood cells (oxygen transport)
- **Network Size**: 19,398 regulatory edges
- **Cell Type Markers**: HBA1, HBA2, HBB (hemoglobin genes)
- **Research Applications**: Anemia, hemoglobinopathies, erythropoiesis
- **Typical Source Tissue**: Peripheral blood, bone marrow

#### 7. **NK Cells** (`nk_cells`)
- **Description**: Natural killer cells (innate lymphocytes)
- **Network Size**: 404 regulatory edges
- **Cell Type Markers**: NCAM1 (CD56), NCR1
- **Research Applications**: Innate immunity, cancer immunosurveillance, viral defense
- **Typical Source Tissue**: Peripheral blood, spleen, lymph nodes

#### 8. **NKT Cells** (`nkt_cells`)
- **Description**: Natural killer T cells (bridge innate and adaptive immunity)
- **Network Size**: 2,509 regulatory edges
- **Cell Type Markers**: CD3E, KLRB1
- **Research Applications**: Immune regulation, autoimmunity, tumor immunity
- **Typical Source Tissue**: Liver, thymus, peripheral blood

#### 9. **Monocyte-Derived Dendritic Cells** (`monocyte-derived_dendritic_cells`)
- **Description**: Antigen-presenting cells derived from monocytes
- **Network Size**: 5,317 regulatory edges
- **Cell Type Markers**: CD1C, ITGAX (CD11c)
- **Research Applications**: Antigen presentation, T cell activation, vaccine development
- **Typical Source Tissue**: Peripheral blood (in vitro differentiated), tissues

---

### Epithelial Cell Type (1 type)

#### 10. **Epithelial Cells** (`epithelial_cell`)
- **Description**: Barrier tissue cells from various epithelial tissues
- **Network Size**: 183,247 regulatory edges (largest network)
- **Cell Type Markers**: EPCAM, KRT8, KRT18
- **Research Applications**: Cancer biology (most carcinomas), barrier function, tissue regeneration
- **Typical Source Tissues**: Lung, intestine, skin, breast, prostate
- **Note**: Likely aggregated from multiple epithelial tissue sources

---

## How to Download the Cell Type Networks

**📥 RECOMMENDED METHOD: Use the RegNetAgents Quickstart Tutorial**

This is how the RegNetAgents networks were obtained. Follow these steps to get the same 10 pre-computed cell type networks.

### Quick Start (TL;DR)

```bash
# 1. Install gdown
pip install gdown

# 2. Download the RegNetAgents tutorial folder
gdown --folder https://drive.google.com/drive/folders/1cMR9HoAC22i6sKSWgfQUEQRf0UP_w3_m?usp=sharing

# 3. Navigate into the downloaded folder and copy networks
cd GREmLN_tutorial
cp -r networks/* /path/to/RegNetAgents/models/networks/
```

### Step-by-Step Download Instructions

#### Prerequisites
Install `gdown` (Google Drive downloader):
```bash
pip install gdown
```

#### Step 1: Visit the Tutorial Page

Go to: [GREmLN Quickstart Tutorial](https://virtualcellmodels.cziscience.com/quickstart/gremln-quickstart)

#### Step 2: Locate the Download Links

On the tutorial page, find the Google Drive download links for:
- RegNetAgents model weights (`model.ckpt`)
- Tutorial data package (includes networks and datasets)

#### Step 3: Download Using gdown

The tutorial provides a Google Drive folder link. Use this exact command:

```bash
# Download the entire RegNetAgents tutorial folder from Google Drive
# This creates a folder named "GREmLN_tutorial" in your current directory
gdown --folder https://drive.google.com/drive/folders/1cMR9HoAC22i6sKSWgfQUEQRf0UP_w3_m?usp=sharing
```

**What happens:**
- `gdown` creates a folder named `GREmLN_tutorial` in your current directory
- Downloads all contents into that folder
- No extraction needed - files are organized and ready to use

The downloaded folder contains:
- `data/` - H5AD datasets (human_immune_cells.h5ad, epithelial_cells.h5ad)
- `networks/` - Pre-computed ARACNe networks for all 10 cell types (network.tsv files)
- `model.ckpt` - RegNetAgents model weights
- `.DS_Store` - System file (can be ignored)

#### Step 4: Verify What You Downloaded

After the download completes, verify the folder structure:

```bash
# Check that the folder was created
ls -la GREmLN_tutorial

# You should see:
# data/
# networks/
# model.ckpt
# .DS_Store (system file, can ignore)
```

The folder structure looks like:
```
GREmLN_tutorial/
  data/
    human_immune_cells.h5ad
    epithelial_cells.h5ad
  networks/
    cd14_monocytes/network.tsv
    cd16_monocytes/network.tsv
    cd20_b_cells/network.tsv
    cd4_t_cells/network.tsv
    cd8_t_cells/network.tsv
    erythrocytes/network.tsv
    nk_cells/network.tsv
    nkt_cells/network.tsv
    monocyte-derived_dendritic_cells/network.tsv
    epithelial_cell/network.tsv
  model.ckpt              # RegNetAgents trained model weights
  .DS_Store               # System file (ignore)
```

#### Step 5: Copy Files to RegNetAgents

Once downloaded, copy the network files to your RegNetAgents project:

```bash
# Navigate into the downloaded folder
cd GREmLN_tutorial

# Copy ONLY the network files to RegNetAgents (this is all you need!)
cp -r networks/* /path/to/RegNetAgents/models/networks/
```

**Important Note**: RegNetAgents only uses the `networks/` folder (network.tsv files). The download also includes:
- `data/` - H5AD expression files (NOT used by RegNetAgents)
- `model.ckpt` - GREmLN model weights (NOT used by RegNetAgents)
- `vocab.csv` - Gene vocabulary (NOT used by RegNetAgents)

These files are for users who want to generate RegNetAgents embeddings or test the foundation model. **You can safely delete them after copying networks** to save ~2-3 GB of disk space:

```bash
# Optional: Clean up unused files to save space
cd GREmLN_tutorial
rm -rf data/
rm model.ckpt
rm vocab.csv
```

#### Step 6: Verify Installation

Check that network files exist:
```bash
ls models/networks/cd14_monocytes/
# Should show: network.tsv
```

Each `network.tsv` file contains ARACNe-inferred regulatory edges in tab-separated format.

**✅ You're done!** The RegNetAgents system will automatically recognize these cell types.

---

## Converting Networks to Pickle Cache (Optional Performance Optimization)

RegNetAgents can use the `network.tsv` files directly, but for better performance, convert them to pickle caches:

```bash
python scripts/build_network_cache.py
```

This creates `network_index.pkl` files with pre-computed PageRank for faster queries.

---

## Troubleshooting Download Issues

### Issue: Can't find download links on tutorial page

**Solution:**
- The tutorial page is hosted at https://virtualcellmodels.cziscience.com/quickstart/gremln-quickstart
- Look for sections titled "Download", "Data", or "Getting Started"
- If the page structure has changed, check the [GREmLN GitHub repository](https://github.com/czi-ai/GREmLN) for updated links
- Contact the GREmLN team at opensource@chanzuckerberg.com if download links are unavailable

### Issue: gdown fails with "Permission denied" or "Access denied"

**Solution:**
```bash
# Make sure gdown is up to date:
pip install --upgrade gdown

# Retry the folder download:
gdown --folder https://drive.google.com/drive/folders/1cMR9HoAC22i6sKSWgfQUEQRf0UP_w3_m?usp=sharing

# Or manually download:
# 1. Visit the Google Drive link in your browser
# 2. Download the entire folder manually
# 3. Extract the zip file if downloaded as archive
```

### Issue: Download is very slow or times out

**Solution:**
```bash
# Try downloading with remaining flag (resumes interrupted downloads):
gdown --folder --remaining-ok https://drive.google.com/drive/folders/1cMR9HoAC22i6sKSWgfQUEQRf0UP_w3_m?usp=sharing

# Or split into smaller downloads if needed:
# Download just the networks folder (if available separately)
```

### Issue: Network files are missing after download

**Solution:**
```bash
# Check that the folder was created
ls -la GREmLN_tutorial

# Navigate into the downloaded folder
cd GREmLN_tutorial

# Verify directory structure:
find . -name "network.tsv"

# Expected locations:
# ./networks/cd14_monocytes/network.tsv
# ./networks/epithelial_cell/network.tsv
# etc.

# Check that networks folder exists:
ls -la networks/

# If files are missing, delete and re-download:
cd ..
rm -rf GREmLN_tutorial
gdown --folder https://drive.google.com/drive/folders/1cMR9HoAC22i6sKSWgfQUEQRf0UP_w3_m?usp=sharing
```

### Issue: Tutorial page or GitHub repo has moved

**Solution:**
- Check the Virtual Cells Platform main page: https://virtualcellmodels.cziscience.com/
- Search for "RegNetAgents" on the Chan Zuckerberg Initiative website
- Check the RegNetAgents publication (bioRxiv 2025.07.03.663009) for updated links
- Contact: opensource@chanzuckerberg.com

### Need Help? Contact CZI

If download issues persist:
- **Email**: opensource@chanzuckerberg.com
- **Subject**: "RegNetAgents Tutorial Data Access Request"
- **Include**: Your research affiliation and intended use case
- **Note**: Data is publicly available; the team can provide alternative download methods

---

## ARACNe Network File Format

**Understanding the Downloaded Network Files**

The `network.tsv` files you downloaded from the RegNetAgents tutorial contain ARACNe-inferred regulatory relationships.

### TSV File Structure (network.tsv)

The ARACNe algorithm outputs tab-separated value files with the following columns:

```
regulator.values    target.values    mi.values    scc.values    count.values    log.p.values
ENSG00000213626    ENSG00000233927    0.120847    0.123325    1    -0.674163
ENSG00000213626    ENSG00000213741    0.11843     0.116792    2    -1.84561
ENSG00000213626    ENSG00000177954    0.167013    0.11157     1    -0.674163
```

### Column Descriptions

| Column | Description | Typical Range |
|--------|-------------|---------------|
| **regulator.values** | Ensembl gene ID of the transcription factor/regulator | ENSG00000* |
| **target.values** | Ensembl gene ID of the target gene | ENSG00000* |
| **mi.values** | Mutual information score (strength of regulatory relationship) | 0.0 - 0.5 |
| **scc.values** | Spearman correlation coefficient (direction: + = activation, - = repression) | -1.0 to +1.0 |
| **count.values** | Number of bootstrap iterations where edge appeared (robustness) | 1 - 100 |
| **log.p.values** | Log-transformed p-value (statistical significance) | Negative values |

### ARACNe-AP Processing Parameters

The networks were generated with the following parameters:

```bash
./aracne3_app_release \
    --input metacells_{cell_type}.txt \
    --output {cell_type}_network.tsv \
    --pvalue 1e-8 \
    --dpi 1.0 \
    --threads 16 \
    --bootstrap 100
```

**Key Parameters**:
- **Mutual Information**: Adaptive partitioning for robust MI estimation
- **DPI (Data Processing Inequality)**: Tolerance = 1.0 to remove indirect edges
- **P-value threshold**: 1e-8 for statistical significance
- **Bootstrap iterations**: 100 for edge robustness assessment
- **Processing time**: 12-14 hours per cell type on HPC (64-128 GB RAM)

### Data Pre-Processing (Before ARACNe)

1. **Quality Control**:
   ```python
   import scanpy as sc

   # Filter cells
   sc.pp.filter_cells(adata, min_genes=200)
   sc.pp.filter_cells(adata, max_genes=8000)  # Doublet removal

   # Filter genes
   sc.pp.filter_genes(adata, min_cells=3)

   # Mitochondrial content
   adata = adata[adata.obs['pct_counts_mt'] < 25]
   ```

2. **Normalization**:
   ```python
   # CPM normalization + log1p
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)
   ```

3. **Highly Variable Genes**:
   ```python
   # Select top 1,024 HVGs
   sc.pp.highly_variable_genes(adata, n_top_genes=1024)
   adata = adata[:, adata.var['highly_variable']]
   ```

4. **Metacell Aggregation**:
   ```python
   # Aggregate every 5 cells into metacells (reduces noise)
   # Creates ~500-2000 metacells from original single cells
   metacells = aggregate_cells(adata, n_cells_per_metacell=5)
   ```

### Conversion to Pickle Format

For fast loading in RegNetAgents, TSV files are converted to NetworkX pickle caches:

```python
import networkx as nx
import pandas as pd
import pickle

# Load ARACNe TSV
edges = pd.read_csv('network.tsv', sep='\t')

# Build directed graph
G = nx.DiGraph()
for _, row in edges.iterrows():
    G.add_edge(
        row['regulator.values'],
        row['target.values'],
        mi=row['mi.values'],
        scc=row['scc.values'],
        count=row['count.values'],
        log_p=row['log.p.values']
    )

# Pre-compute PageRank for fast perturbation analysis
pagerank = nx.pagerank(G, alpha=0.85)
for node in G.nodes():
    G.nodes[node]['pagerank'] = pagerank.get(node, 0)

# Save optimized cache
with open('network_index.pkl', 'wb') as f:
    pickle.dump({
        'graph': G,
        'pagerank': pagerank,
        'metadata': {
            'cell_type': 'cd14_monocytes',
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'source': 'RegNetAgents ARACNe networks'
        }
    }, f)
```

**Loading time**: <0.1 seconds (compared to ~5 seconds for TSV parsing)

---

# APPENDIX: For Advanced Users

**⚠️ The following sections are for advanced users who want to generate networks for additional cell types not included in the RegNetAgents tutorial.**

**Most users should stop here** - the tutorial download provides all necessary networks.

---

## Generating Custom Networks from CELLxGENE Data

If you need cell types beyond the 10 provided in the tutorial, you can generate your own ARACNe networks from CELLxGENE data. **This requires:**
- High-performance computing resources (64-128 GB RAM)
- 12-14 hours processing time per cell type
- Bioinformatics expertise

### Step 1: Download Data from CELLxGENE

Visit: https://cellxgene.cziscience.com/

### Step 2: Quality Control

```python
import scanpy as sc

# Load and validate data
adata = sc.read_h5ad('human_immune_cells.h5ad')

# Quality metrics applied:
# - Min genes per cell: 200
# - Max genes per cell: 8,000 (doublet filtering)
# - Max mitochondrial content: 25%
# - Min cells per gene: 3
```

### Step 3: Cell Type Selection

From each H5AD file, cells were filtered by cell type annotation:

```python
# Example for CD14 monocytes
cd14_monocytes = adata[adata.obs['cell_type'] == 'CD14+ monocyte']

# Metacell aggregation (5 cells per metacell)
# Top 1,024 highly variable genes selected
```

### Step 4: ARACNe Network Generation

```bash
# ARACNe processing (Califano Lab pipeline)
# ~12-14 hours per cell type on HPC
./ARACNe3_app_release \
    --input metacells_cd14_monocytes.txt \
    --output cd14_monocytes_network.tsv \
    --pvalue 1e-8 \
    --threads 16
```

### Step 5: Network Cache Generation

```bash
# Convert ARACNe TSV to optimized pickle cache
python scripts/build_network_cache.py cd14_monocytes
```

---

## Finding Similar Datasets on CellxGene

### For Immune Cell Types

**Recommended Collections:**
1. **Tabula Sapiens** - Pan-tissue human cell atlas
   - URL: Search "Tabula Sapiens" on CellxGene
   - Contains: All major immune cell types from blood and tissues
   - Size: 500,000+ cells

2. **Human Cell Atlas - Blood & Immune**
   - Search: "PBMC human" or "immune cells human"
   - Contains: Comprehensive immune cell populations
   - Typical size: 50,000-200,000 cells per study

**Search Strategy:**
```
# On CellxGene portal:
1. Search: "PBMC human" OR "immune cells human"
2. Filter by:
   - Organism: Homo sapiens
   - Cell count: >50,000 cells
   - Assay: 10x scRNA-seq
3. Check for cell type annotations: CD14+ monocytes, CD4 T cells, etc.
4. Download as H5AD
```

### For Epithelial Cells

**Recommended Collections:**
1. **Tabula Sapiens - Epithelial Tissues**
   - Lung, intestine, kidney epithelial cells
   - Well-annotated cell types

2. **Cancer Atlas Studies**
   - Search: "epithelial human lung/breast/colon"
   - Contains: Normal and tumor epithelial cells
   - Good for cancer research applications

**Search Strategy:**
```
# On CellxGene portal:
1. Search: "epithelial human" + tissue of interest
2. Filter by:
   - Tissue: lung, breast, colon, etc.
   - Cell type contains: "epithelial"
   - Disease: normal (unless cancer research)
3. Download H5AD files
```

---

## Data Quality Requirements for New Cell Types

### Minimum Requirements
- **Cell count**: ≥3,000 cells of target cell type (optimal: 10,000+)
- **Gene count**: ≥15,000 genes detected
- **Quality metrics**:
  - Mean genes/cell: >800
  - Mean UMIs/cell: >2,000
  - Mitochondrial content: <25%

### Cell Type Annotation Requirements
- **Annotation field**: Must have `cell_type` or `celltype` column in `adata.obs`
- **Specificity**: Cell type labels should be specific (e.g., "CD14+ classical monocyte" not just "monocyte")
- **Validation**: Cell type markers should be expressed in filtered cells

### Data Format Requirements
- **File format**: H5AD (AnnData) required
- **Gene IDs**: Ensembl gene IDs or gene symbols (will be converted)
- **Count matrix**: Raw counts or normalized counts (will be re-normalized)

---

## Step-by-Step: Adding New Cell Types

See detailed guides:
- **ACTIVATE_NEW_CELL_TYPES.md** - Complete activation workflow
- **NEW_CELL_TYPES_DATA_GUIDE.md** - Data acquisition for 5 planned cell types
- **PREPROCESSING_PIPELINE_GUIDE.md** - Technical preprocessing details

### Quick Summary

1. **Find Data**: CellxGene Portal → Filter → Download H5AD
2. **Validate Quality**: Check cell count, gene count, annotations
3. **Preprocess**: QC filtering → HVG selection → Metacell aggregation
4. **ARACNe Processing**: Network inference (requires HPC, 12-14 hours)
5. **Cache Generation**: Convert TSV to pickle cache (`build_network_cache.py`)
6. **Integration**: System automatically recognizes new cell types

## Recommended Cell Types for Expansion

Based on research applications and data availability:

### High Priority (Good data availability)
1. **Hepatocytes** - Drug metabolism research
2. **Cardiomyocytes** - Cardiotoxicity studies
3. **Neurons** - Neurological disorders

### Medium Priority
4. **Fibroblasts** - Cancer stroma, wound healing
5. **Endothelial cells** - Vascular biology

### Future Consideration
- Astrocytes, oligodendrocytes (CNS research)
- Smooth muscle cells (cardiovascular)
- Adipocytes (metabolism research)
- Pancreatic beta cells (diabetes research)

---

## Data Citation and Attribution

### Primary Data Sources

When using RegNetAgents, please acknowledge:

1. **RegNetAgents Foundation Model**
   - Citation: Zhang, M., Swamy, V., Cassius, R., Dupire, L., Karaletsos, T., & Califano, A. (2025). "RegNetAgents: A Cellular Regulatory Network-Aware Transcriptomics Foundation Model." *bioRxiv*. doi:10.1101/2025.07.03.663009
   - GitHub: https://github.com/czi-ai/RegNetAgents
   - Virtual Cells Platform: https://virtualcellmodels.cziscience.com/model/regnetagents

2. **CellxGene Data Portal**
   - Citation: Megill, C., et al. (2021). "cellxgene: a performant, scalable exploration platform for high dimensional sparse matrices." *bioRxiv*. doi:10.1101/2021.04.05.438318
   - Portal: https://cellxgene.cziscience.com/

3. **ARACNe Algorithm**
   - Citation: Lachmann, A., et al. (2016). "ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information." *Bioinformatics*, 32(14), 2233-2235.

### RegNetAgents Software Citation

```
RegNetAgents: LLM-Powered Multi-Agent Framework for Gene Regulatory Network Analysis
Built on RegNetAgents pre-computed regulatory networks
```

---

## Troubleshooting Data Acquisition

### Issue: Cannot find specific cell type on CellxGene

**Solution:**
1. Try broader search terms (e.g., "monocyte" instead of "CD14 monocyte")
2. Look at multiple tissues (cell types appear in different tissues)
3. Check alternative atlases:
   - Human Cell Atlas: https://www.humancellatlas.org/
   - Single Cell Portal: https://singlecell.broadinstitute.org/

### Issue: Downloaded dataset has wrong cell type annotations

**Solution:**
1. Check `adata.obs` column names (might be `cell_type`, `celltype`, `cell_type_ontology`)
2. Use marker genes to re-annotate if necessary
3. Validate with known cell type markers using `sc.pl.dotplot()`

### Issue: Dataset is too large/small

**Solution:**
- **Too large**: Sample random subset of cells (aim for 10,000-50,000)
- **Too small**: Combine multiple datasets of same cell type using Scanpy integration

---

## Contact and Support

For questions about data sources or adding new cell types:
- **CellxGene Support**: https://cellxgene.cziscience.com/docs/
- **ARACNe/Network Processing**: Califano Lab documentation
- **Contact**: jbird@birdaisolutions.com

---

**Last Updated**: 2025-01-27
**Maintained by**: RegNetAgents Development Team
