# GREmLN Preprocessing Pipeline - Complete Guide

This document provides comprehensive documentation for the complete preprocessing pipeline used to generate gene regulatory networks for GREmLN analysis.

## Pipeline Overview

```
Raw scRNA-seq Data → Quality Control → Metacells → ARACNe Networks → Optimized Caches → GREmLN Analysis
      (H5AD)            (~2h)           (~2h)        (~12h HPC)        (~1h)           (Real-time)
```

## Stage 1: Data Acquisition and Preparation

### Input Requirements
- **Format**: H5AD (AnnData) or H5 files from single-cell RNA-seq experiments
- **Cell Count**: Minimum 3,000 cells of target cell type, optimal 10,000+
- **Gene Count**: Minimum 15,000 genes detected, optimal 20,000+
- **Quality**: Well-annotated with cell type labels and quality metrics

### Data Sources
1. **CellxGene Portal** (https://cellxgene.cziscience.com/)
   - Primary source for curated datasets
   - Search by tissue, cell type, and organism
   - Download as H5AD format

2. **Human Cell Atlas** (https://www.humancellatlas.org/)
   - Comprehensive reference datasets
   - Standardized processing pipelines
   - High-quality annotations

3. **GEO/SRA Databases**
   - Raw data from published studies
   - Requires additional processing
   - Use when specific studies needed

### Data Validation
```python
import scanpy as sc
import pandas as pd

# Load and inspect data
adata = sc.read_h5ad('dataset.h5ad')
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Cell types: {adata.obs['cell_type'].value_counts()}")

# Basic quality checks
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
print(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"Mean UMIs per cell: {adata.obs['total_counts'].mean():.0f}")
```

## Stage 2: Quality Control and Preprocessing

### Location
Processing scripts are located in `gremln_source/scGraphLLM/scGraphLLM/preprocess.py`

### Quality Control Parameters
Parameters are defined in `gremln_source/scGraphLLM/scGraphLLM/_globals.py`:

```python
# QC Thresholds (customizable)
MIN_GENES_PER_CELL = 200        # Minimum genes detected per cell
MAX_GENES_PER_CELL = 8000       # Maximum genes (filter potential doublets)
MIN_CELLS_PER_GENE = 3          # Minimum cells expressing each gene
MAX_MITOCHONDRIAL_PCT = 25      # Maximum mitochondrial gene percentage
MAX_RIBOSOMAL_PCT = 50          # Maximum ribosomal gene percentage

# Highly Variable Genes
N_TOP_GENES = 1024              # Top highly variable genes for ARACNe
```

### QC Processing Steps

#### Step 1: Cell Filtering
```python
# Filter cells based on QC metrics
sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
sc.pp.filter_cells(adata, max_genes=MAX_GENES_PER_CELL)

# Calculate mitochondrial and ribosomal gene percentages
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

# Filter cells with high mitochondrial content
adata = adata[adata.obs.pct_counts_mt < MAX_MITOCHONDRIAL_PCT, :]
```

#### Step 2: Gene Filtering
```python
# Filter genes expressed in minimum number of cells
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)

# Remove mitochondrial and ribosomal genes for network inference
adata = adata[:, ~adata.var['mt']]
adata = adata[:, ~adata.var['ribo']]
```

#### Step 3: Normalization
```python
# Save raw counts
adata.raw = adata

# Normalize to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
sc.pp.log1p(adata)
```

#### Step 4: Highly Variable Gene Selection
```python
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Keep only top N highly variable genes for ARACNe
hvg_genes = adata.var.highly_variable.sum()
if hvg_genes > N_TOP_GENES:
    # Select top N_TOP_GENES by dispersion
    top_genes = adata.var[adata.var.highly_variable].nlargest(N_TOP_GENES, 'dispersions')
    adata.var['highly_variable'] = False
    adata.var.loc[top_genes.index, 'highly_variable'] = True

# Filter to highly variable genes
adata = adata[:, adata.var.highly_variable]
```

### Expected Outputs
- **Filtered cell count**: 70-90% of original cells retained
- **Filtered gene count**: ~1,024 highly variable genes
- **Quality metrics**: Mean genes/cell >800, mean UMIs/cell >2000

## Stage 3: Metacell Generation

### Purpose
ARACNe requires sufficient statistical power for network inference. Metacells aggregate similar cells to:
- Increase effective sample size
- Reduce noise and technical variation
- Improve network inference quality

### Clustering and Metacell Creation
```python
# Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')

# Neighborhood graph construction
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Leiden clustering for metacell groups
sc.tl.leiden(adata, resolution=0.5)

# Create metacells (aggregate 5 cells per metacell)
CELLS_PER_METACELL = 5
metacells = []

for cluster in adata.obs['leiden'].unique():
    cluster_cells = adata[adata.obs['leiden'] == cluster]
    n_metacells = len(cluster_cells) // CELLS_PER_METACELL

    for i in range(n_metacells):
        start_idx = i * CELLS_PER_METACELL
        end_idx = start_idx + CELLS_PER_METACELL
        metacell_data = cluster_cells[start_idx:end_idx].X.mean(axis=0)
        metacells.append(metacell_data)
```

### Output Format
Metacells are saved in tab-delimited format for ARACNe input:
```
Gene1    Gene2    Gene3    ...    GeneN
0.45     1.23     0.78     ...    2.10
1.67     0.34     1.45     ...    0.89
...
```

## Stage 4: ARACNe Network Inference

### ARACNe Algorithm Details
- **Developer**: Califano Lab, Columbia University
- **Version**: ARACNe3 (modified for large datasets)
- **Method**: Mutual information-based network inference
- **Output**: Regulatory gene pairs with confidence scores

### ARACNe Processing

#### Hardware Requirements
- **CPU**: 16-32 cores recommended
- **Memory**: 64-128 GB RAM
- **Storage**: 50-100 GB temporary space
- **Time**: 12-14 hours per cell type

#### Command Structure
```bash
# ARACNe3 execution
./gremln_source/scripts/ARACNe3_app_release \
    --input metacells_expression.txt \
    --output network_output.tsv \
    --pvalue 1e-8 \
    --threads 16 \
    --mi_threshold 0.1 \
    --seed 12345
```

#### Key Parameters
- **P-value threshold**: 1e-8 (statistical significance)
- **MI threshold**: 0.1 (mutual information cutoff)
- **Bootstrap samples**: 100 (for confidence estimation)
- **Top regulators**: 1000 (maximum regulators to consider)

### ARACNe Output Format
```
regulator.values    target.values    mi.values    scc.values    count.values    log.p.values
ENSG00000134982     ENSG00000157110  0.245        0.156         85              -12.34
ENSG00000134982     ENSG00000111642  0.198        0.123         78              -10.67
...
```

### Network Consolidation
Multiple bootstrap runs are consolidated into a single network:
```python
import pandas as pd

# Load ARACNe bootstrap results
bootstrap_files = ['bootstrap_1.tsv', 'bootstrap_2.tsv', ..., 'bootstrap_100.tsv']
all_edges = []

for file in bootstrap_files:
    df = pd.read_csv(file, sep='\t')
    all_edges.append(df)

# Combine and filter edges by frequency
combined = pd.concat(all_edges)
edge_counts = combined.groupby(['regulator.values', 'target.values']).size()

# Keep edges appearing in >50% of bootstrap runs
BOOTSTRAP_THRESHOLD = 50
final_edges = edge_counts[edge_counts >= BOOTSTRAP_THRESHOLD]
```

## Stage 5: Network Cache Generation

### Purpose
Convert ARACNe TSV output to optimized pickle format for fast analysis.

### Cache Structure
```python
cache_data = {
    'regulator_targets': {
        'ENSG00000134982': ['ENSG00000157110', 'ENSG00000111642', ...],
        'ENSG00000168610': ['ENSG00000183742', 'ENSG00000172493', ...],
        ...
    },
    'target_regulators': {
        'ENSG00000157110': ['ENSG00000134982', 'ENSG00000168610', ...],
        'ENSG00000111642': ['ENSG00000134982', 'ENSG00000067369', ...],
        ...
    },
    'all_genes': ['ENSG00000134982', 'ENSG00000157110', ...],  # Sorted list
    'num_edges': 183247,
    'num_genes': 14628,
    'num_regulons': 1926,
    'created': '2024-09-13 10:07:00'
}
```

### Cache Generation Script
```bash
# Generate cache for specific cell type
python build_network_cache.py epithelial_cell

# Generate cache for all cell types
python build_network_cache.py --all

# Custom input/output directories
python build_network_cache.py hepatocytes \
    --input-dir /path/to/aracne/output \
    --output-dir /path/to/cache/output
```

### Performance Optimization
- **Lookup Time**: O(1) gene lookup using dictionaries
- **Memory Usage**: ~2-10 MB per cell type
- **Loading Time**: <1 second per cache file

## Stage 6: Integration and Validation

### System Integration
1. **Enum Updates**: Add cell types to `CellType` enum
2. **Cache Loading**: Automatic loading during workflow initialization
3. **Tool Registration**: MCP server tool definitions updated
4. **Documentation**: User-facing documentation updated

### Validation Checks
```python
# Network quality validation
def validate_network_quality(cache_data):
    checks = []

    # Minimum network size
    checks.append(cache_data['num_edges'] >= 100)

    # Gene coverage
    checks.append(cache_data['num_genes'] >= 500)

    # Regulatory complexity
    avg_targets = cache_data['num_edges'] / cache_data['num_regulons']
    checks.append(avg_targets >= 2.0)

    # Network connectivity
    isolated_genes = len([g for g, regs in cache_data['target_regulators'].items() if len(regs) == 0])
    connectivity = 1 - (isolated_genes / cache_data['num_genes'])
    checks.append(connectivity >= 0.7)

    return all(checks)
```

### Biological Validation
- **Pathway Enrichment**: Known pathways should be well-represented
- **Cell Type Markers**: Cell-type specific genes should be regulatory hubs
- **Literature Validation**: Compare with known regulatory relationships

## Computational Resource Planning

### Per Cell Type Requirements
- **QC Processing**: 2-4 hours, 16 GB RAM, 4 cores
- **Metacell Generation**: 1-2 hours, 8 GB RAM, 2 cores
- **ARACNe Inference**: 12-14 hours, 64 GB RAM, 16 cores
- **Cache Generation**: 30 minutes, 4 GB RAM, 1 core

### Batch Processing Strategy
```bash
# Parallel processing of multiple cell types
# Submit to SLURM cluster
for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    sbatch --job-name=${cell_type}_aracne \
           --cpus-per-task=16 \
           --mem=64G \
           --time=16:00:00 \
           process_cell_type.sh $cell_type
done
```

## Quality Metrics and Benchmarks

### Expected Network Characteristics
| Cell Type | Genes | Edges | Regulators | Avg Targets/Regulator |
|-----------|-------|-------|------------|----------------------|
| Hepatocytes | 800-1200 | 2000-8000 | 150-400 | 5-20 |
| Cardiomyocytes | 600-1000 | 1500-6000 | 100-300 | 8-25 |
| Neurons | 1000-1500 | 5000-15000 | 200-500 | 10-30 |
| Fibroblasts | 800-1200 | 3000-10000 | 150-400 | 8-25 |
| Endothelial | 700-1000 | 2000-7000 | 120-350 | 6-22 |

### Success Criteria
- **Technical**: All processing steps complete without errors
- **Statistical**: Network statistics within expected ranges
- **Biological**: Known pathways and relationships represented
- **Performance**: Analysis response time <2 seconds

## Troubleshooting Common Issues

### Data Quality Issues
**Problem**: Low cell count after QC filtering
- **Solution**: Relax QC thresholds, combine multiple datasets
- **Prevention**: Validate data quality before processing

**Problem**: Few highly variable genes identified
- **Solution**: Adjust HVG selection parameters, check normalization
- **Prevention**: Use well-processed datasets

### ARACNe Processing Issues
**Problem**: ARACNe runs out of memory
- **Solution**: Increase allocated RAM, reduce gene count
- **Prevention**: Monitor memory usage, use appropriate instance size

**Problem**: Very sparse networks generated
- **Solution**: Lower p-value threshold, increase bootstrap samples
- **Prevention**: Validate metacell quality

### Integration Issues
**Problem**: Cache loading fails
- **Solution**: Check file permissions, validate pickle format
- **Prevention**: Use provided cache generation scripts

## Performance Monitoring

### Processing Time Tracking
```python
import time
import logging

# Log processing stages
logger = logging.getLogger('preprocessing')

start_time = time.time()
# ... QC processing ...
qc_time = time.time() - start_time
logger.info(f"QC processing completed in {qc_time:.2f} seconds")

start_time = time.time()
# ... Metacell generation ...
metacell_time = time.time() - start_time
logger.info(f"Metacell generation completed in {metacell_time:.2f} seconds")
```

### Resource Usage Monitoring
```bash
# Monitor CPU and memory usage during ARACNe
htop -p $(pgrep ARACNe3_app_release)

# Track disk usage
du -sh /tmp/aracne_processing/

# Monitor job status in SLURM
squeue -u $USER --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
```

## Best Practices

### Data Management
1. **Version Control**: Track data sources and processing parameters
2. **Backup Strategy**: Keep copies of intermediate results
3. **Documentation**: Record all processing decisions and parameters
4. **Reproducibility**: Use fixed seeds and documented software versions

### Processing Optimization
1. **Batch Processing**: Process multiple cell types in parallel
2. **Resource Allocation**: Match compute resources to processing stage
3. **Checkpoint Strategy**: Save intermediate results for restart capability
4. **Monitoring**: Track progress and resource usage

### Quality Assurance
1. **Validation Scripts**: Automated quality checks at each stage
2. **Biological Validation**: Compare results with known biology
3. **Cross-Validation**: Compare with other network inference methods
4. **Documentation**: Record all quality metrics and decisions

---

**Next Steps**: Use this pipeline documentation along with the NEW_CELL_TYPES_DATA_GUIDE.md to process data for the 5 new cell types.