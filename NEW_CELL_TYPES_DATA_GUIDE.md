# Data Acquisition Guide for New Cell Types

This guide provides detailed instructions for obtaining and processing single-cell RNA-seq data for the 5 new priority cell types added to GREmLN.

## Overview

**New Cell Types Added:**
1. **Hepatocytes** - Liver parenchymal cells
2. **Cardiomyocytes** - Heart muscle cells
3. **Neurons** - Brain and nervous system cells
4. **Fibroblasts** - Connective tissue cells
5. **Endothelial Cells** - Blood vessel lining cells

**Processing Pipeline:** CellxGene Data → QC/Metacells → ARACNe Networks → TSV Files → System Integration

## Data Sources and Recommendations

### 1. Hepatocytes
**Primary Applications:** Drug metabolism, liver toxicity, hepatic disease research

**Recommended Datasets:**
- **CellxGene Portal**: Search "hepatocytes human liver"
- **Human Liver Cell Atlas**: Comprehensive liver cell type data
- **Dataset ID Example**: `liver_macparland_2018` (MacParland et al., 2018)
- **Expected Cells**: 5,000-50,000 hepatocytes
- **Quality Metrics**: >1,000 genes/cell, <20% mitochondrial genes

**Research Value:**
- Drug metabolism pathways (CYP450 enzymes)
- Hepatotoxicity screening
- Liver disease mechanisms
- Pharmacokinetics modeling

### 2. Cardiomyocytes
**Primary Applications:** Cardiotoxicity, heart disease, cardiac drug development

**Recommended Datasets:**
- **CellxGene Portal**: Search "cardiomyocytes human heart"
- **Human Heart Cell Atlas**: Cardiac cell type reference
- **Dataset ID Example**: `heart_tucker_2020` (Tucker et al., 2020)
- **Expected Cells**: 3,000-30,000 cardiomyocytes
- **Quality Metrics**: >800 genes/cell, cardiac markers (TNNT2, MYH6)

**Research Value:**
- Cardiotoxicity assessment
- Heart failure mechanisms
- Arrhythmia pathways
- Cardiac drug safety

### 3. Neurons
**Primary Applications:** Neurological disorders, CNS drug targets, neurotoxicity

**Recommended Datasets:**
- **CellxGene Portal**: Search "neurons human brain cortex"
- **Allen Brain Cell Types Database**: Comprehensive neuronal data
- **Dataset ID Example**: `brain_hodge_2019` (Hodge et al., 2019)
- **Expected Cells**: 10,000-100,000 neurons
- **Quality Metrics**: >1,500 genes/cell, neuronal markers (MAP2, RBFOX3)

**Research Value:**
- Alzheimer's disease pathways
- Depression/anxiety mechanisms
- Neurotoxicity screening
- CNS drug development

### 4. Fibroblasts
**Primary Applications:** Wound healing, cancer stroma, tissue repair

**Recommended Datasets:**
- **CellxGene Portal**: Search "fibroblasts human skin lung"
- **Tabula Sapiens**: Multi-organ fibroblast data
- **Dataset ID Example**: `skin_reynolds_2021` (Reynolds et al., 2021)
- **Expected Cells**: 5,000-40,000 fibroblasts
- **Quality Metrics**: >1,200 genes/cell, fibroblast markers (COL1A1, VIM)

**Research Value:**
- Cancer-associated fibroblasts (CAFs)
- Wound healing mechanisms
- Fibrosis pathways
- Tissue engineering

### 5. Endothelial Cells
**Primary Applications:** Vascular biology, angiogenesis, vascular toxicity

**Recommended Datasets:**
- **CellxGene Portal**: Search "endothelial human vascular"
- **Human Cell Atlas Vascular**: Vascular cell reference
- **Dataset ID Example**: `lung_travaglini_2020` (Travaglini et al., 2020)
- **Expected Cells**: 3,000-25,000 endothelial cells
- **Quality Metrics**: >1,000 genes/cell, endothelial markers (PECAM1, VWF)

**Research Value:**
- Angiogenesis pathways
- Vascular toxicity
- Atherosclerosis mechanisms
- Anti-angiogenic drug development

## Step-by-Step Data Acquisition Process

### Phase 1: Data Discovery and Download
```bash
# 1. Visit CellxGene Portal (https://cellxgene.cziscience.com/)
# 2. Search for specific cell type: "hepatocytes human liver"
# 3. Filter by:
#    - Organism: Homo sapiens
#    - Cell count: >5,000 cells
#    - Tissue: Relevant organ
# 4. Download as H5AD format
# 5. Verify data quality and cell annotations
```

### Phase 2: Data Preparation
```bash
# Create directories for new cell types
mkdir -p models/networks/hepatocytes
mkdir -p models/networks/cardiomyocytes
mkdir -p models/networks/neurons
mkdir -p models/networks/fibroblasts
mkdir -p models/networks/endothelial_cells

# Place downloaded H5AD files in appropriate directories
cp hepatocytes_data.h5ad models/networks/hepatocytes/
cp cardiomyocytes_data.h5ad models/networks/cardiomyocytes/
# etc.
```

### Phase 3: ARACNe Processing
```bash
# Follow existing pipeline from gremln_source/scripts/README.md
# For each cell type:

# 1. Preprocessing (Quality Control + Metacells)
python gremln_source/scripts/preprocess_cellxgene.py \
    --cell-type hepatocytes \
    --input-data models/networks/hepatocytes/hepatocytes_data.h5ad \
    --output-dir models/networks/hepatocytes/

# 2. ARACNe Network Generation (~12 hours per cell type)
./gremln_source/scripts/ARACNe3_app_release \
    --input models/networks/hepatocytes/metacells.txt \
    --output models/networks/hepatocytes/network.tsv \
    --threads 16

# 3. Cache Generation (automated)
python build_network_cache.py hepatocytes
```

### Phase 4: Validation
```bash
# Test each new cell type
python -c "
from gremln_langgraph_workflow import GREmLNWorkflow, CellType
import asyncio

async def test_cell_type():
    workflow = GREmLNWorkflow()
    result = await workflow.run_analysis('TP53', CellType.HEPATOCYTES, 'basic')
    print(f'Hepatocytes network loaded: {result is not None}')

asyncio.run(test_cell_type())
"
```

## Quality Control Criteria

### Minimum Requirements Per Cell Type
- **Cell Count**: >3,000 cells of target type
- **Gene Count**: >800 genes per cell average
- **Mitochondrial %**: <25% mitochondrial genes
- **Cell Type Markers**: Express known markers for cell type
- **Batch Effects**: Minimal batch effects across samples

### ARACNe Network Quality
- **Network Size**: >100 regulatory edges
- **Gene Coverage**: >500 genes in final network
- **Connectivity**: Well-connected network (not fragmented)
- **Biological Relevance**: Known pathways represented

## Computational Requirements

### HPC Cluster Specifications
- **CPU Cores**: 16-32 cores per ARACNe job
- **Memory**: 64-128 GB RAM per job
- **Storage**: 100 GB temporary space per cell type
- **Time**: ~12-14 hours per cell type
- **Queue**: Submit all 5 cell types in parallel

### Local Requirements
- **Cache Generation**: 4-8 GB RAM, 1-2 hours total
- **Final Storage**: ~50 MB per cell type
- **Validation**: 8 GB RAM, 30 minutes per cell type

## Expected Timeline

### Sequential Processing (Conservative)
- **Week 1**: Data discovery and download (hepatocytes, cardiomyocytes)
- **Week 2**: ARACNe processing for 2 cell types (~24 hours HPC time)
- **Week 3**: Data acquisition for remaining 3 cell types
- **Week 4**: ARACNe processing for remaining 3 cell types (~36 hours HPC time)
- **Week 5**: Cache generation, validation, and integration testing

### Parallel Processing (Optimal)
- **Day 1-2**: Data discovery and download for all 5 cell types
- **Day 3-4**: Submit all 5 ARACNe jobs in parallel (~14 hours HPC time)
- **Day 5**: Cache generation and validation
- **Total**: 1 week with HPC cluster access

## Research Applications by Cell Type

### Hepatocytes
- **Drug Development**: ADME-Tox profiling, drug metabolism
- **Disease Research**: NAFLD, hepatitis, liver cancer
- **Toxicology**: Hepatotoxicity screening, safety assessment

### Cardiomyocytes
- **Drug Safety**: Cardiotoxicity, QT prolongation, arrhythmias
- **Disease Research**: Heart failure, cardiomyopathy, ischemia
- **Therapeutics**: Cardiac regeneration, heart disease drugs

### Neurons
- **Neurological Disorders**: Alzheimer's, Parkinson's, ALS
- **Psychiatric Disorders**: Depression, schizophrenia, autism
- **CNS Drugs**: Neuropharmacology, blood-brain barrier

### Fibroblasts
- **Cancer Research**: Tumor microenvironment, CAF biology
- **Wound Healing**: Tissue repair, regenerative medicine
- **Fibrotic Diseases**: Pulmonary fibrosis, liver cirrhosis

### Endothelial Cells
- **Vascular Biology**: Angiogenesis, endothelial dysfunction
- **Cardiovascular Disease**: Atherosclerosis, hypertension
- **Cancer**: Tumor vasculature, anti-angiogenic therapy

## Troubleshooting Common Issues

### Data Quality Issues
- **Low Cell Count**: Combine multiple datasets from same tissue
- **Poor Gene Expression**: Filter cells more stringently
- **Batch Effects**: Use integration methods (Harmony, Seurat)

### ARACNe Processing Issues
- **Memory Errors**: Increase allocated RAM or reduce gene count
- **Long Runtime**: Use more CPU cores, reduce cell count if needed
- **Empty Networks**: Check cell type annotations and gene filtering

### Integration Issues
- **Cache Loading Errors**: Verify pickle file format and permissions
- **Network Validation Fails**: Check gene ID formats and network structure
- **Missing Cell Types**: Ensure enum and supported types lists match

## Success Metrics

### Technical Success
- [ ] All 5 cell types successfully processed through ARACNe
- [ ] Cache files generated and validated for each cell type
- [ ] System accepts all new cell types in analysis
- [ ] No errors in workflow execution with new cell types

### Scientific Success
- [ ] Networks contain biologically relevant pathways
- [ ] Cell type-specific markers are well-represented
- [ ] Known disease pathways are captured in networks
- [ ] Cross-cell-type comparisons show expected differences

## Next Steps After Implementation

1. **Validation Studies**: Compare networks with known biology
2. **Benchmark Analysis**: Test with well-characterized genes
3. **User Documentation**: Create usage examples for each cell type
4. **Research Applications**: Identify high-impact use cases
5. **Publication Preparation**: Document novel insights from multi-cell-type analysis

---

**Contact Information**: For questions about data acquisition or processing issues, refer to the original GREmLN paper methods or contact the development team.