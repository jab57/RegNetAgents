# Complete Guide: Adding New Cell Types to RegNetAgents

> **Note:** This guide is for users who have received access to the RegNetAgents software framework. The code is available from the corresponding author upon reasonable request for academic research purposes.

**⚠️ STATUS: FUTURE IMPLEMENTATION GUIDE**

This comprehensive guide documents the complete process for adding new cell types to RegNetAgents, from data acquisition through system integration. Currently focused on 5 specific cell types (hepatocytes, cardiomyocytes, neurons, fibroblasts, endothelial cells), but the process applies to any cell type.

**Currently Available Cell Types (10):** cd14_monocytes, cd16_monocytes, cd20_b_cells, cd4_t_cells, cd8_t_cells, erythrocytes, nk_cells, nkt_cells, epithelial_cell, monocyte-derived_dendritic_cells

**See [DATA_SOURCES.md](DATA_SOURCES.md) for complete documentation of current cell types.**

---

## Overview

### Future Cell Types (Examples)

1. **Hepatocytes** - Liver parenchymal cells (drug metabolism)
2. **Cardiomyocytes** - Heart muscle cells (cardiotoxicity)
3. **Neurons** - Brain/nervous system cells (neurological disorders)
4. **Fibroblasts** - Connective tissue cells (cancer stroma, wound healing)
5. **Endothelial Cells** - Vascular lining cells (angiogenesis)

### Complete Pipeline

```
CellxGene Download → QC & Filtering → Metacells → ARACNe → Cache Generation → Integration
    (1-2 days)         (2-4 hours)     (1-2 hours)  (12-14h)     (30 min)       (1 day)
```

### Prerequisites ✅

- [ ] HPC cluster access (64-128 GB RAM, 16+ cores)
- [ ] CellxGene Portal access
- [ ] 500GB+ storage space
- [ ] ARACNe software (`regnetagents_source/scripts/ARACNe3_app_release`)
- [ ] Python environment with scanpy, networkx

### Estimated Timeline

- **Sequential (1 cell type)**: ~2 weeks per cell type
- **Parallel (5 cell types on HPC)**: ~2-3 weeks total with queue time
- **Active work time**: ~3-5 days (rest is compute time)

---

## Phase 1: Environment Setup (Day 1, 2 hours)

### Step 1.1: Create Directory Structure

```bash
# Navigate to RegNetAgents directory
cd /path/to/RegNetAgents

# Create directory structure for all new cell types
for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    mkdir -p models/networks/$cell_type/{raw_data,processed,aracne_output}
done

# Create processing logs directory
mkdir -p logs/{preprocessing,aracne,validation}
```

### Step 1.2: Verify System Integration

```bash
# Test that new cell types are recognized
python -c "
from regnetagents_langgraph_workflow import CellType
print('New cell types available:')
for ct in [CellType.HEPATOCYTES, CellType.CARDIOMYOCYTES, CellType.NEURONS,
           CellType.FIBROBLASTS, CellType.ENDOTHELIAL_CELLS]:
    print(f'  {ct.name} -> {ct.value}')
"

# Verify build script supports new types
python scripts/build_network_cache.py --list-types
```

---

## Phase 2: Data Acquisition (Days 1-3)

### General Data Requirements

**Minimum Requirements:**
- Cell count: ≥3,000 cells of target type (optimal: 10,000+)
- Gene count: ≥15,000 genes detected
- Quality: Mean genes/cell >800, mean UMIs/cell >2,000
- Format: H5AD (AnnData) required
- Annotations: Must have `cell_type` or `celltype` column

**Where to Find Data:**
- **Primary source**: [CellxGene Data Portal](https://cellxgene.cziscience.com/)
- **Alternative**: [Human Cell Atlas](https://www.humancellatlas.org/)
- **Specific studies**: GEO/SRA databases

### Cell Type-Specific Data Sources

#### 1. Hepatocytes (Liver Cells)

**Primary Applications:** Drug metabolism, liver toxicity, hepatic disease research

**Recommended Datasets:**
- **CellxGene Portal**: Search "hepatocytes human liver"
- **Human Liver Cell Atlas**: Comprehensive liver cell type data
- **Dataset Example**: MacParland et al. 2018 liver atlas
- **Expected Cells**: 5,000-50,000 hepatocytes
- **Quality Metrics**: >1,000 genes/cell, <20% mitochondrial genes
- **Key Markers**: ALB, CYP3A4, CYP2E1, APOB

**Research Value:**
- Drug metabolism pathways (CYP450 enzymes)
- Hepatotoxicity screening
- Liver disease mechanisms (NAFLD, hepatitis, cirrhosis)
- Pharmacokinetics modeling

**Download Instructions:**
```bash
# On CellxGene portal:
# 1. Search: "hepatocytes human liver"
# 2. Filter:
#    - Organism: Homo sapiens
#    - Cell count: >10,000 cells
#    - Tissue: liver
#    - Cell type: hepatocyte
# 3. Download as H5AD
# 4. Save to: models/networks/hepatocytes/raw_data/
```

#### 2. Cardiomyocytes (Heart Muscle Cells)

**Primary Applications:** Cardiotoxicity, heart disease, cardiac drug development

**Recommended Datasets:**
- **CellxGene Portal**: Search "cardiomyocytes human heart"
- **Human Heart Cell Atlas**: Cardiac cell type reference
- **Dataset Example**: Tucker et al. 2020 heart atlas
- **Expected Cells**: 3,000-30,000 cardiomyocytes
- **Quality Metrics**: >800 genes/cell, cardiac markers present
- **Key Markers**: TNNT2, MYH6, MYH7, ACTC1

**Research Value:**
- Cardiotoxicity assessment (drug safety)
- Heart failure mechanisms
- Arrhythmia pathways
- Cardiac drug development

#### 3. Neurons (Brain/Nervous System Cells)

**Primary Applications:** Neurological disorders, CNS drug targets, neurotoxicity

**Recommended Datasets:**
- **CellxGene Portal**: Search "neurons human brain cortex"
- **Allen Brain Cell Types Database**: Comprehensive neuronal data
- **Dataset Example**: Hodge et al. 2019 brain atlas
- **Expected Cells**: 10,000-100,000 neurons
- **Quality Metrics**: >1,500 genes/cell, neuronal markers present
- **Key Markers**: MAP2, RBFOX3, SYP, NEUN

**Research Value:**
- Alzheimer's disease pathways
- Depression/anxiety mechanisms
- Neurotoxicity screening
- CNS drug development

#### 4. Fibroblasts (Connective Tissue Cells)

**Primary Applications:** Wound healing, cancer stroma, tissue repair

**Recommended Datasets:**
- **CellxGene Portal**: Search "fibroblasts human skin lung"
- **Tabula Sapiens**: Multi-organ fibroblast data
- **Dataset Example**: Reynolds et al. 2021 skin atlas
- **Expected Cells**: 5,000-40,000 fibroblasts
- **Quality Metrics**: >1,200 genes/cell, fibroblast markers present
- **Key Markers**: COL1A1, VIM, FN1, THY1

**Research Value:**
- Cancer-associated fibroblasts (CAFs) in tumor microenvironment
- Wound healing mechanisms
- Fibrosis pathways (pulmonary, liver)
- Tissue engineering applications

#### 5. Endothelial Cells (Vascular Lining Cells)

**Primary Applications:** Vascular biology, angiogenesis, vascular toxicity

**Recommended Datasets:**
- **CellxGene Portal**: Search "endothelial human vascular"
- **Human Cell Atlas Vascular**: Vascular cell reference
- **Dataset Example**: Travaglini et al. 2020 lung atlas
- **Expected Cells**: 3,000-25,000 endothelial cells
- **Quality Metrics**: >1,000 genes/cell, endothelial markers present
- **Key Markers**: PECAM1, VWF, CDH5, KDR

**Research Value:**
- Angiogenesis pathways
- Vascular toxicity assessment
- Atherosclerosis mechanisms
- Anti-angiogenic drug development

### Data Validation After Download

```python
import scanpy as sc

# Load and inspect data
adata = sc.read_h5ad('models/networks/hepatocytes/raw_data/hepatocytes_dataset.h5ad')
print(f'Cells: {adata.n_obs}, Genes: {adata.n_vars}')
print(f'Cell types: {adata.obs["cell_type"].value_counts()}')

# Check for cell type markers (example: hepatocytes)
hepatocyte_markers = ['ALB', 'CYP3A4', 'CYP2E1', 'APOB']
available_markers = [m for m in hepatocyte_markers if m in adata.var_names]
print(f'Hepatocyte markers available: {available_markers}')

# Basic quality check
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
print(f'Mean genes per cell: {adata.obs["n_genes_by_counts"].mean():.0f}')
print(f'Mean UMIs per cell: {adata.obs["total_counts"].mean():.0f}')
```

---

## Phase 3: Data Processing (Days 4-14)

**See [END_TO_END_DATA_PIPELINE.md](END_TO_END_DATA_PIPELINE.md) for complete technical details.**

### Step 3.1: Quality Control & Preprocessing

```bash
# QC parameters (from scGraphLLM/preprocess.py):
# - MIN_GENES_PER_CELL = 200
# - MAX_GENES_PER_CELL = 8000
# - MAX_MITOCHONDRIAL_PCT = 25
# - TOP_HVG = 1,024

# Process each cell type (example for hepatocytes)
python gremln_source/scGraphLLM/preprocess.py \
    --cell-type hepatocytes \
    --input-file models/networks/hepatocytes/raw_data/hepatocytes_dataset.h5ad \
    --output-dir models/networks/hepatocytes/processed/ \
    --min-genes 200 \
    --max-genes 8000 \
    --max-mito-pct 25 \
    --n-hvg 1024 \
    --cells-per-metacell 5
```

### Step 3.2: ARACNe Network Inference (HPC Required)

```bash
# Create HPC job script
cat > process_cell_type.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=CELLTYPE_aracne
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/aracne/CELLTYPE_%j.out

CELL_TYPE=$1
echo "Processing cell type: $CELL_TYPE"

# ARACNe Network Inference (12-14 hours)
./regnetagents_source/scripts/ARACNe3_app_release \
    --input models/networks/$CELL_TYPE/processed/metacells_expression.txt \
    --output models/networks/$CELL_TYPE/aracne_output/network.tsv \
    --pvalue 1e-8 \
    --threads 16 \
    --mi_threshold 0.1 \
    --bootstrap 100 \
    --seed 12345

echo "Processing completed for $CELL_TYPE"
EOF

chmod +x process_cell_type.sh

# Submit all jobs in parallel
for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    sbatch process_cell_type.sh $cell_type
    sleep 5
done

# Monitor progress
squeue -u $USER
```

### Step 3.3: Monitor Processing

```bash
# Create monitoring script
cat > monitor_progress.sh << 'EOF'
#!/bin/bash
echo "=== RegNetAgents Processing Status ==="
echo "Timestamp: $(date)"
echo

for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    echo "--- $cell_type ---"

    # Check stages
    [ -f "models/networks/$cell_type/raw_data/${cell_type}_dataset.h5ad" ] && echo "  ✓ Raw data" || echo "  ✗ Raw data missing"
    [ -f "models/networks/$cell_type/processed/metacells_expression.txt" ] && echo "  ✓ Preprocessing" || echo "  ⏳ Preprocessing pending"
    [ -f "models/networks/$cell_type/aracne_output/network.tsv" ] && echo "  ✓ ARACNe" || echo "  ⏳ ARACNe pending"
    [ -f "models/networks/$cell_type/network_index.pkl" ] && echo "  ✓ Cache - READY" || echo "  ⏳ Cache pending"
    echo
done
EOF

chmod +x monitor_progress.sh
./monitor_progress.sh
```

---

## Phase 4: Cache Generation & Integration (Day 15)

### Step 4.1: Generate Network Caches

```bash
# Generate cache for each completed cell type
for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    if [ -f "models/networks/$cell_type/aracne_output/network.tsv" ]; then
        echo "Generating cache for $cell_type..."

        # Copy ARACNe output to expected location
        cp models/networks/$cell_type/aracne_output/network.tsv models/networks/$cell_type/

        # Generate optimized cache with pre-computed PageRank
        python scripts/build_network_cache.py $cell_type

        echo "Cache generated for $cell_type"
    fi
done
```

### Step 4.2: Validate Network Quality

```python
#!/usr/bin/env python3
# validate_networks.py

import pickle
import os

def validate_cell_type(cell_type):
    cache_file = f"models/networks/{cell_type}/network_index.pkl"

    if not os.path.exists(cache_file):
        return False, "Cache file missing"

    try:
        with open(cache_file, 'rb') as f:
            data = pickle.load(f)

        # Quality checks
        checks = []
        checks.append(data['num_edges'] >= 100)  # Minimum edges
        checks.append(data['num_genes'] >= 500)  # Minimum genes
        checks.append(data['num_regulons'] >= 50)  # Minimum regulators

        # Network connectivity
        avg_targets = data['num_edges'] / data['num_regulons'] if data['num_regulons'] > 0 else 0
        checks.append(avg_targets >= 2.0)

        if all(checks):
            return True, f"✓ {data['num_genes']} genes, {data['num_edges']} edges, {data['num_regulons']} regulators"
        else:
            return False, f"✗ Quality issues"

    except Exception as e:
        return False, f"Error: {e}"

# Validate all new cell types
cell_types = ['hepatocytes', 'cardiomyocytes', 'neurons', 'fibroblasts', 'endothelial_cells']

print("=== Network Quality Validation ===")
for cell_type in cell_types:
    valid, message = validate_cell_type(cell_type)
    status = "PASS" if valid else "FAIL"
    print(f"{cell_type:<20} {status:<6} {message}")
```

### Step 4.3: Integration Testing

```python
#!/usr/bin/env python3
# test_integration.py

import asyncio
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow, CellType

async def test_cell_type_analysis():
    workflow = RegNetAgentsWorkflow()

    test_cases = [
        (CellType.HEPATOCYTES, "hepatocytes"),
        (CellType.CARDIOMYOCYTES, "cardiomyocytes"),
        (CellType.NEURONS, "neurons"),
        (CellType.FIBROBLASTS, "fibroblasts"),
        (CellType.ENDOTHELIAL_CELLS, "endothelial_cells")
    ]

    print("=== Testing Gene Analysis with New Cell Types ===")

    for cell_type_enum, cell_type_name in test_cases:
        try:
            print(f"\nTesting {cell_type_name}...")
            result = await workflow.run_analysis(
                gene="TP53",
                cell_type=cell_type_enum,
                analysis_depth="basic"
            )

            if result and 'network_analysis' in result:
                analysis = result['network_analysis']
                regulators = len(analysis.get('regulators', []))
                targets = len(analysis.get('targets', []))

                print(f"  ✓ Success: {regulators} regulators, {targets} targets")
            else:
                print(f"  ✗ Failed: No valid result")

        except Exception as e:
            print(f"  ✗ Error: {e}")

asyncio.run(test_cell_type_analysis())
```

---

## Phase 5: Verification & Documentation

### Step 5.1: Final System Test

```bash
# Comprehensive check
python -c "
from regnetagents_langgraph_workflow import RegNetAgentsWorkflow, CellType, RegNetAgentsCache

cache = RegNetAgentsCache()
print(f'Total cell types loaded: {len(cache.network_indices)}')

new_types = ['hepatocytes', 'cardiomyocytes', 'neurons', 'fibroblasts', 'endothelial_cells']
new_loaded = sum(1 for ct in new_types if ct in cache.network_indices)
print(f'New cell types activated: {new_loaded}/5')

print(f'Total cell types in enum: {len(CellType)}')
"
```

### Step 5.2: Update Documentation

```bash
# Record successful activation
echo "
## ✅ NEW CELL TYPES ACTIVATED ($(date +%Y-%m-%d))

Successfully activated 5 additional cell types:
- **Hepatocytes**: Liver drug metabolism research
- **Cardiomyocytes**: Cardiac toxicity studies
- **Neurons**: Neurological disorder research
- **Fibroblasts**: Cancer stroma analysis
- **Endothelial Cells**: Vascular biology research

Total cell types now supported: 15

Processing completed on: $(date)
" >> CELL_TYPE_ACTIVATION_LOG.md

# Update DATA_SOURCES.md to include new cell types
```

---

## Success Criteria Checklist

### Technical Success ✅
- [ ] All 5 cell types processed through ARACNe pipeline
- [ ] Network cache files generated and validated
- [ ] System recognizes all new cell types
- [ ] Analysis workflow works with new cell types
- [ ] Performance benchmarks within acceptable range (<2s per analysis)

### Data Quality ✅
- [ ] Each network has >500 genes and >100 regulatory edges
- [ ] Cell-type specific markers are well-represented in networks
- [ ] Known biological pathways are captured
- [ ] Network connectivity is appropriate (>70% connected genes)

### Integration Success ✅
- [ ] MCP server tools include new cell types in dropdowns
- [ ] Documentation updated with new cell types
- [ ] No errors in gene analysis with new cell types
- [ ] Cross-cell-type comparisons work correctly

---

## Troubleshooting Guide

### Data Quality Issues

**Problem**: Low cell count after QC filtering
- **Solution**: Relax QC thresholds (MAX_MITOCHONDRIAL_PCT = 30 instead of 25)
- **Solution**: Combine multiple datasets from same tissue

**Problem**: Few highly variable genes identified
- **Solution**: Adjust HVG selection parameters
- **Solution**: Check normalization was applied correctly

### ARACNe Processing Issues

**Problem**: ARACNe runs out of memory
- **Solution**: Increase allocated RAM (128 GB instead of 64 GB)
- **Solution**: Reduce gene count or metacell count

**Problem**: Very sparse networks generated (<100 edges)
- **Solution**: Lower p-value threshold: `--pvalue 1e-6` instead of `1e-8`
- **Solution**: Increase bootstrap samples: `--bootstrap 200`

**Problem**: ARACNe fails to start
- **Solution**: Check input file format (genes as rows, metacells as columns)
- **Solution**: Verify no empty rows or columns in expression matrix

### Cache Generation Issues

**Problem**: `build_network_cache.py` fails
- **Solution**: Verify TSV has correct header: `regulator.values  target.values  ...`
- **Solution**: Ensure all gene IDs are Ensembl format (ENSG...)

**Problem**: Cache validation fails
- **Solution**: Check pickle file isn't corrupted
- **Solution**: Re-generate cache from TSV file

**Problem**: PageRank calculation fails
- **Solution**: Check network has sufficient connectivity
- **Solution**: Fall back to on-demand PageRank calculation

### Integration Issues

**Problem**: New cell type not recognized by system
- **Solution**: Verify cell type added to `CellType` enum in workflow
- **Solution**: Check cache file is in correct location: `models/networks/{cell_type}/network_index.pkl`

**Problem**: Analysis times out for new cell type
- **Solution**: Check network size isn't too large (>500K edges may be slow)
- **Solution**: Verify PageRank was pre-computed in cache

---

## Performance Benchmarks

### Expected Timeline (Per Cell Type)

| Stage | Time | Resources |
|-------|------|-----------|
| Data Download | 10-30 min | Internet bandwidth |
| QC & Preprocessing | 2-4 hours | 16 GB RAM, 4 cores |
| Metacell Generation | 1-2 hours | 8 GB RAM, 2 cores |
| ARACNe Inference | 12-14 hours | 64 GB RAM, 16 cores |
| Cache Generation | 10-30 min | 4 GB RAM, 1 core |
| Validation & Testing | 1-2 hours | 8 GB RAM |
| **Total** | **~18-24 hours** | **HPC cluster** |

### Parallel Processing (All 5 Cell Types)

- **Sequential**: ~5-6 days
- **Parallel on HPC**: ~1-2 days + queue time
- **Recommended**: Submit all ARACNe jobs in parallel

---

## Expected Network Characteristics

| Cell Type | Genes (HVGs) | Metacells | Edges (Expected) | Regulators (Expected) |
|-----------|--------------|-----------|------------------|----------------------|
| Hepatocytes | 1,024 | 4,000-10,000 | 2,000-8,000 | 200-500 |
| Cardiomyocytes | 1,024 | 2,000-6,000 | 1,500-6,000 | 150-400 |
| Neurons | 1,024 | 8,000-20,000 | 5,000-15,000 | 300-700 |
| Fibroblasts | 1,024 | 5,000-12,000 | 3,000-10,000 | 200-500 |
| Endothelial | 1,024 | 3,000-8,000 | 2,000-7,000 | 180-450 |

---

## Next Steps After Activation

1. **Research Applications**: Begin using new cell types for specific research questions
2. **Validation Studies**: Compare results with literature and known biology
3. **User Documentation**: Create examples and tutorials for each cell type
4. **Performance Monitoring**: Track analysis speed and optimize if needed
5. **Publication Preparation**: Document methods for papers using new cell types
6. **Community Sharing**: Consider making networks publicly available

---

## Additional Resources

- **[DATA_SOURCES.md](DATA_SOURCES.md)** - Current cell type documentation
- **[END_TO_END_DATA_PIPELINE.md](END_TO_END_DATA_PIPELINE.md)** - Complete technical pipeline
- **[REGNETAGENTS_Analysis_Pipeline.md](REGNETAGENTS_Analysis_Pipeline.md)** - System architecture
- **CellxGene Portal**: https://cellxgene.cziscience.com/
- **ARACNe Documentation**: Califano Lab resources
- **Scanpy Documentation**: https://scanpy.readthedocs.io/

---

**Document Version**: 2.0 (Consolidated)
**Last Updated**: 2025-01-27
**Maintained by**: RegNetAgents Development Team
