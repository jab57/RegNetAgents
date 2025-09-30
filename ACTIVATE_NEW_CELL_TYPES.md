# Step-by-Step Guide: Activating New Cell Types

This guide provides detailed, executable steps to activate the 5 new cell types in GREmLN. Follow these steps in order to complete the activation process.

## Prerequisites ✅

- [x] System updated with new cell types (completed)
- [x] Documentation and guides created (completed)
- [ ] HPC cluster access available
- [ ] CellxGene account/access
- [ ] 500GB+ storage space available

## Activation Roadmap

```
Phase 1: Setup     → Phase 2: Data      → Phase 3: Processing → Phase 4: Integration
  (2 hours)          (1-2 days)           (1-2 weeks)           (1 day)
```

---

## Phase 1: Environment Setup (Day 1, 2 hours)

### Step 1.1: Create Directory Structure
```bash
# Navigate to GREmLN directory
cd C:/Users/josea/OneDrive/Desktop/GREmLN

# Create directory structure for all new cell types
mkdir -p models/networks/hepatocytes/{raw_data,processed,aracne_output}
mkdir -p models/networks/cardiomyocytes/{raw_data,processed,aracne_output}
mkdir -p models/networks/neurons/{raw_data,processed,aracne_output}
mkdir -p models/networks/fibroblasts/{raw_data,processed,aracne_output}
mkdir -p models/networks/endothelial_cells/{raw_data,processed,aracne_output}

# Create processing logs directory
mkdir -p logs/preprocessing
mkdir -p logs/aracne
mkdir -p logs/validation
```

### Step 1.2: Verify System Integration
```bash
# Test that new cell types are recognized
python -c "
from gremln_langgraph_workflow import CellType
print('New cell types available:')
for ct in [CellType.HEPATOCYTES, CellType.CARDIOMYOCYTES, CellType.NEURONS, CellType.FIBROBLASTS, CellType.ENDOTHELIAL_CELLS]:
    print(f'  {ct.name} -> {ct.value}')
"

# Verify build script supports new types
python build_network_cache.py --help
```

### Step 1.3: Prepare Processing Environment
```bash
# Check Python environment
python -c "import scanpy, pandas, numpy, pickle; print('Required packages available')"

# Verify ARACNe availability
ls gremln_source/scripts/ARACNe3_app_release

# Test file permissions
touch models/networks/hepatocytes/test_file.txt && rm models/networks/hepatocytes/test_file.txt
```

---

## Phase 2: Data Acquisition (Days 1-3)

### Step 2.1: Hepatocytes Data
**Target: Liver parenchymal cells for drug metabolism research**

1. **Visit CellxGene Portal**: https://cellxgene.cziscience.com/
2. **Search**: "hepatocytes human liver"
3. **Recommended Dataset**: MacParland et al. 2018 liver atlas
4. **Filters**:
   - Organism: Homo sapiens
   - Cell count: >10,000 cells
   - Tissue: liver
   - Cell type: hepatocyte

```bash
# Download and place data
# After downloading hepatocytes_dataset.h5ad from CellxGene:
mv ~/Downloads/hepatocytes_dataset.h5ad models/networks/hepatocytes/raw_data/

# Verify data quality
python -c "
import scanpy as sc
adata = sc.read_h5ad('models/networks/hepatocytes/raw_data/hepatocytes_dataset.h5ad')
print(f'Hepatocytes dataset: {adata.n_obs} cells, {adata.n_vars} genes')
print(f'Cell types: {adata.obs.get(\"cell_type\", adata.obs.get(\"celltype\", \"Unknown\")).value_counts().head()}')
"
```

### Step 2.2: Cardiomyocytes Data
**Target: Heart muscle cells for cardiotoxicity studies**

1. **Search**: "cardiomyocytes human heart"
2. **Recommended Dataset**: Tucker et al. 2020 heart atlas
3. **Requirements**: >5,000 cardiomyocytes, cardiac markers present

```bash
# Download and verify
mv ~/Downloads/cardiomyocytes_dataset.h5ad models/networks/cardiomyocytes/raw_data/

python -c "
import scanpy as sc
adata = sc.read_h5ad('models/networks/cardiomyocytes/raw_data/cardiomyocytes_dataset.h5ad')
print(f'Cardiomyocytes dataset: {adata.n_obs} cells, {adata.n_vars} genes')
# Check for cardiac markers
cardiac_markers = ['TNNT2', 'MYH6', 'MYH7', 'ACTC1']
available_markers = [m for m in cardiac_markers if m in adata.var_names]
print(f'Cardiac markers available: {available_markers}')
"
```

### Step 2.3: Neurons Data
**Target: Brain cells for neurological disorder research**

1. **Search**: "neurons human brain cortex"
2. **Recommended Dataset**: Hodge et al. 2019 brain atlas
3. **Requirements**: >10,000 neurons, neuronal markers present

```bash
# Download and verify
mv ~/Downloads/neurons_dataset.h5ad models/networks/neurons/raw_data/

python -c "
import scanpy as sc
adata = sc.read_h5ad('models/networks/neurons/raw_data/neurons_dataset.h5ad')
print(f'Neurons dataset: {adata.n_obs} cells, {adata.n_vars} genes')
# Check for neuronal markers
neuronal_markers = ['MAP2', 'RBFOX3', 'SYP', 'NEUN']
available_markers = [m for m in neuronal_markers if m in adata.var_names]
print(f'Neuronal markers available: {available_markers}')
"
```

### Step 2.4: Fibroblasts Data
**Target: Connective tissue cells for cancer stroma research**

1. **Search**: "fibroblasts human skin lung"
2. **Recommended Dataset**: Reynolds et al. 2021 skin atlas
3. **Requirements**: >8,000 fibroblasts, fibroblast markers present

```bash
# Download and verify
mv ~/Downloads/fibroblasts_dataset.h5ad models/networks/fibroblasts/raw_data/

python -c "
import scanpy as sc
adata = sc.read_h5ad('models/networks/fibroblasts/raw_data/fibroblasts_dataset.h5ad')
print(f'Fibroblasts dataset: {adata.n_obs} cells, {adata.n_vars} genes')
# Check for fibroblast markers
fibroblast_markers = ['COL1A1', 'VIM', 'FN1', 'THY1']
available_markers = [m for m in fibroblast_markers if m in adata.var_names]
print(f'Fibroblast markers available: {available_markers}')
"
```

### Step 2.5: Endothelial Cells Data
**Target: Vascular cells for angiogenesis research**

1. **Search**: "endothelial human vascular"
2. **Recommended Dataset**: Travaglini et al. 2020 lung atlas
3. **Requirements**: >5,000 endothelial cells, endothelial markers present

```bash
# Download and verify
mv ~/Downloads/endothelial_dataset.h5ad models/networks/endothelial_cells/raw_data/

python -c "
import scanpy as sc
adata = sc.read_h5ad('models/networks/endothelial_cells/raw_data/endothelial_dataset.h5ad')
print(f'Endothelial dataset: {adata.n_obs} cells, {adata.n_vars} genes')
# Check for endothelial markers
endothelial_markers = ['PECAM1', 'VWF', 'CDH5', 'KDR']
available_markers = [m for m in endothelial_markers if m in adata.var_names]
print(f'Endothelial markers available: {available_markers}')
"
```

---

## Phase 3: Data Processing (Days 4-14)

### Step 3.1: Prepare HPC Processing Scripts

Create processing script for each cell type:

```bash
# Create processing script template
cat > process_cell_type.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=CELLTYPE_aracne
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/aracne/CELLTYPE_%j.out
#SBATCH --error=logs/aracne/CELLTYPE_%j.err

CELL_TYPE=$1
echo "Processing cell type: $CELL_TYPE"
echo "Start time: $(date)"

# Stage 1: Quality Control and Preprocessing (2-3 hours)
echo "Stage 1: QC and Preprocessing"
python gremln_source/scGraphLLM/scGraphLLM/preprocess.py \
    --cell-type $CELL_TYPE \
    --input-file models/networks/$CELL_TYPE/raw_data/${CELL_TYPE}_dataset.h5ad \
    --output-dir models/networks/$CELL_TYPE/processed/ \
    --min-genes 200 \
    --max-genes 8000 \
    --max-mito-pct 25 \
    --n-hvg 1024 \
    --cells-per-metacell 5

# Stage 2: ARACNe Network Inference (12-14 hours)
echo "Stage 2: ARACNe Network Inference"
./gremln_source/scripts/ARACNe3_app_release \
    --input models/networks/$CELL_TYPE/processed/metacells_expression.txt \
    --output models/networks/$CELL_TYPE/aracne_output/network.tsv \
    --pvalue 1e-8 \
    --threads 16 \
    --mi_threshold 0.1 \
    --bootstrap 100

echo "End time: $(date)"
echo "Processing completed for $CELL_TYPE"
EOF

chmod +x process_cell_type.sh
```

### Step 3.2: Submit Processing Jobs

```bash
# Submit all processing jobs in parallel
for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    echo "Submitting $cell_type processing job..."
    sbatch process_cell_type.sh $cell_type
    sleep 5  # Stagger submissions
done

# Monitor job progress
squeue -u $USER --format="%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R"
```

### Step 3.3: Monitor Processing Progress

```bash
# Create monitoring script
cat > monitor_progress.sh << 'EOF'
#!/bin/bash
echo "=== GREmLN Cell Type Processing Status ==="
echo "Timestamp: $(date)"
echo

for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    echo "--- $cell_type ---"

    # Check if raw data exists
    if [ -f "models/networks/$cell_type/raw_data/${cell_type}_dataset.h5ad" ]; then
        echo "  ✓ Raw data available"
    else
        echo "  ✗ Raw data missing"
    fi

    # Check if preprocessing completed
    if [ -f "models/networks/$cell_type/processed/metacells_expression.txt" ]; then
        echo "  ✓ Preprocessing completed"
    else
        echo "  ⏳ Preprocessing in progress/pending"
    fi

    # Check if ARACNe completed
    if [ -f "models/networks/$cell_type/aracne_output/network.tsv" ]; then
        echo "  ✓ ARACNe network generated"
    else
        echo "  ⏳ ARACNe in progress/pending"
    fi

    # Check if cache generated
    if [ -f "models/networks/$cell_type/network_index.pkl" ]; then
        echo "  ✓ Cache generated - READY FOR USE"
    else
        echo "  ⏳ Cache pending"
    fi

    echo
done
EOF

chmod +x monitor_progress.sh

# Run monitoring
./monitor_progress.sh
```

---

## Phase 4: Cache Generation and Integration (Day 15)

### Step 4.1: Generate Network Caches

Once ARACNe processing is complete for each cell type:

```bash
# Generate cache for each completed cell type
for cell_type in hepatocytes cardiomyocytes neurons fibroblasts endothelial_cells; do
    if [ -f "models/networks/$cell_type/aracne_output/network.tsv" ]; then
        echo "Generating cache for $cell_type..."

        # Copy ARACNe output to expected location
        cp models/networks/$cell_type/aracne_output/network.tsv models/networks/$cell_type/

        # Generate optimized cache
        python build_network_cache.py $cell_type

        echo "Cache generated for $cell_type"
    else
        echo "ARACNe output not ready for $cell_type"
    fi
done
```

### Step 4.2: Validate Network Quality

```bash
# Create validation script
cat > validate_networks.py << 'EOF'
#!/usr/bin/env python3
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
        checks.append(avg_targets >= 2.0)  # Average targets per regulator

        if all(checks):
            return True, f"✓ Quality: {data['num_genes']} genes, {data['num_edges']} edges, {data['num_regulons']} regulators"
        else:
            return False, f"✗ Quality issues: {data['num_genes']} genes, {data['num_edges']} edges, {data['num_regulons']} regulators"

    except Exception as e:
        return False, f"Error loading cache: {e}"

# Validate all cell types
cell_types = ['hepatocytes', 'cardiomyocytes', 'neurons', 'fibroblasts', 'endothelial_cells']

print("=== Network Quality Validation ===")
for cell_type in cell_types:
    valid, message = validate_cell_type(cell_type)
    status = "PASS" if valid else "FAIL"
    print(f"{cell_type:<20} {status:<6} {message}")
EOF

python validate_networks.py
```

### Step 4.3: Integration Testing

```bash
# Test analysis with new cell types
cat > test_integration.py << 'EOF'
#!/usr/bin/env python3
import asyncio
from gremln_langgraph_workflow import GREmLNWorkflow, CellType

async def test_cell_type_analysis():
    workflow = GREmLNWorkflow()

    # Test each new cell type with TP53
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
                in_network = analysis.get('in_network', False)

                print(f"  ✓ Analysis successful: {regulators} regulators, {targets} targets, in_network={in_network}")
            else:
                print(f"  ✗ Analysis failed: No valid result")

        except Exception as e:
            print(f"  ✗ Error: {e}")

asyncio.run(test_cell_type_analysis())
EOF

python test_integration.py
```

---

## Phase 5: Verification and Documentation

### Step 5.1: Final System Test

```bash
# Comprehensive system test
python -c "
from gremln_langgraph_workflow import GREmLNWorkflow, CellType, GREmLNCache

# Test cache loading
cache = GREmLNCache()
print(f'Total cell types loaded: {len(cache.network_indices)}')

# Count new vs existing
new_types = ['hepatocytes', 'cardiomyocytes', 'neurons', 'fibroblasts', 'endothelial_cells']
new_loaded = sum(1 for ct in new_types if ct in cache.network_indices)
print(f'New cell types activated: {new_loaded}/5')

# Test enumeration
print(f'Total cell types in enum: {len(CellType)}')
"
```

### Step 5.2: Performance Benchmarking

```bash
# Create performance test
cat > benchmark_performance.py << 'EOF'
#!/usr/bin/env python3
import time
import asyncio
from gremln_langgraph_workflow import GREmLNWorkflow, CellType

async def benchmark_analysis():
    workflow = GREmLNWorkflow()

    # Test genes
    test_genes = ["TP53", "BRCA1", "MYC", "EGFR", "KRAS"]

    # Test new cell types
    new_cell_types = [
        CellType.HEPATOCYTES,
        CellType.CARDIOMYOCYTES,
        CellType.NEURONS,
        CellType.FIBROBLASTS,
        CellType.ENDOTHELIAL_CELLS
    ]

    print("=== Performance Benchmark ===")

    for cell_type in new_cell_types:
        if cell_type.value in workflow.cache.network_indices:
            print(f"\nTesting {cell_type.name}:")

            for gene in test_genes:
                start_time = time.time()
                try:
                    result = await workflow.run_analysis(gene, cell_type, "basic")
                    end_time = time.time()

                    if result:
                        print(f"  {gene}: {(end_time - start_time)*1000:.1f}ms ✓")
                    else:
                        print(f"  {gene}: Failed ✗")

                except Exception as e:
                    print(f"  {gene}: Error - {e}")
        else:
            print(f"\n{cell_type.name}: Not available (data not processed)")

asyncio.run(benchmark_analysis())
EOF

python benchmark_performance.py
```

### Step 5.3: Update System Documentation

```bash
# Update the main pipeline documentation
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
```

---

## Success Metrics Checklist

### Technical Success ✅
- [ ] All 5 cell types processed through ARACNe pipeline
- [ ] Network cache files generated and validated
- [ ] System recognizes all new cell types
- [ ] Analysis workflow works with new cell types
- [ ] Performance benchmarks within acceptable range (<2s per analysis)

### Data Quality ✅
- [ ] Each network has >500 genes and >100 regulatory edges
- [ ] Cell-type specific markers are well-represented
- [ ] Known biological pathways are captured
- [ ] Network connectivity is appropriate (>70% connected genes)

### Integration Success ✅
- [ ] MCP server tools include new cell types in dropdowns
- [ ] Documentation updated with new cell types
- [ ] No errors in gene analysis with new cell types
- [ ] Cross-cell-type comparisons work correctly

---

## Troubleshooting Guide

### Common Issues and Solutions

**Issue**: ARACNe processing fails with memory error
```bash
# Solution: Increase memory allocation
sbatch --mem=128G process_cell_type.sh CELLTYPE
```

**Issue**: Network cache generation fails
```bash
# Solution: Check ARACNe output format
head -5 models/networks/CELLTYPE/aracne_output/network.tsv
# Should show header: regulator.values  target.values  mi.values  scc.values  count.values  log.p.values
```

**Issue**: Analysis fails with new cell type
```bash
# Solution: Verify cache was generated correctly
python -c "
import pickle
with open('models/networks/CELLTYPE/network_index.pkl', 'rb') as f:
    data = pickle.load(f)
print('Required keys:', list(data.keys()))
"
```

**Issue**: Poor network quality (too few edges)
```bash
# Solution: Lower ARACNe p-value threshold
./gremln_source/scripts/ARACNe3_app_release --pvalue 1e-6  # Instead of 1e-8
```

---

## Next Steps After Activation

1. **Research Applications**: Begin using new cell types for specific research questions
2. **Validation Studies**: Compare results with literature and known biology
3. **User Training**: Create examples and tutorials for each cell type
4. **Performance Optimization**: Monitor and optimize analysis speed
5. **Additional Cell Types**: Consider adding more specialized cell types based on research needs

---

**Estimated Timeline**: 2-3 weeks total (depending on HPC queue times)
**Success Rate**: Expected 4-5/5 cell types to activate successfully
**Support**: Refer to PREPROCESSING_PIPELINE_GUIDE.md for detailed technical information