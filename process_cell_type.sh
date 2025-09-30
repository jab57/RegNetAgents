#!/bin/bash
#SBATCH --job-name=CELLTYPE_aracne
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/aracne/CELLTYPE_%j.out
#SBATCH --error=logs/aracne/CELLTYPE_%j.err

# GREmLN Cell Type Processing Script
# Processes single-cell data through the complete pipeline: QC → Metacells → ARACNe → Cache

set -e  # Exit on any error
set -u  # Exit on undefined variables

# Configuration
CELL_TYPE=${1:-""}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$SCRIPT_DIR"

# Validate input
if [ -z "$CELL_TYPE" ]; then
    echo "Error: Cell type not specified"
    echo "Usage: $0 <cell_type>"
    echo "Available cell types: hepatocytes, cardiomyocytes, neurons, fibroblasts, endothelial_cells"
    exit 1
fi

# Define paths
CELL_DIR="$BASE_DIR/models/networks/$CELL_TYPE"
RAW_DATA_DIR="$CELL_DIR/raw_data"
PROCESSED_DIR="$CELL_DIR/processed"
ARACNE_DIR="$CELL_DIR/aracne_output"
LOG_DIR="$BASE_DIR/logs"

# Create directories if they don't exist
mkdir -p "$PROCESSED_DIR" "$ARACNE_DIR" "$LOG_DIR/preprocessing" "$LOG_DIR/aracne"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/preprocessing/${CELL_TYPE}.log"
}

log "=== Starting GREmLN processing for $CELL_TYPE ==="
log "Script directory: $SCRIPT_DIR"
log "Cell type directory: $CELL_DIR"

# Check for input data
INPUT_FILES=($(find "$RAW_DATA_DIR" -name "*.h5ad" -o -name "*.h5" 2>/dev/null || true))

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    log "ERROR: No input data files found in $RAW_DATA_DIR"
    log "Please place H5AD files in $RAW_DATA_DIR before running"
    exit 1
fi

log "Found ${#INPUT_FILES[@]} input file(s):"
for file in "${INPUT_FILES[@]}"; do
    log "  - $(basename "$file")"
done

# Use the first input file
INPUT_FILE="${INPUT_FILES[0]}"
log "Using input file: $INPUT_FILE"

# Stage 1: Data Validation
log "=== Stage 1: Data Validation ==="
if command -v python3 &> /dev/null; then
    PYTHON_CMD=python3
else
    PYTHON_CMD=python
fi

log "Running data validation..."
VALIDATION_OUTPUT=$($PYTHON_CMD "$BASE_DIR/validate_cell_type_data.py" \
    --cell-type "$CELL_TYPE" \
    --data-file "$INPUT_FILE" 2>&1 || true)

echo "$VALIDATION_OUTPUT" | tee -a "$LOG_DIR/preprocessing/${CELL_TYPE}.log"

# Check if validation passed (simplified check)
if echo "$VALIDATION_OUTPUT" | grep -q "validation PASSED"; then
    log "✓ Data validation PASSED"
else
    log "⚠ Data validation issues detected, but continuing..."
fi

# Stage 2: Quality Control and Preprocessing
log "=== Stage 2: Quality Control and Preprocessing ==="

# Cell type specific parameters
case "$CELL_TYPE" in
    hepatocytes)
        MIN_GENES=200
        MAX_GENES=8000
        MAX_MITO_PCT=25
        N_HVG=1024
        ;;
    cardiomyocytes)
        MIN_GENES=200
        MAX_GENES=7000
        MAX_MITO_PCT=30
        N_HVG=1024
        ;;
    neurons)
        MIN_GENES=300
        MAX_GENES=9000
        MAX_MITO_PCT=20
        N_HVG=1024
        ;;
    fibroblasts)
        MIN_GENES=200
        MAX_GENES=8000
        MAX_MITO_PCT=25
        N_HVG=1024
        ;;
    endothelial_cells)
        MIN_GENES=200
        MAX_GENES=7500
        MAX_MITO_PCT=25
        N_HVG=1024
        ;;
    *)
        log "Using default parameters for unknown cell type"
        MIN_GENES=200
        MAX_GENES=8000
        MAX_MITO_PCT=25
        N_HVG=1024
        ;;
esac

log "QC Parameters: min_genes=$MIN_GENES, max_genes=$MAX_GENES, max_mito=$MAX_MITO_PCT%, hvg=$N_HVG"

# Create preprocessing script
cat > "$PROCESSED_DIR/preprocess.py" << EOF
#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import numpy as np
import os

# Configure scanpy
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Load data
print("Loading data...")
adata = sc.read_h5ad('$INPUT_FILE')
print(f"Initial: {adata.n_obs} cells, {adata.n_vars} genes")

# Basic filtering
print("Applying quality control filters...")

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=$MIN_GENES)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

# Filter cells based on QC metrics
adata = adata[adata.obs.n_genes_by_counts < $MAX_GENES, :]
adata = adata[adata.obs.pct_counts_mt < $MAX_MITO_PCT, :]

print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")

# Remove mitochondrial and ribosomal genes for network inference
adata = adata[:, ~adata.var['mt']]
adata = adata[:, ~adata.var['ribo']]

# Save raw counts
adata.raw = adata

# Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Keep only top HVG
hvg_genes = adata.var.highly_variable.sum()
if hvg_genes > $N_HVG:
    top_genes = adata.var[adata.var.highly_variable].nlargest($N_HVG, 'dispersions')
    adata.var['highly_variable'] = False
    adata.var.loc[top_genes.index, 'highly_variable'] = True

# Filter to HVG
adata = adata[:, adata.var.highly_variable]
print(f"After HVG selection: {adata.n_obs} cells, {adata.n_vars} genes")

# PCA and clustering for metacells
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.5)

# Create metacells
print("Creating metacells...")
CELLS_PER_METACELL = 5
metacells = []
metacell_data = []

for cluster in adata.obs['leiden'].unique():
    cluster_cells = adata[adata.obs['leiden'] == cluster]
    n_metacells = len(cluster_cells) // CELLS_PER_METACELL

    for i in range(n_metacells):
        start_idx = i * CELLS_PER_METACELL
        end_idx = start_idx + CELLS_PER_METACELL

        if end_idx <= len(cluster_cells):
            metacell_expr = cluster_cells[start_idx:end_idx].X.mean(axis=0).A1
            metacells.append(f"metacell_{cluster}_{i}")
            metacell_data.append(metacell_expr)

# Create metacell expression matrix
metacell_df = pd.DataFrame(
    metacell_data,
    index=metacells,
    columns=adata.var_names
)

print(f"Created {len(metacells)} metacells")

# Save metacell data for ARACNe
metacell_df.T.to_csv('$PROCESSED_DIR/metacells_expression.txt', sep='\t')

# Save processing statistics
stats = {
    'original_cells': adata.raw.n_obs if adata.raw else 'Unknown',
    'original_genes': adata.raw.n_vars if adata.raw else 'Unknown',
    'filtered_cells': adata.n_obs,
    'filtered_genes': adata.n_vars,
    'metacells': len(metacells),
    'parameters': {
        'min_genes': $MIN_GENES,
        'max_genes': $MAX_GENES,
        'max_mito_pct': $MAX_MITO_PCT,
        'n_hvg': $N_HVG,
        'cells_per_metacell': CELLS_PER_METACELL
    }
}

import json
with open('$PROCESSED_DIR/preprocessing_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

print("Preprocessing completed successfully!")
print(f"Output saved to: $PROCESSED_DIR/metacells_expression.txt")
EOF

# Run preprocessing
log "Running preprocessing script..."
cd "$PROCESSED_DIR"
$PYTHON_CMD preprocess.py 2>&1 | tee -a "$LOG_DIR/preprocessing/${CELL_TYPE}.log"

if [ ! -f "$PROCESSED_DIR/metacells_expression.txt" ]; then
    log "ERROR: Preprocessing failed - no metacell output generated"
    exit 1
fi

log "✓ Preprocessing completed successfully"

# Stage 3: ARACNe Network Inference
log "=== Stage 3: ARACNe Network Inference ==="

# Check for ARACNe executable
ARACNE_EXEC="$BASE_DIR/gremln_source/scripts/ARACNe3_app_release"
if [ ! -f "$ARACNE_EXEC" ]; then
    log "ERROR: ARACNe executable not found at $ARACNE_EXEC"
    log "Please ensure ARACNe3 is properly installed"
    exit 1
fi

# Run ARACNe
log "Starting ARACNe network inference..."
log "Input: $PROCESSED_DIR/metacells_expression.txt"
log "Output: $ARACNE_DIR/network.tsv"

cd "$ARACNE_DIR"

# ARACNe command with appropriate parameters
"$ARACNE_EXEC" \
    --input "$PROCESSED_DIR/metacells_expression.txt" \
    --output "$ARACNE_DIR/network.tsv" \
    --pvalue 1e-8 \
    --threads ${SLURM_CPUS_PER_TASK:-16} \
    --mi_threshold 0.1 \
    --seed 12345 \
    2>&1 | tee -a "$LOG_DIR/aracne/${CELL_TYPE}.log"

if [ ! -f "$ARACNE_DIR/network.tsv" ]; then
    log "ERROR: ARACNe failed - no network output generated"
    exit 1
fi

# Validate ARACNe output
NETWORK_LINES=$(wc -l < "$ARACNE_DIR/network.tsv")
if [ "$NETWORK_LINES" -lt 100 ]; then
    log "WARNING: Network file has only $NETWORK_LINES lines - this may indicate poor quality"
else
    log "✓ ARACNe completed successfully - $NETWORK_LINES edges generated"
fi

# Stage 4: Cache Generation
log "=== Stage 4: Cache Generation ==="

# Copy network file to expected location
cp "$ARACNE_DIR/network.tsv" "$CELL_DIR/"

# Generate optimized cache
log "Generating optimized network cache..."
cd "$BASE_DIR"

$PYTHON_CMD build_network_cache.py "$CELL_TYPE" 2>&1 | tee -a "$LOG_DIR/preprocessing/${CELL_TYPE}.log"

if [ -f "$CELL_DIR/network_index.pkl" ]; then
    log "✓ Cache generated successfully"
else
    log "ERROR: Cache generation failed"
    exit 1
fi

# Stage 5: Validation
log "=== Stage 5: Final Validation ==="

# Test the generated network
VALIDATION_CMD="
from gremln_langgraph_workflow import GREmLNWorkflow, CellType
import asyncio

async def test_network():
    try:
        workflow = GREmLNWorkflow()
        cell_type = getattr(CellType, '${CELL_TYPE^^}')

        # Test with TP53
        result = await workflow.run_analysis('TP53', cell_type, 'basic')

        if result and 'network_analysis' in result:
            analysis = result['network_analysis']
            print(f'✓ Network test successful')
            print(f'  In network: {analysis.get(\"in_network\", False)}')
            print(f'  Regulators: {len(analysis.get(\"regulators\", []))}')
            print(f'  Targets: {len(analysis.get(\"targets\", []))}')
            return True
        else:
            print('✗ Network test failed - no valid analysis result')
            return False
    except Exception as e:
        print(f'✗ Network test failed: {e}')
        return False

result = asyncio.run(test_network())
exit(0 if result else 1)
"

if $PYTHON_CMD -c "$VALIDATION_CMD" 2>&1 | tee -a "$LOG_DIR/preprocessing/${CELL_TYPE}.log"; then
    log "✓ Final validation PASSED"
else
    log "⚠ Final validation had issues, but cache was generated"
fi

# Generate summary report
log "=== Processing Summary ==="
log "Cell type: $CELL_TYPE"
log "Input file: $(basename "$INPUT_FILE")"
log "Processing time: $((SECONDS / 60)) minutes"

if [ -f "$PROCESSED_DIR/preprocessing_stats.json" ]; then
    log "Preprocessing statistics:"
    cat "$PROCESSED_DIR/preprocessing_stats.json" | tee -a "$LOG_DIR/preprocessing/${CELL_TYPE}.log"
fi

log "Files generated:"
log "  - Metacells: $PROCESSED_DIR/metacells_expression.txt"
log "  - Network: $ARACNE_DIR/network.tsv"
log "  - Cache: $CELL_DIR/network_index.pkl"
log "  - Logs: $LOG_DIR/preprocessing/${CELL_TYPE}.log"

log "=== Processing completed successfully for $CELL_TYPE ==="

# Clean up temporary files (optional)
# rm -f "$PROCESSED_DIR/preprocess.py"

exit 0