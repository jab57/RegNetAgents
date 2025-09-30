#!/bin/bash

# GREmLN Cell Type Processing Monitor
# Monitors the progress of new cell type activation

set -e

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELL_TYPES=("hepatocytes" "cardiomyocytes" "neurons" "fibroblasts" "endothelial_cells")

# Colors for output (if terminal supports it)
if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    NC='\033[0m' # No Color
else
    RED=''
    GREEN=''
    YELLOW=''
    BLUE=''
    NC=''
fi

# Status symbols
CHECKMARK="✓"
CROSS="✗"
PENDING="⏳"
WARNING="⚠"

print_header() {
    echo -e "${BLUE}=== GREmLN Cell Type Processing Status ===${NC}"
    echo "Timestamp: $(date)"
    echo "Base directory: $BASE_DIR"
    echo
}

check_file_exists() {
    local file="$1"
    if [ -f "$file" ]; then
        return 0
    else
        return 1
    fi
}

get_file_size() {
    local file="$1"
    if [ -f "$file" ]; then
        if command -v du &> /dev/null; then
            du -h "$file" | cut -f1
        else
            echo "$(wc -c < "$file") bytes"
        fi
    else
        echo "N/A"
    fi
}

get_file_age() {
    local file="$1"
    if [ -f "$file" ]; then
        if command -v stat &> /dev/null; then
            stat -c %Y "$file" 2>/dev/null || stat -f %m "$file" 2>/dev/null || echo "0"
        else
            echo "0"
        fi
    else
        echo "0"
    fi
}

format_time_ago() {
    local timestamp="$1"
    local now=$(date +%s)
    local diff=$((now - timestamp))

    if [ $diff -lt 60 ]; then
        echo "${diff}s ago"
    elif [ $diff -lt 3600 ]; then
        echo "$((diff / 60))m ago"
    elif [ $diff -lt 86400 ]; then
        echo "$((diff / 3600))h ago"
    else
        echo "$((diff / 86400))d ago"
    fi
}

check_slurm_jobs() {
    if command -v squeue &> /dev/null; then
        echo -e "${BLUE}=== SLURM Job Status ===${NC}"

        # Get running jobs related to cell types
        for cell_type in "${CELL_TYPES[@]}"; do
            jobs=$(squeue -u $USER --name="${cell_type}_aracne" --format="%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R" --noheader 2>/dev/null || true)
            if [ -n "$jobs" ]; then
                echo -e "${YELLOW}$cell_type jobs:${NC}"
                echo "$jobs"
            fi
        done
        echo
    fi
}

check_cell_type_status() {
    local cell_type="$1"
    local cell_dir="$BASE_DIR/models/networks/$cell_type"

    echo -e "${BLUE}--- $cell_type ---${NC}"

    # Check raw data
    local raw_data_dir="$cell_dir/raw_data"
    local raw_files=($(find "$raw_data_dir" -name "*.h5ad" -o -name "*.h5" 2>/dev/null || true))

    if [ ${#raw_files[@]} -gt 0 ]; then
        echo -e "  ${GREEN}${CHECKMARK}${NC} Raw data available (${#raw_files[@]} files)"
        for file in "${raw_files[@]}"; do
            echo "    - $(basename "$file") ($(get_file_size "$file"))"
        done
    else
        echo -e "  ${RED}${CROSS}${NC} Raw data missing"
        echo -e "    ${YELLOW}→ Place H5AD files in $raw_data_dir${NC}"
    fi

    # Check preprocessing
    local metacells_file="$cell_dir/processed/metacells_expression.txt"
    if check_file_exists "$metacells_file"; then
        local file_age=$(get_file_age "$metacells_file")
        local time_ago=$(format_time_ago $file_age)
        echo -e "  ${GREEN}${CHECKMARK}${NC} Preprocessing completed ($(get_file_size "$metacells_file"), $time_ago)"

        # Check preprocessing stats
        local stats_file="$cell_dir/processed/preprocessing_stats.json"
        if check_file_exists "$stats_file"; then
            if command -v python3 &> /dev/null; then
                local stats=$(python3 -c "
import json
try:
    with open('$stats_file') as f:
        data = json.load(f)
    print(f\"{data.get('metacells', 'Unknown')} metacells, {data.get('filtered_genes', 'Unknown')} genes\")
except:
    print('Stats unavailable')
" 2>/dev/null || echo "Stats unavailable")
                echo "    - $stats"
            fi
        fi
    else
        # Check if preprocessing is in progress
        local log_file="$BASE_DIR/logs/preprocessing/${cell_type}.log"
        if check_file_exists "$log_file"; then
            local file_age=$(get_file_age "$log_file")
            local time_ago=$(format_time_ago $file_age)
            echo -e "  ${YELLOW}${PENDING}${NC} Preprocessing in progress ($time_ago)"
        else
            echo -e "  ${RED}${CROSS}${NC} Preprocessing not started"
            echo -e "    ${YELLOW}→ Run: sbatch process_cell_type.sh $cell_type${NC}"
        fi
    fi

    # Check ARACNe output
    local aracne_file="$cell_dir/aracne_output/network.tsv"
    if check_file_exists "$aracne_file"; then
        local file_age=$(get_file_age "$aracne_file")
        local time_ago=$(format_time_ago $file_age)
        local line_count=$(wc -l < "$aracne_file" 2>/dev/null || echo "0")
        echo -e "  ${GREEN}${CHECKMARK}${NC} ARACNe network generated ($(get_file_size "$aracne_file"), $line_count edges, $time_ago)"
    else
        # Check if ARACNe is in progress
        local aracne_log="$BASE_DIR/logs/aracne/${cell_type}.log"
        if check_file_exists "$aracne_log"; then
            local file_age=$(get_file_age "$aracne_log")
            local time_ago=$(format_time_ago $file_age)
            echo -e "  ${YELLOW}${PENDING}${NC} ARACNe in progress ($time_ago)"
        else
            echo -e "  ${RED}${CROSS}${NC} ARACNe not started"
        fi
    fi

    # Check cache generation
    local cache_file="$cell_dir/network_index.pkl"
    if check_file_exists "$cache_file"; then
        local file_age=$(get_file_age "$cache_file")
        local time_ago=$(format_time_ago $file_age)
        echo -e "  ${GREEN}${CHECKMARK}${NC} Cache generated ($(get_file_size "$cache_file"), $time_ago)"
        echo -e "  ${GREEN}${CHECKMARK} READY FOR USE${NC}"

        # Test cache quality
        if command -v python3 &> /dev/null; then
            local cache_stats=$(python3 -c "
import pickle
try:
    with open('$cache_file', 'rb') as f:
        data = pickle.load(f)
    print(f\"{data.get('num_genes', 'Unknown')} genes, {data.get('num_edges', 'Unknown')} edges, {data.get('num_regulons', 'Unknown')} regulators\")
except:
    print('Cache stats unavailable')
" 2>/dev/null || echo "Cache stats unavailable")
            echo "    - Network: $cache_stats"
        fi
    else
        echo -e "  ${RED}${CROSS}${NC} Cache not generated"
        if check_file_exists "$aracne_file"; then
            echo -e "    ${YELLOW}→ Run: python build_network_cache.py $cell_type${NC}"
        fi
    fi

    echo
}

show_overall_summary() {
    echo -e "${BLUE}=== Overall Summary ===${NC}"

    local total_cell_types=${#CELL_TYPES[@]}
    local ready_count=0
    local in_progress_count=0
    local not_started_count=0

    for cell_type in "${CELL_TYPES[@]}"; do
        local cache_file="$BASE_DIR/models/networks/$cell_type/network_index.pkl"
        local log_file="$BASE_DIR/logs/preprocessing/${cell_type}.log"

        if check_file_exists "$cache_file"; then
            ready_count=$((ready_count + 1))
        elif check_file_exists "$log_file"; then
            in_progress_count=$((in_progress_count + 1))
        else
            not_started_count=$((not_started_count + 1))
        fi
    done

    echo "Total new cell types: $total_cell_types"
    echo -e "Ready for use: ${GREEN}$ready_count${NC}"
    echo -e "In progress: ${YELLOW}$in_progress_count${NC}"
    echo -e "Not started: ${RED}$not_started_count${NC}"

    if [ $ready_count -eq $total_cell_types ]; then
        echo -e "\n${GREEN}${CHECKMARK} All cell types are ready!${NC}"
    elif [ $in_progress_count -gt 0 ]; then
        echo -e "\n${YELLOW}${PENDING} Processing in progress...${NC}"
    else
        echo -e "\n${RED}${CROSS} No processing started yet${NC}"
    fi
}

show_next_steps() {
    echo -e "\n${BLUE}=== Next Steps ===${NC}"

    local any_ready=false
    local any_missing_data=false

    for cell_type in "${CELL_TYPES[@]}"; do
        local cache_file="$BASE_DIR/models/networks/$cell_type/network_index.pkl"
        local raw_data_dir="$BASE_DIR/models/networks/$cell_type/raw_data"
        local raw_files=($(find "$raw_data_dir" -name "*.h5ad" -o -name "*.h5" 2>/dev/null || true))

        if check_file_exists "$cache_file"; then
            any_ready=true
        elif [ ${#raw_files[@]} -eq 0 ]; then
            any_missing_data=true
        fi
    done

    if [ "$any_missing_data" = true ]; then
        echo "1. Acquire data for missing cell types:"
        echo "   - Follow NEW_CELL_TYPES_DATA_GUIDE.md"
        echo "   - Download H5AD files from CellxGene portal"
        echo "   - Place files in models/networks/CELLTYPE/raw_data/"
    fi

    echo "2. Start processing for cell types with data:"
    echo "   sbatch process_cell_type.sh CELLTYPE"

    echo "3. Monitor progress:"
    echo "   ./monitor_progress.sh"

    if [ "$any_ready" = true ]; then
        echo "4. Test activated cell types:"
        echo "   python -c \"from gremln_langgraph_workflow import GREmLNWorkflow, CellType; import asyncio; print('Testing...')\""
    fi

    echo "5. For detailed guidance:"
    echo "   - Read ACTIVATE_NEW_CELL_TYPES.md"
    echo "   - Read PREPROCESSING_PIPELINE_GUIDE.md"
}

# Main execution
main() {
    print_header

    # Check SLURM jobs if available
    check_slurm_jobs

    # Check each cell type
    for cell_type in "${CELL_TYPES[@]}"; do
        check_cell_type_status "$cell_type"
    done

    # Overall summary
    show_overall_summary

    # Next steps
    show_next_steps
}

# Run with watch option if requested
if [ "$1" = "--watch" ] || [ "$1" = "-w" ]; then
    if command -v watch &> /dev/null; then
        watch -n 30 "$0"
    else
        echo "Watch command not available. Running once..."
        main
    fi
else
    main
fi