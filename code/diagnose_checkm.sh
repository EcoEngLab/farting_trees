#!/bin/bash

# CheckM Diagnostic Script
# This script helps diagnose CheckM issues step by step

echo "=== CheckM Diagnostic Script ==="
echo "Date: $(date)"
echo ""

# Check working directory
echo "Current directory: $(pwd)"
echo ""

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
RESULTS_DIR="${WORK_DIR}/results/04_binning"

echo "=== 1. Checking file paths ==="
echo "Work directory: $WORK_DIR"
echo "Results directory: $RESULTS_DIR"
echo "Results directory exists: $([ -d "$RESULTS_DIR" ] && echo "YES" || echo "NO")"
echo ""

# Initialize conda
echo "=== 2. Checking conda setup ==="
source $(conda info --base)/etc/profile.d/conda.sh
echo "Conda initialized successfully"
echo ""

# Check conda environments
echo "=== 3. Checking conda environments ==="
echo "Available environments:"
conda env list | grep -E "(checkm|metagen|binning)" || echo "No relevant environments found"
echo ""

# Test CheckM environment
echo "=== 4. Testing CheckM environment ==="
if conda activate checkm_env 2>/dev/null; then
    echo "✓ checkm_env activated successfully"
    echo "Python version: $(python --version 2>/dev/null || echo 'Python not available')"
    echo ""
    
    # Check CheckM availability
    echo "=== 5. Checking CheckM availability ==="
    if command -v checkm &>/dev/null; then
        echo "✓ CheckM command found: $(which checkm)"
        
        # Check CheckM version
        echo "CheckM version:"
        checkm --version 2>/dev/null || echo "Version check failed"
        echo ""
        
        # Check CheckM data
        echo "=== 6. Checking CheckM database ==="
        echo "Checking CheckM data location..."
        checkm data setRoot --help 2>/dev/null | head -5 || echo "Data command failed"
        
        # Try to get current data root
        echo "Current CheckM data root:"
        python -c "
import checkm
try:
    from checkm.defaultValues import DefaultValues
    print('Data root:', DefaultValues.DEFAULT_DATA_ROOT)
except:
    print('Could not determine data root')
" 2>/dev/null || echo "Python check failed"
        
    else
        echo "✗ CheckM command not found"
    fi
    echo ""
    
    # Check dependencies
    echo "=== 7. Checking dependencies ==="
    echo "Checking for prodigal..."
    if command -v prodigal &>/dev/null; then
        echo "✓ prodigal found: $(which prodigal)"
    else
        echo "✗ prodigal not found"
    fi
    
    echo "Checking for pplacer..."
    if command -v pplacer &>/dev/null; then
        echo "✓ pplacer found: $(which pplacer)"
    else
        echo "✗ pplacer not found"
    fi
    
    echo "Checking for hmmer..."
    if command -v hmmsearch &>/dev/null; then
        echo "✓ hmmer found: $(which hmmsearch)"
    else
        echo "✗ hmmer not found"
    fi
    echo ""
    
    # Check input files
    echo "=== 8. Checking input files ==="
    echo "Checking for bins to analyze..."
    
    if [ -d "${RESULTS_DIR}/checkm/all_bins" ]; then
        bin_count=$(ls ${RESULTS_DIR}/checkm/all_bins/*.fa 2>/dev/null | wc -l)
        echo "Bins in all_bins directory: $bin_count"
        
        if [ $bin_count -gt 0 ]; then
            echo "Sample bin files:"
            ls ${RESULTS_DIR}/checkm/all_bins/*.fa 2>/dev/null | head -3
            
            # Check bin file sizes
            echo "Bin file sizes:"
            ls -lh ${RESULTS_DIR}/checkm/all_bins/*.fa 2>/dev/null | head -3 | awk '{print $5, $9}'
        fi
    else
        echo "all_bins directory not found, checking source bins..."
        echo "MetaBAT2 bins: $(ls ${RESULTS_DIR}/metabat2/bin*.fa 2>/dev/null | wc -l || echo 0)"
        echo "MaxBin2 bins: $(ls ${RESULTS_DIR}/maxbin2/bin*.fasta 2>/dev/null | wc -l || echo 0)"
        echo "VAMB bins: $(ls ${RESULTS_DIR}/vamb/bins/*.fna 2>/dev/null | wc -l || echo 0)"
    fi
    echo ""
    
else
    echo "✗ Failed to activate checkm_env"
    echo "Available environments:"
    conda env list
fi

echo "=== 9. System resources ==="
echo "Memory usage:"
free -h | grep -E "(Mem|Swap)" || echo 'Memory info not available'
echo "Available disk space:"
df -h . || echo 'Disk info not available'
echo ""

echo "=== 10. Suggested actions ==="
echo "Common CheckM issues and solutions:"
echo "1. Missing database: checkm data setRoot <path>"
echo "2. Missing dependencies: conda install -c bioconda prodigal pplacer hmmer"
echo "3. Memory issues: reduce thread count or process fewer bins"
echo "4. Permission issues: check file/directory permissions"
echo ""

echo "=== Diagnostic complete ==="
