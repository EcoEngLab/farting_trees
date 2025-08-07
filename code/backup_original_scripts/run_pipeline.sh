#!/bin/bash

# Master script to run the complete metagenomic analysis pipeline
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "========================================"
echo "Farting Trees Metagenomic Analysis Pipeline"
echo "========================================"
echo "Project directory: $PROJECT_DIR"
echo "Script directory: $SCRIPT_DIR"
echo ""

# Function to run a script with error handling
run_script() {
    local script_name=$1
    local script_path="${SCRIPT_DIR}/${script_name}"
    
    echo "----------------------------------------"
    echo "Running: $script_name"
    echo "----------------------------------------"
    
    if [ -f "$script_path" ]; then
        cd "$SCRIPT_DIR"
        bash "$script_path"
        echo "✓ $script_name completed successfully"
    else
        echo "✗ Script not found: $script_path"
        exit 1
    fi
    echo ""
}

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first"
    exit 1
fi

# Parse command line arguments
STEP=${1:-"all"}

case $STEP in
    "setup"|"1")
        echo "Setting up conda environments..."
        run_script "01_setup_environment.sh"
        ;;
    "qc"|"2")
        echo "Running quality control..."
        run_script "02_quality_control.sh"
        ;;
    "host"|"2b")
        echo "Running host removal..."
        run_script "02b_host_removal.sh"
        ;;
    "assembly"|"3")
        echo "Running assembly..."
        run_script "03_assembly.sh"
        ;;
    "binning"|"4")
        echo "Running binning..."
        run_script "04_binning.sh"
        ;;
    "annotation"|"5")
        echo "Running annotation..."
        run_script "05_annotation.sh"
        ;;
    "target"|"6")
        echo "Running targeted analysis..."
        run_script "06_target_analysis.sh"
        ;;
    "all")
        echo "Running complete pipeline..."
        run_script "01_setup_environment.sh"
        run_script "02_quality_control.sh"
        run_script "02b_host_removal.sh"
        run_script "03_assembly.sh"
        run_script "04_binning.sh"
        run_script "05_annotation.sh"
        run_script "06_target_analysis.sh"
        ;;
    *)
        echo "Usage: $0 [step]"
        echo ""
        echo "Available steps:"
        echo "  setup|1    - Setup conda environments"
        echo "  qc|2       - Quality control and preprocessing"
        echo "  host|2b    - Host (plant) genome removal"
        echo "  assembly|3 - Metagenomic assembly"
        echo "  binning|4  - Metagenomic binning"
        echo "  annotation|5 - Annotation and taxonomic classification"
        echo "  target|6   - Targeted analysis for methanogens/methanotrophs"
        echo "  all        - Run complete pipeline (default)"
        echo ""
        echo "Example:"
        echo "  $0 all          # Run complete pipeline"
        echo "  $0 qc           # Run only quality control"
        echo "  $0 host         # Run only host removal"
        echo "  $0 assembly     # Run only assembly"
        exit 1
        ;;
esac

echo "========================================"
echo "Pipeline step(s) completed successfully!"
echo "========================================"

# Generate final summary if running complete pipeline
if [ "$STEP" = "all" ]; then
    echo ""
    echo "📊 Analysis Summary:"
    echo "  - Raw data: data/"
    echo "  - Quality control: results/01_qc/"
    echo "  - Host removal: results/02_host_removal/"
    echo "  - Assembly: results/03_assembly/"
    echo "  - Binning: results/04_binning/"
    echo "  - Annotation: results/05_annotation/"
    echo "  - Targeted analysis: results/06_taxonomy/"
    echo ""
    echo "📋 Key reports to check:"
    echo "  - MultiQC report: results/01_qc/multiqc/multiqc_report.html"
    echo "  - Host removal stats: results/02_host_removal/host_removal_summary.txt"
    echo "  - Assembly stats: results/03_assembly/quast/"
    echo "  - Bin quality: results/04_binning/checkm/"
    echo "  - Taxonomic classification: results/05_annotation/gtdbtk/"
    echo "  - Methanogen/methanotroph analysis: results/06_taxonomy/"
fi
