#!/bin/bash

# Individual sample annotation for host-specific microbiome analysis
# Author: Generated for farting_trees project  
# Date: August 2025
# Purpose: Annotate bins from individual sample binning (05b_individual_binning.sh)

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
INDIVIDUAL_BINS_DIR="${WORK_DIR}/results/04b_individual_binning/binning"
RESULTS_DIR="${WORK_DIR}/results/05b_annotation_individual"
THREADS=8

# Activate required environments
source $(conda info --base)/etc/profile.d/conda.sh

echo "========================================"
echo "INDIVIDUAL SAMPLE ANNOTATION PIPELINE"
echo "========================================"
echo "This script annotates bins from individual sample binning"
echo "Complementary to co-assembly annotation (05_annotation.sh)"
echo ""

# Create output directories
mkdir -p ${RESULTS_DIR}/{prokka,gtdbtk,eggnog,summary}

# Only process the large samples that were individually binned
SAMPLES=("53395_A16" "53398_A16S" "53399_B10S")

echo "Processing individual bins from samples: ${SAMPLES[*]}"
echo ""

# Check if individual binning was completed
missing_samples=()
for sample in "${SAMPLES[@]}"; do
    if [ ! -d "${INDIVIDUAL_BINS_DIR}/${sample}/metabat2" ]; then
        missing_samples+=("$sample")
    fi
done

if [ ${#missing_samples[@]} -gt 0 ]; then
    echo "⚠ Warning: Individual binning not found for: ${missing_samples[*]}"
    echo "Please run 05b_individual_binning.sh first"
    echo "Continuing with available samples..."
fi

# Function to check tool availability
check_tool() {
    if ! command -v $1 &> /dev/null; then
        echo "❌ Error: $1 not found. Please install or activate the appropriate conda environment."
        return 1
    else
        echo "✓ $1 found"
        return 0
    fi
}

# ==========================================
# STEP 1: PROKKA ANNOTATION
# ==========================================
echo "=========================================="
echo "STEP 1: PROKKA ANNOTATION (Gene Prediction)"
echo "=========================================="

conda activate metagenome_annotation

echo "Checking Prokka availability..."
if check_tool prokka; then
    echo "Starting Prokka annotation for individual bins..."
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        sample_bins_dir="${INDIVIDUAL_BINS_DIR}/${sample}/metabat2"
        
        if [ ! -d "$sample_bins_dir" ]; then
            echo "⚠ Skipping $sample - no bins found"
            continue
        fi
        
        # Count bins for this sample
        bin_count=$(ls ${sample_bins_dir}/${sample}_bin*.fa 2>/dev/null | wc -l)
        if [ $bin_count -eq 0 ]; then
            echo "⚠ No bins found for sample $sample"
            continue
        fi
        
        echo "Processing $bin_count bins for sample $sample..."
        
        # Create sample-specific output directory
        mkdir -p ${RESULTS_DIR}/prokka/${sample}
        
        # Process each bin
        for bin_file in ${sample_bins_dir}/${sample}_bin*.fa; do
            if [ -f "$bin_file" ]; then
                bin_name=$(basename "$bin_file" .fa)
                output_dir="${RESULTS_DIR}/prokka/${sample}/${bin_name}"
                
                if [ ! -d "$output_dir" ] || [ ! -f "${output_dir}/${bin_name}.gff" ]; then
                    echo "  Annotating ${bin_name}..."
                    
                    prokka --outdir "$output_dir" \
                           --prefix "$bin_name" \
                           --kingdom Bacteria \
                           --metagenome \
                           --cpus $THREADS \
                           --force \
                           "$bin_file"
                    
                    echo "  ✓ ${bin_name} completed"
                else
                    echo "  ✓ ${bin_name} already annotated"
                fi
            fi
        done
        
        echo "✓ Sample $sample Prokka annotation completed"
        echo ""
    done
    
    echo "✓ All Prokka annotations completed!"
else
    echo "❌ Prokka not available - skipping gene prediction"
fi

echo ""

# ==========================================
# STEP 2: GTDB-TK TAXONOMIC CLASSIFICATION  
# ==========================================
echo "=========================================="
echo "STEP 2: GTDB-TK TAXONOMIC CLASSIFICATION"
echo "=========================================="

echo "Checking GTDB-Tk availability..."
if check_tool gtdbtk; then
    echo "Starting GTDB-Tk taxonomic classification..."
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        sample_bins_dir="${INDIVIDUAL_BINS_DIR}/${sample}/metabat2"
        
        if [ ! -d "$sample_bins_dir" ]; then
            echo "⚠ Skipping $sample - no bins found"
            continue
        fi
        
        bin_count=$(ls ${sample_bins_dir}/${sample}_bin*.fa 2>/dev/null | wc -l)
        if [ $bin_count -eq 0 ]; then
            echo "⚠ No bins found for sample $sample"
            continue
        fi
        
        echo "Running GTDB-Tk for $bin_count bins from sample $sample..."
        
        sample_output_dir="${RESULTS_DIR}/gtdbtk/${sample}"
        mkdir -p "$sample_output_dir"
        
        if [ ! -f "${sample_output_dir}/gtdbtk.bac120.summary.tsv" ] && \
           [ ! -f "${sample_output_dir}/gtdbtk.ar53.summary.tsv" ]; then
            
            # Set GTDB-Tk database path (adjust if needed)
            export GTDBTK_DATA_PATH="/home/jiayi-chen/miniconda3/envs/metagenome_annotation/share/gtdbtk-2.4.0/db/"
            
            gtdbtk classify_wf --genome_dir "$sample_bins_dir" \
                              --out_dir "$sample_output_dir" \
                              --cpus $THREADS \
                              --extension fa
            
            echo "✓ GTDB-Tk completed for sample $sample"
        else
            echo "✓ GTDB-Tk results already exist for sample $sample"
        fi
        echo ""
    done
    
    echo "✓ All GTDB-Tk classifications completed!"
else
    echo "❌ GTDB-Tk not available - skipping taxonomic classification"
    echo "Note: You may need to install GTDB-Tk database"
fi

echo ""

# ==========================================
# STEP 3: EGGNOG-MAPPER FUNCTIONAL ANNOTATION
# ==========================================
echo "=========================================="
echo "STEP 3: EGGNOG-MAPPER FUNCTIONAL ANNOTATION"
echo "=========================================="

echo "Checking eggNOG-mapper availability..."
if check_tool emapper.py; then
    echo "Starting eggNOG-mapper functional annotation..."
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        sample_bins_dir="${INDIVIDUAL_BINS_DIR}/${sample}/metabat2"
        
        if [ ! -d "$sample_bins_dir" ]; then
            echo "⚠ Skipping $sample - no bins found"
            continue
        fi
        
        bin_count=$(ls ${sample_bins_dir}/${sample}_bin*.fa 2>/dev/null | wc -l)
        if [ $bin_count -eq 0 ]; then
            echo "⚠ No bins found for sample $sample"
            continue
        fi
        
        echo "Running eggNOG-mapper for sample $sample..."
        
        # Create sample-specific output directory
        mkdir -p ${RESULTS_DIR}/eggnog/${sample}
        
        # Process each bin
        for bin_file in ${sample_bins_dir}/${sample}_bin*.fa; do
            if [ -f "$bin_file" ]; then
                bin_name=$(basename "$bin_file" .fa)
                output_prefix="${RESULTS_DIR}/eggnog/${sample}/${bin_name}"
                
                if [ ! -f "${output_prefix}.emapper.annotations" ]; then
                    echo "  Annotating ${bin_name}..."
                    
                    emapper.py -i "$bin_file" \
                               -o "$bin_name" \
                               --output_dir "${RESULTS_DIR}/eggnog/${sample}" \
                               --data_dir /home/jiayi-chen/miniconda3/envs/metagenome_annotation/lib/python3.9/site-packages/data \
                               --cpu $THREADS \
                               --override
                    
                    echo "  ✓ ${bin_name} eggNOG annotation completed"
                else
                    echo "  ✓ ${bin_name} eggNOG annotation already exists"
                fi
            fi
        done
        
        echo "✓ Sample $sample eggNOG annotation completed"
        echo ""
    done
    
    echo "✓ All eggNOG annotations completed!"
else
    echo "❌ eggNOG-mapper not available - skipping functional annotation"
fi

echo ""

# ==========================================
# STEP 4: GENERATE SUMMARY REPORT
# ==========================================
echo "=========================================="
echo "STEP 4: GENERATING SUMMARY REPORT"
echo "=========================================="

cat > ${RESULTS_DIR}/summary/individual_annotation_summary.md << 'EOF'
# Individual Sample Annotation Summary

## Overview
This analysis provides taxonomic classification and functional annotation for bins generated from individual sample binning (complementary to co-assembly annotation).

## Processed Samples
Only large samples with sufficient data for individual binning:

EOF

# Add sample-specific results
for sample in "${SAMPLES[@]}"; do
    echo "### Sample: $sample" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    
    sample_bins_dir="${INDIVIDUAL_BINS_DIR}/${sample}/metabat2"
    if [ -d "$sample_bins_dir" ]; then
        bin_count=$(ls ${sample_bins_dir}/${sample}_bin*.fa 2>/dev/null | wc -l)
        echo "- **Bins processed**: $bin_count" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
        
        # Check annotation completeness
        prokka_done=0
        gtdbtk_done=0
        eggnog_done=0
        
        if [ -d "${RESULTS_DIR}/prokka/${sample}" ]; then
            prokka_done=$(find ${RESULTS_DIR}/prokka/${sample} -name "*.gff" | wc -l)
        fi
        
        if [ -f "${RESULTS_DIR}/gtdbtk/${sample}/gtdbtk.bac120.summary.tsv" ] || \
           [ -f "${RESULTS_DIR}/gtdbtk/${sample}/gtdbtk.ar53.summary.tsv" ]; then
            gtdbtk_done=1
        fi
        
        if [ -d "${RESULTS_DIR}/eggnog/${sample}" ]; then
            eggnog_done=$(find ${RESULTS_DIR}/eggnog/${sample} -name "*.emapper.annotations" | wc -l)
        fi
        
        echo "- **Prokka annotations**: $prokka_done/$bin_count completed" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
        echo "- **GTDB-Tk classification**: $([ $gtdbtk_done -eq 1 ] && echo "✓ Completed" || echo "❌ Not completed")" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
        echo "- **eggNOG annotations**: $eggnog_done/$bin_count completed" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    else
        echo "- **Status**: No individual bins found (run 05b_individual_binning.sh first)" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    fi
    echo "" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
done

# Add comparison with co-assembly
echo "## Comparison with Co-assembly Annotation" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md

if [ -d "${WORK_DIR}/results/06_annotation" ]; then
    echo "- **Co-assembly annotation**: Available in results/06_annotation/" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    echo "- **Individual annotation**: Available in results/06b_annotation_individual/" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    echo "" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    echo "### Analysis Strategy" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    echo "1. **Co-assembly bins** → Shared microorganisms and comparative analysis" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    echo "2. **Individual bins** → Host-specific microorganisms and detailed characterization" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
    echo "3. **Combined analysis** → Comprehensive understanding of tree microbiome diversity" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
else
    echo "- **Co-assembly annotation**: Not found (run 06_annotation.sh for co-assembly bins)" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
fi

echo "" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "## Next Steps" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "1. **Compare taxonomies**: Identify host-specific vs shared taxa" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "2. **Functional analysis**: Focus on methane metabolism pathways" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "3. **Host comparison**: Analyze differences between A-group vs B-group trees" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md
echo "4. **Integration**: Combine with co-assembly results for comprehensive analysis" >> ${RESULTS_DIR}/summary/individual_annotation_summary.md

# Summary statistics
echo ""
echo "=========================================="
echo "ANNOTATION SUMMARY STATISTICS"
echo "=========================================="

total_individual_bins=0
total_prokka=0
total_eggnog=0
total_gtdbtk=0

for sample in "${SAMPLES[@]}"; do
    sample_bins_dir="${INDIVIDUAL_BINS_DIR}/${sample}/metabat2"
    if [ -d "$sample_bins_dir" ]; then
        sample_bins=$(ls ${sample_bins_dir}/${sample}_bin*.fa 2>/dev/null | wc -l)
        total_individual_bins=$((total_individual_bins + sample_bins))
        
        if [ -d "${RESULTS_DIR}/prokka/${sample}" ]; then
            sample_prokka=$(find ${RESULTS_DIR}/prokka/${sample} -name "*.gff" | wc -l)
            total_prokka=$((total_prokka + sample_prokka))
        fi
        
        if [ -d "${RESULTS_DIR}/eggnog/${sample}" ]; then
            sample_eggnog=$(find ${RESULTS_DIR}/eggnog/${sample} -name "*.emapper.annotations" | wc -l)
            total_eggnog=$((total_eggnog + sample_eggnog))
        fi
        
        if [ -f "${RESULTS_DIR}/gtdbtk/${sample}/gtdbtk.bac120.summary.tsv" ] || \
           [ -f "${RESULTS_DIR}/gtdbtk/${sample}/gtdbtk.ar53.summary.tsv" ]; then
            total_gtdbtk=$((total_gtdbtk + 1))
        fi
    fi
done

echo "Total individual bins processed: $total_individual_bins"
echo "Prokka annotations completed: $total_prokka/$total_individual_bins"
echo "GTDB-Tk classifications completed: $total_gtdbtk/${#SAMPLES[@]} samples"
echo "eggNOG annotations completed: $total_eggnog/$total_individual_bins"

echo ""
echo "=========================================="
echo "INDIVIDUAL ANNOTATION COMPLETED!"
echo "=========================================="
echo ""
echo "Results saved in: ${RESULTS_DIR}"
echo "Summary report: ${RESULTS_DIR}/summary/individual_annotation_summary.md"
echo ""
echo "This complements your co-assembly annotation for comprehensive analysis!"
echo ""
echo "Next: Compare individual vs co-assembly results to identify:"
echo "  - Host-specific microorganisms (individual bins)"  
echo "  - Shared microorganisms (co-assembly bins)"
echo "  - Methane metabolism differences between tree hosts"
