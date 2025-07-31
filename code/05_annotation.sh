#!/bin/bash

# Taxonomic classification and functional annotation of MAGs
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
BINNING_DIR="../results/03_binning"
RESULTS_DIR="../results/04_annotation"
THREADS=8

# Activate annotation environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_annotation

echo "Starting MAG annotation and taxonomic classification..."

# Create output directories
mkdir -p ${RESULTS_DIR}/{prokka,gtdbtk,eggnog,summary}

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

# Step 1: Gene prediction and annotation with Prokka
echo "Running Prokka for gene prediction and annotation..."
for sample in "${SAMPLES[@]}"; do
    echo "Prokka annotation for sample: $sample"
    
    # Process DAS Tool bins (high quality)
    if [ -d "${BINNING_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins" ]; then
        for bin in ${BINNING_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins/*.fa; do
            if [ -f "$bin" ]; then
                bin_name=$(basename "$bin" .fa)
                echo "  Annotating bin: $bin_name"
                
                prokka --outdir ${RESULTS_DIR}/prokka/${sample}_${bin_name} \
                       --prefix ${sample}_${bin_name} \
                       --kingdom Bacteria \
                       --cpus ${THREADS} \
                       --force \
                       $bin
            fi
        done
    fi
    
    # Process MetaBAT2 bins if DAS Tool bins not available
    if [ ! -d "${BINNING_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins" ]; then
        for bin in ${BINNING_DIR}/metabat2/${sample}_bin*.fa; do
            if [ -f "$bin" ]; then
                bin_name=$(basename "$bin" .fa)
                echo "  Annotating MetaBAT2 bin: $bin_name"
                
                prokka --outdir ${RESULTS_DIR}/prokka/${sample}_${bin_name} \
                       --prefix ${sample}_${bin_name} \
                       --kingdom Bacteria \
                       --cpus ${THREADS} \
                       --force \
                       $bin
            fi
        done
    fi
done

# Step 2: Taxonomic classification with GTDB-Tk
echo "Running GTDB-Tk for taxonomic classification..."

# Collect all high-quality bins for batch processing
mkdir -p ${RESULTS_DIR}/gtdbtk/input_bins

for sample in "${SAMPLES[@]}"; do
    # Use DAS Tool bins preferentially
    if [ -d "${BINNING_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins" ]; then
        for bin in ${BINNING_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins/*.fa; do
            if [ -f "$bin" ]; then
                bin_name="${sample}_$(basename "$bin")"
                cp "$bin" ${RESULTS_DIR}/gtdbtk/input_bins/$bin_name
            fi
        done
    else
        # Use MetaBAT2 bins as fallback
        for bin in ${BINNING_DIR}/metabat2/${sample}_bin*.fa; do
            if [ -f "$bin" ]; then
                bin_name="${sample}_$(basename "$bin")"
                cp "$bin" ${RESULTS_DIR}/gtdbtk/input_bins/$bin_name
            fi
        done
    fi
done

# Run GTDB-Tk classify workflow
if [ "$(ls -A ${RESULTS_DIR}/gtdbtk/input_bins)" ]; then
    echo "Running GTDB-Tk classification..."
    gtdbtk classify_wf --genome_dir ${RESULTS_DIR}/gtdbtk/input_bins \
                       --out_dir ${RESULTS_DIR}/gtdbtk/output \
                       --cpus ${THREADS} \
                       --extension fa
else
    echo "No bins found for GTDB-Tk classification"
fi

# Step 3: Functional annotation with eggNOG-mapper
echo "Running eggNOG-mapper for functional annotation..."
for sample in "${SAMPLES[@]}"; do
    echo "eggNOG annotation for sample: $sample"
    
    # Process protein files from Prokka
    for prokka_dir in ${RESULTS_DIR}/prokka/${sample}_*; do
        if [ -d "$prokka_dir" ]; then
            bin_name=$(basename "$prokka_dir")
            protein_file="$prokka_dir/${bin_name}.faa"
            
            if [ -f "$protein_file" ]; then
                echo "  Processing proteins for: $bin_name"
                
                emapper.py -i "$protein_file" \
                          --output ${bin_name}_eggnog \
                          --output_dir ${RESULTS_DIR}/eggnog \
                          --cpu ${THREADS} \
                          --override
            fi
        fi
    done
done

# Step 4: Generate summary report
echo "Generating summary report..."
cat > ${RESULTS_DIR}/summary/analysis_summary.md << 'EOF'
# MAG Analysis Summary

## Overview
This report summarizes the results of metagenomic assembly, binning, and annotation.

## Quality Metrics
- **High-quality MAGs**: >90% completeness, <5% contamination
- **Medium-quality MAGs**: >50% completeness, <10% contamination
- **Low-quality MAGs**: <50% completeness or >10% contamination

## Key Findings

### Sample Summary
EOF

# Add CheckM results to summary
for sample in "${SAMPLES[@]}"; do
    echo "### Sample: $sample" >> ${RESULTS_DIR}/summary/analysis_summary.md
    
    # Add DAS Tool quality if available
    if [ -f "${BINNING_DIR}/checkm/${sample}_dastool_quality.txt" ]; then
        echo "#### DAS Tool Bins Quality:" >> ${RESULTS_DIR}/summary/analysis_summary.md
        echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
        head -20 "${BINNING_DIR}/checkm/${sample}_dastool_quality.txt" >> ${RESULTS_DIR}/summary/analysis_summary.md
        echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
    fi
    
    # Add MetaBAT2 quality if available
    if [ -f "${BINNING_DIR}/checkm/${sample}_metabat2_quality.txt" ]; then
        echo "#### MetaBAT2 Bins Quality:" >> ${RESULTS_DIR}/summary/analysis_summary.md
        echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
        head -20 "${BINNING_DIR}/checkm/${sample}_metabat2_quality.txt" >> ${RESULTS_DIR}/summary/analysis_summary.md
        echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
    fi
    
    echo "" >> ${RESULTS_DIR}/summary/analysis_summary.md
done

# Add GTDB-Tk results to summary
if [ -f "${RESULTS_DIR}/gtdbtk/output/gtdbtk.bac120.summary.tsv" ]; then
    echo "## Taxonomic Classification (GTDB-Tk)" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
    cat "${RESULTS_DIR}/gtdbtk/output/gtdbtk.bac120.summary.tsv" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

echo "Annotation and taxonomic classification completed!"
echo "Results saved in: ${RESULTS_DIR}"
echo "Summary report: ${RESULTS_DIR}/summary/analysis_summary.md"
