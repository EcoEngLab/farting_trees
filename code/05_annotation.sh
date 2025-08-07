#!/bin/bash

# Taxonomic classification and functional annotation of MAGs
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
BINNING_DIR="${WORK_DIR}/results/04_binning"
RESULTS_DIR="${WORK_DIR}/results/05_annotation"
THREADS=8

# Activate annotation environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_annotation

echo "Starting MAG annotation and taxonomic classification..."

# Create output directories
mkdir -p ${RESULTS_DIR}/{prokka,gtdbtk,eggnog,summary}

echo "Working with co-assembly bins from MetaBAT2..."

# Step 1: Gene prediction and annotation with Prokka
echo "Running Prokka for gene prediction and annotation..."

# Process MetaBAT2 co-assembly bins
if [ -d "${BINNING_DIR}/metabat2" ]; then
    echo "Processing MetaBAT2 co-assembly bins..."
    for bin in ${BINNING_DIR}/metabat2/coassembly_bin.*.fa; do
        if [ -f "$bin" ]; then
            bin_name=$(basename "$bin" .fa)
            echo "  Annotating bin: $bin_name"
            
            # Check if Prokka is available
            if command -v prokka &> /dev/null; then
                prokka --outdir ${RESULTS_DIR}/prokka/${bin_name} \
                       --prefix ${bin_name} \
                       --kingdom Bacteria \
                       --cpus ${THREADS} \
                       --force \
                       $bin
            else
                echo "⚠ Warning: Prokka not found, skipping gene annotation"
                echo "You can install Prokka with: conda install -c bioconda prokka"
                mkdir -p ${RESULTS_DIR}/prokka
                echo "Skipping Prokka due to missing installation" > ${RESULTS_DIR}/prokka/SKIPPED.txt
                break
            fi
        fi
    done
else
    echo "⚠ No MetaBAT2 bins found in ${BINNING_DIR}/metabat2"
fi

# Process DAS Tool bins if available
if [ -d "${BINNING_DIR}/das_tool" ] && [ -d "${BINNING_DIR}/das_tool/coassembly_DASToolBins_DASTool_bins" ]; then
    echo "Processing DAS Tool co-assembly bins..."
    for bin in ${BINNING_DIR}/das_tool/coassembly_DASToolBins_DASTool_bins/*.fa; do
        if [ -f "$bin" ]; then
            bin_name="dastool_$(basename "$bin" .fa)"
            echo "  Annotating DAS Tool bin: $bin_name"
            
            if command -v prokka &> /dev/null; then
                prokka --outdir ${RESULTS_DIR}/prokka/${bin_name} \
                       --prefix ${bin_name} \
                       --kingdom Bacteria \
                       --cpus ${THREADS} \
                       --force \
                       $bin
            fi
        fi
    done
fi

# Step 2: Taxonomic classification with GTDB-Tk
echo "Running GTDB-Tk for taxonomic classification..."

# Collect all high-quality bins for batch processing
mkdir -p ${RESULTS_DIR}/gtdbtk/input_bins

# Copy MetaBAT2 co-assembly bins
if [ -d "${BINNING_DIR}/metabat2" ]; then
    echo "Copying MetaBAT2 co-assembly bins for GTDB-Tk..."
    for bin in ${BINNING_DIR}/metabat2/coassembly_bin.*.fa; do
        if [ -f "$bin" ]; then
            bin_name=$(basename "$bin")
            cp "$bin" ${RESULTS_DIR}/gtdbtk/input_bins/$bin_name
            echo "  Copied: $bin_name"
        fi
    done
fi

# Copy DAS Tool bins if available
if [ -d "${BINNING_DIR}/das_tool/coassembly_DASToolBins_DASTool_bins" ]; then
    echo "Copying DAS Tool co-assembly bins for GTDB-Tk..."
    for bin in ${BINNING_DIR}/das_tool/coassembly_DASToolBins_DASTool_bins/*.fa; do
        if [ -f "$bin" ]; then
            bin_name="dastool_$(basename "$bin")"
            cp "$bin" ${RESULTS_DIR}/gtdbtk/input_bins/$bin_name
            echo "  Copied: $bin_name"
        fi
    done
fi

# Run GTDB-Tk classify workflow
if [ "$(ls -A ${RESULTS_DIR}/gtdbtk/input_bins 2>/dev/null)" ]; then
    echo "Running GTDB-Tk classification..."
    
    # Check if GTDB-Tk is available
    if command -v gtdbtk &> /dev/null; then
        gtdbtk classify_wf --genome_dir ${RESULTS_DIR}/gtdbtk/input_bins \
                           --out_dir ${RESULTS_DIR}/gtdbtk/output \
                           --cpus ${THREADS} \
                           --extension fa
    else
        echo "⚠ Warning: GTDB-Tk not found, skipping taxonomic classification"
        echo "You can install GTDB-Tk with: conda install -c bioconda gtdbtk"
        mkdir -p ${RESULTS_DIR}/gtdbtk
        echo "Skipping GTDB-Tk due to missing installation" > ${RESULTS_DIR}/gtdbtk/SKIPPED.txt
    fi
else
    echo "No bins found for GTDB-Tk classification"
fi

# Step 3: Functional annotation with eggNOG-mapper
echo "Running eggNOG-mapper for functional annotation..."

# Process protein files from Prokka for all bins
if [ -d "${RESULTS_DIR}/prokka" ]; then
    for prokka_dir in ${RESULTS_DIR}/prokka/*; do
        if [ -d "$prokka_dir" ]; then
            bin_name=$(basename "$prokka_dir")
            protein_file="$prokka_dir/${bin_name}.faa"
            
            if [ -f "$protein_file" ]; then
                echo "  Processing proteins for: $bin_name"
                
                # Check if eggNOG-mapper is available
                if command -v emapper.py &> /dev/null; then
                    emapper.py -i "$protein_file" \
                              --output ${bin_name}_eggnog \
                              --output_dir ${RESULTS_DIR}/eggnog \
                              --cpu ${THREADS} \
                              --override
                else
                    echo "⚠ Warning: eggNOG-mapper not found, skipping functional annotation"
                    echo "You can install eggNOG-mapper with: conda install -c bioconda eggnog-mapper"
                    mkdir -p ${RESULTS_DIR}/eggnog
                    echo "Skipping eggNOG due to missing installation" > ${RESULTS_DIR}/eggnog/SKIPPED.txt
                    break
                fi
            else
                echo "  No protein file found for: $bin_name"
            fi
        fi
    done
else
    echo "No Prokka results found, skipping eggNOG annotation"
fi

# Step 4: Generate summary report
echo "Generating summary report..."
cat > ${RESULTS_DIR}/summary/analysis_summary.md << 'EOF'
# MAG Analysis Summary - Co-Assembly Approach

## Overview
This report summarizes the results of metagenomic co-assembly, binning, and annotation.
The analysis used a co-assembly approach with differential coverage binning.

## Quality Metrics
- **High-quality MAGs**: >90% completeness, <5% contamination
- **Medium-quality MAGs**: >50% completeness, <10% contamination  
- **Low-quality MAGs**: <50% completeness or >10% contamination

## Co-Assembly Results

### Binning Summary
EOF

# Add MetaBAT2 bin count
if [ -d "${BINNING_DIR}/metabat2" ]; then
    bin_count=$(ls ${BINNING_DIR}/metabat2/coassembly_bin.*.fa 2>/dev/null | wc -l)
    echo "- **MetaBAT2 bins**: $bin_count" >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

# Add DAS Tool bin count if available
if [ -d "${BINNING_DIR}/das_tool/coassembly_DASToolBins_DASTool_bins" ]; then
    dastool_bin_count=$(ls ${BINNING_DIR}/das_tool/coassembly_DASToolBins_DASTool_bins/*.fa 2>/dev/null | wc -l)
    echo "- **DAS Tool integrated bins**: $dastool_bin_count" >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

# Add CheckM results to summary if available
echo "" >> ${RESULTS_DIR}/summary/analysis_summary.md
echo "### Quality Assessment" >> ${RESULTS_DIR}/summary/analysis_summary.md

if [ -f "${BINNING_DIR}/checkm/metabat2_quality.txt" ]; then
    echo "#### MetaBAT2 Bins Quality:" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
    cat "${BINNING_DIR}/checkm/metabat2_quality.txt" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

if [ -f "${BINNING_DIR}/checkm/dastool_quality.txt" ]; then
    echo "#### DAS Tool Bins Quality:" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
    cat "${BINNING_DIR}/checkm/dastool_quality.txt" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

# Add basic bin statistics
echo "" >> ${RESULTS_DIR}/summary/analysis_summary.md
echo "### Bin Statistics" >> ${RESULTS_DIR}/summary/analysis_summary.md

if [ -f "${WORK_DIR}/code/analyze_bins.py" ]; then
    echo "Running bin analysis for summary..."
    python3 "${WORK_DIR}/code/analyze_bins.py" >> ${RESULTS_DIR}/summary/analysis_summary.md 2>/dev/null || echo "Basic bin analysis not available" >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

# Add GTDB-Tk results to summary
if [ -f "${RESULTS_DIR}/gtdbtk/output/gtdbtk.bac120.summary.tsv" ]; then
    echo "" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo "## Taxonomic Classification (GTDB-Tk)" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
    cat "${RESULTS_DIR}/gtdbtk/output/gtdbtk.bac120.summary.tsv" >> ${RESULTS_DIR}/summary/analysis_summary.md
    echo '```' >> ${RESULTS_DIR}/summary/analysis_summary.md
fi

echo "Annotation and taxonomic classification completed!"
echo "Results saved in: ${RESULTS_DIR}"
echo "Summary report: ${RESULTS_DIR}/summary/analysis_summary.md"
