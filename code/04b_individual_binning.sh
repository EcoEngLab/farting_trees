#!/bin/bash

# Individual sample binning for host-specific microbiome analysis
# Author: Generated for farting_trees project  
# Date: July 2025

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/03_assembly"
RESULTS_DIR="${WORK_DIR}/results/05b_individual_binning"
THREADS=8

# Activate required environments
source $(conda info --base)/etc/profile.d/conda.sh

echo "Starting individual sample binning analysis..."
echo "This will complement the co-assembly approach for host-specific analysis"

# Create output directories
mkdir -p ${RESULTS_DIR}/{assemblies,binning,summary}

# Sample list - ONLY large samples with sufficient data for individual binning
SAMPLES=("53395_A16" "53398_A16S" "53399_B10S")

# Group samples by host type
A_GROUP=("53395_A16" "53398_A16S")  # Large A samples
B_GROUP=("53399_B10S")              # Large B sample

# Small samples (insufficient for individual binning): 53394_A15 (146MB), 53396_B6 (170MB), 53397_B10 (201MB)
# These benefit more from co-assembly approach

echo "Sample grouping (large samples only):"
echo "A group (large): ${A_GROUP[*]}"
echo "B group (large): ${B_GROUP[*]}"
echo "Note: Small samples (53394_A15, 53396_B6, 53397_B10) use co-assembly approach"

# Function to run individual assembly and binning
run_individual_binning() {
    local sample=$1
    local sample_dir="${RESULTS_DIR}/assemblies/${sample}"
    local binning_dir="${RESULTS_DIR}/binning/${sample}"
    
    echo "========================================"
    echo "Processing sample: $sample"
    echo "========================================"
    
    # Create sample directories
    mkdir -p ${sample_dir} ${binning_dir}/{metabat2,checkm}
    
    # Input files
    R1="${WORK_DIR}/results/02_host_removal/clean_reads/${sample}_R1_clean.fastq.gz"
    R2="${WORK_DIR}/results/02_host_removal/clean_reads/${sample}_R2_clean.fastq.gz"
    
    # Check if files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "⚠ Warning: Missing files for $sample, skipping..."
        return
    fi
    
    # Step 1: Individual assembly with MEGAHIT
    echo "Running individual assembly for $sample..."
    conda activate metagenome_assembly
    
    if [ ! -f "${sample_dir}/${sample}_contigs.fa" ]; then
        megahit -1 $R1 -2 $R2 \
                -o ${sample_dir}/megahit_output \
                --min-contig-len 1000 \
                --k-min 21 \
                --k-max 141 \
                --k-step 12 \
                -t ${THREADS}
        
        # Copy and rename contigs
        cp ${sample_dir}/megahit_output/final.contigs.fa ${sample_dir}/${sample}_contigs.fa
        
        echo "✓ Assembly completed for $sample"
    else
        echo "✓ Assembly already exists for $sample"
    fi
    
    # Step 2: Self-mapping for coverage
    echo "Mapping $sample to its own assembly..."
    
    # BWA index
    if [ ! -f "${sample_dir}/${sample}_contigs.fa.bwt" ]; then
        bwa index ${sample_dir}/${sample}_contigs.fa
    fi
    
    # Mapping
    if [ ! -f "${binning_dir}/${sample}_sorted.bam" ]; then
        bwa mem -t ${THREADS} ${sample_dir}/${sample}_contigs.fa $R1 $R2 | \
        samtools view -@ ${THREADS} -bS -F 4 - | \
        samtools sort -@ ${THREADS} -o ${binning_dir}/${sample}_sorted.bam
        
        samtools index ${binning_dir}/${sample}_sorted.bam
        echo "✓ Self-mapping completed for $sample"
    else
        echo "✓ Self-mapping already exists for $sample"
    fi
    
    # Step 3: Generate depth file
    echo "Generating depth file for $sample..."
    conda activate metagenome_binning
    
    if [ ! -f "${binning_dir}/${sample}_depth.txt" ]; then
        jgi_summarize_bam_contig_depths --outputDepth ${binning_dir}/${sample}_depth.txt \
                                       ${binning_dir}/${sample}_sorted.bam
        echo "✓ Depth file generated for $sample"
    else
        echo "✓ Depth file already exists for $sample"
    fi
    
    # Step 4: MetaBAT2 binning
    echo "Running MetaBAT2 binning for $sample..."
    
    if [ ! -d "${binning_dir}/metabat2" ] || [ -z "$(ls -A ${binning_dir}/metabat2 2>/dev/null)" ]; then
        metabat2 -i ${sample_dir}/${sample}_contigs.fa \
                 -a ${binning_dir}/${sample}_depth.txt \
                 -o ${binning_dir}/metabat2/${sample}_bin \
                 --minContig 2500 \
                 --numThreads ${THREADS}
        
        echo "✓ MetaBAT2 binning completed for $sample"
    else
        echo "✓ MetaBAT2 binning already exists for $sample"
    fi
    
    # Step 5: Basic bin statistics
    echo "Analyzing bins for $sample..."
    
    bin_count=$(ls ${binning_dir}/metabat2/${sample}_bin*.fa 2>/dev/null | wc -l)
    echo "Generated $bin_count bins for sample $sample"
    
    # Calculate total assembly size
    if [ -f "${sample_dir}/${sample}_contigs.fa" ]; then
        assembly_size=$(grep -v "^>" ${sample_dir}/${sample}_contigs.fa | wc -c)
        assembly_size_mb=$((assembly_size / 1000000))
        echo "Assembly size: ${assembly_size_mb} MB"
    fi
    
    echo "Sample $sample processing completed!"
    echo ""
}

# Process all samples
for sample in "${SAMPLES[@]}"; do
    run_individual_binning $sample
done

# Generate comparative summary
echo "========================================"
echo "GENERATING COMPARATIVE SUMMARY"
echo "========================================"

# Create summary report
cat > ${RESULTS_DIR}/summary/individual_binning_summary.md << 'EOF'
# Individual Sample Binning Summary

## Overview
This analysis complements the co-assembly approach by binning each sample individually to identify host-specific microorganisms.

## Sample Groups
- **A Group (A15/A16)**: 53394_A15, 53395_A16, 53398_A16S
- **B Group (B6/B10)**: 53396_B6, 53397_B10, 53399_B10S

## Individual Sample Results

EOF

# Add results for each sample
for sample in "${SAMPLES[@]}"; do
    echo "### Sample: $sample" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
    
    if [ -d "${RESULTS_DIR}/binning/${sample}/metabat2" ]; then
        bin_count=$(ls ${RESULTS_DIR}/binning/${sample}/metabat2/${sample}_bin*.fa 2>/dev/null | wc -l)
        echo "- **Bins generated**: $bin_count" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
        
        if [ -f "${RESULTS_DIR}/assemblies/${sample}/${sample}_contigs.fa" ]; then
            assembly_size=$(grep -v "^>" ${RESULTS_DIR}/assemblies/${sample}/${sample}_contigs.fa | wc -c)
            assembly_size_mb=$((assembly_size / 1000000))
            echo "- **Assembly size**: ${assembly_size_mb} MB" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
        fi
    else
        echo "- **Status**: Not processed" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
    fi
    echo "" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
done

# Compare with co-assembly results
echo "## Comparison with Co-assembly" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
echo "" >> ${RESULTS_DIR}/summary/individual_binning_summary.md

if [ -d "${WORK_DIR}/results/05_binning/metabat2" ]; then
    coassembly_bins=$(ls ${WORK_DIR}/results/05_binning/metabat2/coassembly_bin*.fa 2>/dev/null | wc -l)
    echo "- **Co-assembly bins**: $coassembly_bins" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
    
    total_individual_bins=0
    for sample in "${SAMPLES[@]}"; do
        if [ -d "${RESULTS_DIR}/binning/${sample}/metabat2" ]; then
            sample_bins=$(ls ${RESULTS_DIR}/binning/${sample}/metabat2/${sample}_bin*.fa 2>/dev/null | wc -l)
            total_individual_bins=$((total_individual_bins + sample_bins))
        fi
    done
    echo "- **Total individual bins**: $total_individual_bins" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
    echo "" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
fi

echo "## Recommendations" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
echo "" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
echo "1. **Co-assembly bins**: Use for comparative analysis and shared microorganisms" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
echo "2. **Individual bins**: Use for host-specific microorganisms and detailed characterization" >> ${RESULTS_DIR}/summary/individual_binning_summary.md
echo "3. **Combined approach**: Merge unique bins from both strategies for comprehensive analysis" >> ${RESULTS_DIR}/summary/individual_binning_summary.md

echo ""
echo "=========================================="
echo "INDIVIDUAL BINNING COMPLETED!"
echo "=========================================="
echo ""
echo "Results saved in: ${RESULTS_DIR}"
echo "Summary report: ${RESULTS_DIR}/summary/individual_binning_summary.md"
echo ""
echo "Next steps:"
echo "1. Compare individual vs co-assembly results"
echo "2. Identify host-specific vs shared microorganisms"  
echo "3. Proceed with annotation of unique bins"
