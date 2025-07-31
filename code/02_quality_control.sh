#!/bin/bash

# Quality control and preprocessing of metagenomic reads
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${PROJECT_DIR}/results/01_qc"
THREADS=8

# Activate QC environment
source $(conda info --base)/etc/profile.d/conda.sh

# Try to activate the main QC environment, fall back to alternatives
if conda activate metagenome_qc 2>/dev/null; then
    echo "Using metagenome_qc environment"
elif conda activate trimming 2>/dev/null; then
    echo "Using trimming environment (alternative)"
else
    echo "Warning: No QC environment found. Please run setup script first."
    echo "Attempting to continue with base environment..."
fi

echo "Starting quality control analysis..."

# Create output directories
mkdir -p ${RESULTS_DIR}/{fastqc_raw,fastqc_clean,trimmed,multiqc}

# Sample list - direct file prefixes (the actual filenames without .fastq.gz)
SAMPLES=("53394_R1_A15" "53394_R2_A15" "53395_R1_A16" "53395_R2_A16" "53396_R1_B6" "53396_R2_B6" "53397_R1_B10" "53397_R2_B10" "53398_R1_A16S" "53398_R2_A16S" "53399_R1_B10S" "53399_R2_B10S")

# Sample pairs for processing (R1/R2 pairs)
SAMPLE_PAIRS=("53394_A15:53394_R1_A15:53394_R2_A15" "53395_A16:53395_R1_A16:53395_R2_A16" "53396_B6:53396_R1_B6:53396_R2_B6" "53397_B10:53397_R1_B10:53397_R2_B10" "53398_A16S:53398_R1_A16S:53398_R2_A16S" "53399_B10S:53399_R1_B10S:53399_R2_B10S")

# Step 1: FastQC on raw reads
echo "Running FastQC on raw reads..."
for sample_pair in "${SAMPLE_PAIRS[@]}"; do
    IFS=':' read -r sample_name r1_file r2_file <<< "$sample_pair"
    echo "Processing sample: $sample_name (files: ${r1_file}.fastq.gz, ${r2_file}.fastq.gz)"
    
    fastqc ${DATA_DIR}/${r1_file}.fastq.gz \
           ${DATA_DIR}/${r2_file}.fastq.gz \
           -o ${RESULTS_DIR}/fastqc_raw -t ${THREADS}
done

# Step 2: Quality trimming with fastp
echo "Quality trimming with fastp..."
for sample_pair in "${SAMPLE_PAIRS[@]}"; do
    IFS=':' read -r sample_name r1_file r2_file <<< "$sample_pair"
    echo "Trimming sample: $sample_name"
    
    fastp -i ${DATA_DIR}/${r1_file}.fastq.gz \
          -I ${DATA_DIR}/${r2_file}.fastq.gz \
          -o ${RESULTS_DIR}/trimmed/${sample_name}_R1_trimmed.fastq.gz \
          -O ${RESULTS_DIR}/trimmed/${sample_name}_R2_trimmed.fastq.gz \
          --thread ${THREADS} \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 40 \
          --length_required 50 \
          --html ${RESULTS_DIR}/trimmed/${sample_name}_fastp.html \
          --json ${RESULTS_DIR}/trimmed/${sample_name}_fastp.json
done

# Step 3: FastQC on trimmed reads
echo "Running FastQC on trimmed reads..."
for sample_pair in "${SAMPLE_PAIRS[@]}"; do
    IFS=':' read -r sample_name r1_file r2_file <<< "$sample_pair"
    echo "QC for trimmed sample: $sample_name"
    fastqc ${RESULTS_DIR}/trimmed/${sample_name}_R1_trimmed.fastq.gz \
           ${RESULTS_DIR}/trimmed/${sample_name}_R2_trimmed.fastq.gz \
           -o ${RESULTS_DIR}/fastqc_clean -t ${THREADS}
done

# Step 4: MultiQC report
echo "Generating MultiQC report..."
multiqc ${RESULTS_DIR}/ -o ${RESULTS_DIR}/multiqc

echo "Quality control completed!"
echo "Results saved in: ${RESULTS_DIR}"
echo "Check MultiQC report: ${RESULTS_DIR}/multiqc/multiqc_report.html"
