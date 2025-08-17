#!/bin/bash

# Simple MetaSPAdes Assembly
set -e

# Configuration
INPUT_DIR="/home/jiayi-chen/Documents/farting_trees/results/01_qc/primer_removed"
OUTPUT_DIR="/home/jiayi-chen/Documents/farting_trees/results/02_assembly"
THREADS=8

# Setup - Create new environment with compatible Python version
source $(conda info --base)/etc/profile.d/conda.sh

# Create new environment with Python 3.8 and SPAdes
ENV_NAME="spades_assembly"
echo "Setting up SPAdes environment..."

if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating new conda environment: ${ENV_NAME}"
    conda create -n ${ENV_NAME} python=3.8 spades -c bioconda -y
else
    echo "Environment ${ENV_NAME} already exists"
fi

conda activate ${ENV_NAME}
mkdir -p "$OUTPUT_DIR"

# Verify SPAdes installation
echo "Verifying SPAdes installation..."
spades.py --version

echo "MetaSPAdes Assembly for individual samples..."

# Find and process sample pairs
for r1_file in "$INPUT_DIR"/*_R1_no_primers.fastq.gz; do
    [[ -f "$r1_file" ]] || continue
    r2_file="${r1_file/_R1_no_primers/_R2_no_primers}"
    [[ -f "$r2_file" ]] || continue
    
    sample_name=$(basename "$r1_file" | sed 's/_R1_no_primers.fastq.gz$//')
    sample_output="$OUTPUT_DIR/${sample_name}"
    
    # Skip if already assembled
    [[ -f "$sample_output/contigs.fasta" ]] && { echo "✓ $sample_name (skipped)"; continue; }
    
    echo "Assembling $sample_name..."
    
    # Run metaSPAdes
    spades.py --meta \
        -o "$sample_output" \
        --only-assembler \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -t $THREADS
    
    echo "✓ $sample_name completed"
done

echo "Assembly completed! Output in: $OUTPUT_DIR"
