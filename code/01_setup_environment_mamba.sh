#!/bin/bash

# Alternative setup script using mamba for better dependency resolution
# Author: Generated for farting_trees project
# Date: July 2025

set -e

echo "Setting up conda environments for metagenomic analysis (using mamba)..."

# Check if mamba is available, if not install it
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Installing mamba..."
    conda install mamba -n base -c conda-forge -y
fi

# Create base directory for results
mkdir -p results/{01_qc,02_assembly,03_binning,04_annotation,05_taxonomy}

# Create QC environment with flexible versions
echo "Creating QC environment..."
mamba create -n metagenome_qc -c bioconda -c conda-forge -y \
    fastqc \
    multiqc \
    fastp \
    trimmomatic

# Create assembly environment  
echo "Creating assembly environment..."
mamba create -n metagenome_assembly -c bioconda -c conda-forge -y \
    megahit \
    spades \
    quast \
    bowtie2 \
    samtools \
    bwa

# Create binning environment
echo "Creating binning environment..."
mamba create -n metagenome_binning -c bioconda -c conda-forge -y \
    metabat2 \
    maxbin2 \
    das_tool \
    checkm-genome \
    vamb

# Create annotation environment
echo "Creating annotation environment..."
mamba create -n metagenome_annotation -c bioconda -c conda-forge -y \
    prokka \
    gtdbtk \
    eggnog-mapper \
    diamond

echo "Environment setup complete!"
echo "Activate environments using:"
echo "  conda activate metagenome_qc"
echo "  conda activate metagenome_assembly" 
echo "  conda activate metagenome_binning"
echo "  conda activate metagenome_annotation"
