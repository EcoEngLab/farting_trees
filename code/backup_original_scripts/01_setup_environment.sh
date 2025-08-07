#!/bin/bash

# Setup conda environments for metagenomic analysis
# Author: Generated for farting_trees project
# Date: July 2025

set -e

echo "Setting up conda environments for metagenomic analysis..."

# Create base directory for results
mkdir -p results/{01_qc,02_assembly,03_binning,04_annotation,05_taxonomy}

# Create QC environment
echo "Creating QC environment..."
conda create -n metagenome_qc -c bioconda -c conda-forge -y \
    fastqc=0.12.1 \
    multiqc=1.15 \
    fastp \
    trimmomatic=0.39

# Create assembly environment  
echo "Creating assembly environment..."
conda create -n metagenome_assembly -c bioconda -c conda-forge -y \
    megahit=1.2.9 \
    spades=3.15.5 \
    quast=5.2.0 \
    bowtie2=2.5.1 \
    samtools=1.17 \
    bwa=0.7.17

# Create binning environment
echo "Creating binning environment..."
conda create -n metagenome_binning -c bioconda -c conda-forge -y \
    metabat2=2.15 \
    maxbin2=2.2.7 \
    das_tool=1.1.6 \
    checkm-genome=1.2.2 \
    vamb=4.1.3

# Create annotation environment
echo "Creating annotation environment..."
conda create -n metagenome_annotation -c bioconda -c conda-forge -y \
    prokka=1.14.6 \
    gtdbtk=2.3.2 \
    eggnog-mapper=2.1.9 \
    diamond=2.1.8

echo "Environment setup complete!"
echo "Activate environments using:"
echo "  conda activate metagenome_qc"
echo "  conda activate metagenome_assembly" 
echo "  conda activate metagenome_binning"
echo "  conda activate metagenome_annotation"
