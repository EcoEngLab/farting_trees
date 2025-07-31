#!/bin/bash

# Gradual environment setup to avoid dependency conflicts
# Author: Generated for farting_trees project
# Date: July 2025

set -e

echo "Setting up environments step by step to avoid conflicts..."

# Function to create environment with error handling
create_env_safe() {
    local env_name=$1
    local packages=("${@:2}")
    
    echo "Creating environment: $env_name"
    
    # Remove existing environment if it exists
    conda env remove -n $env_name -y 2>/dev/null || true
    
    # Create base environment
    conda create -n $env_name -c conda-forge -y python=3.9
    
    # Activate and install packages one by one
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $env_name
    
    for package in "${packages[@]}"; do
        echo "Installing $package..."
        conda install -c bioconda -c conda-forge -y $package || {
            echo "Failed to install $package with conda, trying mamba..."
            mamba install -c bioconda -c conda-forge -y $package || {
                echo "Warning: Could not install $package"
            }
        }
    done
    
    conda deactivate
    echo "Environment $env_name created successfully"
}

# Create base directory for results
mkdir -p results/{01_qc,02_assembly,03_binning,04_annotation,05_taxonomy}

# Install mamba if not available
if ! command -v mamba &> /dev/null; then
    echo "Installing mamba..."
    conda install mamba -n base -c conda-forge -y
fi

# Create QC environment
create_env_safe "metagenome_qc" "fastqc" "multiqc" "fastp" "trimmomatic"

# Create assembly environment
create_env_safe "metagenome_assembly" "megahit" "spades" "quast" "bowtie2" "samtools" "bwa"

# Create binning environment (this might take longer)
create_env_safe "metagenome_binning" "metabat2" "maxbin2" "checkm-genome"

# Create a separate environment for DAS Tool and VAMB (they can be problematic)
create_env_safe "metagenome_binning_extra" "das_tool" "vamb"

# Create annotation environment
create_env_safe "metagenome_annotation" "prokka" "diamond"

# Create separate environment for GTDB-Tk (it has many dependencies)
create_env_safe "gtdbtk_env" "gtdbtk"

# Create separate environment for eggNOG-mapper
create_env_safe "eggnog_env" "eggnog-mapper"

echo ""
echo "=========================================="
echo "Environment setup completed!"
echo "=========================================="
echo ""
echo "Available environments:"
echo "  conda activate metagenome_qc"
echo "  conda activate metagenome_assembly"
echo "  conda activate metagenome_binning"
echo "  conda activate metagenome_binning_extra"
echo "  conda activate metagenome_annotation"
echo "  conda activate gtdbtk_env"
echo "  conda activate eggnog_env"
echo ""
echo "Note: Some tools are in separate environments to avoid conflicts."
