#!/bin/bash

# Comprehensive environment setup for metagenomic analysis
# Combines features from all setup approaches with error handling
# Author: Generated for farting_trees project
# Date: July 2025

set -e

echo "=============================================="
echo "Metagenomic Analysis Environment Setup"
echo "=============================================="
echo "Setting up conda environments with intelligent fallbacks..."
echo ""

# Configuration
SETUP_MODE="auto"  # auto, fast, safe
USE_MAMBA=true

# Function to detect best setup approach
detect_setup_mode() {
    echo "Detecting optimal setup approach..."
    
    # Check available resources
    AVAILABLE_MEMORY=$(free -g | awk '/^Mem:/{print $2}')
    AVAILABLE_CORES=$(nproc)
    
    echo "  Available memory: ${AVAILABLE_MEMORY}GB"
    echo "  Available cores: ${AVAILABLE_CORES}"
    
    # Check if mamba is available
    if command -v mamba >/dev/null 2>&1; then
        echo "  Mamba found: Using fast installation"
        USE_MAMBA=true
        SETUP_MODE="fast"
    elif [ "$AVAILABLE_MEMORY" -lt 4 ]; then
        echo "  Low memory detected: Using safe installation"
        SETUP_MODE="safe"
    else
        echo "  Using standard installation"
        SETUP_MODE="standard"
    fi
    
    echo "  Setup mode: $SETUP_MODE"
    echo ""
}

# Function to install mamba if needed
install_mamba() {
    if ! command -v mamba >/dev/null 2>&1; then
        echo "Installing mamba for faster dependency resolution..."
        conda install mamba -n base -c conda-forge -y
        echo "✓ Mamba installed successfully"
    else
        echo "✓ Mamba already available"
    fi
}

# Function to create environment with multiple fallback strategies
create_env_intelligent() {
    local env_name=$1
    local description=$2
    shift 2
    local packages=("$@")
    
    echo "----------------------------------------"
    echo "Creating environment: $env_name ($description)"
    echo "Packages: ${packages[*]}"
    
    # Remove existing environment if it exists
    conda env remove -n $env_name -y 2>/dev/null || true
    
    # Strategy 1: Try mamba with all packages at once (fastest)
    if [ "$USE_MAMBA" = true ] && [ "$SETUP_MODE" = "fast" ]; then
        echo "  Strategy 1: Mamba bulk installation..."
        if mamba create -n $env_name -c bioconda -c conda-forge -y "${packages[@]}" 2>/dev/null; then
            echo "  ✓ Environment $env_name created successfully with mamba"
            return 0
        else
            echo "  ✗ Mamba bulk installation failed, trying conda..."
        fi
    fi
    
    # Strategy 2: Try conda with all packages at once
    if [ "$SETUP_MODE" != "safe" ]; then
        echo "  Strategy 2: Conda bulk installation..."
        if conda create -n $env_name -c bioconda -c conda-forge -y "${packages[@]}" 2>/dev/null; then
            echo "  ✓ Environment $env_name created successfully with conda"
            return 0
        else
            echo "  ✗ Conda bulk installation failed, trying safe mode..."
        fi
    fi
    
    # Strategy 3: Safe mode - create base environment and install packages individually
    echo "  Strategy 3: Safe installation (individual packages)..."
    
    # Create base environment
    conda create -n $env_name -c conda-forge -y python=3.9
    
    # Activate environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $env_name
    
    local failed_packages=()
    local installer="conda"
    
    # Try mamba first if available
    if [ "$USE_MAMBA" = true ]; then
        installer="mamba"
    fi
    
    for package in "${packages[@]}"; do
        echo "    Installing $package with $installer..."
        if $installer install -c bioconda -c conda-forge -y $package 2>/dev/null; then
            echo "    ✓ $package installed successfully"
        else
            echo "    ✗ $package failed with $installer"
            failed_packages+=("$package")
            
            # Try alternative installer
            if [ "$installer" = "mamba" ]; then
                echo "      Trying with conda..."
                if conda install -c bioconda -c conda-forge -y $package 2>/dev/null; then
                    echo "    ✓ $package installed with conda"
                    failed_packages=("${failed_packages[@]/$package}")
                fi
            fi
        fi
    done
    
    conda deactivate
    
    if [ ${#failed_packages[@]} -eq 0 ]; then
        echo "  ✓ Environment $env_name created successfully (safe mode)"
    else
        echo "  ⚠ Environment $env_name created with some missing packages:"
        printf "    %s\n" "${failed_packages[@]}"
        echo "    You can install these manually later if needed"
    fi
    
    return 0
}

# Create directory structure
echo "Creating directory structure..."
mkdir -p results/{01_qc,02_host_removal,03_assembly,04_functional_genes,05_binning,06_annotation,07_taxonomy}
echo "✓ Directory structure created"
echo ""

# Detect optimal setup approach
detect_setup_mode

# Install mamba if beneficial
if [ "$USE_MAMBA" = true ]; then
    install_mamba
fi

echo ""
echo "Creating conda environments..."
echo ""

# Environment 1: Quality Control and Preprocessing
create_env_intelligent "metagenome_qc" "Quality Control & Preprocessing" \
    "fastqc=0.12.1" \
    "multiqc=1.15" \
    "fastp" \
    "trimmomatic=0.39"

# Environment 2: Assembly and Host Removal  
create_env_intelligent "metagenome_assembly_new" "Assembly & Host Removal" \
    "megahit=1.2.9" \
    "spades=3.15.5" \
    "quast=5.2.0" \
    "bowtie2=2.5.1" \
    "samtools=1.17" \
    "bwa=0.7.17" \
    "seqkit" \
    "parallel"

# Environment 3: Functional Gene Analysis
create_env_intelligent "metagenome_functional" "Functional Gene Analysis" \
    "diamond=2.1.8" \
    "hmmer=3.3.2" \
    "blast=2.14.0" \
    "seqkit" \
    "prodigal"

# Environment 4: Binning
create_env_intelligent "metagenome_binning" "Genome Binning" \
    "metabat2=2.15" \
    "maxbin2=2.2.7" \
    "checkm-genome=1.2.2" \
    "drep=3.4.3"

# Environment 5: Advanced Binning (separate to avoid conflicts)
create_env_intelligent "metagenome_binning_advanced" "Advanced Binning Tools" \
    "das_tool=1.1.6" \
    "vamb=4.1.3" \
    "concoct=1.1.0"

# Environment 6: Annotation
create_env_intelligent "metagenome_annotation" "Genome Annotation" \
    "prokka=1.14.6" \
    "diamond=2.1.8" \
    "hmmer=3.3.2"

# Environment 7: Taxonomy (separate due to large databases)
create_env_intelligent "metagenome_taxonomy" "Taxonomic Classification" \
    "gtdbtk=2.3.2" \
    "kraken2=2.1.2" \
    "bracken=2.8"

# Environment 8: Metabolic Analysis
create_env_intelligent "metagenome_metabolism" "Metabolic Analysis" \
    "eggnog-mapper=2.1.9" \
    "kegg-decoder" \
    "metabolic"

echo ""
echo "=============================================="
echo "ENVIRONMENT SETUP COMPLETED!"
echo "=============================================="
echo ""
echo "✓ Created environments:"
echo "  metagenome_qc - Quality control and trimming"
echo "  metagenome_assembly_new - Assembly and host removal"  
echo "  metagenome_functional - Functional gene screening"
echo "  metagenome_binning - Basic genome binning"
echo "  metagenome_binning_advanced - Advanced binning tools"
echo "  metagenome_annotation - Genome annotation"
echo "  metagenome_taxonomy - Taxonomic classification"
echo "  metagenome_metabolism - Metabolic pathway analysis"
echo ""
echo "Usage examples:"
echo "  conda activate metagenome_qc"
echo "  conda activate metagenome_assembly_new"
echo "  conda activate metagenome_functional"
echo ""
echo "Next steps:"
echo "  1. Run quality control: ./code/02_quality_control.sh"
echo "  2. Run host removal: ./code/03_host_removal.sh"  
echo "  3. Run assembly: ./code/04_assembly.sh"
echo ""
echo "Note: Some tools are in separate environments to prevent conflicts."
echo "This ensures maximum compatibility and stability."
echo ""
