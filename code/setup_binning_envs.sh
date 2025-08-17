#!/bin/bash
# Setup script for binning environments
# Creates separate conda environments for each binning tool

set -e

echo "Setting up conda environments for metagenomic binning tools..."
echo "This may take 15-30 minutes depending on your internet connection."

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# 1. Create MetaBAT2 environment
echo "Creating MetaBAT2 environment..."
conda create -n metabat2_env -c bioconda metabat2 -y
echo "MetaBAT2 environment created"

# 2. Create MaxBin2 environment
echo "Creating MaxBin2 environment..."
conda create -n maxbin2_env -c bioconda maxbin2 -y
echo "MaxBin2 environment created"

# 3. Create VAMB environment
echo "Creating VAMB environment..."
conda create -n vamb_env -c bioconda vamb -y
echo "VAMB environment created"

# 4. Create DAS Tool environment
echo "Creating DAS Tool environment..."
# Try comprehensive installation, fallback to step-by-step if it fails
if ! conda create -n dastool_env -c conda-forge -c bioconda r-base r-data.table r-ggplot2 r-magrittr das_tool -y; then
    echo "Comprehensive installation failed, trying step-by-step approach..."
    
    # Create base environment with R
    conda create -n dastool_env -c conda-forge r-base -y
    
    # Activate and install R packages
    conda activate dastool_env
    conda install -c conda-forge r-data.table r-ggplot2 r-magrittr -y
    
    # Install DAS Tool
    if ! conda install -c bioconda das_tool -y; then
        echo "Warning: DAS Tool installation failed. The pipeline will skip DAS Tool optimization."
        echo "You can manually install DAS Tool later or run without optimization."
    fi
else
    echo "DAS Tool environment created successfully"
fi

# 5. Create CheckM2 environment
echo "Creating CheckM2 environment..."
conda create -n checkm2_env -c conda-forge -c bioconda checkm2 -y
echo "CheckM2 environment created"

# Test installations
echo ""
echo "Testing installations..."

# Test MetaBAT2
conda activate metabat2_env
if command -v metabat2 &> /dev/null; then
    echo "✓ MetaBAT2 installed successfully"
    metabat2 --help | head -3
else
    echo "✗ MetaBAT2 installation failed"
fi

# Test MaxBin2
conda activate maxbin2_env
if command -v run_MaxBin.pl &> /dev/null; then
    echo "✓ MaxBin2 installed successfully"
else
    echo "✗ MaxBin2 installation failed"
fi

# Test VAMB
conda activate vamb_env
if command -v vamb &> /dev/null; then
    echo "✓ VAMB installed successfully"
    vamb --version
else
    echo "✗ VAMB installation failed"
fi

# Test DAS Tool
conda activate dastool_env
if command -v DAS_Tool &> /dev/null; then
    echo "✓ DAS Tool installed successfully"
else
    echo "✗ DAS Tool installation failed"
fi

# Test CheckM2
conda activate checkm2_env
if command -v checkm2 &> /dev/null; then
    echo "✓ CheckM2 installed successfully"
    checkm2 --version
else
    echo "✗ CheckM2 installation failed"
fi

# Setup CheckM2 database
echo ""
echo "Setting up CheckM2 database..."
conda activate checkm2_env
CHECKM2_DB="$HOME/.checkm2db"
if [ ! -d "$CHECKM2_DB" ]; then
    echo "Downloading CheckM2 database (much smaller than CheckM1)..."
    checkm2 database --download --path "$CHECKM2_DB"
    echo "CheckM2 database setup completed"
else
    echo "CheckM2 database already configured"
fi

echo ""
echo "===========================================" 
echo "ENVIRONMENT SETUP COMPLETED!"
echo "==========================================="
echo ""
echo "Created environments:"
echo "  - metabat2_env: MetaBAT2"
echo "  - maxbin2_env: MaxBin2" 
echo "  - vamb_env: VAMB"
echo "  - dastool_env: DAS Tool"
echo "  - checkm2_env: CheckM2"
echo ""
echo "You can now run the binning pipeline:"
echo "  bash code/04_binning.sh all"
