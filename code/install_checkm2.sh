#!/bin/bash
# Install CheckM2 only
# CheckM2 is faster and more efficient than CheckM1

set -e

echo "Installing CheckM2..."

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Create CheckM2 environment
echo "Creating CheckM2 environment..."
conda create -n checkm2_env -c conda-forge -c bioconda checkm2 -y
echo "CheckM2 environment created"

# Test installation
echo ""
echo "Testing CheckM2 installation..."
conda activate checkm2_env
if command -v checkm2 &> /dev/null; then
    echo "✓ CheckM2 installed successfully"
    checkm2 --version
else
    echo "✗ CheckM2 installation failed"
    exit 1
fi

# Setup CheckM2 database
echo ""
echo "Setting up CheckM2 database..."
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
echo "CHECKM2 INSTALLATION COMPLETED!"
echo "==========================================="
echo ""
echo "Environment created: checkm2_env"
echo "Database location: $CHECKM2_DB"
echo ""
echo "To use CheckM2:"
echo "  conda activate checkm2_env"
echo "  checkm2 predict --input <bin_directory> --output-directory <output_dir>"
