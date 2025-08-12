#!/bin/bash

# Fix Differential Coverage for Binning
# Maps each sample individually to coassembly to get true coverage differences
# Author: Generated for farting_trees project

set -e

WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/03_assembly"
BINNING_DIR="${WORK_DIR}/results/04_binning"
DATA_DIR="${WORK_DIR}/data"
THREADS=8

# Assembly file
COASSEMBLY_FASTA="${ASSEMBLY_DIR}/megahit/coassembly_contigs.fa"

echo "=============================================="
echo "FIXING DIFFERENTIAL COVERAGE FOR BINNING"
echo "=============================================="
echo "Assembly: $COASSEMBLY_FASTA"
echo "Data: $DATA_DIR"
echo "Output: $BINNING_DIR"
echo ""

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_assembly

# Create mapping directory
mkdir -p ${BINNING_DIR}/individual_mapping

# Check if BWA index exists
if [ ! -f "${COASSEMBLY_FASTA}.bwt" ]; then
    echo "Creating BWA index for coassembly..."
    bwa index ${COASSEMBLY_FASTA}
fi

echo "Step 1: Individual sample mapping to coassembly..."

# Define sample pairs
declare -A SAMPLES=(
    ["53394"]="53394_R1_A15.fastq.gz 53394_R2_A15.fastq.gz"
    ["53395"]="53395_R1_A16.fastq.gz 53395_R2_A16.fastq.gz"
    ["53396"]="53396_R1_B6.fastq.gz 53396_R2_B6.fastq.gz"
    ["53397"]="53397_R1_B10.fastq.gz 53397_R2_B10.fastq.gz"
    ["53398"]="53398_R1_A16S.fastq.gz 53398_R2_A16S.fastq.gz"
    ["53399"]="53399_R1_B10S.fastq.gz 53399_R2_B10S.fastq.gz"
)

# Map each sample individually
for SAMPLE in "${!SAMPLES[@]}"; do
    echo "  Processing sample $SAMPLE..."
    
    # Get R1 and R2 files
    READ_FILES=(${SAMPLES[$SAMPLE]})
    R1_FILE="${DATA_DIR}/${READ_FILES[0]}"
    R2_FILE="${DATA_DIR}/${READ_FILES[1]}"
    
    # Check if files exist
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "    Warning: Files not found for $SAMPLE"
        echo "    Expected: $R1_FILE, $R2_FILE"
        continue
    fi
    
    echo "    R1: $(basename $R1_FILE)"
    echo "    R2: $(basename $R2_FILE)"
    
    # Map reads to coassembly
    echo "    Mapping reads..."
    bwa mem -t ${THREADS} \
            ${COASSEMBLY_FASTA} \
            ${R1_FILE} \
            ${R2_FILE} | \
    samtools view -bS - | \
    samtools sort -@ ${THREADS} -o ${BINNING_DIR}/individual_mapping/${SAMPLE}_coassembly_sorted.bam
    
    # Index BAM file
    echo "    Indexing BAM..."
    samtools index ${BINNING_DIR}/individual_mapping/${SAMPLE}_coassembly_sorted.bam
    
    echo "    ✓ Sample $SAMPLE mapping completed"
done

echo ""
echo "Step 2: Generate proper differential coverage matrix..."

# Remove old identical BAM files
echo "  Removing old identical BAM files..."
rm -f ${BINNING_DIR}/5339*_coassembly_sorted.bam*

# Copy new individual BAM files to binning directory
echo "  Installing new individual BAM files..."
cp ${BINNING_DIR}/individual_mapping/*_coassembly_sorted.bam* ${BINNING_DIR}/

# Generate new depth file with true differential coverage
echo "  Generating new depth matrix..."
conda activate metagenome_binning
jgi_summarize_bam_contig_depths --outputDepth ${BINNING_DIR}/coassembly_depth_corrected.txt \
                               ${BINNING_DIR}/individual_mapping/*_coassembly_sorted.bam

# Backup old depth file and replace
mv ${BINNING_DIR}/coassembly_depth.txt ${BINNING_DIR}/coassembly_depth_old_identical.txt
mv ${BINNING_DIR}/coassembly_depth_corrected.txt ${BINNING_DIR}/coassembly_depth.txt

echo ""
echo "Step 3: Verify differential coverage..."
echo "  Checking coverage diversity across samples..."

# Show coverage statistics
python3 -c "
import pandas as pd
import numpy as np

# Read depth file
try:
    df = pd.read_csv('${BINNING_DIR}/coassembly_depth.txt', sep='\t')
    
    # Extract sample columns (every 2nd column starting from 4th)
    sample_cols = [col for col in df.columns if col.endswith('_coassembly_sorted.bam') and not col.endswith('-var')]
    
    print('Coverage Statistics:')
    print('=' * 50)
    for col in sample_cols:
        sample_name = col.replace('_coassembly_sorted.bam', '')
        mean_cov = df[col].mean()
        std_cov = df[col].std()
        print(f'{sample_name}: Mean={mean_cov:.2f}, Std={std_cov:.2f}')
    
    # Check for diversity
    print('\nCoverage Diversity Check:')
    print('=' * 30)
    sample_means = [df[col].mean() for col in sample_cols]
    overall_std = np.std(sample_means)
    
    if overall_std < 0.1:
        print('⚠ WARNING: Low coverage diversity detected!')
        print('  Samples may still be too similar for effective binning.')
    else:
        print('✓ Good coverage diversity detected!')
        print(f'  Standard deviation between samples: {overall_std:.2f}')
        
except Exception as e:
    print(f'Error analyzing coverage: {e}')
    print('Please check the depth file manually.')
"

echo ""
echo "=============================================="
echo "DIFFERENTIAL COVERAGE FIX COMPLETED!"
echo "=============================================="
echo ""
echo "Results:"
echo "  Individual BAM files: ${BINNING_DIR}/individual_mapping/"
echo "  New depth matrix: ${BINNING_DIR}/coassembly_depth.txt"
echo "  Old (identical) depth matrix: ${BINNING_DIR}/coassembly_depth_old_identical.txt"
echo ""
echo "Next steps:"
echo "  1. Run binning pipeline again: ./code/04_binning.sh"
echo "  2. Check if you get more bins with proper differential coverage"
echo "  3. If still getting 0 bins, consider relaxed parameters"
echo ""
echo "Note: If coverage is still too similar across samples,"
echo "this might indicate that your samples have very similar"
echo "microbial communities, which is biologically possible."
