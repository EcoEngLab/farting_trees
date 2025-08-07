#!/bin/bash

# Standalone Coverage Calculation for Assembled Contigs
# Author: Generated for farting_trees project
# Date: August 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HOST_REMOVAL_DIR="${PROJECT_DIR}/results/02_host_removal/clean_reads"
QC_DIR="${PROJECT_DIR}/results/01_qc/trimmed"
RESULTS_DIR="${PROJECT_DIR}/results/03_assembly"

# Use all available cores for speed (but limit to prevent memory issues)
AVAILABLE_CORES=$(nproc)
# Limit threads to prevent memory issues with BWA
if [ $AVAILABLE_CORES -gt 8 ]; then
    THREADS=8
else
    THREADS=$AVAILABLE_CORES
fi

echo "=========================================="
echo "STANDALONE COVERAGE CALCULATION"
echo "=========================================="
echo "Using $THREADS CPU cores (limited for stability)"

# Check available memory
AVAILABLE_MEMORY=$(free -g | awk '/^Mem:/{print $2}')
echo "Available memory: ${AVAILABLE_MEMORY}GB"

if [ $AVAILABLE_MEMORY -lt 8 ]; then
    echo "⚠ Warning: Low memory detected. BWA may fail with large datasets."
    echo "  Consider running with fewer threads or on a machine with more RAM."
fi
echo ""

# Activate assembly environment
source $(conda info --base)/etc/profile.d/conda.sh
if conda activate metagenome_assembly_new 2>/dev/null; then
    echo "✓ Using metagenome_assembly_new environment"
elif conda activate metagenome_assembly 2>/dev/null; then
    echo "✓ Using metagenome_assembly environment"
else
    echo "✗ No suitable environment found. Please run setup script first."
    exit 1
fi

# Check required tools
echo "Checking tools..."
for tool in bwa samtools; do
    if command -v $tool >/dev/null 2>&1; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool missing - installing..."
        conda install -c bioconda $tool -y
    fi
done

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

# Check for input reads (same logic as assembly script)
echo "Locating input reads..."
if [ -f "${HOST_REMOVAL_DIR}/53394_A15_R1_clean.fastq.gz" ]; then
    INPUT_DIR="$HOST_REMOVAL_DIR"
    R1_SUFFIX="_R1_clean.fastq.gz"
    R2_SUFFIX="_R2_clean.fastq.gz"
    echo "✓ Using host-removed reads: $INPUT_DIR"
else
    INPUT_DIR="$QC_DIR"
    R1_SUFFIX="_R1_trimmed.fastq.gz"
    R2_SUFFIX="_R2_trimmed.fastq.gz"
    echo "✓ Using trimmed reads: $INPUT_DIR"
fi

# Prepare input file lists
R1_FILES=()
R2_FILES=()
total_reads=0

for sample in "${SAMPLES[@]}"; do
    R1_file="${INPUT_DIR}/${sample}${R1_SUFFIX}"
    R2_file="${INPUT_DIR}/${sample}${R2_SUFFIX}"
    
    if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then
        R1_FILES+=("$R1_file")
        R2_FILES+=("$R2_file")
        reads=$(( $(zcat "$R1_file" | wc -l) / 4 ))
        total_reads=$((total_reads + reads))
        echo "  ✓ $sample: $reads reads"
    else
        echo "  ✗ $sample: files not found"
        exit 1
    fi
done

echo "Total input reads: $total_reads"
echo ""

# Verify assembly and indexes exist
assembly_file="${RESULTS_DIR}/megahit/coassembly_contigs.fa"

echo "Verifying assembly and indexes..."
if [ ! -f "$assembly_file" ]; then
    echo "✗ Assembly file not found: $assembly_file"
    echo "Please run the assembly script first!"
    exit 1
fi

# Check if BWA indexes exist
if [ ! -f "${assembly_file}.bwt" ]; then
    echo "BWA indexes not found. Creating them..."
    bwa index "$assembly_file"
    echo "✓ BWA index created"
else
    echo "✓ BWA indexes already exist"
fi

# Check assembly file size
assembly_size=$(stat -c%s "$assembly_file" 2>/dev/null || echo 0)
if [ $assembly_size -lt 1000 ]; then
    echo "✗ Assembly file too small (${assembly_size} bytes)"
    exit 1
fi

echo "✓ Assembly file verified: $(( assembly_size / 1024 / 1024 )) MB"
echo ""

# Create mapping directory
mkdir -p "${RESULTS_DIR}/mapping"

# Start coverage calculation
start_time=$(date +%s)

echo "=========================================="
echo "MAPPING READS TO ASSEMBLY"
echo "=========================================="

# Create temporary SAM file first to check for errors
temp_sam="${RESULTS_DIR}/mapping/temp_alignment.sam"

# Run BWA mem with error checking
echo "Running BWA alignment..."
echo "  Mapping ${total_reads} reads to assembly..."
echo "  Assembly file: $assembly_file"
echo "  R1 files: ${R1_FILES[*]}"
echo "  R2 files: ${R2_FILES[*]}"
echo "  Threads: ${THREADS}"
echo ""

# Test BWA command first - process each sample pair individually
echo "Processing samples individually and merging results..."
echo "Assembly file: $assembly_file"
echo "Number of sample pairs: ${#R1_FILES[@]}"
echo ""

# Create temporary directory for individual SAM files
temp_dir="${RESULTS_DIR}/mapping/temp_sams"
mkdir -p "$temp_dir"

# Process each sample pair individually
sam_files=()
for i in "${!R1_FILES[@]}"; do
    r1_file="${R1_FILES[$i]}"
    r2_file="${R2_FILES[$i]}"
    sample_name=$(basename "$r1_file" | sed 's/_R1_clean\.fastq\.gz//')
    individual_sam="${temp_dir}/${sample_name}.sam"
    
    echo "Processing sample $((i+1))/${#R1_FILES[@]}: $sample_name"
    echo "  R1: $(basename "$r1_file")"
    echo "  R2: $(basename "$r2_file")"
    
    if ! bwa mem -t ${THREADS} -M "$assembly_file" "$r1_file" "$r2_file" > "$individual_sam" 2>"${temp_dir}/${sample_name}_error.log"; then
        echo "✗ BWA alignment failed for $sample_name"
        echo "Error log: ${temp_dir}/${sample_name}_error.log"
        cat "${temp_dir}/${sample_name}_error.log"
        exit 1
    fi
    
    echo "  ✓ Alignment completed for $sample_name"
    sam_files+=("$individual_sam")
done

echo ""
echo "Merging SAM files..."
# Merge all SAM files - take header from first file, then all alignments
head -n 100 "${sam_files[0]}" | grep "^@" > "$temp_sam"
for sam_file in "${sam_files[@]}"; do
    grep -v "^@" "$sam_file" >> "$temp_sam"
done

# Clean up individual SAM files
rm -rf "$temp_dir"

echo "✓ BWA alignment completed"

# Check if SAM file has valid header
if head -1 "$temp_sam" | grep -q "^@"; then
    echo "✓ Valid SAM header detected"
    
    # Convert to BAM and sort
    echo "Converting to sorted BAM..."
    samtools view -@ ${THREADS} -bS -F 4 "$temp_sam" 2>/dev/null | \
    samtools sort -@ ${THREADS} -o "${RESULTS_DIR}/mapping/coassembly_sorted.bam" 2>/dev/null
    
    # Clean up temporary file
    rm -f "$temp_sam"
    echo "✓ BAM file created"
    
else
    echo "✗ Invalid SAM header, alignment may have failed"
    rm -f "$temp_sam"
    exit 1
fi

# Index BAM file
if [ -f "${RESULTS_DIR}/mapping/coassembly_sorted.bam" ]; then
    echo "Indexing BAM file..."
    samtools index "${RESULTS_DIR}/mapping/coassembly_sorted.bam" 2>/dev/null && {
        echo "✓ BAM indexing completed"
    } || {
        echo "✗ BAM indexing failed"
        exit 1
    }
else
    echo "✗ BAM file not created"
    exit 1
fi

# Calculate basic coverage statistics
echo ""
echo "=========================================="
echo "COVERAGE STATISTICS"
echo "=========================================="

bam_file="${RESULTS_DIR}/mapping/coassembly_sorted.bam"

# Basic mapping stats
total_mapped=$(samtools view -c -F 4 "$bam_file")
total_reads_mapped=$((total_mapped / 2))  # Paired reads
mapping_rate=$(echo "scale=2; $total_reads_mapped * 100 / $total_reads" | bc -l)

echo "Mapping summary:"
echo "  Total input reads: $total_reads"
echo "  Mapped reads: $total_reads_mapped"
echo "  Mapping rate: ${mapping_rate}%"

# Coverage depth calculation
echo ""
echo "Calculating coverage depth..."
samtools depth "$bam_file" > "${RESULTS_DIR}/mapping/coverage_depth.txt"

# Basic coverage stats
if [ -f "${RESULTS_DIR}/mapping/coverage_depth.txt" ]; then
    avg_coverage=$(awk '{sum+=$3; count++} END {if(count>0) print sum/count; else print 0}' "${RESULTS_DIR}/mapping/coverage_depth.txt")
    max_coverage=$(awk 'BEGIN{max=0} {if($3>max) max=$3} END{print max}' "${RESULTS_DIR}/mapping/coverage_depth.txt")
    
    echo "  Average coverage: ${avg_coverage}x"
    echo "  Maximum coverage: ${max_coverage}x"
    
    # Coverage distribution
    echo ""
    echo "Coverage distribution:"
    awk '{
        if($3 == 0) zero++
        else if($3 < 5) low++
        else if($3 < 20) medium++
        else high++
        total++
    } END {
        printf "  0x coverage: %.1f%% (%d positions)\n", zero*100/total, zero
        printf "  1-4x coverage: %.1f%% (%d positions)\n", low*100/total, low
        printf "  5-19x coverage: %.1f%% (%d positions)\n", medium*100/total, medium
        printf "  ≥20x coverage: %.1f%% (%d positions)\n", high*100/total, high
    }' "${RESULTS_DIR}/mapping/coverage_depth.txt"
fi

# Final summary
end_time=$(date +%s)
runtime=$((end_time - start_time))

echo ""
echo "=========================================="
echo "COVERAGE CALCULATION COMPLETED!"
echo "=========================================="
echo "Runtime: ${runtime} seconds"
echo ""
echo "Output files:"
echo "  Sorted BAM: ${RESULTS_DIR}/mapping/coassembly_sorted.bam"
echo "  BAM index: ${RESULTS_DIR}/mapping/coassembly_sorted.bam.bai"
echo "  Coverage depth: ${RESULTS_DIR}/mapping/coverage_depth.txt"
echo ""
echo "Next steps:"
echo "  - Use coverage data for genome binning"
echo "  - Analyze coverage patterns for abundance estimation"
echo "  - Proceed with: ./code/06_binning.sh"
echo ""
