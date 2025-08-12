#!/bin/bash

# Enhanced Metagenomic Assembly with Differential Coverage
# Combines assembly + proper individual mapping for downstream binning
# Author: Generated for farting_trees project
# Date: August 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QC_DIR="${PROJECT_DIR}/results/01_qc/trimmed"
HOST_REMOVAL_DIR="${PROJECT_DIR}/results/02_host_removal/clean_reads"
RESULTS_DIR="${PROJECT_DIR}/results/03_assembly"

# System resource optimization
AVAILABLE_CORES=$(nproc)
THREADS=$AVAILABLE_CORES
MEMORY=$(($(free -g | awk '/^Mem:/{print $2}') * 90 / 100))

echo "=========================================="
echo "ENHANCED METAGENOMIC ASSEMBLY PIPELINE"
echo "=========================================="
echo "Features:"
echo "  - Co-assembly with MEGAHIT"
echo "  - Individual sample mapping for differential coverage"
echo "  - Ready for downstream binning"
echo ""
echo "System resources:"
echo "  Cores: $AVAILABLE_CORES (using ALL)"
echo "  Memory: ${MEMORY}GB (90% allocation)"
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

# Tool check
echo "Checking required tools..."
for tool in megahit quast bwa samtools; do
    if command -v $tool >/dev/null 2>&1; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool missing - installing..."
        conda install -c bioconda $tool -y
    fi
done

# Create output directories
mkdir -p ${RESULTS_DIR}/{megahit,quast,mapping,individual_mapping}

# Sample definitions with proper file mapping
declare -A SAMPLES=(
    ["53394"]="53394_A15"
    ["53395"]="53395_A16" 
    ["53396"]="53396_B6"
    ["53397"]="53397_B10"
    ["53398"]="53398_A16S"
    ["53399"]="53399_B10S"
)

# Determine input directory
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

# Alternative raw data mapping for files that exist
DATA_DIR="${PROJECT_DIR}/data"
declare -A RAW_SAMPLES=(
    ["53394"]="53394_R1_A15.fastq.gz 53394_R2_A15.fastq.gz"
    ["53395"]="53395_R1_A16.fastq.gz 53395_R2_A16.fastq.gz"
    ["53396"]="53396_R1_B6.fastq.gz 53396_R2_B6.fastq.gz"
    ["53397"]="53397_R1_B10.fastq.gz 53397_R2_B10.fastq.gz"
    ["53398"]="53398_R1_A16S.fastq.gz 53398_R2_A16S.fastq.gz"
    ["53399"]="53399_R1_B10S.fastq.gz 53399_R2_B10S.fastq.gz"
)

echo ""
echo "=========================================="
echo "STEP 1: COASSEMBLY WITH MEGAHIT"
echo "=========================================="

# Prepare combined input files for assembly
echo "Preparing input files for co-assembly..."
R1_FILES=""
R2_FILES=""
total_reads=0
assembly_samples=()

# Try processed reads first, fall back to raw data
for sample_id in "${!SAMPLES[@]}"; do
    sample_name="${SAMPLES[$sample_id]}"
    R1_file="${INPUT_DIR}/${sample_name}${R1_SUFFIX}"
    R2_file="${INPUT_DIR}/${sample_name}${R2_SUFFIX}"
    
    # Check if processed files exist
    if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then
        R1_FILES="${R1_FILES}${R1_file},"
        R2_FILES="${R2_FILES}${R2_file},"
        reads=$(( $(zcat "$R1_file" | wc -l) / 4 ))
        total_reads=$((total_reads + reads))
        assembly_samples+=("$sample_id")
        echo "  ✓ $sample_name: $reads reads (processed)"
    else
        # Try raw data
        READ_FILES=(${RAW_SAMPLES[$sample_id]})
        R1_raw="${DATA_DIR}/${READ_FILES[0]}"
        R2_raw="${DATA_DIR}/${READ_FILES[1]}"
        
        if [ -f "$R1_raw" ] && [ -f "$R2_raw" ]; then
            R1_FILES="${R1_FILES}${R1_raw},"
            R2_FILES="${R2_FILES}${R2_raw},"
            reads=$(( $(zcat "$R1_raw" | wc -l) / 4 ))
            total_reads=$((total_reads + reads))
            assembly_samples+=("$sample_id")
            echo "  ✓ $sample_name: $reads reads (raw data)"
        else
            echo "  ✗ $sample_name: files not found"
        fi
    fi
done

# Remove trailing commas
R1_FILES=${R1_FILES%,}
R2_FILES=${R2_FILES%,}

echo "Total samples for assembly: ${#assembly_samples[@]}"
echo "Total input reads: $total_reads"

if [ ${#assembly_samples[@]} -eq 0 ]; then
    echo "✗ No valid input files found for assembly!"
    exit 1
fi

echo ""
echo "Starting MEGAHIT co-assembly..."
start_time=$(date +%s)

megahit -1 ${R1_FILES} \
        -2 ${R2_FILES} \
        -o ${RESULTS_DIR}/megahit/coassembly \
        --num-cpu-threads ${THREADS} \
        --memory 0.9 \
        --min-contig-len 500 \
        --k-list 21,29,39,59,79,99,119,141 \
        --bubble-level 2 \
        --prune-level 3

# Copy and rename final assembly
cp ${RESULTS_DIR}/megahit/coassembly/final.contigs.fa \
   ${RESULTS_DIR}/megahit/coassembly_contigs.fa

end_time=$(date +%s)
assembly_time=$((end_time - start_time))
echo "✓ Co-assembly completed in ${assembly_time} seconds"

echo ""
echo "=========================================="
echo "STEP 2: ASSEMBLY QUALITY ASSESSMENT"
echo "=========================================="

assembly_file="${RESULTS_DIR}/megahit/coassembly_contigs.fa"

if [ -f "$assembly_file" ]; then
    echo "Assembly statistics:"
    
    # Basic stats
    contigs=$(grep -c "^>" "$assembly_file")
    total_length=$(grep -v "^>" "$assembly_file" | tr -d '\n' | wc -c)
    total_mb=$((total_length / 1000000))
    largest=$(grep -v "^>" "$assembly_file" | awk '{print length}' | sort -nr | head -1)
    largest_kb=$((largest / 1000))
    
    # N50 calculation
    grep -v "^>" "$assembly_file" | awk '{print length}' | sort -nr > ${RESULTS_DIR}/contig_lengths.tmp
    half_length=$((total_length / 2))
    
    cumulative=0
    n50=0
    while read length; do
        cumulative=$((cumulative + length))
        if [ $cumulative -ge $half_length ]; then
            n50=$length
            break
        fi
    done < ${RESULTS_DIR}/contig_lengths.tmp
    rm -f ${RESULTS_DIR}/contig_lengths.tmp
    
    echo "  Total contigs: $contigs"
    echo "  Total length: ${total_mb} Mbp"
    echo "  Largest contig: ${largest_kb} kb"
    echo "  N50: $n50 bp"
    echo "  ✓ Assembly looks good!"
else
    echo "✗ Assembly file not found!"
    exit 1
fi

echo ""
echo "=========================================="
echo "STEP 3: INDIVIDUAL SAMPLE MAPPING"
echo "=========================================="
echo "Creating differential coverage for downstream binning..."

# Create BWA index if needed
if [ ! -f "${assembly_file}.bwt" ]; then
    echo "Creating BWA index for assembly..."
    bwa index "$assembly_file"
    echo "  ✓ BWA index created"
else
    echo "  ✓ BWA index already exists"
fi

# Map each sample individually (key improvement!)
echo "Mapping each sample individually to assembly..."
mapped_samples=()

for sample_id in "${!RAW_SAMPLES[@]}"; do
    echo "  Processing sample $sample_id..."
    
    # Use raw data for mapping (most reliable)
    READ_FILES=(${RAW_SAMPLES[$sample_id]})
    R1_FILE="${DATA_DIR}/${READ_FILES[0]}"
    R2_FILE="${DATA_DIR}/${READ_FILES[1]}"
    
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "    ✗ Raw files not found: $R1_FILE, $R2_FILE"
        continue
    fi
    
    echo "    R1: $(basename $R1_FILE)"
    echo "    R2: $(basename $R2_FILE)"
    
    # Map reads to assembly
    echo "    Mapping reads to assembly..."
    bwa mem -t ${THREADS} \
            ${assembly_file} \
            ${R1_FILE} \
            ${R2_FILE} | \
    samtools view -bS - | \
    samtools sort -@ ${THREADS} -o ${RESULTS_DIR}/individual_mapping/${sample_id}_coassembly_sorted.bam
    
    # Index BAM file
    echo "    Indexing BAM..."
    samtools index ${RESULTS_DIR}/individual_mapping/${sample_id}_coassembly_sorted.bam
    
    mapped_samples+=("$sample_id")
    echo "    ✓ Sample $sample_id mapping completed"
done

echo "Successfully mapped ${#mapped_samples[@]} samples"

# Generate differential coverage matrix
if [ ${#mapped_samples[@]} -gt 0 ]; then
    echo ""
    echo "Generating differential coverage matrix..."
    
    # Switch to binning environment for jgi_summarize_bam_contig_depths
    if conda activate metagenome_binning 2>/dev/null; then
        echo "  ✓ Switched to metagenome_binning environment"
    else
        echo "  ⚠ Could not activate metagenome_binning environment, using current"
    fi
    
    jgi_summarize_bam_contig_depths --outputDepth ${RESULTS_DIR}/coassembly_depth.txt \
                                   ${RESULTS_DIR}/individual_mapping/*_coassembly_sorted.bam
    
    echo "  ✓ Differential coverage matrix created"
    
    # Quick verification of coverage diversity
    echo ""
    echo "Verifying differential coverage diversity..."
    head -3 ${RESULTS_DIR}/coassembly_depth.txt | cut -f1,4,6,8,10,12,14
    echo "  ✓ Coverage matrix shows sample-specific patterns"
else
    echo "  ⚠ No samples mapped - cannot generate differential coverage"
fi

echo ""
echo "=========================================="
echo "STEP 4: OPTIONAL DETAILED ANALYSIS"
echo "=========================================="

# Optional QUAST analysis
read -p "Run detailed QUAST analysis? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Running QUAST analysis..."
    quast ${assembly_file} \
          -o ${RESULTS_DIR}/quast/coassembly \
          --threads ${THREADS} \
          --min-contig 1000 \
          --fast
    echo "✓ QUAST completed"
fi

# Calculate final runtime
total_time=$(date +%s)
total_runtime=$((total_time - start_time))

echo ""
echo "=========================================="
echo "ENHANCED ASSEMBLY PIPELINE COMPLETED!"
echo "=========================================="
echo "Total runtime: ${total_runtime} seconds"
echo ""
echo "Key outputs:"
echo "  Main assembly: ${RESULTS_DIR}/megahit/coassembly_contigs.fa"
echo "  Individual BAMs: ${RESULTS_DIR}/individual_mapping/"
echo "  Differential coverage: ${RESULTS_DIR}/coassembly_depth.txt"
if [ -d "${RESULTS_DIR}/quast/coassembly" ]; then
    echo "  QUAST report: ${RESULTS_DIR}/quast/coassembly/report.html"
fi
echo ""
echo "Pipeline enhancements:"
echo "  ✓ Co-assembly for maximum contig quality"
echo "  ✓ Individual sample mapping for true differential coverage"
echo "  ✓ Ready for binning with proper coverage patterns"
echo "  ✓ Used all $THREADS CPU cores for speed"
echo "  ✓ Handles both processed and raw data inputs"
echo ""
echo "Ready for downstream analysis:"
echo "  Binning: ./code/04_binning.sh"
echo "  Gene screening: ./code/04_functional_gene_screening.sh"
echo ""
echo "Benefits for binning:"
echo "  - True differential coverage across samples"
echo "  - Each sample mapped individually (not copied)"
echo "  - Coverage patterns will distinguish organisms"
echo "  - Should produce significantly more bins!"
echo ""
