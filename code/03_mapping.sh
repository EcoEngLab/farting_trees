#!/bin/bash
# Read Mapping for Binning Pipeline
# Maps reads back to assemblies for coverage calculation

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/02_assembly"
READS_DIR="${WORK_DIR}/data"
RESULTS_DIR="${WORK_DIR}/results/04_binning"
THREADS=8

# Initialize conda
eval "$(conda shell.bash hook)"

echo "READ MAPPING PIPELINE"
echo "Assembly: $ASSEMBLY_DIR"
echo "Reads: $READS_DIR"
echo "Results: $RESULTS_DIR"

# Create results directory
mkdir -p "$RESULTS_DIR"

# Find read files
SAMPLES=(53394 53395 53396 53397 53398 53399)

# Function to map reads for each sample
map_reads() {
    local sample="$1"
    local assembly_file="${ASSEMBLY_DIR}/${sample}/contigs.fasta"
    local r1_file="${READS_DIR}/${sample}_R1_*.fastq.gz"
    local r2_file="${READS_DIR}/${sample}_R2_*.fastq.gz"
    
    # Expand wildcards
    r1_file=$(ls $r1_file 2>/dev/null | head -1)
    r2_file=$(ls $r2_file 2>/dev/null | head -1)
    
    if [ ! -f "$assembly_file" ]; then
        echo "Warning: Assembly not found for $sample: $assembly_file"
        return
    fi
    
    if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
        echo "Warning: Read files not found for $sample"
        echo "  R1: $r1_file"
        echo "  R2: $r2_file"
        return
    fi
    
    echo "Processing sample $sample..."
    echo "  Assembly: $assembly_file"
    echo "  R1: $r1_file"
    echo "  R2: $r2_file"
    
    # Activate BWA environment (assuming it exists or create it)
    if ! conda info --envs | grep -q bwa_env; then
        echo "Creating BWA environment..."
        conda create -n bwa_env -c bioconda bwa samtools -y
    fi
    
    conda activate bwa_env
    
    # Index assembly if not already done
    if [ ! -f "${assembly_file}.bwt" ]; then
        echo "Indexing assembly for $sample..."
        bwa index "$assembly_file"
    fi
    
    # Map reads
    echo "Mapping reads for $sample..."
    bwa mem -t $THREADS "$assembly_file" "$r1_file" "$r2_file" | \
        samtools view -@ $THREADS -bS - | \
        samtools sort -@ $THREADS -o "${RESULTS_DIR}/${sample}_sorted.bam" -
    
    # Index BAM file
    samtools index "${RESULTS_DIR}/${sample}_sorted.bam"
    
    echo "Completed mapping for $sample"
}

# Map reads for all samples
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "${RESULTS_DIR}/${sample}_sorted.bam" ]; then
        map_reads "$sample"
    else
        echo "BAM file already exists for $sample, skipping..."
    fi
done

echo "Read mapping completed!"
echo "BAM files created in: $RESULTS_DIR"
