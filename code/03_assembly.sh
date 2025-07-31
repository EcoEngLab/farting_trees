#!/bin/bash

# Metagenomic assembly using MEGAHIT and metaSPAdes
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QC_DIR="${PROJECT_DIR}/results/01_qc/trimmed"
HOST_REMOVAL_DIR="${PROJECT_DIR}/results/02_host_removal/clean_reads"
RESULTS_DIR="${PROJECT_DIR}/results/03_assembly"

# Optimize CPU usage - detect available cores
AVAILABLE_CORES=$(nproc)
THREADS=$((AVAILABLE_CORES > 16 ? 16 : AVAILABLE_CORES))  # Cap at 16 for stability
MEMORY=$(($(free -g | awk '/^Mem:/{print $2}') * 80 / 100))  # Use 80% of available RAM

echo "System resources detected:"
echo "  Available cores: $AVAILABLE_CORES"
echo "  Using threads: $THREADS"
echo "  Available memory: ${MEMORY}GB"

# Activate assembly environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_assembly_new

# Install GNU parallel if not available (for faster processing)
if ! command -v parallel >/dev/null 2>&1; then
    echo "Installing GNU parallel for faster processing..."
    conda install -c conda-forge parallel -y
fi

echo "Starting metagenomic assembly..."

# Create output directories
mkdir -p ${RESULTS_DIR}/{megahit,spades,quast,mapping}

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

# Determine which reads to use (host-removed if available, otherwise trimmed)
echo "Checking for host-removed reads..."
USE_HOST_REMOVED=true
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "${HOST_REMOVAL_DIR}/${sample}_R1_clean.fastq.gz" ] || [ ! -f "${HOST_REMOVAL_DIR}/${sample}_R2_clean.fastq.gz" ]; then
        echo "Host-removed reads not found for $sample, will use trimmed reads"
        USE_HOST_REMOVED=false
        break
    fi
done

if [ "$USE_HOST_REMOVED" = true ]; then
    INPUT_DIR="$HOST_REMOVAL_DIR"
    R1_SUFFIX="_R1_clean.fastq.gz"
    R2_SUFFIX="_R2_clean.fastq.gz"
    echo "Using host-removed reads from: $INPUT_DIR"
else
    INPUT_DIR="$QC_DIR"
    R1_SUFFIX="_R1_trimmed.fastq.gz"
    R2_SUFFIX="_R2_trimmed.fastq.gz"
    echo "Using trimmed reads from: $INPUT_DIR"
fi

# Step 1: Individual sample assembly with MEGAHIT (parallel processing)
echo "Running MEGAHIT assembly for individual samples (parallel)..."

# Function to run MEGAHIT for a single sample
run_megahit() {
    local sample=$1
    local input_dir=$2
    local r1_suffix=$3
    local r2_suffix=$4
    local results_dir=$5
    local threads=$6
    
    echo "Starting MEGAHIT assembly for sample: $sample"
    
    R1_file="${input_dir}/${sample}${r1_suffix}"
    R2_file="${input_dir}/${sample}${r2_suffix}"
    
    if [ ! -f "$R1_file" ] || [ ! -f "$R2_file" ]; then
        echo "Warning: Input files not found for $sample, skipping..."
        return
    fi
    
    # Run MEGAHIT with optimized settings
    megahit -1 "$R1_file" \
            -2 "$R2_file" \
            -o "${results_dir}/megahit/${sample}" \
            --num-cpu-threads "$threads" \
            --memory 0.9 \
            --min-contig-len 1000 \
            --k-list 21,29,39,59,79,99,119 \
            --verbose
    
    # Rename final contigs
    cp "${results_dir}/megahit/${sample}/final.contigs.fa" \
       "${results_dir}/megahit/${sample}_contigs.fa"
    
    echo "Completed MEGAHIT assembly for sample: $sample"
}

# Export function for parallel execution
export -f run_megahit

# Calculate threads per job for parallel execution
PARALLEL_JOBS=$((THREADS > 8 ? THREADS / 4 : 2))  # 4 cores per job minimum
THREADS_PER_JOB=$((THREADS / PARALLEL_JOBS))

echo "Parallel execution settings:"
echo "  Parallel jobs: $PARALLEL_JOBS"
echo "  Threads per job: $THREADS_PER_JOB"

# Run assemblies in parallel using GNU parallel or xargs
if command -v parallel >/dev/null 2>&1; then
    echo "Using GNU parallel for faster processing..."
    printf '%s\n' "${SAMPLES[@]}" | \
    parallel -j $PARALLEL_JOBS run_megahit {} "$INPUT_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$RESULTS_DIR" "$THREADS_PER_JOB"
else
    echo "GNU parallel not found, using sequential processing..."
    for sample in "${SAMPLES[@]}"; do
        run_megahit "$sample" "$INPUT_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$RESULTS_DIR" "$THREADS"
    done
fi

# Step 2: Co-assembly with MEGAHIT (optimized settings)
echo "Running co-assembly with MEGAHIT (using all available cores)..."
R1_FILES=""
R2_FILES=""
for sample in "${SAMPLES[@]}"; do
    R1_file="${INPUT_DIR}/${sample}${R1_SUFFIX}"
    R2_file="${INPUT_DIR}/${sample}${R2_SUFFIX}"
    
    if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then
        R1_FILES="${R1_FILES}${R1_file},"
        R2_FILES="${R2_FILES}${R2_file},"
    fi
done
# Remove trailing comma
R1_FILES=${R1_FILES%,}
R2_FILES=${R2_FILES%,}

# Use all cores for co-assembly (most intensive step)
megahit -1 ${R1_FILES} \
        -2 ${R2_FILES} \
        -o ${RESULTS_DIR}/megahit/coassembly \
        --num-cpu-threads ${AVAILABLE_CORES} \
        --memory 0.9 \
        --min-contig-len 1000 \
        --k-list 21,29,39,59,79,99,119,141 \
        --verbose

cp ${RESULTS_DIR}/megahit/coassembly/final.contigs.fa \
   ${RESULTS_DIR}/megahit/coassembly_contigs.fa

# Step 3: Alternative assembly with metaSPAdes (optional, memory-intensive)
echo "Running metaSPAdes assembly for comparison (optional - can be skipped for speed)..."
read -p "Run metaSPAdes assemblies? This is slower but higher quality (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    # Function for parallel SPAdes execution
    run_spades() {
        local sample=$1
        local input_dir=$2
        local r1_suffix=$3
        local r2_suffix=$4
        local results_dir=$5
        local threads=$6
        local memory=$7
        
        echo "Starting metaSPAdes assembly for sample: $sample"
        
        R1_file="${input_dir}/${sample}${r1_suffix}"
        R2_file="${input_dir}/${sample}${r2_suffix}"
        
        if [ ! -f "$R1_file" ] || [ ! -f "$R2_file" ]; then
            echo "Warning: Input files not found for $sample, skipping..."
            return
        fi
        
        spades.py --meta \
                  -1 "$R1_file" \
                  -2 "$R2_file" \
                  -o "${results_dir}/spades/${sample}" \
                  --threads "$threads" \
                  --memory "$memory" \
                  --only-assembler
        
        # Copy scaffolds (contigs)
        cp "${results_dir}/spades/${sample}/scaffolds.fasta" \
           "${results_dir}/spades/${sample}_scaffolds.fa"
        
        echo "Completed metaSPAdes assembly for sample: $sample"
    }
    
    export -f run_spades
    
    # Run SPAdes with fewer parallel jobs due to memory requirements
    SPADES_JOBS=$((MEMORY > 64 ? 2 : 1))  # Limit based on memory
    MEMORY_PER_JOB=$((MEMORY / SPADES_JOBS))
    
    echo "SPAdes execution settings:"
    echo "  Parallel SPAdes jobs: $SPADES_JOBS"
    echo "  Memory per job: ${MEMORY_PER_JOB}GB"
    
    if command -v parallel >/dev/null 2>&1 && [ $SPADES_JOBS -gt 1 ]; then
        printf '%s\n' "${SAMPLES[@]}" | \
        parallel -j $SPADES_JOBS run_spades {} "$INPUT_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$RESULTS_DIR" "$THREADS_PER_JOB" "$MEMORY_PER_JOB"
    else
        for sample in "${SAMPLES[@]}"; do
            run_spades "$sample" "$INPUT_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$RESULTS_DIR" "$THREADS" "$MEMORY"
        done
    fi
else
    echo "Skipping metaSPAdes assemblies for faster processing..."
fi

# Step 4: Assembly quality assessment with QUAST
echo "Running QUAST for assembly quality assessment..."
# For MEGAHIT assemblies
quast ${RESULTS_DIR}/megahit/*_contigs.fa \
      -o ${RESULTS_DIR}/quast/megahit \
      --threads ${THREADS} \
      --min-contig 1000

# For metaSPAdes assemblies
quast ${RESULTS_DIR}/spades/*_scaffolds.fa \
      -o ${RESULTS_DIR}/quast/spades \
      --threads ${THREADS} \
      --min-contig 1000

# Step 5: Read mapping back to assemblies for coverage calculation (parallel)
echo "Mapping reads back to assemblies for coverage calculation (parallel)..."

# Function for parallel mapping
run_mapping() {
    local sample=$1
    local input_dir=$2
    local r1_suffix=$3
    local r2_suffix=$4
    local results_dir=$5
    local threads=$6
    
    echo "Starting mapping for sample: $sample"
    
    R1_file="${input_dir}/${sample}${r1_suffix}"
    R2_file="${input_dir}/${sample}${r2_suffix}"
    
    if [ ! -f "$R1_file" ] || [ ! -f "$R2_file" ]; then
        echo "Warning: Input files not found for $sample, skipping..."
        return
    fi
    
    # Check if assembly exists
    if [ ! -f "${results_dir}/megahit/${sample}_contigs.fa" ]; then
        echo "Warning: Assembly not found for $sample, skipping mapping..."
        return
    fi
    
    # Index assembly
    bwa index "${results_dir}/megahit/${sample}_contigs.fa"
    
    # Map reads with optimized settings
    bwa mem -t "$threads" \
            -M \
            "${results_dir}/megahit/${sample}_contigs.fa" \
            "$R1_file" "$R2_file" | \
    samtools view -@ $((threads/2)) -bS -F 4 - | \
    samtools sort -@ $((threads/2)) -o "${results_dir}/mapping/${sample}_sorted.bam"
    
    # Index BAM
    samtools index "${results_dir}/mapping/${sample}_sorted.bam"
    
    echo "Completed mapping for sample: $sample"
}

export -f run_mapping

# Run mapping in parallel
if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${SAMPLES[@]}" | \
    parallel -j $PARALLEL_JOBS run_mapping {} "$INPUT_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$RESULTS_DIR" "$THREADS_PER_JOB"
else
    for sample in "${SAMPLES[@]}"; do
        run_mapping "$sample" "$INPUT_DIR" "$R1_SUFFIX" "$R2_SUFFIX" "$RESULTS_DIR" "$THREADS"
    done
fi

echo "Assembly completed!"
echo "Results saved in: ${RESULTS_DIR}"
echo "Check QUAST reports for assembly statistics"

# Display final summary
echo ""
echo "=========================================="
echo "ASSEMBLY SUMMARY"
echo "=========================================="
echo "System resources used:"
echo "  Total cores: $AVAILABLE_CORES"
echo "  Threads used: $THREADS"
echo "  Memory used: ${MEMORY}GB"
echo "  Parallel jobs: $PARALLEL_JOBS"
echo ""
echo "Output files:"
echo "  MEGAHIT assemblies: ${RESULTS_DIR}/megahit/"
echo "  Co-assembly: ${RESULTS_DIR}/megahit/coassembly_contigs.fa"
if [ -d "${RESULTS_DIR}/spades" ]; then
    echo "  SPAdes assemblies: ${RESULTS_DIR}/spades/"
fi
echo "  QUAST reports: ${RESULTS_DIR}/quast/"
echo "  Read mappings: ${RESULTS_DIR}/mapping/"
echo ""
echo "Next step: Run binning"
echo "  ./code/04_binning.sh"
