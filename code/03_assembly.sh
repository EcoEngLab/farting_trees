#!/bin/bash

# Fast Metagenomic Assembly - Optimized for Speed
# Author: Generated for farting_trees project
# Date: August 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QC_DIR="${PROJECT_DIR}/results/01_qc/trimmed"
HOST_REMOVAL_DIR="${PROJECT_DIR}/results/02_host_removal/clean_reads"
RESULTS_DIR="${PROJECT_DIR}/results/03_assembly"

# Aggressive CPU optimization - use all available cores
AVAILABLE_CORES=$(nproc)
THREADS=$AVAILABLE_CORES  # Use ALL cores for speed
MEMORY=$(($(free -g | awk '/^Mem:/{print $2}') * 90 / 100))  # Use 90% of available RAM

echo "=========================================="
echo "FAST METAGENOMIC ASSEMBLY"
echo "=========================================="
echo "System resources:"
echo "  Cores: $AVAILABLE_CORES (using ALL)"
echo "  Memory: ${MEMORY}GB (90% allocation)"
echo "  Strategy: Co-assembly ONLY for maximum speed"
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

# Quick tool check
echo "Checking tools..."
for tool in megahit quast bwa samtools; do
    if command -v $tool >/dev/null 2>&1; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool missing - installing..."
        conda install -c bioconda $tool -y
    fi
done

# Create output directories
mkdir -p ${RESULTS_DIR}/{megahit,quast,mapping}

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

# Check for host-removed reads (fast check)
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

echo ""
echo "STRATEGY: Fast co-assembly only (skip individual assemblies for speed)"
echo "This gives you the best assembly in minimum time!"
echo ""

# FAST STEP 1: Co-assembly ONLY with MEGAHIT (maximum speed)
echo "=========================================="
echo "RUNNING FAST CO-ASSEMBLY WITH MEGAHIT"
echo "=========================================="

# Prepare combined input files
echo "Preparing input files..."
R1_FILES=""
R2_FILES=""
total_reads=0

for sample in "${SAMPLES[@]}"; do
    R1_file="${INPUT_DIR}/${sample}${R1_SUFFIX}"
    R2_file="${INPUT_DIR}/${sample}${R2_SUFFIX}"
    
    if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then
        R1_FILES="${R1_FILES}${R1_file},"
        R2_FILES="${R2_FILES}${R2_file},"
        reads=$(( $(zcat "$R1_file" | wc -l) / 4 ))
        total_reads=$((total_reads + reads))
        echo "  ✓ $sample: $reads reads"
    else
        echo "  ✗ $sample: files not found"
    fi
done

# Remove trailing commas
R1_FILES=${R1_FILES%,}
R2_FILES=${R2_FILES%,}

echo "Total input reads: $total_reads"
echo ""

# Run MEGAHIT co-assembly with ALL cores and optimized settings
echo "Starting MEGAHIT co-assembly..."
echo "Using $THREADS cores and ${MEMORY}GB memory"

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

end_time=$(date +%s)
assembly_time=$((end_time - start_time))

# Copy and rename final assembly
cp ${RESULTS_DIR}/megahit/coassembly/final.contigs.fa \
   ${RESULTS_DIR}/megahit/coassembly_contigs.fa

echo "✓ Co-assembly completed in ${assembly_time} seconds"

# FAST STEP 2: Quick quality assessment
echo ""
echo "=========================================="
echo "QUICK ASSEMBLY QUALITY CHECK"
echo "=========================================="

# Basic assembly stats (much faster than full QUAST)
assembly_file="${RESULTS_DIR}/megahit/coassembly_contigs.fa"

if [ -f "$assembly_file" ]; then
    echo "Co-assembly statistics:"
    
    # Count contigs
    contigs=$(grep -c "^>" "$assembly_file")
    
    # Total length
    total_length=$(grep -v "^>" "$assembly_file" | tr -d '\n' | wc -c)
    total_mb=$((total_length / 1000000))
    
    # Largest contig
    largest=$(grep -v "^>" "$assembly_file" | awk '{print length}' | sort -nr | head -1)
    largest_kb=$((largest / 1000))
    
    # N50 calculation (simplified)
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

# OPTIONAL: Run full QUAST only if requested
echo ""
read -p "Run detailed QUAST analysis? (adds ~1-2 minutes) (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Running detailed QUAST analysis..."
    quast ${RESULTS_DIR}/megahit/coassembly_contigs.fa \
          -o ${RESULTS_DIR}/quast/coassembly \
          --threads ${THREADS} \
          --min-contig 1000 \
          --fast
    echo "✓ QUAST completed"
else
    echo "Skipping detailed QUAST for speed"
fi

# OPTIONAL STEP 3: Coverage mapping (enhanced with working implementation)
echo ""
read -p "Calculate coverage by mapping reads back? (adds ~10-15 minutes) (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "=========================================="
    echo "CALCULATING COVERAGE"
    echo "=========================================="
    
    assembly_file="${RESULTS_DIR}/megahit/coassembly_contigs.fa"
    
    # Verify assembly file exists and is valid
    if [ ! -f "$assembly_file" ]; then
        echo "  ✗ Assembly file not found: $assembly_file"
        echo "  Cannot proceed with coverage calculation"
    else
        # Check assembly file size
        assembly_size=$(stat -c%s "$assembly_file" 2>/dev/null || echo 0)
        if [ $assembly_size -lt 1000 ]; then
            echo "  ✗ Assembly file too small (${assembly_size} bytes)"
            echo "  Cannot proceed with coverage calculation"
        else
            echo "  Assembly file verified: $(( assembly_size / 1024 / 1024 )) MB"
            
            # Index assembly (BWA indexes may already exist from assembly step)
            if [ ! -f "${assembly_file}.bwt" ]; then
                echo "Creating BWA indexes..."
                if bwa index "$assembly_file" 2>/dev/null; then
                    echo "  ✓ BWA index created successfully"
                else
                    echo "  ✗ BWA indexing failed"
                    echo "  Cannot proceed with coverage calculation"
                fi
            else
                echo "  ✓ BWA indexes already exist"
            fi
            
            # Only proceed if BWA indexes exist
            if [ -f "${assembly_file}.bwt" ]; then
                echo "Processing coverage calculation..."
                
                # Create temporary directory for individual SAM files
                temp_dir="${RESULTS_DIR}/mapping/temp_sams"
                mkdir -p "$temp_dir"
                
                # Process each sample pair individually (fixes the BWA command issue)
                echo "Mapping reads from ${#SAMPLES[@]} samples individually..."
                sam_files=()
                
                for sample in "${SAMPLES[@]}"; do
                    R1_file="${INPUT_DIR}/${sample}${R1_SUFFIX}"
                    R2_file="${INPUT_DIR}/${sample}${R2_SUFFIX}"
                    
                    if [ -f "$R1_file" ] && [ -f "$R2_file" ]; then
                        individual_sam="${temp_dir}/${sample}.sam"
                        
                        echo "  Mapping sample: $sample"
                        if bwa mem -t ${THREADS} -M "$assembly_file" "$R1_file" "$R2_file" > "$individual_sam" 2>"${temp_dir}/${sample}_error.log"; then
                            echo "    ✓ Alignment completed for $sample"
                            sam_files+=("$individual_sam")
                        else
                            echo "    ✗ Alignment failed for $sample"
                            echo "    Check error log: ${temp_dir}/${sample}_error.log"
                        fi
                    else
                        echo "    ✗ Input files not found for $sample"
                    fi
                done
                
                # Merge SAM files if we have successful alignments
                if [ ${#sam_files[@]} -gt 0 ]; then
                    echo "Merging ${#sam_files[@]} successful alignments..."
                    temp_sam="${RESULTS_DIR}/mapping/temp_alignment.sam"
                    
                    # Take header from first file
                    head -n 100 "${sam_files[0]}" | grep "^@" > "$temp_sam"
                    
                    # Add all alignments
                    for sam_file in "${sam_files[@]}"; do
                        grep -v "^@" "$sam_file" >> "$temp_sam"
                    done
                    
                    # Convert to BAM and sort
                    echo "Converting to sorted BAM..."
                    if samtools view -@ ${THREADS} -bS -F 4 "$temp_sam" 2>/dev/null | \
                       samtools sort -@ ${THREADS} -o "${RESULTS_DIR}/mapping/coassembly_sorted.bam" 2>/dev/null; then
                        echo "  ✓ BAM file created successfully"
                        
                        # Index BAM file
                        echo "Indexing BAM file..."
                        if samtools index "${RESULTS_DIR}/mapping/coassembly_sorted.bam" 2>/dev/null; then
                            echo "  ✓ BAM indexing completed"
                            
                            # Calculate basic coverage statistics
                            bam_file="${RESULTS_DIR}/mapping/coassembly_sorted.bam"
                            total_mapped=$(samtools view -c -F 4 "$bam_file")
                            
                            echo "Coverage summary:"
                            echo "  Total mapped alignments: $total_mapped"
                            echo "  BAM file: ${RESULTS_DIR}/mapping/coassembly_sorted.bam"
                            echo "  ✓ Coverage calculation completed successfully"
                            
                        else
                            echo "  ✗ BAM indexing failed"
                        fi
                    else
                        echo "  ✗ BAM conversion failed"
                    fi
                    
                    # Clean up temporary files
                    rm -f "$temp_sam"
                    rm -rf "$temp_dir"
                    
                else
                    echo "  ✗ No successful alignments to process"
                fi
            fi
        fi
    fi
else
    echo "Skipping coverage calculation for maximum speed"
fi

# Final summary
total_time=$(date +%s)
total_runtime=$((total_time - start_time))

echo ""
echo "=========================================="
echo "FAST ASSEMBLY COMPLETED!"
echo "=========================================="
echo "Total runtime: ${total_runtime} seconds"
echo ""
echo "Key output:"
echo "  Main assembly: ${RESULTS_DIR}/megahit/coassembly_contigs.fa"
echo "  Basic stats: displayed above"
if [ -d "${RESULTS_DIR}/quast/coassembly" ]; then
    echo " QUAST report: ${RESULTS_DIR}/quast/coassembly/report.html"
fi
if [ -f "${RESULTS_DIR}/mapping/coassembly_sorted.bam" ]; then
    echo " Coverage: ${RESULTS_DIR}/mapping/coassembly_sorted.bam"
fi
echo ""
echo "SPEED OPTIMIZATIONS USED:"
echo "  ✓ Co-assembly only (skipped individual assemblies)"
echo "  ✓ Used all $THREADS CPU cores"
echo "  ✓ Optimized k-mer list for speed"
echo "  ✓ Disabled expensive mercy kmers"
echo "  ✓ Optional detailed analysis steps"
echo ""
echo "Next steps:"
echo " Functional gene screening: ./code/04_functional_gene_screening.sh"
echo " Genome binning: ./code/04_binning.sh"
echo ""
