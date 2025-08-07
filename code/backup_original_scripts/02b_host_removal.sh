#!/bin/bash

# Complete host removal script for tree microbiome samples
# Handles tree genome setup and host removal using Bowtie2
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QC_DIR="${PROJECT_DIR}/results/01_qc/trimmed"
RESULTS_DIR="${PROJECT_DIR}/results/02_host_removal"
GENOMES_DIR="${RESULTS_DIR}/tree_genomes"
THREADS=8

echo "============================================"
echo "Complete Host Removal for Tree Microbiomes"
echo "============================================"
echo "Tree species: Shorea balangeran, Xylopia fusca, Alnus glutinosa, Betula pubescens"
echo ""

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
if conda activate metagenome_assembly_new 2>/dev/null; then
    echo "✓ Using metagenome_assembly_new environment"
elif conda activate metagenome_assembly 2>/dev/null; then
    echo "✓ Using metagenome_assembly environment"
else
    echo "✗ No suitable environment found. Please run setup script first."
    exit 1
fi

# Check and install required tools
echo "Checking required tools..."
MISSING_TOOLS=()

if ! command -v bowtie2-build >/dev/null 2>&1; then
    MISSING_TOOLS+=("bowtie2")
fi

if ! command -v samtools >/dev/null 2>&1; then
    MISSING_TOOLS+=("samtools")
fi

if ! command -v bc >/dev/null 2>&1; then
    MISSING_TOOLS+=("bc")
fi

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo "Installing missing tools: ${MISSING_TOOLS[*]}"
    if [[ " ${MISSING_TOOLS[@]} " =~ " bowtie2 " ]] || [[ " ${MISSING_TOOLS[@]} " =~ " samtools " ]]; then
        conda install -c bioconda bowtie2 samtools bc -y
    else
        conda install bc -y
    fi
    echo "✓ Tools installed successfully"
else
    echo "✓ All required tools are available"
fi

# Create output directories
mkdir -p ${RESULTS_DIR}/{clean_reads,stats,mapping}
mkdir -p ${GENOMES_DIR}

# Sample classification
TREE_SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10")
SOIL_SAMPLES=("53398_A16S" "53399_B10S")
ALL_SAMPLES=("${TREE_SAMPLES[@]}" "${SOIL_SAMPLES[@]}")

echo "Sample classification:"
echo "  Tree samples (host removal): ${TREE_SAMPLES[@]}"
echo "  Soil samples (copy only): ${SOIL_SAMPLES[@]}"
echo ""

# STEP 1: Download and setup tree reference genomes
echo "STEP 1: Setting up tree reference genomes..."
cd ${GENOMES_DIR}

# Check if we need to download genomes
GENOMES_EXIST=false
if [ -f "combined_tree_genomes.fna" ] && [ -s "combined_tree_genomes.fna" ]; then
    GENOMES_EXIST=true
    echo "✓ Combined tree genomes already exist"
fi

if [ "$GENOMES_EXIST" = false ]; then
    echo "Downloading tree reference genomes..."
    
    # Install ncbi-datasets-cli if needed
    if ! command -v datasets >/dev/null 2>&1; then
        echo "Installing NCBI datasets CLI..."
        conda install -c conda-forge -c bioconda ncbi-datasets-cli -y
    fi
    
    # 1. Betula pubescens (Downy birch) - TaxID: 3506
    echo "Downloading Betula pubescens genome..."
    if [ ! -f "betula_pubescens.fna" ]; then
        datasets download genome taxon 3506 --reference --filename betula.zip 2>/dev/null && {
            unzip -q betula.zip
            find ncbi_dataset/data/ -name "*.fna" -exec cp {} betula_pubescens.fna \; 2>/dev/null
            rm -rf ncbi_dataset betula.zip
            echo "  ✓ Betula pubescens downloaded"
        } || {
            echo "  Datasets failed, trying direct download..."
            wget -q -O betula_pubescens.fna.gz \
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/067/625/GCA_905067625.1_Bpub_v1/GCA_905067625.1_Bpub_v1_genomic.fna.gz" && \
            gunzip betula_pubescens.fna.gz && echo "  ✓ Betula pubescens downloaded" || echo "  ✗ Betula download failed"
        }
    fi
    
    # 2. Alnus glutinosa (European alder) - TaxID: 3517
    echo "Downloading Alnus glutinosa genome..."
    if [ ! -f "alnus_glutinosa.fna" ]; then
        datasets download genome taxon 3517 --reference --filename alnus.zip 2>/dev/null && {
            unzip -q alnus.zip
            find ncbi_dataset/data/ -name "*.fna" -exec cp {} alnus_glutinosa.fna \; 2>/dev/null
            rm -rf ncbi_dataset alnus.zip
            echo "  ✓ Alnus glutinosa downloaded"
        } || {
            echo "  Datasets failed, trying direct download..."
            wget -q -O alnus_glutinosa.fna.gz \
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/418/595/GCA_030418595.1_Agl_v1.0/GCA_030418595.1_Agl_v1.0_genomic.fna.gz" && \
            gunzip alnus_glutinosa.fna.gz && echo "  ✓ Alnus glutinosa downloaded" || echo "  ✗ Alnus download failed"
        }
    fi
    
    # 3. Related genome for Shorea balangeran (Dipterocarpaceae family)
    echo "Downloading related genome for Shorea balangeran..."
    if [ ! -f "dipterocarp_related.fna" ]; then
        wget -q -O dipterocarp_related.fna.gz \
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/519/115/GCA_900519115.1_Dryobalanops_aromatica_v1.0/GCA_900519115.1_Dryobalanops_aromatica_v1.0_genomic.fna.gz" && \
        gunzip dipterocarp_related.fna.gz && echo "  ✓ Dipterocarp relative downloaded" || echo "  ✗ Dipterocarp download failed"
    fi
    
    # 4. Related genome for Xylopia fusca (Annonaceae family)
    echo "Downloading related genome for Xylopia fusca..."
    if [ ! -f "annonaceae_related.fna" ]; then
        datasets download genome taxon 49314 --reference --filename annona.zip 2>/dev/null && {
            unzip -q annona.zip
            find ncbi_dataset/data/ -name "*.fna" -exec cp {} annonaceae_related.fna \; 2>/dev/null
            rm -rf ncbi_dataset annona.zip
            echo "  ✓ Annonaceae relative downloaded"
        } || echo "  ✗ No Annonaceae genome available"
    fi
    
    # Check downloaded genomes and combine
    echo ""
    echo "Checking downloaded genomes:"
    total_size=0
    for file in *.fna; do
        if [ -f "$file" ] && [ -s "$file" ]; then
            size=$(stat --format=%s "$file")
            size_mb=$((size / 1024 / 1024))
            contigs=$(grep -c "^>" "$file" 2>/dev/null || echo "0")
            echo "  $file: ${size_mb}MB, $contigs contigs"
            total_size=$((total_size + size))
        fi
    done
    
    # Combine all available genomes
    if ls *.fna 1> /dev/null 2>&1; then
        echo ""
        echo "Combining tree genomes..."
        cat *.fna > combined_tree_genomes.fna
        
        combined_size=$(($(stat --format=%s combined_tree_genomes.fna) / 1024 / 1024))
        combined_contigs=$(grep -c "^>" combined_tree_genomes.fna)
        echo "  ✓ Combined genome: ${combined_size}MB, $combined_contigs total contigs"
    else
        echo "  ✗ No tree genomes successfully downloaded!"
        echo "  Will skip host removal and copy trimmed reads directly"
        SKIP_HOST_REMOVAL=true
    fi
fi

# STEP 2: Build Bowtie2 index
TREE_INDEX=""
if [ -f "combined_tree_genomes.fna" ] && [ -s "combined_tree_genomes.fna" ]; then
    if [ ! -f "tree_index.1.bt2" ]; then
        echo ""
        echo "STEP 2: Building Bowtie2 index..."
        if command -v bowtie2-build >/dev/null 2>&1; then
            echo "  Building index with $(nproc) threads..."
            bowtie2-build --threads ${THREADS} combined_tree_genomes.fna tree_index 2>/dev/null && {
                echo "  ✓ Bowtie2 index built successfully"
            } || {
                echo "  ✗ Bowtie2 index building failed, will skip host removal"
                SKIP_HOST_REMOVAL=true
            }
        else
            echo "  ✗ bowtie2-build not found! Installing..."
            conda install -c bioconda bowtie2 -y
            if command -v bowtie2-build >/dev/null 2>&1; then
                bowtie2-build --threads ${THREADS} combined_tree_genomes.fna tree_index 2>/dev/null && {
                    echo "  ✓ Bowtie2 index built successfully"
                } || {
                    echo "  ✗ Bowtie2 index building failed"
                    SKIP_HOST_REMOVAL=true
                }
            else
                echo "  ✗ Failed to install bowtie2-build"
                SKIP_HOST_REMOVAL=true
            fi
        fi
    else
        echo "STEP 2: ✓ Bowtie2 index already exists"
    fi
    
    if [ -f "tree_index.1.bt2" ]; then
        TREE_INDEX="${GENOMES_DIR}/tree_index"
        echo "  Using index: $TREE_INDEX"
        # Verify all index files are present
        index_files=(tree_index.1.bt2 tree_index.2.bt2 tree_index.3.bt2 tree_index.4.bt2 tree_index.rev.1.bt2 tree_index.rev.2.bt2)
        for idx_file in "${index_files[@]}"; do
            if [ ! -f "$idx_file" ]; then
                echo "  ✗ Index file $idx_file missing, rebuilding..."
                rm -f tree_index.*.bt2
                bowtie2-build --threads ${THREADS} combined_tree_genomes.fna tree_index
                break
            fi
        done
    fi
else
    echo "STEP 2: ✗ No combined genomes found, skipping index building"
    SKIP_HOST_REMOVAL=true
fi

echo ""
echo ""

# STEP 3: Process samples with host removal
echo "STEP 3: Processing samples..."
cd ${PROJECT_DIR}

for sample in "${ALL_SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "Processing sample: $sample"
    
    R1_trimmed="${QC_DIR}/${sample}_R1_trimmed.fastq.gz"
    R2_trimmed="${QC_DIR}/${sample}_R2_trimmed.fastq.gz"
    
    if [ ! -f "$R1_trimmed" ] || [ ! -f "$R2_trimmed" ]; then
        echo "  ✗ Trimmed reads not found for $sample, skipping..."
        continue
    fi
    
    # Count original reads
    echo "  Counting original reads..."
    original_reads=$(( $(zcat $R1_trimmed | wc -l) / 4 ))
    echo "  Original reads: $original_reads"
    
    # Check if this is a soil sample or if we should skip host removal
    if [[ " ${SOIL_SAMPLES[@]} " =~ " ${sample} " ]]; then
        echo "  Sample type: SOIL - copying trimmed reads (no host removal)"
        cp $R1_trimmed ${RESULTS_DIR}/clean_reads/${sample}_R1_clean.fastq.gz
        cp $R2_trimmed ${RESULTS_DIR}/clean_reads/${sample}_R2_clean.fastq.gz
        
        # Stats for soil samples
        clean_reads=$original_reads
        host_reads=0
        host_percent="0.00"
        clean_percent="100.00"
        
    elif [ -z "$TREE_INDEX" ] || [ "$SKIP_HOST_REMOVAL" = true ]; then
        echo "  Sample type: TREE - but no tree index available, copying trimmed reads"
        cp $R1_trimmed ${RESULTS_DIR}/clean_reads/${sample}_R1_clean.fastq.gz
        cp $R2_trimmed ${RESULTS_DIR}/clean_reads/${sample}_R2_clean.fastq.gz
        
        # Stats for copied tree samples
        clean_reads=$original_reads
        host_reads=0
        host_percent="0.00"
        clean_percent="100.00"
        
    else
        echo "  Sample type: TREE - performing host removal with custom tree genomes"
        
        # Map to tree genomes with Bowtie2
        echo "  Mapping reads to tree genomes..."
        bowtie2 -x $TREE_INDEX \
                -1 $R1_trimmed \
                -2 $R2_trimmed \
                --threads ${THREADS} \
                --very-fast \
                --no-mixed \
                --no-discordant \
                -S ${RESULTS_DIR}/mapping/${sample}_host_mapped.sam 2>/dev/null
        
        echo "  Extracting non-host reads..."
        # Extract properly paired unmapped reads (non-host)
        # Both reads must be unmapped (flags 77 and 141)
        samtools view -f 77 ${RESULTS_DIR}/mapping/${sample}_host_mapped.sam | \
        awk '{print "@"$1"/1\n"$10"\n+\n"$11}' > ${RESULTS_DIR}/clean_reads/${sample}_R1_clean.fastq
        
        samtools view -f 141 ${RESULTS_DIR}/mapping/${sample}_host_mapped.sam | \
        awk '{print "@"$1"/2\n"$10"\n+\n"$11}' > ${RESULTS_DIR}/clean_reads/${sample}_R2_clean.fastq
        
        # Compress clean reads
        gzip ${RESULTS_DIR}/clean_reads/${sample}_R1_clean.fastq
        gzip ${RESULTS_DIR}/clean_reads/${sample}_R2_clean.fastq
        
        # Clean up SAM file to save space
        rm ${RESULTS_DIR}/mapping/${sample}_host_mapped.sam
        
        # Count clean reads
        clean_reads=$(( $(zcat ${RESULTS_DIR}/clean_reads/${sample}_R1_clean.fastq.gz | wc -l) / 4 ))
        host_reads=$((original_reads - clean_reads))
        
        # Calculate percentages using bc for precision
        if [ $original_reads -gt 0 ]; then
            host_percent=$(echo "scale=2; $host_reads * 100 / $original_reads" | bc -l)
            clean_percent=$(echo "scale=2; $clean_reads * 100 / $original_reads" | bc -l)
        else
            host_percent="0.00"
            clean_percent="0.00"
        fi
    fi
    
    # Write statistics
    echo -e "Sample\tOriginal_Reads\tHost_Reads\tClean_Reads\tHost_Percent\tClean_Percent" > ${RESULTS_DIR}/stats/${sample}_host_removal_stats.txt
    echo -e "$sample\t$original_reads\t$host_reads\t$clean_reads\t$host_percent\t$clean_percent" >> ${RESULTS_DIR}/stats/${sample}_host_removal_stats.txt
    
    echo "  Result: $original_reads → $clean_reads reads (${clean_percent}% retained)"
done

# STEP 4: Generate summary statistics
echo ""
echo "STEP 4: Generating summary statistics..."
echo -e "Sample\tOriginal_Reads\tHost_Reads\tClean_Reads\tHost_Percent\tClean_Percent" > ${RESULTS_DIR}/host_removal_summary.txt

for sample in "${ALL_SAMPLES[@]}"; do
    if [ -f "${RESULTS_DIR}/stats/${sample}_host_removal_stats.txt" ]; then
        tail -n +2 ${RESULTS_DIR}/stats/${sample}_host_removal_stats.txt >> ${RESULTS_DIR}/host_removal_summary.txt
    fi
done

echo ""
echo "============================================"
echo "HOST REMOVAL COMPLETED SUCCESSFULLY!"
echo "============================================"
echo ""
echo "Summary:"
echo "  Tree samples processed: ${TREE_SAMPLES[@]}"
echo "  Soil samples processed: ${SOIL_SAMPLES[@]}"
if [ ! -z "$TREE_INDEX" ] && [ "$SKIP_HOST_REMOVAL" != true ]; then
    echo "  Used custom tree genome index: $TREE_INDEX"
else
    echo "  No host removal performed (copied trimmed reads)"
fi
echo ""
echo "Output files:"
echo "  Clean reads: ${RESULTS_DIR}/clean_reads/"
echo "  Statistics: ${RESULTS_DIR}/stats/"
echo "  Summary: ${RESULTS_DIR}/host_removal_summary.txt"
echo ""
echo "Next step: Run assembly using clean reads"
echo "  ./code/03_assembly.sh"
echo ""

# Display summary table
if [ -f "${RESULTS_DIR}/host_removal_summary.txt" ]; then
    echo "SUMMARY TABLE:"
    column -t -s $'\t' ${RESULTS_DIR}/host_removal_summary.txt
fi

