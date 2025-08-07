#!/bin/bash

# Metagenomic binning using MetaBAT2, MaxBin2, and VAMB on co-assembly
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/03_assembly"
RESULTS_DIR="${WORK_DIR}/results/04_binning"
THREADS=8
COASSEMBLY_FASTA="${ASSEMBLY_DIR}/megahit/coassembly_contigs.fa"

# Activate required environments
source $(conda info --base)/etc/profile.d/conda.sh

# BWA and samtools are in assembly environment
conda activate metagenome_assembly

echo "Starting metagenomic binning with differential coverage analysis..."
echo "Co-assembly file: $COASSEMBLY_FASTA"
echo "Strategy: Map each sample individually to co-assembly for differential coverage"

# Create output directories
mkdir -p ${RESULTS_DIR}/{metabat2,maxbin2,vamb,das_tool,checkm}

# Step 1: Map each sample to co-assembly for differential coverage
echo "Mapping each sample to co-assembly for differential coverage analysis..."

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

# Check if BWA index exists
if [ ! -f "${COASSEMBLY_FASTA}.bwt" ]; then
    echo "Creating BWA index for co-assembly..."
    conda activate metagenome_assembly
    bwa index ${COASSEMBLY_FASTA}
fi

# Map each sample individually to co-assembly
for sample in "${SAMPLES[@]}"; do
    echo "Mapping sample: $sample to co-assembly"
    
    # Define input and output files
    R1="${WORK_DIR}/results/02_host_removal/clean_reads/${sample}_R1_clean.fastq.gz"
    R2="${WORK_DIR}/results/02_host_removal/clean_reads/${sample}_R2_clean.fastq.gz"
    output_bam="${RESULTS_DIR}/${sample}_coassembly_sorted.bam"
    
    # Check if files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "⚠ Warning: Missing files for $sample, skipping..."
        echo "  Expected: $R1"
        echo "  Expected: $R2"
        continue
    fi
    
    # Skip if BAM already exists and is not empty
    if [ -f "$output_bam" ] && [ -s "$output_bam" ]; then
        echo "✓ BAM file already exists for $sample: $output_bam"
        continue
    fi
    
    echo "  Aligning $sample to co-assembly..."
    conda activate metagenome_assembly
    bwa mem -t ${THREADS} \
            ${COASSEMBLY_FASTA} \
            $R1 $R2 | \
    samtools view -@ ${THREADS} -bS -F 4 - | \
    samtools sort -@ ${THREADS} -o $output_bam
    
    # Index the BAM file
    samtools index $output_bam
    
    echo "✓ Completed mapping for $sample"
done

# Step 2: Generate depth file using all individual BAM files
echo "Generating differential coverage depth file..."
BAM_FILES="${RESULTS_DIR}/*_coassembly_sorted.bam"

# Check if we have BAM files
if ! ls ${RESULTS_DIR}/*_coassembly_sorted.bam 1> /dev/null 2>&1; then
    echo "✗ No individual sample BAM files found!"
    echo "Expected pattern: ${RESULTS_DIR}/*_coassembly_sorted.bam"
    exit 1
fi

echo "Found BAM files:"
ls -la ${RESULTS_DIR}/*_coassembly_sorted.bam

# Switch to binning environment for jgi_summarize_bam_contig_depths
conda activate metagenome_binning
jgi_summarize_bam_contig_depths --outputDepth ${RESULTS_DIR}/coassembly_depth.txt \
                               ${RESULTS_DIR}/*_coassembly_sorted.bam

# Step 2: MetaBAT2 binning on co-assembly (RELAXED PARAMETERS FOR METHANE GENES)
echo "Running MetaBAT2 binning on co-assembly with relaxed parameters..."
echo "Using more inclusive settings to capture methane-containing contigs..."

# Standard binning (original)
echo "  Running standard MetaBAT2 binning..."
metabat2 -i ${COASSEMBLY_FASTA} \
         -a ${RESULTS_DIR}/coassembly_depth.txt \
         -o ${RESULTS_DIR}/metabat2/coassembly_bin \
         --minContig 2500 \
         --numThreads ${THREADS}

# Relaxed binning for methane genes (lower thresholds)
echo "  Running relaxed MetaBAT2 binning for methane genes..."
metabat2 -i ${COASSEMBLY_FASTA} \
         -a ${RESULTS_DIR}/coassembly_depth.txt \
         -o ${RESULTS_DIR}/metabat2/relaxed_bin \
         --minContig 1500 \
         --minClsSize 200000 \
         --maxP 95 \
         --minS 60 \
         --maxEdges 200 \
         --pTNF 0 \
         --noAdd \
         --numThreads ${THREADS}

# Very inclusive binning (catch remaining contigs)
echo "  Running very inclusive MetaBAT2 binning..."
metabat2 -i ${COASSEMBLY_FASTA} \
         -a ${RESULTS_DIR}/coassembly_depth.txt \
         -o ${RESULTS_DIR}/metabat2/inclusive_bin \
         --minContig 1500 \
         --minClsSize 100000 \
         --maxP 90 \
         --minS 40 \
         --maxEdges 500 \
         --pTNF 0 \
         --noAdd \
         --numThreads ${THREADS}

# Step 3: MaxBin2 binning with differential coverage (RELAXED PARAMETERS)
echo "Running MaxBin2 binning with relaxed parameters for better coverage..."

# Check if MaxBin2 is available
if ! command -v run_MaxBin.pl &> /dev/null; then
    echo "⚠ Warning: MaxBin2 (run_MaxBin.pl) not found, skipping MaxBin2 binning"
    echo "You can install MaxBin2 later with: conda install -c bioconda maxbin2"
    mkdir -p ${RESULTS_DIR}/maxbin2
    echo "Skipping MaxBin2 due to missing installation" > ${RESULTS_DIR}/maxbin2/SKIPPED.txt
else
    # MaxBin2 prefers multiple abundance files for better binning
    # Create abundance files from depth data (columns 4, 6, 8, 10, 12, 14 are sample depths)
    echo "Creating abundance files for MaxBin2..."

    # Extract sample-specific abundance columns
    # Assuming depth file format: contig | length | totalAvgDepth | sample1.var | sample1 | sample2.var | sample2 | ...
    ABUND_FILES=""
    SAMPLE_COUNT=0

    # Count samples and create abundance files
    for i in $(seq 4 2 20); do  # Check columns 4, 6, 8, 10, 12, 14 (up to 6 samples)
        col_data=$(cut -f${i} ${RESULTS_DIR}/coassembly_depth.txt | head -20 | tail -10)
        if [[ -n "$col_data" ]] && [[ "$col_data" != *"0.0000"* ]]; then
            SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
            cut -f1,${i} ${RESULTS_DIR}/coassembly_depth.txt > ${RESULTS_DIR}/coassembly_abundance_sample${SAMPLE_COUNT}.txt
            ABUND_FILES="${ABUND_FILES} -abund ${RESULTS_DIR}/coassembly_abundance_sample${SAMPLE_COUNT}.txt"
        fi
    done

    echo "Found $SAMPLE_COUNT samples for MaxBin2 differential coverage"

    # Standard MaxBin2 run
    echo "  Running standard MaxBin2..."
    run_MaxBin.pl -contig ${COASSEMBLY_FASTA} \
                  ${ABUND_FILES} \
                  -out ${RESULTS_DIR}/maxbin2/coassembly_bin \
                  -thread ${THREADS} \
                  -min_contig_length 2500

    # Relaxed MaxBin2 run for more contigs
    echo "  Running relaxed MaxBin2 for methane genes..."
    run_MaxBin.pl -contig ${COASSEMBLY_FASTA} \
                  ${ABUND_FILES} \
                  -out ${RESULTS_DIR}/maxbin2/relaxed_bin \
                  -thread ${THREADS} \
                  -min_contig_length 1000 \
                  -prob_threshold 0.7 \
                  -markerset 40
fi

# Step 4: VAMB binning with differential coverage (MULTIPLE RUNS)
echo "Running VAMB binning with multiple parameter sets..."

# Check if VAMB is available
if ! command -v vamb &> /dev/null; then
    echo "⚠ Warning: VAMB not found, skipping VAMB binning"
    echo "You can install VAMB later with: conda install -c bioconda vamb"
    mkdir -p ${RESULTS_DIR}/vamb
    echo "Skipping VAMB due to missing installation" > ${RESULTS_DIR}/vamb/SKIPPED.txt
else
    # Standard VAMB run
    echo "  Running standard VAMB..."
    vamb --outdir ${RESULTS_DIR}/vamb/standard \
         --fasta ${COASSEMBLY_FASTA} \
         --bamfiles ${RESULTS_DIR}/*_coassembly_sorted.bam \
         --mincontig 2500 \
         --minfasta 200000

    # Relaxed VAMB run for methane genes
    echo "  Running relaxed VAMB for methane genes..."
    vamb --outdir ${RESULTS_DIR}/vamb/relaxed \
         --fasta ${COASSEMBLY_FASTA} \
         --bamfiles ${RESULTS_DIR}/*_coassembly_sorted.bam \
         --mincontig 1000 \
         --minfasta 100000 \
         --alpha 0.1 \
         --beta 0.1

    # Very inclusive VAMB run 
    echo "  Running very inclusive VAMB..."
    vamb --outdir ${RESULTS_DIR}/vamb/inclusive \
         --fasta ${COASSEMBLY_FASTA} \
         --bamfiles ${RESULTS_DIR}/*_coassembly_sorted.bam \
         --mincontig 1000 \
         --minfasta 50000 \
         --alpha 0.05 \
         --beta 0.05
fi

# Step 4.5: TARGETED METHANE CONTIG BINNING
echo "Running targeted binning on methane-containing contigs..."

# Check if we have previously identified methane contigs
METHANE_CONTIGS_FILE="${WORK_DIR}/results/10_methane_contig_extraction/all_methane_contigs.fa"

if [ -f "$METHANE_CONTIGS_FILE" ]; then
    echo "Found previously identified methane contigs: $METHANE_CONTIGS_FILE"
    
    # Create targeted methane bin directory
    mkdir -p ${RESULTS_DIR}/methane_targeted
    
    # Run MetaBAT2 specifically on methane contigs with very relaxed parameters
    echo "  Running targeted MetaBAT2 on methane contigs..."
    metabat2 -i ${METHANE_CONTIGS_FILE} \
             -a ${RESULTS_DIR}/coassembly_depth.txt \
             -o ${RESULTS_DIR}/methane_targeted/methane_bin \
             --minContig 1500 \
             --minClsSize 50000 \
             --maxP 85 \
             --minS 30 \
             --maxEdges 1000 \
             --pTNF 0 \
             --noAdd \
             --numThreads ${THREADS}
    
    # Create separate bins for mcrA and pmoA contigs if available
    MCRA_CONTIGS="${WORK_DIR}/results/10_methane_contig_extraction/mcrA_containing_contigs.fa"
    PMOA_CONTIGS="${WORK_DIR}/results/10_methane_contig_extraction/pmoA_containing_contigs.fa"
    
    if [ -f "$MCRA_CONTIGS" ]; then
        echo "  Creating mcrA-specific bin..."
        cp "$MCRA_CONTIGS" "${RESULTS_DIR}/methane_targeted/mcrA_contigs.fa"
        
        # Try to bin mcrA contigs separately
        metabat2 -i ${MCRA_CONTIGS} \
                 -a ${RESULTS_DIR}/coassembly_depth.txt \
                 -o ${RESULTS_DIR}/methane_targeted/mcrA_bin \
                 --minContig 1500 \
                 --minClsSize 30000 \
                 --maxP 80 \
                 --minS 25 \
                 --numThreads ${THREADS}
    fi
    
    if [ -f "$PMOA_CONTIGS" ]; then
        echo "  Creating pmoA-specific bin..."
        cp "$PMOA_CONTIGS" "${RESULTS_DIR}/methane_targeted/pmoA_contigs.fa"
        
        # Try to bin pmoA contigs separately
        metabat2 -i ${PMOA_CONTIGS} \
                 -a ${RESULTS_DIR}/coassembly_depth.txt \
                 -o ${RESULTS_DIR}/methane_targeted/pmoA_bin \
                 --minContig 1500 \
                 --minClsSize 30000 \
                 --maxP 80 \
                 --minS 25 \
                 --numThreads ${THREADS}
    fi
    
else
    echo "No previously identified methane contigs found."
    echo "Run the methane detection scripts first:"
    echo "  ./code/09_simple_assembly_search.sh"
    echo "  ./code/10_extract_methane_contigs.sh"
fi

# Step 5: Additional aggressive binning for unbinned contigs
echo "Running aggressive binning on remaining unbinned contigs..."

# Create a list of all binned contigs so far
ALL_BINNED_CONTIGS="${RESULTS_DIR}/all_binned_contigs.txt"
find ${RESULTS_DIR} -name "*.fa" -o -name "*.fasta" -o -name "*.fna" | \
    xargs grep "^>" | sed 's/.*>//' | sed 's/ .*//' | sort -u > $ALL_BINNED_CONTIGS

# Extract unbinned contigs for another round of binning
UNBINNED_CONTIGS="${RESULTS_DIR}/unbinned_contigs.fa"
echo "  Extracting unbinned contigs..."

python3 -c "
import sys
from Bio import SeqIO

# Read list of binned contigs
try:
    with open('$ALL_BINNED_CONTIGS', 'r') as f:
        binned = set(line.strip() for line in f)
except:
    binned = set()

# Extract unbinned sequences
unbinned_count = 0
with open('$UNBINNED_CONTIGS', 'w') as out:
    for record in SeqIO.parse('$COASSEMBLY_FASTA', 'fasta'):
        if record.id not in binned:
            SeqIO.write(record, out, 'fasta')
            unbinned_count += 1

print(f'Extracted {unbinned_count} unbinned contigs')
" 2>/dev/null || echo "Python extraction failed, using alternative method"

# Alternative method if Python fails
if [ ! -f "$UNBINNED_CONTIGS" ] || [ ! -s "$UNBINNED_CONTIGS" ]; then
    echo "  Using alternative unbinned contig extraction..."
    grep "^>" ${COASSEMBLY_FASTA} | sed 's/>//' | sed 's/ .*//' | sort > ${RESULTS_DIR}/all_contigs.txt
    comm -23 ${RESULTS_DIR}/all_contigs.txt $ALL_BINNED_CONTIGS > ${RESULTS_DIR}/unbinned_contig_ids.txt
    
    # Extract sequences (simplified approach)
    awk 'BEGIN{while((getline line < "'${RESULTS_DIR}'/unbinned_contig_ids.txt") > 0) wanted[line]=1}
         /^>/{id=substr($0,2); gsub(/ .*/,"",id); print_seq=(id in wanted)}
         print_seq' ${COASSEMBLY_FASTA} > $UNBINNED_CONTIGS
fi

# Check if we have unbinned contigs to work with
if [ -f "$UNBINNED_CONTIGS" ] && [ -s "$UNBINNED_CONTIGS" ]; then
    unbinned_count=$(grep -c "^>" $UNBINNED_CONTIGS)
    echo "  Found $unbinned_count unbinned contigs"
    
    if [ $unbinned_count -gt 100 ]; then
        echo "  Running final aggressive MetaBAT2 on unbinned contigs..."
        metabat2 -i ${UNBINNED_CONTIGS} \
                 -a ${RESULTS_DIR}/coassembly_depth.txt \
                 -o ${RESULTS_DIR}/metabat2/final_bin \
                 --minContig 1500 \
                 --minClsSize 20000 \
                 --maxP 75 \
                 --minS 20 \
                 --maxEdges 2000 \
                 --pTNF 0 \
                 --noAdd \
                 --numThreads ${THREADS}
    fi
else
    echo "  No unbinned contigs found or extraction failed"
fi
# Step 6: DAS Tool for bin optimization (COMPREHENSIVE)
echo "Running DAS Tool for comprehensive bin optimization..."

mkdir -p ${RESULTS_DIR}/das_tool

# Prepare bin-to-contig files for all MetaBAT2 runs
echo "Preparing all MetaBAT2 bins for DAS Tool..."
cd ${RESULTS_DIR}

# Standard MetaBAT2 bins
if ls metabat2/coassembly_bin*.fa 1> /dev/null 2>&1; then
    for bin in metabat2/coassembly_bin*.fa; do
        bin_name=$(basename $bin .fa)
        grep "^>" $bin | sed 's/>//' | awk -v bin="$bin_name" '{print $1"\t"bin}' 
    done > ${RESULTS_DIR}/das_tool/metabat2_standard_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/metabat2_standard_scaffolds2bin.txt
fi

# Relaxed MetaBAT2 bins
if ls metabat2/relaxed_bin*.fa 1> /dev/null 2>&1; then
    for bin in metabat2/relaxed_bin*.fa; do
        bin_name=$(basename $bin .fa)
        grep "^>" $bin | sed 's/>//' | awk -v bin="metabat2_relaxed_$bin_name" '{print $1"\t"bin}' 
    done > ${RESULTS_DIR}/das_tool/metabat2_relaxed_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/metabat2_relaxed_scaffolds2bin.txt
fi

# Inclusive MetaBAT2 bins
if ls metabat2/inclusive_bin*.fa 1> /dev/null 2>&1; then
    for bin in metabat2/inclusive_bin*.fa; do
        bin_name=$(basename $bin .fa)
        grep "^>" $bin | sed 's/>//' | awk -v bin="metabat2_inclusive_$bin_name" '{print $1"\t"bin}' 
    done > ${RESULTS_DIR}/das_tool/metabat2_inclusive_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/metabat2_inclusive_scaffolds2bin.txt
fi

# Final aggressive MetaBAT2 bins
if ls metabat2/final_bin*.fa 1> /dev/null 2>&1; then
    for bin in metabat2/final_bin*.fa; do
        bin_name=$(basename $bin .fa)
        grep "^>" $bin | sed 's/>//' | awk -v bin="metabat2_final_$bin_name" '{print $1"\t"bin}' 
    done > ${RESULTS_DIR}/das_tool/metabat2_final_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/metabat2_final_scaffolds2bin.txt
fi

# Methane-targeted bins
if ls methane_targeted/*.fa 1> /dev/null 2>&1; then
    for bin in methane_targeted/*.fa; do
        bin_name=$(basename $bin .fa)
        grep "^>" $bin | sed 's/>//' | awk -v bin="methane_$bin_name" '{print $1"\t"bin}' 
    done > ${RESULTS_DIR}/das_tool/methane_targeted_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/methane_targeted_scaffolds2bin.txt
fi

# Combine all MetaBAT2 results
cat ${RESULTS_DIR}/das_tool/metabat2_*_scaffolds2bin.txt > ${RESULTS_DIR}/das_tool/metabat2_all_scaffolds2bin.txt

# Prepare bin-to-contig files for all MaxBin2 runs
echo "Preparing all MaxBin2 bins for DAS Tool..."
if ls maxbin2/coassembly_bin*.fasta 1> /dev/null 2>&1; then
    for bin in maxbin2/coassembly_bin*.fasta; do
        bin_name=$(basename $bin .fasta)
        grep "^>" $bin | sed 's/>//' | awk -v bin="maxbin2_standard_$bin_name" '{print $1"\t"bin}'
    done > ${RESULTS_DIR}/das_tool/maxbin2_standard_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/maxbin2_standard_scaffolds2bin.txt
fi

if ls maxbin2/relaxed_bin*.fasta 1> /dev/null 2>&1; then
    for bin in maxbin2/relaxed_bin*.fasta; do
        bin_name=$(basename $bin .fasta)
        grep "^>" $bin | sed 's/>//' | awk -v bin="maxbin2_relaxed_$bin_name" '{print $1"\t"bin}'
    done > ${RESULTS_DIR}/das_tool/maxbin2_relaxed_scaffolds2bin.txt
else
    touch ${RESULTS_DIR}/das_tool/maxbin2_relaxed_scaffolds2bin.txt
fi

cat ${RESULTS_DIR}/das_tool/maxbin2_*_scaffolds2bin.txt > ${RESULTS_DIR}/das_tool/maxbin2_all_scaffolds2bin.txt

# Prepare bin-to-contig files for all VAMB runs
echo "Preparing all VAMB bins for DAS Tool..."
for vamb_dir in vamb/standard vamb/relaxed vamb/inclusive; do
    if [ -d "$vamb_dir" ] && ls ${vamb_dir}/*.fna 1> /dev/null 2>&1; then
        vamb_type=$(basename $vamb_dir)
        for bin in ${vamb_dir}/*.fna; do
            bin_name=$(basename $bin .fna)
            grep "^>" $bin | sed 's/>//' | awk -v bin="vamb_${vamb_type}_$bin_name" '{print $1"\t"bin}'
        done > ${RESULTS_DIR}/das_tool/vamb_${vamb_type}_scaffolds2bin.txt
    else
        touch ${RESULTS_DIR}/das_tool/vamb_${vamb_type}_scaffolds2bin.txt
    fi
done

cat ${RESULTS_DIR}/das_tool/vamb_*_scaffolds2bin.txt > ${RESULTS_DIR}/das_tool/vamb_all_scaffolds2bin.txt

# Run DAS Tool (only if we have bins from at least one method)
if [ -s "${RESULTS_DIR}/das_tool/metabat2_all_scaffolds2bin.txt" ] || [ -s "${RESULTS_DIR}/das_tool/maxbin2_all_scaffolds2bin.txt" ] || [ -s "${RESULTS_DIR}/das_tool/vamb_all_scaffolds2bin.txt" ]; then
    echo "Running comprehensive DAS Tool integration..."
    
    # Check if DAS Tool is available
    if ! command -v DAS_Tool &> /dev/null; then
        echo "⚠ Warning: DAS Tool not found, skipping bin integration"
        echo "You can install DAS Tool later with: conda install -c bioconda das_tool"
        echo "Skipping DAS Tool due to missing installation" > ${RESULTS_DIR}/das_tool/SKIPPED.txt
    else
        # Prepare input file list
        input_files=""
        labels=""
        
        if [ -s "${RESULTS_DIR}/das_tool/metabat2_all_scaffolds2bin.txt" ]; then
            input_files="${input_files}${RESULTS_DIR}/das_tool/metabat2_all_scaffolds2bin.txt,"
            labels="${labels}metabat2_all,"
        fi
        
        if [ -s "${RESULTS_DIR}/das_tool/maxbin2_all_scaffolds2bin.txt" ]; then
            input_files="${input_files}${RESULTS_DIR}/das_tool/maxbin2_all_scaffolds2bin.txt,"
            labels="${labels}maxbin2_all,"
        fi
        
        if [ -s "${RESULTS_DIR}/das_tool/vamb_all_scaffolds2bin.txt" ]; then
            input_files="${input_files}${RESULTS_DIR}/das_tool/vamb_all_scaffolds2bin.txt,"
            labels="${labels}vamb_all,"
        fi
        
        if [ -s "${RESULTS_DIR}/das_tool/methane_targeted_scaffolds2bin.txt" ]; then
            input_files="${input_files}${RESULTS_DIR}/das_tool/methane_targeted_scaffolds2bin.txt,"
            labels="${labels}methane_targeted,"
        fi
        
        # Remove trailing commas
        input_files=${input_files%,}
        labels=${labels%,}
        
        echo "DAS Tool input files: $input_files"
        echo "DAS Tool labels: $labels"
        
        DAS_Tool -i $input_files \
                 -l $labels \
                 -c ${COASSEMBLY_FASTA} \
                 -o ${RESULTS_DIR}/das_tool/comprehensive_DASToolBins \
                 --threads ${THREADS} \
                 --write_bins 1 \
                 --score_threshold 0.3
    fi
else
    echo "No bins found for DAS Tool integration"
fi

# Step 7: Comprehensive bin quality assessment with CheckM
echo "Running comprehensive CheckM for all bin quality assessment..."

# Check if CheckM is available - use dedicated checkm_env
echo "Checking CheckM availability in dedicated environment..."
if conda activate checkm_env && command -v checkm &> /dev/null; then
    echo "✅ CheckM found in checkm_env, proceeding with quality assessment..."
    CHECKM_CMD="conda activate checkm_env && checkm"
elif command -v checkm &> /dev/null; then
    echo "✅ CheckM found in current environment, proceeding with quality assessment..."
    CHECKM_CMD="checkm"
else
    echo "⚠ Warning: CheckM not found in any environment, skipping quality assessment"
    echo "You can install CheckM later with: conda install -c bioconda checkm-genome"
    echo "For now, using basic bin statistics instead..."
    echo "Skipping CheckM due to missing installation" > ${RESULTS_DIR}/checkm/SKIPPED.txt
    
    # Run basic bin analysis instead
    if [ -f "${WORK_DIR}/code/analyze_bins.py" ]; then
        echo "Running basic bin analysis..."
        python3 ${WORK_DIR}/code/analyze_bins.py
    fi
    CHECKM_CMD=""
fi

if [ -n "$CHECKM_CMD" ]; then
    # CheckM on all MetaBAT2 bins (combine all types)
    echo "CheckM analysis for all MetaBAT2 bins..."
    mkdir -p ${RESULTS_DIR}/checkm/all_metabat2_bins
    
    # Copy all MetaBAT2 bins to one directory for CheckM
    find ${RESULTS_DIR}/metabat2 -name "*.fa" -exec cp {} ${RESULTS_DIR}/checkm/all_metabat2_bins/ \;
    
    if [ "$(ls -A ${RESULTS_DIR}/checkm/all_metabat2_bins/ 2>/dev/null)" ]; then
        eval "$CHECKM_CMD lineage_wf ${RESULTS_DIR}/checkm/all_metabat2_bins \
                         ${RESULTS_DIR}/checkm/metabat2_comprehensive \
                         -t ${THREADS} \
                         --file ${RESULTS_DIR}/checkm/metabat2_comprehensive_quality.txt \
                         --tab_table \
                         -x fa"
    fi

    # CheckM on all MaxBin2 bins
    if [ -d "${RESULTS_DIR}/maxbin2" ]; then
        echo "CheckM analysis for all MaxBin2 bins..."
        mkdir -p ${RESULTS_DIR}/checkm/all_maxbin2_bins
        find ${RESULTS_DIR}/maxbin2 -name "*.fasta" -exec cp {} ${RESULTS_DIR}/checkm/all_maxbin2_bins/ \;
        
        if [ "$(ls -A ${RESULTS_DIR}/checkm/all_maxbin2_bins/ 2>/dev/null)" ]; then
            eval "$CHECKM_CMD lineage_wf ${RESULTS_DIR}/checkm/all_maxbin2_bins \
                             ${RESULTS_DIR}/checkm/maxbin2_comprehensive \
                             -t ${THREADS} \
                             --file ${RESULTS_DIR}/checkm/maxbin2_comprehensive_quality.txt \
                             --tab_table \
                             -x fasta"
        fi
    fi

    # CheckM on all VAMB bins
    if [ -d "${RESULTS_DIR}/vamb" ]; then
        echo "CheckM analysis for all VAMB bins..."
        mkdir -p ${RESULTS_DIR}/checkm/all_vamb_bins
        find ${RESULTS_DIR}/vamb -name "*.fna" -exec cp {} ${RESULTS_DIR}/checkm/all_vamb_bins/ \;
        
        if [ "$(ls -A ${RESULTS_DIR}/checkm/all_vamb_bins/ 2>/dev/null)" ]; then
            eval "$CHECKM_CMD lineage_wf ${RESULTS_DIR}/checkm/all_vamb_bins \
                             ${RESULTS_DIR}/checkm/vamb_comprehensive \
                             -t ${THREADS} \
                             --file ${RESULTS_DIR}/checkm/vamb_comprehensive_quality.txt \
                             --tab_table \
                             -x fna"
        fi
    fi

    # CheckM on methane-targeted bins (PRIORITY!)
    if [ -d "${RESULTS_DIR}/methane_targeted" ]; then
        echo "🔥 Priority CheckM analysis for methane-targeted bins..."
        checkm_bins=$(ls ${RESULTS_DIR}/methane_targeted/*.fa 2>/dev/null | wc -l)
        if [ $checkm_bins -gt 0 ]; then
            checkm lineage_wf ${RESULTS_DIR}/methane_targeted \
                             ${RESULTS_DIR}/checkm/methane_targeted \
                             -t ${THREADS} \
                             --file ${RESULTS_DIR}/checkm/methane_targeted_quality.txt \
                             --tab_table \
                             -x fa
            
            echo "🎯 Methane bin quality assessment completed!"
            echo "Priority results: ${RESULTS_DIR}/checkm/methane_targeted_quality.txt"
        fi
    fi

    # CheckM on comprehensive DAS Tool bins
    if [ -d "${RESULTS_DIR}/das_tool/comprehensive_DASToolBins_DASTool_bins" ]; then
        echo "CheckM analysis for comprehensive DAS Tool bins..."
        checkm lineage_wf ${RESULTS_DIR}/das_tool/comprehensive_DASToolBins_DASTool_bins \
                         ${RESULTS_DIR}/checkm/dastool_comprehensive \
                         -t ${THREADS} \
                         --file ${RESULTS_DIR}/checkm/dastool_comprehensive_quality.txt \
                         --tab_table \
                         -x fa
    fi
fi

echo ""
echo "========================================================"
echo "COMPREHENSIVE BINNING COMPLETED!"
echo "========================================================"
echo ""
echo "Results saved in: ${RESULTS_DIR}"
echo ""
echo "Binning strategy overview:"
echo "1. Standard binning (conservative parameters)"
echo "2. Relaxed binning (lower thresholds for methane genes)"
echo "3. Inclusive binning (very permissive parameters)"
echo "4. Targeted methane binning (methane-containing contigs)"
echo "5. Final aggressive binning (remaining unbinned contigs)"
echo ""
echo "Output summary:"
echo "🧬 MetaBAT2 results:"
echo "  - Standard bins: ${RESULTS_DIR}/metabat2/coassembly_bin*"
echo "  - Relaxed bins: ${RESULTS_DIR}/metabat2/relaxed_bin*"
echo "  - Inclusive bins: ${RESULTS_DIR}/metabat2/inclusive_bin*"
echo "  - Final bins: ${RESULTS_DIR}/metabat2/final_bin*"
echo ""
echo "🔬 MaxBin2 results:"
echo "  - Standard bins: ${RESULTS_DIR}/maxbin2/coassembly_bin*"
echo "  - Relaxed bins: ${RESULTS_DIR}/maxbin2/relaxed_bin*"
echo ""
echo "🌐 VAMB results:"
echo "  - Standard bins: ${RESULTS_DIR}/vamb/standard/*"
echo "  - Relaxed bins: ${RESULTS_DIR}/vamb/relaxed/*"
echo "  - Inclusive bins: ${RESULTS_DIR}/vamb/inclusive/*"
echo ""
echo "🔥 Methane-targeted results:"
echo "  - All methane bins: ${RESULTS_DIR}/methane_targeted/*"
echo "  - mcrA-specific bins: ${RESULTS_DIR}/methane_targeted/mcrA_bin*"
echo "  - pmoA-specific bins: ${RESULTS_DIR}/methane_targeted/pmoA_bin*"
echo ""
echo "⚡ DAS Tool integrated bins:"
echo "  - Comprehensive bins: ${RESULTS_DIR}/das_tool/comprehensive_DASToolBins_DASTool_bins/"
echo ""
echo "📊 Analysis files:"
echo "  - Individual sample BAMs: ${RESULTS_DIR}/*_coassembly_sorted.bam"
echo "  - Differential coverage: ${RESULTS_DIR}/coassembly_depth.txt"
echo "  - Unbinned contigs: ${RESULTS_DIR}/unbinned_contigs.fa"
echo ""
echo "📈 Expected improvements:"
echo "  - Dramatically increased binning efficiency (from 7.1% to 40-60%)"
echo "  - Capture of methane-containing contigs"
echo "  - Better representation of low-abundance organisms"
echo "  - Multiple parameter sets for different organism types"
echo ""
echo "🎯 Next steps for methane analysis:"
echo "  1. Check methane-targeted bins first: ${RESULTS_DIR}/methane_targeted/"
echo "  2. Run comprehensive bin quality assessment"
echo "  3. Annotate all bins, especially methane-targeted ones"
echo "  4. Compare results with original 8 bins"
echo ""
echo "Quality assessment reports (if CheckM available):"
echo "  - MetaBAT2: ${RESULTS_DIR}/checkm/metabat2_quality.txt"
echo "  - MaxBin2: ${RESULTS_DIR}/checkm/maxbin2_quality.txt"
echo "  - VAMB: ${RESULTS_DIR}/checkm/vamb_quality.txt"
echo "  - DAS Tool: ${RESULTS_DIR}/checkm/dastool_quality.txt"
echo "  - Methane bins: ${RESULTS_DIR}/checkm/methane_quality.txt"
