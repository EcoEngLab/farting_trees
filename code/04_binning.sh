#!/bin/bash

# Resumable Metagenomic Binning Pipeline
# This script allows you to run individual steps or resume from a specific step
# Author: Generated for farting_trees project
# Date: August 2025

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/03_assembly"
RESULTS_DIR="${WORK_DIR}/results/04_binning"
THREADS=8
COASSEMBLY_FASTA="${ASSEMBLY_DIR}/megahit/coassembly_contigs.fa"

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Usage function
usage() {
    echo "Usage: $0 [STEP]"
    echo ""
    echo "Available steps:"
    echo "  prep     - Prepare input data (depth file)"
    echo "  metabat  - Run MetaBAT2 binning"
    echo "  maxbin   - Run MaxBin2 binning"
    echo "  vamb     - Run VAMB binning"
    echo "  dastool  - Run DAS Tool optimization"
    echo "  checkm   - Run CheckM quality assessment"
    echo "  extract  - Extract high-quality bins"
    echo "  all      - Run all steps (default)"
    echo ""
    echo "Example:"
    echo "  $0 vamb     # Run only VAMB step"
    echo "  $0 all      # Run all steps"
    echo ""
    exit 1
}

# Get step from command line argument
STEP=${1:-all}

echo "==============================================================="
echo "RESUMABLE METAGENOMIC BINNING PIPELINE"
echo "==============================================================="
echo "Running step: $STEP"
echo "Co-assembly: $COASSEMBLY_FASTA"
echo "Results: $RESULTS_DIR"
echo ""

# Create output directories
mkdir -p ${RESULTS_DIR}/{metabat2,maxbin2,vamb,dastool,checkm,final_bins}

# Step functions
prep_data() {
    echo "STEP 1: Preparing input data..."
    
    # Check if depth file exists
    DEPTH_FILE="${RESULTS_DIR}/coassembly_depth.txt"
    if [ ! -f "$DEPTH_FILE" ]; then
        echo "  Generating differential coverage depth file..."
        conda activate metagenome_assembly
        
        # Check for BAM files
        if ! ls ${RESULTS_DIR}/*_coassembly_sorted.bam 1> /dev/null 2>&1; then
            echo "  ERROR: No BAM files found! Run read mapping first."
            echo "  Expected: ${RESULTS_DIR}/*_coassembly_sorted.bam"
            exit 1
        fi
        
        echo "  Found BAM files:"
        ls -la ${RESULTS_DIR}/*_coassembly_sorted.bam
        
        conda activate metagenome_binning
        jgi_summarize_bam_contig_depths --outputDepth ${DEPTH_FILE} \
                                       ${RESULTS_DIR}/*_coassembly_sorted.bam
    else
        echo "  ✓ Depth file exists: $DEPTH_FILE"
    fi
    
    echo "  ✓ Input data ready: assembly.fa + depth.txt"
}

run_metabat2() {
    echo "STEP 2.1: Running MetaBAT2 binning..."
    conda activate metagenome_binning
    
    # Check if MetaBAT2 already completed successfully
    if [ -d "${RESULTS_DIR}/metabat2" ] && [ $(ls ${RESULTS_DIR}/metabat2/bin*.fa 2>/dev/null | wc -l) -gt 0 ]; then
        METABAT2_BINS=$(ls ${RESULTS_DIR}/metabat2/bin*.fa 2>/dev/null | wc -l)
        echo "  ✓ MetaBAT2 already completed: $METABAT2_BINS bins found, skipping..."
    else
        echo "  Running MetaBAT2 with standard parameters..."
        metabat2 -i ${COASSEMBLY_FASTA} \
                 -a ${RESULTS_DIR}/coassembly_depth.txt \
                 -o ${RESULTS_DIR}/metabat2/bin \
                 --minContig 2500 \
                 --numThreads ${THREADS}

        METABAT2_BINS=$(ls ${RESULTS_DIR}/metabat2/bin*.fa 2>/dev/null | wc -l)
        echo "  ✓ MetaBAT2 completed: $METABAT2_BINS bins generated"
    fi
}

run_maxbin2() {
    echo "STEP 2.2: Running MaxBin2 binning..."
    conda activate maxbin2_env
    
    # Check if MaxBin2 already completed successfully
    if [ -d "${RESULTS_DIR}/maxbin2" ] && [ $(ls ${RESULTS_DIR}/maxbin2/bin*.fasta 2>/dev/null | wc -l) -gt 0 ]; then
        MAXBIN2_BINS=$(ls ${RESULTS_DIR}/maxbin2/bin*.fasta 2>/dev/null | wc -l)
        echo "  ✓ MaxBin2 already completed: $MAXBIN2_BINS bins found, skipping..."
    elif ! command -v run_MaxBin.pl &> /dev/null; then
        echo "  ⚠ MaxBin2 not found, skipping..."
        mkdir -p ${RESULTS_DIR}/maxbin2
        echo "MaxBin2 not available" > ${RESULTS_DIR}/maxbin2/SKIPPED.txt
        MAXBIN2_BINS=0
    else
        # Create abundance files for MaxBin2
        echo "  Preparing abundance files..."
        ABUND_FILES=""
        SAMPLE_COUNT=0
        
        for i in $(seq 4 2 20); do  # Sample depth columns
            col_data=$(cut -f${i} ${RESULTS_DIR}/coassembly_depth.txt | head -20 | tail -10)
            if [[ -n "$col_data" ]] && [[ "$col_data" != *"0.0000"* ]]; then
                SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
                cut -f1,${i} ${RESULTS_DIR}/coassembly_depth.txt > ${RESULTS_DIR}/maxbin2_abundance_${SAMPLE_COUNT}.txt
                if [ $SAMPLE_COUNT -eq 1 ]; then
                    ABUND_FILES="${ABUND_FILES} -abund ${RESULTS_DIR}/maxbin2_abundance_${SAMPLE_COUNT}.txt"
                else
                    ABUND_FILES="${ABUND_FILES} -abund${SAMPLE_COUNT} ${RESULTS_DIR}/maxbin2_abundance_${SAMPLE_COUNT}.txt"
                fi
            fi
        done
        
        echo "  Running MaxBin2 with $SAMPLE_COUNT abundance files..."
        
        # Check if this is the old MaxBin that requires a seed file
        if MaxBin 2>&1 | grep -q "seed file"; then
            echo "  Detected old MaxBin version requiring seed file, creating proper seed..."
            
            # Create seed file by extracting only contig names (headers without '>')
            echo "  Extracting first 10 contig names as seeds..."
            grep "^>" ${COASSEMBLY_FASTA} | head -10 | sed 's/^>//' | sed 's/ .*//' > ${RESULTS_DIR}/maxbin2_seed.txt
            
            # Verify seed file contains contig names
            seed_count=$(wc -l < ${RESULTS_DIR}/maxbin2_seed.txt 2>/dev/null || echo 0)
            echo "  Created seed file with $seed_count contig names"
            
            if [ $seed_count -lt 1 ]; then
                echo "  ⚠ Seed file creation failed, skipping MaxBin2..."
                mkdir -p ${RESULTS_DIR}/maxbin2
                echo "Seed file creation failed" > ${RESULTS_DIR}/maxbin2/SKIPPED.txt
                MAXBIN2_BINS=0
            else
                # Use direct MaxBin command with proper seed file (contig names only)
                echo "  Running MaxBin with seed file containing $seed_count contig names..."
                MaxBin -fasta ${COASSEMBLY_FASTA} \
                       -seed ${RESULTS_DIR}/maxbin2_seed.txt \
                       ${ABUND_FILES} \
                       -out ${RESULTS_DIR}/maxbin2/bin \
                       -thread ${THREADS}
            fi
        else
            # Try the wrapper for newer versions
            run_MaxBin.pl -fasta ${COASSEMBLY_FASTA} \
                          ${ABUND_FILES} \
                          -out ${RESULTS_DIR}/maxbin2/bin \
                          -thread ${THREADS} \
                          -min_contig_length 2500
        fi
        
        # Count MaxBin2 bins
        MAXBIN2_BINS=$(ls ${RESULTS_DIR}/maxbin2/bin*.fasta 2>/dev/null | wc -l)
        echo "  ✓ MaxBin2 completed: $MAXBIN2_BINS bins generated"
    fi
}

run_vamb() {
    echo "STEP 2.3: Running VAMB binning..."
    conda activate metagenome_binning
    
    # Check if VAMB already completed successfully
    if [ -d "${RESULTS_DIR}/vamb/bins" ] && [ $(ls ${RESULTS_DIR}/vamb/bins/*.fna 2>/dev/null | wc -l) -gt 0 ]; then
        VAMB_BINS=$(ls ${RESULTS_DIR}/vamb/bins/*.fna 2>/dev/null | wc -l)
        echo "  ✓ VAMB already completed: $VAMB_BINS bins found, skipping..."
    elif ! command -v vamb &> /dev/null; then
        echo "  ⚠ VAMB not found, skipping..."
        mkdir -p ${RESULTS_DIR}/vamb
        echo "VAMB not available" > ${RESULTS_DIR}/vamb/SKIPPED.txt
        VAMB_BINS=0
    else
        echo "  Running VAMB..."
        # Remove existing VAMB directory to avoid FileExistsError
        rm -rf ${RESULTS_DIR}/vamb
        
        vamb --outdir ${RESULTS_DIR}/vamb \
             --fasta ${COASSEMBLY_FASTA} \
             --bamfiles ${RESULTS_DIR}/*_coassembly_sorted.bam \
             -m 2500 \
             --minfasta 200000
        
        # Count VAMB bins
        VAMB_BINS=$(ls ${RESULTS_DIR}/vamb/bins/*.fna 2>/dev/null | wc -l)
        echo "  ✓ VAMB completed: $VAMB_BINS bins generated"
    fi
}

run_dastool() {
    echo "STEP 3: Running DAS Tool for bin optimization..."
    conda activate dastool_env
    
    # Check if DAS Tool already completed successfully
    if [ -d "${RESULTS_DIR}/dastool/optimized_DASTool_bins" ] && [ $(ls ${RESULTS_DIR}/dastool/optimized_DASTool_bins/*.fa 2>/dev/null | wc -l) -gt 0 ]; then
        echo "  ✓ DAS Tool already completed, skipping..."
    elif ! command -v DAS_Tool &> /dev/null; then
        echo "  ⚠ DAS Tool not found, skipping optimization..."
        mkdir -p ${RESULTS_DIR}/dastool
        echo "DAS Tool not available" > ${RESULTS_DIR}/dastool/SKIPPED.txt
    else
        # Try to fix DAS Tool R dependency issues
        echo "  Checking DAS Tool R dependencies..."
        
        # Test if DAS Tool works with a simple help command
        if ! DAS_Tool --help &>/dev/null; then
            echo "  ⚠ DAS Tool has R dependency issues, attempting to fix..."
            
            # Try to reinstall docopt in R
            R --slave -e "if (!require('docopt', quietly=TRUE)) { install.packages('docopt', repos='https://cran.r-project.org/', dependencies=TRUE) }" 2>/dev/null || true
            
            # If still not working, skip DAS Tool
            if ! DAS_Tool --help &>/dev/null; then
                echo "  ⚠ Could not fix DAS Tool R dependencies, skipping..."
                echo "  You can manually fix this later with:"
                echo "    conda activate dastool_env"
                echo "    R --slave -e \"install.packages('docopt', repos='https://cran.r-project.org/')\""
                mkdir -p ${RESULTS_DIR}/dastool
                echo "DAS Tool R dependency error" > ${RESULTS_DIR}/dastool/SKIPPED.txt
                return
            fi
        fi
        
        # Find DAS Tool script directory
        DASTOOL_SCRIPTS=$(dirname $(which DAS_Tool))
        
        # Prepare bin-to-contig tables for DAS Tool
        echo "  Preparing bin-to-contig mapping files..."
        
        # Initialize variables
        BIN_SETS=""
        BIN_LABELS=""
        
        # MetaBAT2 mapping
        if [ -d "${RESULTS_DIR}/metabat2" ] && [ $(ls ${RESULTS_DIR}/metabat2/bin*.fa 2>/dev/null | wc -l) -gt 0 ]; then
            ${DASTOOL_SCRIPTS}/Fasta_to_Contig2Bin.sh -i ${RESULTS_DIR}/metabat2 \
                                                       -e fa > ${RESULTS_DIR}/metabat2_contigs2bins.tsv
            BIN_SETS="${BIN_SETS},${RESULTS_DIR}/metabat2_contigs2bins.tsv"
            BIN_LABELS="${BIN_LABELS},MetaBAT2"
        fi
        
        # MaxBin2 mapping
        if [ -d "${RESULTS_DIR}/maxbin2" ] && [ $(ls ${RESULTS_DIR}/maxbin2/bin*.fasta 2>/dev/null | wc -l) -gt 0 ]; then
            ${DASTOOL_SCRIPTS}/Fasta_to_Contig2Bin.sh -i ${RESULTS_DIR}/maxbin2 \
                                                       -e fasta > ${RESULTS_DIR}/maxbin2_contigs2bins.tsv
            BIN_SETS="${BIN_SETS},${RESULTS_DIR}/maxbin2_contigs2bins.tsv"
            BIN_LABELS="${BIN_LABELS},MaxBin2"
        fi
        
        # VAMB mapping
        if [ -d "${RESULTS_DIR}/vamb/bins" ] && [ $(ls ${RESULTS_DIR}/vamb/bins/*.fna 2>/dev/null | wc -l) -gt 0 ]; then
            ${DASTOOL_SCRIPTS}/Fasta_to_Contig2Bin.sh -i ${RESULTS_DIR}/vamb/bins \
                                                       -e fna > ${RESULTS_DIR}/vamb_contigs2bins.tsv
            BIN_SETS="${BIN_SETS},${RESULTS_DIR}/vamb_contigs2bins.tsv"
            BIN_LABELS="${BIN_LABELS},VAMB"
        fi
        
        # Remove leading commas
        BIN_SETS=${BIN_SETS#,}
        BIN_LABELS=${BIN_LABELS#,}
        
        if [ -n "$BIN_SETS" ]; then
            echo "  Running DAS Tool with bin sets: $BIN_LABELS"
            
            # Try running DAS Tool with error handling
            echo "  Command: DAS_Tool -i \"$BIN_SETS\" -l \"$BIN_LABELS\" -c \"${COASSEMBLY_FASTA}\" -o \"${RESULTS_DIR}/dastool/optimized\" --threads ${THREADS} --write_bins"
            
            if DAS_Tool -i "$BIN_SETS" \
                        -l "$BIN_LABELS" \
                        -c "${COASSEMBLY_FASTA}" \
                        -o "${RESULTS_DIR}/dastool/optimized" \
                        --threads ${THREADS} \
                        --write_bins; then
                echo "  ✓ DAS Tool optimization completed"
            else
                echo "  ⚠ DAS Tool failed, but continuing pipeline..."
                echo "DAS Tool execution failed" > ${RESULTS_DIR}/dastool/FAILED.txt
            fi
        else
            echo "  ⚠ No valid bin sets found for optimization"
        fi
    fi
}

run_checkm() {
    echo "STEP 4: Running CheckM quality assessment..."
    conda activate checkm_env
    
    # Check if CheckM already completed successfully
    if [ -f "${RESULTS_DIR}/checkm/results/storage/bin_stats_ext.tsv" ]; then
        echo "  ✓ CheckM already completed, skipping..."
    elif ! command -v checkm &> /dev/null; then
        echo "  ⚠ CheckM not found, skipping quality assessment..."
        mkdir -p ${RESULTS_DIR}/checkm
        echo "CheckM not available" > ${RESULTS_DIR}/checkm/SKIPPED.txt
    else
        # Check CheckM setup first
        echo "  Checking CheckM setup..."
        
        # Test CheckM basic functionality
        if ! checkm --version &>/dev/null; then
            echo "  ⚠ CheckM version check failed, there may be setup issues"
        fi
        
        # Check for required dependencies
        missing_deps=""
        if ! command -v prodigal &>/dev/null; then
            missing_deps="${missing_deps} prodigal"
        fi
        if ! command -v pplacer &>/dev/null; then
            missing_deps="${missing_deps} pplacer"
        fi
        if ! command -v hmmsearch &>/dev/null; then
            missing_deps="${missing_deps} hmmer"
        fi
        
        if [ -n "$missing_deps" ]; then
            echo "  ⚠ Missing dependencies:$missing_deps"
            echo "  Install with: conda install -c bioconda$missing_deps"
        fi
        
        # Assess quality of all bins
        echo "  Assessing bin quality with CheckM..."
        
        # Create directory for all bins
        mkdir -p ${RESULTS_DIR}/checkm/all_bins
        
        # Copy all bins to one directory for CheckM
        echo "  Collecting bins from different binning methods..."
        bin_count=0
        
        if [ -d "${RESULTS_DIR}/metabat2" ]; then
            metabat2_bins=$(ls ${RESULTS_DIR}/metabat2/bin*.fa 2>/dev/null | wc -l)
            if [ $metabat2_bins -gt 0 ]; then
                cp ${RESULTS_DIR}/metabat2/bin*.fa ${RESULTS_DIR}/checkm/all_bins/ 2>/dev/null || true
                bin_count=$((bin_count + metabat2_bins))
                echo "    MetaBAT2: $metabat2_bins bins"
            fi
        fi
        
        if [ -d "${RESULTS_DIR}/maxbin2" ]; then
            maxbin2_bins=$(ls ${RESULTS_DIR}/maxbin2/bin*.fasta 2>/dev/null | wc -l)
            if [ $maxbin2_bins -gt 0 ]; then
                cp ${RESULTS_DIR}/maxbin2/bin*.fasta ${RESULTS_DIR}/checkm/all_bins/ 2>/dev/null || true
                bin_count=$((bin_count + maxbin2_bins))
                echo "    MaxBin2: $maxbin2_bins bins"
            fi
        fi
        
        if [ -d "${RESULTS_DIR}/vamb/bins" ]; then
            vamb_bins=$(ls ${RESULTS_DIR}/vamb/bins/*.fna 2>/dev/null | wc -l)
            if [ $vamb_bins -gt 0 ]; then
                cp ${RESULTS_DIR}/vamb/bins/*.fna ${RESULTS_DIR}/checkm/all_bins/ 2>/dev/null || true
                bin_count=$((bin_count + vamb_bins))
                echo "    VAMB: $vamb_bins bins"
            fi
        fi
        
        # If DAS Tool was successful, also include optimized bins
        if [ -d "${RESULTS_DIR}/dastool/optimized_DASTool_bins" ]; then
            dastool_bins=$(ls ${RESULTS_DIR}/dastool/optimized_DASTool_bins/*.fa 2>/dev/null | wc -l)
            if [ $dastool_bins -gt 0 ]; then
                cp ${RESULTS_DIR}/dastool/optimized_DASTool_bins/*.fa ${RESULTS_DIR}/checkm/all_bins/ 2>/dev/null || true
                bin_count=$((bin_count + dastool_bins))
                echo "    DAS Tool: $dastool_bins bins"
            fi
        fi
        
        # Convert all bin files to .fa extension for consistency
        echo "  Standardizing file extensions..."
        cd ${RESULTS_DIR}/checkm/all_bins
        for file in *.fasta; do
            [ -f "$file" ] && mv "$file" "${file%.fasta}.fa"
        done
        for file in *.fna; do
            [ -f "$file" ] && mv "$file" "${file%.fna}.fa"
        done
        cd - > /dev/null
        
        # Final count after processing
        final_bin_count=$(ls ${RESULTS_DIR}/checkm/all_bins/*.fa 2>/dev/null | wc -l)
        echo "  Total bins collected: $final_bin_count"
        
        if [ $final_bin_count -eq 0 ]; then
            echo "  ⚠ No bins found for CheckM analysis"
            echo "  Make sure binning steps completed successfully first"
            mkdir -p ${RESULTS_DIR}/checkm
            echo "No bins found for analysis" > ${RESULTS_DIR}/checkm/SKIPPED.txt
            return
        fi
        
        # Check bin file sizes
        echo "  Checking bin sizes..."
        small_bins=$(find ${RESULTS_DIR}/checkm/all_bins -name "*.fa" -size -1k | wc -l)
        if [ $small_bins -gt 0 ]; then
            echo "  ⚠ Warning: $small_bins bins are smaller than 1KB"
        fi
        
        # Run CheckM lineage workflow with error handling
        echo "  Running CheckM lineage workflow..."
        echo "  This may take a while depending on the number of bins..."
        
        # Create results directory
        mkdir -p ${RESULTS_DIR}/checkm/results
        
        # Run CheckM with timeout and error handling
        if timeout 3600 checkm lineage_wf -t ${THREADS} \
                                         -x fa \
                                         ${RESULTS_DIR}/checkm/all_bins \
                                         ${RESULTS_DIR}/checkm/results 2>&1 | tee ${RESULTS_DIR}/checkm/checkm.log; then
            
            # Check if results were generated
            if [ -f "${RESULTS_DIR}/checkm/results/storage/bin_stats_ext.tsv" ]; then
                echo "  ✓ CheckM quality assessment completed successfully"
                
                # Show summary
                result_lines=$(wc -l < ${RESULTS_DIR}/checkm/results/storage/bin_stats_ext.tsv)
                echo "  Results: $((result_lines - 1)) bins analyzed"
                
                # Show high-quality bins preview
                echo "  High-quality bins (>90% complete, <5% contamination):"
                awk 'NR>1 && $12>=90 && $13<=5 {print "    " $1 ": " $12 "% complete, " $13 "% contamination"}' \
                    ${RESULTS_DIR}/checkm/results/storage/bin_stats_ext.tsv | head -5
                
            else
                echo "  ⚠ CheckM completed but results file not found"
                echo "  Check log: ${RESULTS_DIR}/checkm/checkm.log"
            fi
        else
            echo "  ⚠ CheckM failed or timed out"
            echo "  Check log: ${RESULTS_DIR}/checkm/checkm.log"
            echo "  You can try running with fewer threads or bins"
            echo "CheckM execution failed" > ${RESULTS_DIR}/checkm/FAILED.txt
        fi
    fi
}

extract_hq_bins() {
    echo "STEP 5: Extracting high-quality bins..."
    
    mkdir -p ${RESULTS_DIR}/high_quality_bins
    
    if [ -f "${RESULTS_DIR}/checkm/results/storage/bin_stats_ext.tsv" ]; then
        echo "  Filtering bins by quality criteria:"
        echo "    - Completeness ≥ 90%"
        echo "    - Contamination ≤ 5%"
        
        # Extract high-quality bins based on CheckM results
        awk 'NR>1 && $12>=90 && $13<=5 {print $1}' ${RESULTS_DIR}/checkm/results/storage/bin_stats_ext.tsv > ${RESULTS_DIR}/high_quality_bin_list.txt
        
        HQ_COUNT=0
        while read -r bin_name; do
            if [ -f "${RESULTS_DIR}/checkm/all_bins/${bin_name}.fa" ]; then
                cp "${RESULTS_DIR}/checkm/all_bins/${bin_name}.fa" "${RESULTS_DIR}/high_quality_bins/"
                HQ_COUNT=$((HQ_COUNT + 1))
            fi
        done < ${RESULTS_DIR}/high_quality_bin_list.txt
        
        echo "  ✓ Extracted $HQ_COUNT high-quality bins"
    else
        echo "  ⚠ CheckM results not found, extracting by size (>1MB) as proxy..."
        
        # Fallback: extract bins larger than 1MB as potential high-quality
        find ${RESULTS_DIR}/checkm/all_bins -name "*.fa" -size +1M -exec cp {} ${RESULTS_DIR}/high_quality_bins/ \;
        
        HQ_COUNT=$(ls -1 ${RESULTS_DIR}/high_quality_bins/*.fa 2>/dev/null | wc -l)
        echo "  ✓ Extracted $HQ_COUNT large bins (>1MB) as potential high-quality"
    fi
}

# Main execution logic
case $STEP in
    prep)
        prep_data
        ;;
    metabat)
        prep_data
        run_metabat2
        ;;
    maxbin)
        prep_data
        run_maxbin2
        ;;
    vamb)
        prep_data
        run_vamb
        ;;
    dastool)
        run_dastool
        ;;
    checkm)
        run_checkm
        ;;
    extract)
        extract_hq_bins
        ;;
    all)
        prep_data
        run_metabat2
        run_maxbin2
        run_vamb
        run_dastool
        run_checkm
        extract_hq_bins
        ;;
    *)
        echo "Error: Unknown step '$STEP'"
        usage
        ;;
esac

echo ""
echo "========================================================"
echo "STEP '$STEP' COMPLETED!"
echo "========================================================"
echo "Results location: ${RESULTS_DIR}"
echo ""
