#!/bin/bash
# Metagenomic Binning Pipeline
# Author: Generated for farting_trees project

set -e

# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/02_assembly"
MAPPING_DIR="${WORK_DIR}/results/04_mapping"
RESULTS_DIR="${WORK_DIR}/results/05_binning"
THREADS=8
SAMPLES=(53394 53395 53396 53397 53398 53399)

# Initialize conda
eval "$(conda shell.bash hook)"

STEP=${1:-all}

echo "METAGENOMIC BINNING PIPELINE"
echo "Step: $STEP"
echo "Results: $RESULTS_DIR"

mkdir -p ${RESULTS_DIR}/{metabat2,maxbin2,dastool,checkm2,final_bins}

# Create sample-specific directories
for sample in "${SAMPLES[@]}"; do
    mkdir -p ${RESULTS_DIR}/metabat2/${sample}
    mkdir -p ${RESULTS_DIR}/maxbin2/${sample}
done

prep_data() {
    echo "STEP 1: Preparing input data for individual samples..."
    
    # Check that mapping results exist
    for sample in "${SAMPLES[@]}"; do
        bam_file="${MAPPING_DIR}/${sample}_sorted.bam"
        if [ ! -f "$bam_file" ]; then
            echo "ERROR: BAM file not found for $sample: $bam_file"
            echo "Run read mapping (04_mapping.sh) first!"
            exit 1
        fi
    done
    
    echo "All BAM files found, proceeding with binning..."
}

run_metabat2() {
    echo "STEP 2.1: Running MetaBAT2 on individual samples..."
    conda activate metabat218
    
    total_bins=0
    
    for sample in "${SAMPLES[@]}"; do
        echo "Processing sample $sample with MetaBAT2..."
        
        assembly_file="${ASSEMBLY_DIR}/${sample}/contigs.fasta"
        bam_file="${MAPPING_DIR}/${sample}_sorted.bam"
        depth_file="${RESULTS_DIR}/${sample}_depth.txt"
        sample_output_dir="${RESULTS_DIR}/metabat2/${sample}"
        
        mkdir -p "$sample_output_dir"
        
        if [ $(ls ${sample_output_dir}/bin*.fa 2>/dev/null | wc -l) -gt 0 ]; then
            echo "MetaBAT2 already completed for $sample, skipping..."
            continue
        fi
        
        # Generate depth file for this sample
        if [ ! -f "$depth_file" ]; then
            jgi_summarize_bam_contig_depths --outputDepth "$depth_file" "$bam_file"
        fi
        
        # Run MetaBAT2
        metabat2 -i "$assembly_file" \
                 -a "$depth_file" \
                 -o "${sample_output_dir}/bin" \
                 --minContig 1500 \
                 --numThreads ${THREADS}
        
        sample_bins=$(ls ${sample_output_dir}/bin*.fa 2>/dev/null | wc -l)
        total_bins=$((total_bins + sample_bins))
        echo "MetaBAT2 completed for $sample: $sample_bins bins"
    done
    
    echo "MetaBAT2 completed for all samples: $total_bins total bins"
}

run_maxbin2() {
    echo "STEP 2.2: Running MaxBin2 on individual samples..."
    conda activate maxbin2_env
    
    if ! command -v run_MaxBin.pl &> /dev/null; then
        echo "MaxBin2 not found, skipping..."
        return
    fi
    
    # Ensure HMMER3 tools are in PATH
    export PATH="${CONDA_PREFIX}/bin:$PATH"
    
    echo "Testing HMMER3 availability:"
    which hmmsearch || echo "hmmsearch not found in PATH"
    hmmsearch -h 2>/dev/null | head -2 || echo "hmmsearch test failed"
    
    # Try without setting file first (relying on PATH)
    total_bins=0
    
    for sample in "${SAMPLES[@]}"; do
        echo "Processing sample $sample with MaxBin2..."
        
        assembly_file="${ASSEMBLY_DIR}/${sample}/contigs.fasta"
        depth_file="${RESULTS_DIR}/${sample}_depth.txt"
        sample_output_dir="${RESULTS_DIR}/maxbin2/${sample}"
        
        mkdir -p "$sample_output_dir"
        
        if [ $(ls ${sample_output_dir}/bin*.fasta 2>/dev/null | wc -l) -gt 0 ]; then
            echo "MaxBin2 already completed for $sample, skipping..."
            continue
        fi
        
        # Check if depth file exists (created by MetaBAT2)
        if [ ! -f "$depth_file" ]; then
            echo "Depth file not found for $sample, skipping MaxBin2..."
            continue
        fi
        
        # Create abundance file for MaxBin2 (using column 4 from depth file)
        abund_file="${sample_output_dir}/abundance.txt"
        cut -f1,4 "$depth_file" > "$abund_file"
        
        echo "Running MaxBin2 for sample $sample..."
        echo "Assembly: $assembly_file"
        echo "Abundance: $abund_file"
        echo "Output: ${sample_output_dir}/bin"
        
        # Try MaxBin2 with verbose output to see what's happening
        if run_MaxBin.pl -contig "$assembly_file" \
                         -abund "$abund_file" \
                         -out "${sample_output_dir}/bin" \
                         -thread ${THREADS} \
                         -min_contig_length 500 \
                         -verbose 2>&1 | grep -v "Cannot run HMMER3"; then
            echo "MaxBin2 completed successfully for $sample"
        else
            echo "MaxBin2 failed for $sample, trying with custom setting..."
            
            # Create setting file in the sample directory
            cat > "${sample_output_dir}/setting" << 'EOF'
[FragGeneScan] /home/jiayi-chen/anaconda3/envs/maxbin2_env/bin
[Bowtie2] /home/jiayi-chen/anaconda3/envs/maxbin2_env/bin  
[HMMER3] /home/jiayi-chen/anaconda3/envs/maxbin2_env/bin
EOF
            
            # Run from sample directory where setting file is located
            cd "${sample_output_dir}"
            run_MaxBin.pl -contig "$assembly_file" \
                         -abund "$abund_file" \
                         -out "bin" \
                         -thread ${THREADS} \
                         -min_contig_length 500 \
                         -verbose
            cd - > /dev/null
        fi
        
        sample_bins=$(ls ${sample_output_dir}/bin*.fasta 2>/dev/null | wc -l)
        total_bins=$((total_bins + sample_bins))
        echo "MaxBin2 completed for $sample: $sample_bins bins"
    done
    
    echo "MaxBin2 completed for all samples: $total_bins total bins"
}


run_dastool() {
    echo "STEP 3: Running DAS Tool optimization..."
    conda activate dastool_env
    
    if [ $(ls ${RESULTS_DIR}/dastool/optimized_DASTool_bins/*.fa 2>/dev/null | wc -l) -gt 0 ]; then
        echo "DAS Tool already completed, skipping..."
        return
    fi
    
    if ! command -v DAS_Tool &> /dev/null; then
        echo "DAS Tool not found, skipping..."
        return
    fi
    
    DASTOOL_SCRIPTS=$(dirname $(which DAS_Tool))
    
    # DAS Tool requires a single reference assembly, but we have individual assemblies
    # Let's create a concatenated assembly for DAS Tool if it doesn't exist
    CONCAT_FASTA="${RESULTS_DIR}/concatenated_contigs.fasta"
    if [ ! -f "$CONCAT_FASTA" ]; then
        echo "Creating concatenated assembly for DAS Tool..."
        > "$CONCAT_FASTA"  # Create empty file
        
        for sample in "${SAMPLES[@]}"; do
            assembly_file="${ASSEMBLY_DIR}/${sample}/contigs.fasta"
            if [ -f "$assembly_file" ]; then
                echo "Adding $sample contigs to concatenated assembly..."
                # Add sample prefix to contig names to avoid conflicts
                sed "s/^>/>${sample}_/" "$assembly_file" >> "$CONCAT_FASTA"
            fi
        done
    fi
    
    # Prepare bin-to-contig tables
    BIN_SETS=""
    BIN_LABELS=""
    
    # Collect all MetaBAT2 bins from all samples
    if [ $(find ${RESULTS_DIR}/metabat2 -name "bin*.fa" 2>/dev/null | wc -l) -gt 0 ]; then
        echo "Processing MetaBAT2 bins for DAS Tool..."
        # Create a temporary directory for all MetaBAT2 bins
        mkdir -p ${RESULTS_DIR}/dastool_temp/metabat2_all
        for sample in "${SAMPLES[@]}"; do
            if [ -d "${RESULTS_DIR}/metabat2/${sample}" ]; then
                for bin in "${RESULTS_DIR}/metabat2/${sample}"/bin*.fa; do
                    if [ -f "$bin" ]; then
                        # Copy with sample prefix to avoid name conflicts
                        cp "$bin" "${RESULTS_DIR}/dastool_temp/metabat2_all/${sample}_$(basename "$bin")"
                    fi
                done
            fi
        done
        
        # Create contig2bin mapping for MetaBAT2 (only 2 columns: contig, bin)
        > ${RESULTS_DIR}/metabat2_contigs2bins.tsv  # Create empty file
        for bin in ${RESULTS_DIR}/dastool_temp/metabat2_all/*.fa; do
            if [ -f "$bin" ]; then
                bin_name=$(basename "$bin" .fa)
                # Extract sample prefix from bin name (e.g., 53394 from 53394_bin.1)
                sample_prefix=$(echo "$bin_name" | cut -d'_' -f1)
                # Extract contig names and map to bin (only first field, no depth info)
                grep "^>" "$bin" | sed 's/^>//' | cut -f1 > /tmp/contigs_$$.txt
                while read contig; do
                    # Add sample prefix to contig name to match concatenated assembly
                    echo -e "${sample_prefix}_${contig}\t${bin_name}" >> ${RESULTS_DIR}/metabat2_contigs2bins.tsv
                done < /tmp/contigs_$$.txt
                rm -f /tmp/contigs_$$.txt
            fi
        done
        
        BIN_SETS="${BIN_SETS},${RESULTS_DIR}/metabat2_contigs2bins.tsv"
        BIN_LABELS="${BIN_LABELS},MetaBAT2"
    fi
    
    # Collect all MaxBin2 bins from all samples
    if [ $(find ${RESULTS_DIR}/maxbin2 -name "bin*.fasta" 2>/dev/null | wc -l) -gt 0 ]; then
        echo "Processing MaxBin2 bins for DAS Tool..."
        # Create a temporary directory for all MaxBin2 bins
        mkdir -p ${RESULTS_DIR}/dastool_temp/maxbin2_all
        for sample in "${SAMPLES[@]}"; do
            if [ -d "${RESULTS_DIR}/maxbin2/${sample}" ]; then
                for bin in "${RESULTS_DIR}/maxbin2/${sample}"/bin*.fasta; do
                    if [ -f "$bin" ]; then
                        # Copy with sample prefix to avoid name conflicts
                        cp "$bin" "${RESULTS_DIR}/dastool_temp/maxbin2_all/${sample}_$(basename "$bin")"
                    fi
                done
            fi
        done
        
        # Create contig2bin mapping for MaxBin2 (only 2 columns: contig, bin)
        > ${RESULTS_DIR}/maxbin2_contigs2bins.tsv  # Create empty file
        for bin in ${RESULTS_DIR}/dastool_temp/maxbin2_all/*.fasta; do
            if [ -f "$bin" ]; then
                bin_name=$(basename "$bin" .fasta)
                # Extract sample prefix from bin name (e.g., 53394 from 53394_bin.001.fasta)
                sample_prefix=$(echo "$bin_name" | cut -d'_' -f1)
                # Extract contig names and map to bin
                grep "^>" "$bin" | sed 's/^>//' > /tmp/contigs_$$.txt
                while read contig; do
                    # Add sample prefix to contig name to match concatenated assembly
                    echo -e "${sample_prefix}_${contig}\t${bin_name}" >> ${RESULTS_DIR}/maxbin2_contigs2bins.tsv
                done < /tmp/contigs_$$.txt
                rm -f /tmp/contigs_$$.txt
            fi
        done
        
        BIN_SETS="${BIN_SETS},${RESULTS_DIR}/maxbin2_contigs2bins.tsv"
        BIN_LABELS="${BIN_LABELS},MaxBin2"
    fi
    
    # Remove leading commas
    BIN_SETS=${BIN_SETS#,}
    BIN_LABELS=${BIN_LABELS#,}
    
    if [ -n "$BIN_SETS" ]; then
        echo "Running DAS Tool with bins from: $BIN_LABELS"
        DAS_Tool -i "$BIN_SETS" \
                 -l "$BIN_LABELS" \
                 -c "$CONCAT_FASTA" \
                 -o "${RESULTS_DIR}/dastool/optimized" \
                 --threads ${THREADS} \
                 --write_bins
        echo "DAS Tool optimization completed"
        
        # Clean up temporary directories
        rm -rf ${RESULTS_DIR}/dastool_temp
    else
        echo "No bins found for DAS Tool optimization"
    fi
}

run_checkm2() {
    echo "STEP 4: Running CheckM2 quality assessment..."
    conda activate checkm2_env
    
    if [ -f "${RESULTS_DIR}/checkm2/quality_report.tsv" ]; then
        echo "CheckM2 already completed, skipping..."
        return
    fi
    
    if ! command -v checkm2 &> /dev/null; then
        echo "CheckM2 not found, skipping..."
        return
    fi
    
    # Use DAS Tool optimized bins directly
    INPUT_DIR="${RESULTS_DIR}/checkm2/input"
    OUTPUT_DIR="${RESULTS_DIR}/checkm2"
    
    # Ensure input directory exists and copy DAS Tool bins
    mkdir -p "$INPUT_DIR"
    
    if [ -d "${RESULTS_DIR}/dastool/optimized_DASTool_bins" ]; then
        # Copy DAS Tool bins to input directory
        cp "${RESULTS_DIR}/dastool/optimized_DASTool_bins"/*.fa "$INPUT_DIR/" 2>/dev/null || true
        
        bin_count=$(ls "$INPUT_DIR"/*.fa 2>/dev/null | wc -l)
        echo "Found $bin_count DAS Tool optimized bins for CheckM2 analysis"
        
        if [ $bin_count -eq 0 ]; then
            echo "No bins found for CheckM2 analysis"
            return
        fi
        
        # Run CheckM2
        echo "Running CheckM2 on optimized bins..."
        if checkm2 predict \
            --input "$INPUT_DIR" \
            --output-directory "$OUTPUT_DIR" \
            --threads "$THREADS" \
            --force \
            -x .fa; then
            
            echo "CheckM2 completed successfully"
            
            if [ -f "${OUTPUT_DIR}/quality_report.tsv" ]; then
                # Quality summary
                total=$(tail -n +2 "${OUTPUT_DIR}/quality_report.tsv" | wc -l)
                hq=$(tail -n +2 "${OUTPUT_DIR}/quality_report.tsv" | awk -F'\t' '$2 >= 90 && $3 <= 5' | wc -l)
                mq=$(tail -n +2 "${OUTPUT_DIR}/quality_report.tsv" | awk -F'\t' '$2 >= 50 && $3 <= 10' | wc -l)
                
                echo "Quality Summary:"
                echo "  Total bins: $total"
                echo "  High quality (>=90% complete, <=5% contamination): $hq"
                echo "  Medium quality (>=50% complete, <=10% contamination): $mq"
            fi
        else
            echo "CheckM2 failed"
        fi
    else
        echo "No DAS Tool optimized bins found. Run DAS Tool first."
    fi
}

extract_hq_bins() {
    echo "STEP 5: Extracting high-quality bins..."
    
    mkdir -p ${RESULTS_DIR}/high_quality_bins
    
    CHECKM_RESULTS="${RESULTS_DIR}/checkm2/quality_report.tsv"
    
    if [ -f "$CHECKM_RESULTS" ]; then
        # Extract bins with completeness >= 50% and contamination <= 10%
        awk 'NR>1 && $2>=50 && $3<=10 {print $1}' "$CHECKM_RESULTS" > ${RESULTS_DIR}/hq_bin_list.txt
        
        HQ_COUNT=0
        SOURCE_DIRS=(
            "${RESULTS_DIR}/checkm2/clean_bins"
            "${RESULTS_DIR}/metabat2"
            "${RESULTS_DIR}/maxbin2"
            "${RESULTS_DIR}/dastool/optimized_DASTool_bins"
        )
        
        while read -r bin_name; do
            copied=false
            for source_dir in "${SOURCE_DIRS[@]}"; do
                if [ -d "$source_dir" ]; then
                    for ext in ".fa" ".fasta" ".fna"; do
                        source_file="${source_dir}/${bin_name}${ext}"
                        if [ -f "$source_file" ]; then
                            cp "$source_file" "${RESULTS_DIR}/high_quality_bins/${bin_name}.fa"
                            HQ_COUNT=$((HQ_COUNT + 1))
                            copied=true
                            break 2
                        fi
                    done
                fi
            done
        done < ${RESULTS_DIR}/hq_bin_list.txt
        
        echo "Extracted $HQ_COUNT high-quality bins"
    else
        echo "CheckM results not found, extracting large bins (>1MB)..."
        
        HQ_COUNT=0
        for source_dir in "${RESULTS_DIR}"/metabat2 "${RESULTS_DIR}"/maxbin2; do
            if [ -d "$source_dir" ]; then
                for ext in "*.fa" "*.fasta"; do
                    for bin_file in "$source_dir"/$ext; do
                        if [ -f "$bin_file" ] && [ $(stat -c%s "$bin_file" 2>/dev/null) -gt 1048576 ]; then
                            bin_name=$(basename "$bin_file")
                            cp "$bin_file" "${RESULTS_DIR}/high_quality_bins/${bin_name%.*}.fa"
                            HQ_COUNT=$((HQ_COUNT + 1))
                        fi
                    done
                done
            fi
        done
        
        echo "Extracted $HQ_COUNT large bins as potential high-quality"
    fi
}

# Main execution
case $STEP in
    prep) prep_data ;;
    metabat) prep_data; run_metabat2 ;;
    maxbin) prep_data; run_maxbin2 ;;
    vamb) prep_data; run_vamb ;;
    dastool) run_dastool ;;
    checkm2) run_checkm2 ;;
    extract) extract_hq_bins ;;
    all)
        prep_data
        run_metabat2
        run_maxbin2
        run_vamb
        run_dastool
        run_checkm2
        extract_hq_bins
        ;;
    *)
        echo "Usage: $0 {prep|metabat|maxbin|vamb|dastool|checkm2|extract|all}"
        exit 1
        ;;
esac

echo "STEP '$STEP' COMPLETED!"
echo "Results: ${RESULTS_DIR}"
