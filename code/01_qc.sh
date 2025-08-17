#!/bin/bash

# Fastp Auto Primer Removal - Concise Version
set -e

# Configuration
DATA_DIR="/home/jiayi-chen/Documents/farting_trees/data"
OUTPUT_DIR="/home/jiayi-chen/Documents/farting_trees/results/01_qc/primer_removed"
THREADS=8

# Setup
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_qc
mkdir -p "$OUTPUT_DIR"

echo "Auto-detecting and processing paired FASTQ files..."

# Find and process sample pairs
declare -A sample_pairs
for r1_file in "$DATA_DIR"/*_R1_*.fastq.gz; do
    [[ -f "$r1_file" ]] || continue
    r2_file="${r1_file/_R1_/_R2_}"
    [[ -f "$r2_file" ]] || continue
    
    sample_name=$(basename "$r1_file" | sed 's/_R1_.*$//')
    sample_pairs["$sample_name"]="$r1_file,$r2_file"
    echo "Found: $sample_name"
done

echo "Processing ${#sample_pairs[@]} samples..."

# Process each sample
for sample_name in "${!sample_pairs[@]}"; do
    IFS=',' read -r r1_file r2_file <<< "${sample_pairs[$sample_name]}"
    
    out_r1="$OUTPUT_DIR/${sample_name}_R1_no_primers.fastq.gz"
    out_r2="$OUTPUT_DIR/${sample_name}_R2_no_primers.fastq.gz"
    
    # Skip if already processed
    [[ -f "$out_r1" && -f "$out_r2" ]] && { echo "✓ $sample_name (skipped)"; continue; }
    
    echo "Processing $sample_name..."
    
    # Run fastp with auto-detection
    fastp \
        --in1 "$r1_file" --in2 "$r2_file" \
        --out1 "$out_r1" --out2 "$out_r2" \
        --detect_adapter_for_pe \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread $THREADS \
        --html "$OUTPUT_DIR/${sample_name}.html" \
        --json "$OUTPUT_DIR/${sample_name}.json" \
        --report_title "$sample_name" \
        --quiet
    
    echo "✓ $sample_name completed"
done

echo "Done! Output in: $OUTPUT_DIR"

# Run FastQC on primer-removed files
echo "Running FastQC on primer-removed files..."
mkdir -p "$OUTPUT_DIR/fastqc"

fastqc "$OUTPUT_DIR"/*_no_primers.fastq.gz \
    -o "$OUTPUT_DIR/fastqc" \
    -t $THREADS \
    --quiet

echo "✓ FastQC completed"

# Generate MultiQC report
echo "Generating MultiQC report..."
cd "$OUTPUT_DIR"
multiqc . \
    -o . \
    -n primer_removal_multiqc_report \
    --title "Primer Removal QC Report" \
    --comment "Quality control after primer removal with fastp auto-detection" \
    --force \
    --quiet

echo "✓ MultiQC report generated: $OUTPUT_DIR/primer_removal_multiqc_report.html"
echo "All done! Check the MultiQC report for quality assessment."
    
