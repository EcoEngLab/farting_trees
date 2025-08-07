#!/bin/bash

# Complete the coverage calculation using existing SAM files
# This script finishes what calculate_coverage.sh started

set -e

RESULTS_DIR="/home/jiayi-chen/Documents/farting_trees/results/03_assembly"
temp_dir="${RESULTS_DIR}/mapping/temp_sams"
temp_sam="${RESULTS_DIR}/mapping/temp_alignment.sam"
THREADS=8

echo "=========================================="
echo "COMPLETING COVERAGE CALCULATION"
echo "=========================================="

# Check if individual SAM files exist
if [ ! -d "$temp_dir" ]; then
    echo "✗ Individual SAM files not found. Please run the full script."
    exit 1
fi

# Count existing SAM files
sam_count=$(ls "$temp_dir"/*.sam 2>/dev/null | wc -l)
echo "Found $sam_count individual SAM files to merge"

if [ $sam_count -eq 0 ]; then
    echo "✗ No SAM files found in $temp_dir"
    exit 1
fi

# Activate assembly environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_assembly_new

echo ""
echo "Step 1: Merging SAM files..."
sam_files=("$temp_dir"/*.sam)

# Take header from first file
echo "Extracting header from first SAM file..."
head -n 100 "${sam_files[0]}" | grep "^@" > "$temp_sam"

# Add all alignments (non-header lines)
echo "Merging alignment records..."
for sam_file in "${sam_files[@]}"; do
    sample_name=$(basename "$sam_file" .sam)
    echo "  Adding alignments from $sample_name..."
    grep -v "^@" "$sam_file" >> "$temp_sam"
done

echo "✓ SAM files merged"
echo "Merged file size: $(ls -lh "$temp_sam" | awk '{print $5}')"

echo ""
echo "Step 2: Converting to sorted BAM..."
# Convert to BAM and sort
samtools view -@ ${THREADS} -bS -F 4 "$temp_sam" | \
samtools sort -@ ${THREADS} -o "${RESULTS_DIR}/mapping/coassembly_sorted.bam"

echo "✓ BAM file created"
echo "BAM file size: $(ls -lh "${RESULTS_DIR}/mapping/coassembly_sorted.bam" | awk '{print $5}')"

echo ""
echo "Step 3: Indexing BAM file..."
samtools index "${RESULTS_DIR}/mapping/coassembly_sorted.bam"
echo "✓ BAM indexing completed"

echo ""
echo "Step 4: Calculating coverage statistics..."
bam_file="${RESULTS_DIR}/mapping/coassembly_sorted.bam"

# Basic mapping stats
total_mapped=$(samtools view -c -F 4 "$bam_file")
total_reads_mapped=$((total_mapped / 2))  # Paired reads

echo "Mapping summary:"
echo "  Total mapped reads: $total_reads_mapped"
echo "  Total alignments: $total_mapped"

# Coverage depth calculation
echo ""
echo "Calculating coverage depth (this may take a few minutes)..."
samtools depth "$bam_file" > "${RESULTS_DIR}/mapping/coverage_depth.txt"

# Basic coverage stats
if [ -f "${RESULTS_DIR}/mapping/coverage_depth.txt" ]; then
    avg_coverage=$(awk '{sum+=$3; count++} END {if(count>0) print sum/count; else print 0}' "${RESULTS_DIR}/mapping/coverage_depth.txt")
    max_coverage=$(awk 'BEGIN{max=0} {if($3>max) max=$3} END{print max}' "${RESULTS_DIR}/mapping/coverage_depth.txt")
    
    echo "✓ Coverage statistics:"
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

echo ""
echo "Step 5: Cleaning up temporary files..."
# Clean up temporary files
rm -f "$temp_sam"
rm -rf "$temp_dir"
echo "✓ Temporary files cleaned up"

echo ""
echo "=========================================="
echo "COVERAGE CALCULATION COMPLETED!"
echo "=========================================="
echo ""
echo "✓ Final output files:"
echo "  Sorted BAM: ${RESULTS_DIR}/mapping/coassembly_sorted.bam"
echo "  BAM index: ${RESULTS_DIR}/mapping/coassembly_sorted.bam.bai"
echo "  Coverage depth: ${RESULTS_DIR}/mapping/coverage_depth.txt"
echo ""
echo "Next steps:"
echo "  - Use coverage data for genome binning"
echo "  - Proceed with: ./code/06_binning.sh"
echo ""
