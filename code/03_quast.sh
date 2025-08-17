#!/bin/bash
# QUAST Assembly Quality Assessment for farting_trees project
set -e

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"
conda activate quast

command -v quast.py &> /dev/null || { echo "❌ QUAST not found. Install with: conda install -c bioconda quast"; exit 1; }
 
# Configuration
WORK_DIR="/home/jiayi-chen/Documents/farting_trees"
ASSEMBLY_DIR="${WORK_DIR}/results/02_assembly"
OUTPUT_DIR="${WORK_DIR}/results/03_quast"
THREADS=$(nproc)
MIN_CONTIG_LENGTH=500

echo "Running QUAST analysis..."
mkdir -p "$OUTPUT_DIR"

# Find assembly files
assembly_files=()
sample_names=()
for dir in "$ASSEMBLY_DIR"/*/; do
    contigs_file="${dir}contigs.fasta"
    [[ -f "$contigs_file" ]] && {
        assembly_files+=("$contigs_file")
        sample_names+=($(basename "$dir"))
    }
done

[[ ${#assembly_files[@]} -eq 0 ]] && { echo "❌ No contigs.fasta files found"; exit 1; }

# Run QUAST on individual assemblies
for i in "${!assembly_files[@]}"; do
    sample_name="${sample_names[$i]}"
    sample_output_dir="$OUTPUT_DIR/${sample_name}"
    
    [[ -f "$sample_output_dir/report.html" ]] && continue
    
    quast.py "${assembly_files[$i]}" -o "$sample_output_dir" --threads $THREADS \
        --min-contig $MIN_CONTIG_LENGTH --labels "$sample_name" --no-plots --silent
done

# Run comparative analysis
comparative_output="$OUTPUT_DIR/comparative_analysis"
if [[ ! -f "$comparative_output/report.html" ]]; then
    labels=$(IFS=','; echo "${sample_names[*]}")
    quast.py "${assembly_files[@]}" -o "$comparative_output" --threads $THREADS \
        --min-contig $MIN_CONTIG_LENGTH --labels "$labels" --silent
fi

# Generate summary
summary_file="$OUTPUT_DIR/assembly_summary.txt"
{
    echo "ASSEMBLY QUALITY SUMMARY - $(date)"
    echo "Total assemblies: ${#assembly_files[@]}"
    echo ""
    for sample_name in "${sample_names[@]}"; do
        report_file="$OUTPUT_DIR/${sample_name}/report.txt"
        if [[ -f "$report_file" ]]; then
            echo "=== $sample_name ==="
            grep -E "(# contigs|Total length|N50|GC %)" "$report_file" | head -4
            echo ""
        fi
    done
} > "$summary_file"

echo "✓ Analysis complete!"
echo "Summary: $summary_file"
echo "Comparative report: $comparative_output/report.html"
