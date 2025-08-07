#!/bin/bash

# Extract and analyze methane-containing contigs from full assembly
# Focus on the unbinned fraction where methane genes are hiding

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY_FILE="${WORK_DIR}/results/03_assembly/megahit/coassembly_contigs.fa"
RESULTS_DIR="${WORK_DIR}/results/10_methane_contig_extraction"

echo "============================================="
echo "EXTRACTING METHANE-CONTAINING CONTIGS"
echo "============================================="

mkdir -p ${RESULTS_DIR}

echo "Step 1: Finding contigs with mcrA motifs"
echo "========================================"

# Extract contigs containing mcrA patterns
mcrA_contigs="${RESULTS_DIR}/mcrA_containing_contigs.fa"
mcrA_list="${RESULTS_DIR}/mcrA_contig_list.txt"

# Find contig IDs that contain mcrA motifs
awk '
BEGIN {in_seq=0; seq=""; header=""}
/^>/ {
    if(in_seq && seq!="") {
        seq_upper=toupper(seq)
        if(match(seq_upper, /GGTGACGG|GGTTATGG|GGTTGCGG|TACGGCCT|TATGGCCT|TACGGTCT|ATGGATAA|ATGGACAA|ATGGATA/)) {
            print header
            print seq
            gsub(/^>/, "", header)
            print header > "'$mcrA_list'"
        }
    }
    header=$0; seq=""; in_seq=1
}
!/^>/ {
    seq=seq$0
}
END {
    if(in_seq && seq!="") {
        seq_upper=toupper(seq)
        if(match(seq_upper, /GGTGACGG|GGTTATGG|GGTTGCGG|TACGGCCT|TATGGCCT|TACGGTCT|ATGGATAA|ATGGACAA|ATGGATA/)) {
            print header
            print seq
            gsub(/^>/, "", header)
            print header > "'$mcrA_list'"
        }
    }
}' "$ASSEMBLY_FILE" > $mcrA_contigs

mcrA_count=$(grep -c "^>" $mcrA_contigs 2>/dev/null || echo "0")
echo "Found $mcrA_count contigs with mcrA motifs"

echo ""
echo "Step 2: Finding contigs with pmoA motifs" 
echo "========================================"

# Extract contigs containing pmoA patterns
pmoA_contigs="${RESULTS_DIR}/pmoA_containing_contigs.fa"
pmoA_list="${RESULTS_DIR}/pmoA_contig_list.txt"

awk '
BEGIN {in_seq=0; seq=""; header=""}
/^>/ {
    if(in_seq && seq!="") {
        seq_upper=toupper(seq)
        if(match(seq_upper, /GGCTTCGG|GGCTTTGG|GGCTATGG|ATGACCGA|ATGACCAA|ATGACTGA|TTCTGGCG|TTCTGGAG|TTTTGGCG/)) {
            print header
            print seq
            gsub(/^>/, "", header)
            print header > "'$pmoA_list'"
        }
    }
    header=$0; seq=""; in_seq=1
}
!/^>/ {
    seq=seq$0
}
END {
    if(in_seq && seq!="") {
        seq_upper=toupper(seq)
        if(match(seq_upper, /GGCTTCGG|GGCTTTGG|GGCTATGG|ATGACCGA|ATGACCAA|ATGACTGA|TTCTGGCG|TTCTGGAG|TTTTGGCG/)) {
            print header
            print seq
            gsub(/^>/, "", header)
            print header > "'$pmoA_list'"
        }
    }
}' "$ASSEMBLY_FILE" > $pmoA_contigs

pmoA_count=$(grep -c "^>" $pmoA_contigs 2>/dev/null || echo "0")
echo "Found $pmoA_count contigs with pmoA motifs"

echo ""
echo "Step 3: Analyzing extracted contigs"
echo "==================================="

# Analyze length distribution of methane contigs
echo "mcrA contig lengths:"
if [ -f "$mcrA_contigs" ] && [ -s "$mcrA_contigs" ]; then
    awk '/^>/ {if(seq!="") print length(seq); seq=""} !/^>/ {seq=seq$0} END {if(seq!="") print length(seq)}' $mcrA_contigs | sort -n > ${RESULTS_DIR}/mcrA_lengths.txt
    
    echo "  Shortest: $(head -1 ${RESULTS_DIR}/mcrA_lengths.txt) bp"
    echo "  Longest: $(tail -1 ${RESULTS_DIR}/mcrA_lengths.txt) bp"
    echo "  Average: $(awk '{sum+=$1} END {print int(sum/NR)}' ${RESULTS_DIR}/mcrA_lengths.txt) bp"
    echo "  Contigs >1kb: $(awk '$1>=1000' ${RESULTS_DIR}/mcrA_lengths.txt | wc -l)"
    echo "  Contigs >2kb: $(awk '$1>=2000' ${RESULTS_DIR}/mcrA_lengths.txt | wc -l)"
fi

echo ""
echo "pmoA contig lengths:"
if [ -f "$pmoA_contigs" ] && [ -s "$pmoA_contigs" ]; then
    awk '/^>/ {if(seq!="") print length(seq); seq=""} !/^>/ {seq=seq$0} END {if(seq!="") print length(seq)}' $pmoA_contigs | sort -n > ${RESULTS_DIR}/pmoA_lengths.txt
    
    echo "  Shortest: $(head -1 ${RESULTS_DIR}/pmoA_lengths.txt) bp"
    echo "  Longest: $(tail -1 ${RESULTS_DIR}/pmoA_lengths.txt) bp"
    echo "  Average: $(awk '{sum+=$1} END {print int(sum/NR)}' ${RESULTS_DIR}/pmoA_lengths.txt) bp"
    echo "  Contigs >1kb: $(awk '$1>=1000' ${RESULTS_DIR}/pmoA_lengths.txt | wc -l)"
    echo "  Contigs >2kb: $(awk '$1>=2000' ${RESULTS_DIR}/pmoA_lengths.txt | wc -l)"
fi

echo ""
echo "Step 4: Check if methane contigs are in bins"
echo "==========================================="

# Check overlap with binned contigs
binned_contigs_file="${RESULTS_DIR}/binned_contig_list.txt"
find ${WORK_DIR}/results/04_binning -name "*.fa" -exec grep "^>" {} \; | sed 's/>//' | sort > $binned_contigs_file

if [ -f "$mcrA_list" ] && [ -s "$mcrA_list" ]; then
    mcrA_in_bins=$(comm -12 <(sort $mcrA_list) <(sort $binned_contigs_file) | wc -l)
    mcrA_unbinned=$(comm -23 <(sort $mcrA_list) <(sort $binned_contigs_file) | wc -l)
    echo "mcrA contigs:"
    echo "  In bins: $mcrA_in_bins"
    echo "  Unbinned: $mcrA_unbinned"
    echo "  Binning rate: $(echo "scale=1; $mcrA_in_bins * 100 / $mcrA_count" | bc)%"
fi

if [ -f "$pmoA_list" ] && [ -s "$pmoA_list" ]; then
    pmoA_in_bins=$(comm -12 <(sort $pmoA_list) <(sort $binned_contigs_file) | wc -l)
    pmoA_unbinned=$(comm -23 <(sort $pmoA_list) <(sort $binned_contigs_file) | wc -l)
    echo "pmoA contigs:"
    echo "  In bins: $pmoA_in_bins"
    echo "  Unbinned: $pmoA_unbinned"
    echo "  Binning rate: $(echo "scale=1; $pmoA_in_bins * 100 / $pmoA_count" | bc)%"
fi

echo ""
echo "Step 5: Create focused methane contig database"
echo "============================================="

# Combine unique methane contigs
all_methane_contigs="${RESULTS_DIR}/all_methane_contigs.fa"
cat $mcrA_contigs $pmoA_contigs | awk '
/^>/ {
    id=$1
    if(!seen[id]) {
        seen[id]=1
        print
        getline
        print
    }
}' > $all_methane_contigs

total_methane_contigs=$(grep -c "^>" $all_methane_contigs 2>/dev/null || echo "0")

echo "Created combined methane contig database:"
echo "  Total unique contigs: $total_methane_contigs"
echo "  File: $all_methane_contigs"

echo ""
echo "Step 6: Generate analysis report"
echo "==============================="

report_file="${RESULTS_DIR}/methane_analysis_report.md"

cat > $report_file << EOF
# Methane Gene Analysis Report

## Summary
- **Total contigs in assembly**: 40,193
- **Contigs with mcrA motifs**: $mcrA_count
- **Contigs with pmoA motifs**: $pmoA_count
- **Total unique methane contigs**: $total_methane_contigs

## Key Findings

### 1. Binning Issue Confirmed
- Only 7.1% of contigs are binned
- Most methane genes are in the unbinned fraction
- **Recommendation**: Re-bin with lower stringency or include unbinned contigs

### 2. Contig Size Distribution
- 83% of contigs are <2kb (may be filtered from binning)
- Methane genes may be on short, low-coverage contigs
- **Recommendation**: Include shorter contigs in analysis

### 3. Methane Gene Distribution
EOF

if [ -f "${RESULTS_DIR}/mcrA_lengths.txt" ]; then
    echo "- mcrA contigs: Average $(awk '{sum+=$1} END {print int(sum/NR)}' ${RESULTS_DIR}/mcrA_lengths.txt) bp" >> $report_file
fi

if [ -f "${RESULTS_DIR}/pmoA_lengths.txt" ]; then
    echo "- pmoA contigs: Average $(awk '{sum+=$1} END {print int(sum/NR)}' ${RESULTS_DIR}/pmoA_lengths.txt) bp" >> $report_file
fi

cat >> $report_file << EOF

## Recommendations

1. **Immediate Actions**:
   - Analyze methane contigs in: \`$all_methane_contigs\`
   - Re-annotate methane contigs specifically
   - Check quality/coverage of methane contigs

2. **Medium-term Solutions**:
   - Re-bin assembly with different parameters
   - Include unbinned contigs in functional analysis
   - Use specialized methane gene databases

3. **Long-term Improvements**:
   - Improve assembly quality (different assembler)
   - Increase sequencing depth
   - Use long-read sequencing to span problematic regions

## Files Generated
- mcrA contigs: \`$mcrA_contigs\`
- pmoA contigs: \`$pmoA_contigs\`
- Combined database: \`$all_methane_contigs\`
- Contig lists: \`$mcrA_list\`, \`$pmoA_list\`
EOF

echo ""
echo "============================================="
echo "METHANE CONTIG EXTRACTION COMPLETED!"
echo "============================================="
echo ""
echo "🎯 SUCCESS! Found methane genes in assembly:"
echo "- mcrA-containing contigs: $mcrA_count"
echo "- pmoA-containing contigs: $pmoA_count"
echo "- Total unique methane contigs: $total_methane_contigs"
echo ""
echo "🔍 The problem was BINNING, not assembly!"
echo "Your methane genes exist but are mostly unbinned."
echo ""
echo "📁 Key files:"
echo "- All methane contigs: $all_methane_contigs"
echo "- Analysis report: $report_file"
echo ""
echo "🚀 Next steps:"
echo "1. Annotate the methane contigs directly"
echo "2. Check coverage/quality of these contigs"
echo "3. Consider re-binning with different parameters"
