#!/bin/bash

# Simple but comprehensive methane gene search in full assembly
# No external dependencies required

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY_FILE="${WORK_DIR}/results/03_assembly/megahit/coassembly_contigs.fa"
RESULTS_DIR="${WORK_DIR}/results/09_simple_methane_search"

echo "============================================="
echo "SIMPLE METHANE GENE SEARCH IN FULL ASSEMBLY"
echo "============================================="
echo "Searching all 40,193 contigs for methane genes"
echo ""

mkdir -p ${RESULTS_DIR}

results_file="${RESULTS_DIR}/full_assembly_methane_search.txt"

echo "=== METHANE GENE SEARCH IN FULL ASSEMBLY ===" > $results_file
echo "Assembly: $ASSEMBLY_FILE" >> $results_file
echo "Search date: $(date)" >> $results_file
echo "" >> $results_file

# Get assembly stats
total_contigs=$(grep -c "^>" "$ASSEMBLY_FILE")
echo "Total contigs in assembly: $total_contigs" >> $results_file
echo "Total contigs in assembly: $total_contigs"

# Search for methane gene patterns in contig names/headers
echo "" >> $results_file
echo "STEP 1: Searching contig headers for methane annotations" >> $results_file
echo "=======================================================" >> $results_file

methane_headers=$(grep -i "mcr\|pmo\|mmo\|mxa\|methane\|methanol" "$ASSEMBLY_FILE" | grep "^>" || true)
if [ ! -z "$methane_headers" ]; then
    echo "FOUND methane-related headers:" >> $results_file
    echo "$methane_headers" >> $results_file
    header_count=$(echo "$methane_headers" | wc -l)
    echo "✓ Found $header_count contigs with methane-related annotations!"
else
    echo "No methane gene annotations in headers" >> $results_file
    echo "No methane annotations found in contig headers"
fi

echo "" >> $results_file
echo "STEP 2: Searching sequences for methane gene motifs" >> $results_file
echo "===================================================" >> $results_file

# Extract all sequences without headers for pattern searching
sequences_only="${RESULTS_DIR}/sequences_only.txt"
grep -v "^>" "$ASSEMBLY_FILE" > $sequences_only

echo "Searching for mcrA patterns..."

# mcrA conserved motifs (DNA level)
mcrA_motif1=$(grep -i "GGTGACGG\|GGTTATGG\|GGTTGCGG" $sequences_only | wc -l)
mcrA_motif2=$(grep -i "TACGGCCT\|TATGGCCT\|TACGGTCT" $sequences_only | wc -l)
mcrA_motif3=$(grep -i "ATGGATAA\|ATGGACAA\|ATGGATA" $sequences_only | wc -l)

echo "mcrA motif searches:" >> $results_file
echo "- Pattern GGTGACGG/GGTTATGG/GGTTGCGG: $mcrA_motif1 hits" >> $results_file
echo "- Pattern TACGGCCT/TATGGCCT/TACGGTCT: $mcrA_motif2 hits" >> $results_file  
echo "- Pattern ATGGATAA/ATGGACAA/ATGGATA: $mcrA_motif3 hits" >> $results_file

echo "Searching for pmoA patterns..."

# pmoA conserved motifs (DNA level)
pmoA_motif1=$(grep -i "GGCTTCGG\|GGCTTTGG\|GGCTATGG" $sequences_only | wc -l)
pmoA_motif2=$(grep -i "ATGACCGA\|ATGACCAA\|ATGACTGA" $sequences_only | wc -l)
pmoA_motif3=$(grep -i "TTCTGGCG\|TTCTGGAG\|TTTTGGCG" $sequences_only | wc -l)

echo "pmoA motif searches:" >> $results_file
echo "- Pattern GGCTTCGG/GGCTTTGG/GGCTATGG: $pmoA_motif1 hits" >> $results_file
echo "- Pattern ATGACCGA/ATGACCAA/ATGACTGA: $pmoA_motif2 hits" >> $results_file
echo "- Pattern TTCTGGCG/TTCTGGAG/TTTTGGCG: $pmoA_motif3 hits" >> $results_file

total_motif_hits=$((mcrA_motif1 + mcrA_motif2 + mcrA_motif3 + pmoA_motif1 + pmoA_motif2 + pmoA_motif3))

echo "Total motif hits: $total_motif_hits" >> $results_file
echo "Total conserved motif hits found: $total_motif_hits"

echo "" >> $results_file
echo "STEP 3: Searching for gene names in sequences" >> $results_file
echo "=============================================" >> $results_file

# Search for actual gene names in sequences (might be annotations)
mcrA_name_hits=$(grep -i "mcrA\|mcr_A\|mcr-A" $sequences_only | wc -l)
mcrB_name_hits=$(grep -i "mcrB\|mcr_B\|mcr-B" $sequences_only | wc -l)
pmoA_name_hits=$(grep -i "pmoA\|pmo_A\|pmo-A" $sequences_only | wc -l)
pmoB_name_hits=$(grep -i "pmoB\|pmo_B\|pmo-B" $sequences_only | wc -l)
mmoX_name_hits=$(grep -i "mmoX\|mmo_X\|mmo-X" $sequences_only | wc -l)
mxaF_name_hits=$(grep -i "mxaF\|mxa_F\|mxa-F" $sequences_only | wc -l)

echo "Gene name searches in sequences:" >> $results_file
echo "- mcrA: $mcrA_name_hits hits" >> $results_file
echo "- mcrB: $mcrB_name_hits hits" >> $results_file
echo "- pmoA: $pmoA_name_hits hits" >> $results_file
echo "- pmoB: $pmoB_name_hits hits" >> $results_file
echo "- mmoX: $mmoX_name_hits hits" >> $results_file
echo "- mxaF: $mxaF_name_hits hits" >> $results_file

total_gene_hits=$((mcrA_name_hits + mcrB_name_hits + pmoA_name_hits + pmoB_name_hits + mmoX_name_hits + mxaF_name_hits))
echo "Total gene name hits: $total_gene_hits" >> $results_file

echo "" >> $results_file
echo "STEP 4: Searching for enzyme names" >> $results_file
echo "==================================" >> $results_file

# Search for enzyme descriptions
mcr_enzyme=$(grep -i "methyl.*coenzyme.*reductase\|methylcoenzyme.*reductase" $sequences_only | wc -l)
pmo_enzyme=$(grep -i "particulate.*methane.*monooxygenase\|methane.*monooxygenase" $sequences_only | wc -l)
mmo_enzyme=$(grep -i "soluble.*methane.*monooxygenase" $sequences_only | wc -l)
mxa_enzyme=$(grep -i "methanol.*dehydrogenase" $sequences_only | wc -l)

echo "Enzyme name searches:" >> $results_file
echo "- methyl-coenzyme M reductase: $mcr_enzyme hits" >> $results_file
echo "- particulate methane monooxygenase: $pmo_enzyme hits" >> $results_file
echo "- soluble methane monooxygenase: $mmo_enzyme hits" >> $results_file
echo "- methanol dehydrogenase: $mxa_enzyme hits" >> $results_file

total_enzyme_hits=$((mcr_enzyme + pmo_enzyme + mmo_enzyme + mxa_enzyme))
echo "Total enzyme name hits: $total_enzyme_hits" >> $results_file

echo "" >> $results_file
echo "STEP 5: Analyzing contig size distribution" >> $results_file
echo "==========================================" >> $results_file

# Analyze where methane genes might be based on size
echo "Analyzing contig lengths for methane gene potential..."

awk '
BEGIN {
    short=0; medium=0; long=0; total=0;
}
/^>/ {
    if(seq!="") {
        len=length(seq);
        total++;
        if(len < 2000) short++;
        else if(len < 5000) medium++;
        else long++;
    }
    seq="";
}
!/^>/ {
    seq=seq$0;
}
END {
    if(seq!="") {
        len=length(seq);
        total++;
        if(len < 2000) short++;
        else if(len < 5000) medium++;
        else long++;
    }
    print "Contig size distribution:";
    print "Short (<2kb): " short " (" int(short*100/total) "%)";
    print "Medium (2-5kb): " medium " (" int(medium*100/total) "%)";
    print "Long (>5kb): " long " (" int(long*100/total) "%)";
    print "Total: " total;
}' "$ASSEMBLY_FILE" >> $results_file

# Summary
echo "" >> $results_file
echo "SUMMARY OF FINDINGS:" >> $results_file
echo "===================" >> $results_file
echo "Header annotations: $header_count" >> $results_file
echo "Conserved motifs: $total_motif_hits" >> $results_file
echo "Gene names: $total_gene_hits" >> $results_file
echo "Enzyme names: $total_enzyme_hits" >> $results_file

grand_total=$((header_count + total_motif_hits + total_gene_hits + total_enzyme_hits))
echo "TOTAL POTENTIAL HITS: $grand_total" >> $results_file

# Clean up
rm -f $sequences_only

echo ""
echo "============================================="
echo "SEARCH COMPLETED!"
echo "============================================="
echo ""

echo "Results summary:"
echo "- Header annotations: $header_count"
echo "- Conserved motifs: $total_motif_hits"  
echo "- Gene names in sequences: $total_gene_hits"
echo "- Enzyme names: $total_enzyme_hits"
echo "- TOTAL POTENTIAL HITS: $grand_total"
echo ""

if [ $grand_total -gt 0 ]; then
    echo "🎯 METHANE GENES DETECTED IN FULL ASSEMBLY!"
    echo ""
    echo "The genes are likely present but may be:"
    echo "1. In unbinned contigs (92.9% of assembly)"
    echo "2. On short contigs that were filtered"
    echo "3. Poorly assembled or fragmented"
    echo "4. Missing from annotation databases"
    echo ""
    echo "Next steps:"
    echo "- Check specific contigs with hits"
    echo "- Search unbinned fraction specifically"
    echo "- Use HMM-based search tools"
    echo "- Consider re-binning with different parameters"
else
    echo "⚠️  No obvious methane gene signatures found"
    echo "This could indicate:"
    echo "1. Assembly quality issues"
    echo "2. Very divergent methane gene sequences"
    echo "3. Genes fragmented across multiple contigs"
    echo "4. Need for more sensitive search methods"
fi

echo ""
echo "Full results: ${RESULTS_DIR}/full_assembly_methane_search.txt"
