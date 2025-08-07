#!/bin/bash

# Comprehensive methane gene search in annotated MAGs
# Author: Generated for farting_trees project
# Date: August 2025

set -e

# Configuration
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOTATION_DIR="${WORK_DIR}/results/05_annotation/prokka"
RESULTS_DIR="${WORK_DIR}/results/07_methane_search"

echo "============================================="
echo "COMPREHENSIVE METHANE GENE SEARCH"
echo "============================================="
echo "Searching for methane metabolism genes in your 8 co-assembly bins"
echo ""

# Check if annotation directory exists
if [ ! -d "$ANNOTATION_DIR" ]; then
    echo "Error: Annotation directory not found: $ANNOTATION_DIR"
    echo "Please run the annotation script first"
    exit 1
fi

# Create output directory
mkdir -p ${RESULTS_DIR}

echo "Using annotation directory: $ANNOTATION_DIR"
echo "Results will be saved to: $RESULTS_DIR"
echo ""

# Define methane-specific genes with biological significance
echo "Gene categories being searched:"
echo "1. Core methanogenesis: mcrA, mcrB, mcrG (methyl-coenzyme M reductase)"
echo "2. Core methanotrophy: pmoA, pmoB, pmoC (particulate methane monooxygenase)"
echo "3. Alternative methanotrophy: mmoX, mmoY, mmoZ (soluble methane monooxygenase)"
echo "4. Methanol oxidation: mxaF, mxaI (methanol dehydrogenase)"
echo "5. Auxiliary methanogenesis: fmd, ftr, mtd (C1 metabolism enzymes)"
echo ""

# Initialize summary file
summary_file="${RESULTS_DIR}/methane_gene_summary.txt"
echo "=== COMPREHENSIVE METHANE GENE SEARCH RESULTS ===" > $summary_file
echo "" >> $summary_file
echo "Search Date: $(date)" >> $summary_file
echo "MAGs Analyzed: 8 co-assembly bins" >> $summary_file
echo "" >> $summary_file

# Initialize counters
total_mcrA=0
total_pmoA=0
total_mmoX=0
total_mxaF=0
total_fmd=0
total_ftr=0
total_mtd=0
total_other_methane=0
bins_with_methane=0

echo "Detailed results by MAG:" >> $summary_file
echo "========================" >> $summary_file

# Search each bin
for bin_dir in ${ANNOTATION_DIR}/coassembly_bin.*; do
    if [ -d "$bin_dir" ]; then
        bin_name=$(basename "$bin_dir")
        echo "Analyzing $bin_name..."
        
        gff_file="$bin_dir/${bin_name}.gff"
        tsv_file="$bin_dir/${bin_name}.tsv"
        bin_results="${RESULTS_DIR}/${bin_name}_methane_genes.txt"
        
        echo "=== $bin_name METHANE GENE ANALYSIS ===" > $bin_results
        echo "" >> $bin_results
        
        # Count total genes
        if [ -f "$gff_file" ]; then
            total_genes=$(grep -c "^[^#].*CDS" "$gff_file" 2>/dev/null || echo "0")
            echo "Total CDS features: $total_genes" >> $bin_results
            echo "" >> $bin_results
        else
            echo "GFF file not found for $bin_name" >> $bin_results
            continue
        fi
        
        bin_methane_genes=0
        
        echo "CORE METHANOGENESIS GENES:" >> $bin_results
        echo "=========================" >> $bin_results
        
        # mcrA - KEY methanogenesis marker
        mcrA_hits=$(grep -i "mcrA\|methyl.*coenzyme.*reductase.*alpha" "$gff_file" 2>/dev/null || true)
        mcrA_count=0
        if [ ! -z "$mcrA_hits" ]; then
            mcrA_count=$(echo "$mcrA_hits" | wc -l)
            echo "✓ mcrA (methanogenesis marker): $mcrA_count copies" >> $bin_results
            echo "$mcrA_hits" >> $bin_results
            echo "" >> $bin_results
            total_mcrA=$((total_mcrA + mcrA_count))
            bin_methane_genes=$((bin_methane_genes + mcrA_count))
        fi
        
        # mcrB and mcrG
        other_mcr=$(grep -i "mcrB\|mcrG\|methyl.*coenzyme.*reductase.*beta\|methyl.*coenzyme.*reductase.*gamma" "$gff_file" 2>/dev/null || true)
        other_mcr_count=0
        if [ ! -z "$other_mcr" ]; then
            other_mcr_count=$(echo "$other_mcr" | wc -l)
            echo "✓ mcrB/mcrG: $other_mcr_count copies" >> $bin_results
            echo "$other_mcr" >> $bin_results
            echo "" >> $bin_results
            total_other_methane=$((total_other_methane + other_mcr_count))
            bin_methane_genes=$((bin_methane_genes + other_mcr_count))
        fi
        
        echo "" >> $bin_results
        echo "CORE METHANOTROPHY GENES:" >> $bin_results
        echo "========================" >> $bin_results
        
        # pmoA - KEY methanotrophy marker
        pmoA_hits=$(grep -i "pmoA\|particulate.*methane.*monooxygenase.*alpha" "$gff_file" 2>/dev/null || true)
        pmoA_count=0
        if [ ! -z "$pmoA_hits" ]; then
            pmoA_count=$(echo "$pmoA_hits" | wc -l)
            echo "✓ pmoA (methanotrophy marker): $pmoA_count copies" >> $bin_results
            echo "$pmoA_hits" >> $bin_results
            echo "" >> $bin_results
            total_pmoA=$((total_pmoA + pmoA_count))
            bin_methane_genes=$((bin_methane_genes + pmoA_count))
        fi
        
        # pmoB and pmoC
        other_pmo=$(grep -i "pmoB\|pmoC\|particulate.*methane.*monooxygenase.*beta\|particulate.*methane.*monooxygenase.*gamma" "$gff_file" 2>/dev/null || true)
        other_pmo_count=0
        if [ ! -z "$other_pmo" ]; then
            other_pmo_count=$(echo "$other_pmo" | wc -l)
            echo "✓ pmoB/pmoC: $other_pmo_count copies" >> $bin_results
            echo "$other_pmo" >> $bin_results
            echo "" >> $bin_results
            total_other_methane=$((total_other_methane + other_pmo_count))
            bin_methane_genes=$((bin_methane_genes + other_pmo_count))
        fi
        
        echo "" >> $bin_results
        echo "ALTERNATIVE METHANOTROPHY:" >> $bin_results
        echo "=========================" >> $bin_results
        
        # mmoX - Alternative methanotrophy
        mmoX_hits=$(grep -i "mmoX\|soluble.*methane.*monooxygenase.*alpha" "$gff_file" 2>/dev/null || true)
        mmoX_count=0
        if [ ! -z "$mmoX_hits" ]; then
            mmoX_count=$(echo "$mmoX_hits" | wc -l)
            echo "✓ mmoX (alternative methanotrophy): $mmoX_count copies" >> $bin_results
            echo "$mmoX_hits" >> $bin_results
            echo "" >> $bin_results
            total_mmoX=$((total_mmoX + mmoX_count))
            bin_methane_genes=$((bin_methane_genes + mmoX_count))
        fi
        
        # mmoY and mmoZ
        other_mmo=$(grep -i "mmoY\|mmoZ\|soluble.*methane.*monooxygenase.*beta\|soluble.*methane.*monooxygenase.*gamma" "$gff_file" 2>/dev/null || true)
        other_mmo_count=0
        if [ ! -z "$other_mmo" ]; then
            other_mmo_count=$(echo "$other_mmo" | wc -l)
            echo "✓ mmoY/mmoZ: $other_mmo_count copies" >> $bin_results
            echo "$other_mmo" >> $bin_results
            echo "" >> $bin_results
            total_other_methane=$((total_other_methane + other_mmo_count))
            bin_methane_genes=$((bin_methane_genes + other_mmo_count))
        fi
        
        echo "" >> $bin_results
        echo "METHANOL OXIDATION:" >> $bin_results
        echo "==================" >> $bin_results
        
        # mxaF - Methanol oxidation
        mxaF_hits=$(grep -i "mxaF\|methanol.*dehydrogenase.*large\|methanol.*dehydrogenase.*alpha" "$gff_file" 2>/dev/null || true)
        mxaF_count=0
        if [ ! -z "$mxaF_hits" ]; then
            mxaF_count=$(echo "$mxaF_hits" | wc -l)
            echo "✓ mxaF (methanol oxidation): $mxaF_count copies" >> $bin_results
            echo "$mxaF_hits" >> $bin_results
            echo "" >> $bin_results
            total_mxaF=$((total_mxaF + mxaF_count))
            bin_methane_genes=$((bin_methane_genes + mxaF_count))
        fi
        
        # mxaI
        mxaI_hits=$(grep -i "mxaI\|methanol.*dehydrogenase.*small" "$gff_file" 2>/dev/null || true)
        mxaI_count=0
        if [ ! -z "$mxaI_hits" ]; then
            mxaI_count=$(echo "$mxaI_hits" | wc -l)
            echo "✓ mxaI: $mxaI_count copies" >> $bin_results
            echo "$mxaI_hits" >> $bin_results
            echo "" >> $bin_results
            total_other_methane=$((total_other_methane + mxaI_count))
            bin_methane_genes=$((bin_methane_genes + mxaI_count))
        fi
        
        echo "" >> $bin_results
        echo "AUXILIARY METHANOGENESIS:" >> $bin_results
        echo "========================" >> $bin_results
        
        # fmd - Formylmethanofuran dehydrogenase
        fmd_hits=$(grep -i "fmd\|formylmethanofuran.*dehydrogenase" "$gff_file" 2>/dev/null || true)
        fmd_count=0
        if [ ! -z "$fmd_hits" ]; then
            fmd_count=$(echo "$fmd_hits" | wc -l)
            echo "✓ fmd (formylmethanofuran dehydrogenase): $fmd_count copies" >> $bin_results
            echo "$fmd_hits" >> $bin_results
            echo "" >> $bin_results
            total_fmd=$((total_fmd + fmd_count))
            bin_methane_genes=$((bin_methane_genes + fmd_count))
        fi
        
        # ftr - Formyltetrahydrofolate synthetase
        ftr_hits=$(grep -i "ftr\|formyltetrahydrofolate.*synthetase" "$gff_file" 2>/dev/null || true)
        ftr_count=0
        if [ ! -z "$ftr_hits" ]; then
            ftr_count=$(echo "$ftr_hits" | wc -l)
            echo "✓ ftr (formyltetrahydrofolate synthetase): $ftr_count copies" >> $bin_results
            echo "$ftr_hits" >> $bin_results
            echo "" >> $bin_results
            total_ftr=$((total_ftr + ftr_count))
            bin_methane_genes=$((bin_methane_genes + ftr_count))
        fi
        
        # mtd - Methylenetetrahydrofolate dehydrogenase
        mtd_hits=$(grep -i "mtd\|methylenetetrahydrofolate.*dehydrogenase" "$gff_file" 2>/dev/null || true)
        mtd_count=0
        if [ ! -z "$mtd_hits" ]; then
            mtd_count=$(echo "$mtd_hits" | wc -l)
            echo "✓ mtd (methylenetetrahydrofolate dehydrogenase): $mtd_count copies" >> $bin_results
            echo "$mtd_hits" >> $bin_results
            echo "" >> $bin_results
            total_mtd=$((total_mtd + mtd_count))
            bin_methane_genes=$((bin_methane_genes + mtd_count))
        fi
        
        # Summary for this bin
        if [ $bin_methane_genes -gt 0 ]; then
            bins_with_methane=$((bins_with_methane + 1))
            echo "$bin_name: $bin_methane_genes methane genes (mcrA:$mcrA_count, pmoA:$pmoA_count, mmoX:$mmoX_count, mxaF:$mxaF_count, fmd:$fmd_count, ftr:$ftr_count, mtd:$mtd_count, other:$((other_mcr_count + other_pmo_count + other_mmo_count + mxaI_count)))" >> $summary_file
            echo "✓ $bin_name: FOUND $bin_methane_genes methane-related genes"
        else
            echo "$bin_name: 0 methane genes" >> $summary_file
            echo "No methane-specific genes found" >> $bin_results
            echo "  $bin_name: No methane genes found"
        fi
    fi
done

# Final summary
echo "" >> $summary_file
echo "FINAL SUMMARY:" >> $summary_file
echo "==============" >> $summary_file
echo "MAGs analyzed: 8" >> $summary_file
echo "MAGs with methane genes: $bins_with_methane" >> $summary_file
echo "" >> $summary_file
echo "Gene-specific counts:" >> $summary_file
echo "mcrA (methanogenesis marker): $total_mcrA" >> $summary_file
echo "pmoA (methanotrophy marker): $total_pmoA" >> $summary_file
echo "mmoX (alternative methanotrophy): $total_mmoX" >> $summary_file
echo "mxaF (methanol oxidation): $total_mxaF" >> $summary_file
echo "fmd (auxiliary methanogenesis): $total_fmd" >> $summary_file
echo "ftr (auxiliary methanogenesis): $total_ftr" >> $summary_file
echo "mtd (auxiliary methanogenesis): $total_mtd" >> $summary_file
echo "Other methane-related genes: $total_other_methane" >> $summary_file
echo "TOTAL methane-specific genes: $((total_mcrA + total_pmoA + total_mmoX + total_mxaF + total_fmd + total_ftr + total_mtd + total_other_methane))" >> $summary_file

echo ""
echo "============================================="
echo "METHANE GENE SEARCH COMPLETED!"
echo "============================================="
echo ""
echo "Key findings:"
echo "- mcrA (methanogenesis): $total_mcrA copies"
echo "- pmoA (methanotrophy): $total_pmoA copies" 
echo "- mmoX (alt. methanotrophy): $total_mmoX copies"
echo "- mxaF (methanol oxidation): $total_mxaF copies"
echo "- Auxiliary genes (fmd/ftr/mtd): $((total_fmd + total_ftr + total_mtd)) copies"
echo "- Other methane genes: $total_other_methane copies"
echo ""
echo "Total methane-specific genes: $((total_mcrA + total_pmoA + total_mmoX + total_mxaF + total_fmd + total_ftr + total_mtd + total_other_methane))"
echo "MAGs with methane genes: $bins_with_methane out of 8"
echo ""

# Biological interpretation
if [ $total_mcrA -gt 0 ]; then
    echo "METHANOGENESIS DETECTED!"
    echo "✓ Found $total_mcrA mcrA genes - CONFIRMED methane-producing organisms!"
fi

if [ $total_pmoA -gt 0 ]; then
    echo "METHANOTROPHY DETECTED!"
    echo "✓ Found $total_pmoA pmoA genes - CONFIRMED methane-consuming organisms!"
fi

if [ $total_mmoX -gt 0 ]; then
    echo "ALTERNATIVE METHANOTROPHY DETECTED!"
    echo "✓ Found $total_mmoX mmoX genes - Alternative methane oxidation pathway!"
fi

if [ $total_mxaF -gt 0 ]; then
    echo "METHANOL METABOLISM DETECTED!"
    echo "✓ Found $total_mxaF mxaF genes - Methanol oxidation potential!"
fi

if [ $((total_mcrA + total_pmoA + total_mmoX)) -eq 0 ]; then
    echo "NO CORE METHANE GENES FOUND"
    echo "Your tree microbiomes do not contain the key enzymes for methane production or consumption."
    if [ $((total_fmd + total_ftr + total_mtd)) -gt 0 ]; then
        echo "However, auxiliary C1 metabolism genes were found, indicating general carbon metabolism."
    fi
fi

echo ""
echo "Results saved in: ${RESULTS_DIR}/"
echo "Summary: ${RESULTS_DIR}/methane_gene_summary.txt"
echo "Individual results: ${RESULTS_DIR}/coassembly_bin.*_methane_genes.txt"
echo ""
echo "Biological significance:"
echo "- mcrA = Methane production (methanogenesis)"
echo "- pmoA = Methane consumption (methanotrophy)"  
echo "- mmoX = Alternative methane consumption"
echo "- mxaF = Methanol utilization (related to methane)"
echo "- fmd/ftr/mtd = Supporting C1 metabolism"
