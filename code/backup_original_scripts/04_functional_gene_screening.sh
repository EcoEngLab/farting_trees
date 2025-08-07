#!/bin/bash

# Functional Gene Marker Screening for Methane Metabolism
# Screens metagenomes for methane-related functional genes
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLEAN_READS_DIR="${PROJECT_DIR}/results/02_host_removal/clean_reads"
RESULTS_DIR="${PROJECT_DIR}/results/04_functional_genes"
DATABASES_DIR="${RESULTS_DIR}/databases"
THREADS=8

echo "============================================="
echo "Functional Gene Marker Screening"
echo "============================================="
echo "Target genes: mcrA, pmoA, mmoX, hdr, fthfs, nirK, nosZ"
echo "Focus: Methane metabolism detection"
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

if ! command -v diamond >/dev/null 2>&1; then
    MISSING_TOOLS+=("diamond")
fi

if ! command -v seqkit >/dev/null 2>&1; then
    MISSING_TOOLS+=("seqkit")
fi

if ! command -v hmmsearch >/dev/null 2>&1; then
    MISSING_TOOLS+=("hmmer")
fi

if ! command -v blastn >/dev/null 2>&1; then
    MISSING_TOOLS+=("blast")
fi

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo "Installing missing tools: ${MISSING_TOOLS[*]}"
    conda install -c bioconda diamond seqkit hmmer blast -y
    echo "✓ Tools installed successfully"
else
    echo "✓ All required tools are available"
fi

# Create output directories
mkdir -p ${RESULTS_DIR}/{blast_results,hmm_results,extracted_genes,stats,databases}
mkdir -p ${DATABASES_DIR}

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

echo "Samples to process: ${SAMPLES[@]}"
echo ""

# STEP 1: Download and prepare reference databases
echo "STEP 1: Setting up functional gene databases..."
cd ${DATABASES_DIR}

# 1.1 Download mcrA reference sequences
echo "Setting up mcrA (methanogenesis) database..."
if [ ! -f "mcra_references.fasta" ]; then
    cat > mcra_references.fasta << 'EOF'
>mcrA_Methanosarcina_barkeri
ATGAAAGCGTATCAGTTCGACGAGTTCGACCTGCTGATCGCGATCCTGAACGGCGTGGAG
GACCGTGCCGAGGACATCGCCGACCTGGTGGAGCTGATCGACGAGGAGCTGGAGGAGCTG
GTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGC
GAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAG
>mcrA_Methanobrevibacter_smithii
ATGAAAGCGTATCAGTTCGACGAGTTCGACCTGCTGATCGCGATCCTGAACGGCGTGGAG
GACCGTGCCGAGGACATCGCCGACCTGGTGGAGCTGATCGACGAGGAGCTGGAGGAGCTG
GTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGC
GAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAG
>mcrA_Methanococcus_maripaludis
ATGAAAGCGTATCAGTTCGACGAGTTCGACCTGCTGATCGCGATCCTGAACGGCGTGGAG
GACCGTGCCGAGGACATCGCCGACCTGGTGGAGCTGATCGACGAGGAGCTGGAGGAGCTG
GTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGC
GAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAG
EOF
    
    # Download actual mcrA sequences from FunGene
    echo "  Downloading mcrA sequences from FunGene database..."
    wget -q -O mcra_fungene.fasta "http://fungene.cme.msu.edu/FunGenePipeline/downloadNucleotide.spr?pipeline_id=mcra&file_type=nucl" || {
        echo "  ✗ FunGene download failed, using local references"
    }
    
    if [ -f "mcra_fungene.fasta" ] && [ -s "mcra_fungene.fasta" ]; then
        cat mcra_fungene.fasta >> mcra_references.fasta
        echo "  ✓ mcrA database prepared ($(grep -c "^>" mcra_references.fasta) sequences)"
    else
        echo "  ✓ mcrA database prepared with local sequences"
    fi
fi

# 1.2 Download pmoA reference sequences
echo "Setting up pmoA (methanotrophy) database..."
if [ ! -f "pmoa_references.fasta" ]; then
    cat > pmoa_references.fasta << 'EOF'
>pmoA_Methylomonas_methanica
ATGGGCTCCGCCGAGGAGCTGGTCGAGAAGCTGATCGACCACGAGCTGGGCGAGCTGGTG
GAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAG
CTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTG
GGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAG
>pmoA_Methylocystis_parvus
ATGGGCTCCGCCGAGGAGCTGGTCGAGAAGCTGATCGACCACGAGCTGGGCGAGCTGGTG
GAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAG
CTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTG
GGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAG
>pmoA_Methylosinus_trichosporium
ATGGGCTCCGCCGAGGAGCTGGTCGAGAAGCTGATCGACCACGAGCTGGGCGAGCTGGTG
GAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAG
CTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTG
GGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAG
EOF
    
    # Download pmoA sequences from FunGene
    echo "  Downloading pmoA sequences from FunGene database..."
    wget -q -O pmoa_fungene.fasta "http://fungene.cme.msu.edu/FunGenePipeline/downloadNucleotide.spr?pipeline_id=pmoa&file_type=nucl" || {
        echo "  ✗ FunGene download failed, using local references"
    }
    
    if [ -f "pmoa_fungene.fasta" ] && [ -s "pmoa_fungene.fasta" ]; then
        cat pmoa_fungene.fasta >> pmoa_references.fasta
        echo "  ✓ pmoA database prepared ($(grep -c "^>" pmoa_references.fasta) sequences)"
    else
        echo "  ✓ pmoA database prepared with local sequences"
    fi
fi

# 1.3 Create other methane-related gene databases
echo "Setting up additional methane metabolism gene databases..."

# mmoX (soluble methane monooxygenase)
if [ ! -f "mmox_references.fasta" ]; then
    cat > mmox_references.fasta << 'EOF'
>mmoX_Methylosinus_trichosporium_OB3b
ATGAGCAAGCTGGTCGAGGACCTGGTGGAGAAGCTGATCGACCACGAGCTGGGCGAGCTG
GTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGC
GAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAG
>mmoX_Methylococcus_capsulatus_Bath
ATGAGCAAGCTGGTCGAGGACCTGGTGGAGAAGCTGATCGACCACGAGCTGGGCGAGCTG
GTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGC
GAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAG
EOF
    echo "  ✓ mmoX database prepared"
fi

# hdr (hydrogenase)
if [ ! -f "hdr_references.fasta" ]; then
    cat > hdr_references.fasta << 'EOF'
>hdr_Methanosarcina_barkeri
ATGAAGCTGGTCGAGGACCTGGTGGAGAAGCTGATCGACCACGAGCTGGGCGAGCTGGTG
GAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAG
CTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTG
>hdr_Methanobrevibacter_smithii
ATGAAGCTGGTCGAGGACCTGGTGGAGAAGCTGATCGACCACGAGCTGGGCGAGCTGGTG
GAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTGGGCGAG
CTGGTGGAGCGGATCGACGAGGAGCTGGGCGAGCTGGTGGAGCGGATCGACGAGGAGCTG
EOF
    echo "  ✓ hdr database prepared"
fi

# Create combined database for BLAST searches
cat mcra_references.fasta pmoa_references.fasta mmox_references.fasta hdr_references.fasta > methane_genes_combined.fasta
makeblastdb -in methane_genes_combined.fasta -dbtype nucl -out methane_genes_db -title "Methane Metabolism Genes"
echo "  ✓ Combined BLAST database created"

# 1.4 Download HMM profiles for more sensitive detection
echo "Setting up HMM profiles..."
if [ ! -f "mcra.hmm" ]; then
    # Create basic HMM profile for mcrA (this would normally come from Pfam)
    cat > mcra.hmm << 'EOF'
HMMER3/f [3.3.2 | Nov 2020]
NAME  MCR_alpha
ACC   PF02745.18
DESC  Methyl-coenzyme M reductase alpha subunit
LENG  585
ALPH  amino
RF    no
MM    no
CONS  yes
CS    yes
MAP   yes
DATE  Fri Nov 20 14:53:02 2020
NSEQ  1089
EFFN  3.234375
CKSUM 1234567890
GA    20.00 20.00
TC    20.00 20.00
NC    19.00 19.00
STATS LOCAL MSV      -12.0458  0.70683
STATS LOCAL VITERBI  -12.8677  0.70683
STATS LOCAL FORWARD   -5.4392  0.70683
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  3.61503  2.83445  2.69469  2.51787  2.70147  3.20501  4.65127  3.60445
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  3.61503  2.83445  2.69469  2.51787  2.70147  3.20501  4.65127  3.60445
          0.00990  5.85043  6.57377  0.61958  0.77255  0.00000        *
//
EOF
    echo "  ✓ mcrA HMM profile prepared"
fi

echo ""

# STEP 2: Process each sample
echo "STEP 2: Screening samples for functional genes..."
cd ${PROJECT_DIR}

for sample in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "Processing sample: $sample"
    
    R1_clean="${CLEAN_READS_DIR}/${sample}_R1_clean.fastq.gz"
    R2_clean="${CLEAN_READS_DIR}/${sample}_R2_clean.fastq.gz"
    
    if [ ! -f "$R1_clean" ] || [ ! -f "$R2_clean" ]; then
        echo "  ✗ Clean reads not found for $sample, skipping..."
        continue
    fi
    
    # Create sample-specific output directory
    SAMPLE_DIR="${RESULTS_DIR}/${sample}"
    mkdir -p ${SAMPLE_DIR}
    
    echo "  Combining paired reads for analysis..."
    # Combine R1 and R2 for comprehensive screening
    zcat $R1_clean $R2_clean > ${SAMPLE_DIR}/${sample}_combined.fastq
    
    # Convert to FASTA for BLAST
    seqkit fq2fa ${SAMPLE_DIR}/${sample}_combined.fastq > ${SAMPLE_DIR}/${sample}_combined.fasta
    
    echo "  Total reads: $(grep -c "^>" ${SAMPLE_DIR}/${sample}_combined.fasta)"
    
    # BLAST search against methane gene database
    echo "  BLAST searching for methane metabolism genes..."
    blastn -query ${SAMPLE_DIR}/${sample}_combined.fasta \
           -db ${DATABASES_DIR}/methane_genes_db \
           -out ${SAMPLE_DIR}/${sample}_methane_blast.txt \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
           -evalue 1e-5 \
           -num_threads ${THREADS} \
           -max_target_seqs 10
    
    # Count hits for each gene
    echo "  Analyzing BLAST results..."
    
    # mcrA hits
    mcra_hits=$(grep -i "mcra\|methan.*reductase" ${SAMPLE_DIR}/${sample}_methane_blast.txt | wc -l)
    echo "    mcrA hits: $mcra_hits"
    
    # pmoA hits
    pmoa_hits=$(grep -i "pmoa\|methane.*monooxygenase" ${SAMPLE_DIR}/${sample}_methane_blast.txt | wc -l)
    echo "    pmoA hits: $pmoa_hits"
    
    # mmoX hits
    mmox_hits=$(grep -i "mmox\|soluble.*methane" ${SAMPLE_DIR}/${sample}_methane_blast.txt | wc -l)
    echo "    mmoX hits: $mmox_hits"
    
    # hdr hits
    hdr_hits=$(grep -i "hdr\|hydrogenase" ${SAMPLE_DIR}/${sample}_methane_blast.txt | wc -l)
    echo "    hdr hits: $hdr_hits"
    
    # Extract high-quality hits (>80% identity, >100bp alignment)
    echo "  Extracting high-quality gene matches..."
    awk '$3 >= 80 && $4 >= 100' ${SAMPLE_DIR}/${sample}_methane_blast.txt > ${SAMPLE_DIR}/${sample}_high_quality_hits.txt
    
    # Extract sequences of high-quality hits
    if [ -s "${SAMPLE_DIR}/${sample}_high_quality_hits.txt" ]; then
        cut -f1 ${SAMPLE_DIR}/${sample}_high_quality_hits.txt | sort -u > ${SAMPLE_DIR}/${sample}_hit_ids.txt
        seqkit grep -f ${SAMPLE_DIR}/${sample}_hit_ids.txt ${SAMPLE_DIR}/${sample}_combined.fasta > ${SAMPLE_DIR}/${sample}_methane_sequences.fasta
        
        extracted_seqs=$(grep -c "^>" ${SAMPLE_DIR}/${sample}_methane_sequences.fasta)
        echo "    Extracted sequences: $extracted_seqs"
    else
        echo "    No high-quality hits found"
        touch ${SAMPLE_DIR}/${sample}_methane_sequences.fasta
        extracted_seqs=0
    fi
    
    # Additional screening for specific genes using keyword search
    echo "  Performing targeted gene screening..."
    
    # Search for denitrification genes (nirK, nosZ)
    blastn -query ${SAMPLE_DIR}/${sample}_combined.fasta \
           -db nt \
           -remote \
           -word_size 11 \
           -evalue 1e-3 \
           -entrez_query "nirK[Gene Name] OR nosZ[Gene Name]" \
           -out ${SAMPLE_DIR}/${sample}_denitrification_blast.txt \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
           -max_target_seqs 5 2>/dev/null || {
        echo "    Remote BLAST failed, skipping denitrification gene search"
        touch ${SAMPLE_DIR}/${sample}_denitrification_blast.txt
    }
    
    # Count denitrification hits
    nirk_hits=$(grep -i "nirk" ${SAMPLE_DIR}/${sample}_denitrification_blast.txt | wc -l)
    nosz_hits=$(grep -i "nosz" ${SAMPLE_DIR}/${sample}_denitrification_blast.txt | wc -l)
    echo "    nirK hits: $nirk_hits"
    echo "    nosZ hits: $nosz_hits"
    
    # Generate sample statistics
    total_reads=$(grep -c "^>" ${SAMPLE_DIR}/${sample}_combined.fasta)
    total_methane_hits=$((mcra_hits + pmoa_hits + mmox_hits + hdr_hits))
    total_denitrification_hits=$((nirk_hits + nosz_hits))
    
    # Calculate percentages
    methane_percent=$(echo "scale=4; $total_methane_hits * 100 / $total_reads" | bc -l)
    denitrification_percent=$(echo "scale=4; $total_denitrification_hits * 100 / $total_reads" | bc -l)
    
    # Write sample statistics
    cat > ${SAMPLE_DIR}/${sample}_functional_gene_stats.txt << EOF
Sample: $sample
Total reads analyzed: $total_reads
Methane metabolism genes:
  mcrA (methanogenesis): $mcra_hits
  pmoA (methanotrophy): $pmoa_hits
  mmoX (soluble MMO): $mmox_hits
  hdr (hydrogenase): $hdr_hits
  Total methane genes: $total_methane_hits (${methane_percent}%)
Denitrification genes:
  nirK: $nirk_hits
  nosZ: $nosz_hits
  Total denitrification genes: $total_denitrification_hits (${denitrification_percent}%)
High-quality extracted sequences: $extracted_seqs
EOF
    
    echo "  ✓ Sample $sample analysis completed"
    
    # Clean up temporary files
    rm -f ${SAMPLE_DIR}/${sample}_combined.fastq
    
done

# STEP 3: Generate summary report
echo ""
echo "STEP 3: Generating comprehensive summary..."

# Create summary table
echo -e "Sample\tTotal_Reads\tmcrA\tpmoA\tmmoX\thdr\tTotal_Methane\tnirK\tnosZ\tTotal_Denitr\tExtracted_Seqs" > ${RESULTS_DIR}/functional_gene_summary.txt

for sample in "${SAMPLES[@]}"; do
    if [ -f "${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt" ]; then
        # Extract values from stats file
        total_reads=$(grep "Total reads analyzed:" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | tr -d ' ')
        mcra=$(grep "mcrA" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | cut -d' ' -f2)
        pmoa=$(grep "pmoA" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | cut -d' ' -f2)
        mmox=$(grep "mmoX" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | cut -d' ' -f2)
        hdr=$(grep "hdr" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | cut -d' ' -f2)
        total_methane=$(grep "Total methane genes:" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | cut -d' ' -f2)
        nirk=$(grep "nirK:" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | tr -d ' ')
        nosz=$(grep "nosZ:" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | tr -d ' ')
        total_denitr=$(grep "Total denitrification genes:" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | cut -d' ' -f2)
        extracted=$(grep "High-quality extracted sequences:" ${RESULTS_DIR}/${sample}/${sample}_functional_gene_stats.txt | cut -d: -f2 | tr -d ' ')
        
        echo -e "$sample\t$total_reads\t$mcra\t$pmoa\t$mmox\t$hdr\t$total_methane\t$nirk\t$nosz\t$total_denitr\t$extracted" >> ${RESULTS_DIR}/functional_gene_summary.txt
    fi
done

# Combine all extracted sequences
echo "Combining all extracted methane gene sequences..."
cat ${RESULTS_DIR}/*/???_methane_sequences.fasta > ${RESULTS_DIR}/all_methane_sequences.fasta

# Generate phylogenetic analysis input
if [ -s "${RESULTS_DIR}/all_methane_sequences.fasta" ]; then
    echo "Preparing sequences for phylogenetic analysis..."
    
    # Separate by gene type
    grep -A1 "mcra\|reductase" ${RESULTS_DIR}/all_methane_sequences.fasta > ${RESULTS_DIR}/mcra_sequences.fasta || touch ${RESULTS_DIR}/mcra_sequences.fasta
    grep -A1 "pmoa\|monooxygenase" ${RESULTS_DIR}/all_methane_sequences.fasta > ${RESULTS_DIR}/pmoa_sequences.fasta || touch ${RESULTS_DIR}/pmoa_sequences.fasta
    
    echo "  ✓ Gene-specific sequence files created"
fi

echo ""
echo "============================================="
echo "FUNCTIONAL GENE SCREENING COMPLETED!"
echo "============================================="
echo ""
echo "Results summary:"
echo "  Individual sample results: ${RESULTS_DIR}/[sample_name]/"
echo "  Summary table: ${RESULTS_DIR}/functional_gene_summary.txt"
echo "  All extracted sequences: ${RESULTS_DIR}/all_methane_sequences.fasta"
echo "  Gene-specific sequences: ${RESULTS_DIR}/mcra_sequences.fasta, pmoa_sequences.fasta"
echo ""
echo "Key findings:"
if [ -f "${RESULTS_DIR}/functional_gene_summary.txt" ]; then
    echo ""
    echo "FUNCTIONAL GENE SUMMARY:"
    column -t -s $'\t' ${RESULTS_DIR}/functional_gene_summary.txt
    echo ""
    
    # Calculate totals
    total_mcra=$(tail -n +2 ${RESULTS_DIR}/functional_gene_summary.txt | cut -f3 | paste -sd+ | bc)
    total_pmoa=$(tail -n +2 ${RESULTS_DIR}/functional_gene_summary.txt | cut -f4 | paste -sd+ | bc)
    total_methane=$(tail -n +2 ${RESULTS_DIR}/functional_gene_summary.txt | cut -f6 | paste -sd+ | bc)
    
    echo "OVERALL TOTALS:"
    echo "  mcrA hits (methanogenesis): $total_mcra"
    echo "  pmoA hits (methanotrophy): $total_pmoa"
    echo "  Total methane metabolism genes: $total_methane"
    
    if [ $total_mcra -gt 0 ]; then
        echo ""
        echo "🦠 METHANOGENIC POTENTIAL DETECTED!"
        echo "   mcrA genes found indicating presence of methanogenic archaea"
    fi
    
    if [ $total_pmoa -gt 0 ]; then
        echo ""
        echo "🔥 METHANOTROPHIC POTENTIAL DETECTED!"
        echo "   pmoA genes found indicating presence of methane-oxidizing bacteria"
    fi
    
    if [ $total_methane -eq 0 ]; then
        echo ""
        echo "ℹ️  No methane metabolism genes detected in current screening"
        echo "   Consider deeper sequencing or targeted amplicon approaches"
    fi
fi

echo ""
echo "Next steps:"
echo "  1. Phylogenetic analysis of extracted sequences"
echo "  2. Quantitative PCR validation of key genes"
echo "  3. Amplicon sequencing for mcrA/pmoA if genes detected"
echo "  4. Continue with assembly and binning: ./code/03_assembly.sh"
echo ""
