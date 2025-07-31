#!/bin/bash

# Metagenomic binning using MetaBAT2, MaxBin2, and VAMB
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
ASSEMBLY_DIR="../results/02_assembly"
MAPPING_DIR="../results/02_assembly/mapping"
RESULTS_DIR="../results/03_binning"
THREADS=8

# Activate binning environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_binning

echo "Starting metagenomic binning..."

# Create output directories
mkdir -p ${RESULTS_DIR}/{metabat2,maxbin2,vamb,das_tool,checkm}

# Sample list
SAMPLES=("53394_A15" "53395_A16" "53396_B6" "53397_B10" "53398_A16S" "53399_B10S")

# Step 1: Prepare depth files for MetaBAT2
echo "Preparing depth files for MetaBAT2..."
for sample in "${SAMPLES[@]}"; do
    echo "Calculating depth for sample: $sample"
    jgi_summarize_bam_contig_depths --outputDepth ${RESULTS_DIR}/${sample}_depth.txt \
                                   ${MAPPING_DIR}/${sample}_sorted.bam
done

# Step 2: MetaBAT2 binning
echo "Running MetaBAT2 binning..."
for sample in "${SAMPLES[@]}"; do
    echo "MetaBAT2 binning for sample: $sample"
    metabat2 -i ${ASSEMBLY_DIR}/megahit/${sample}_contigs.fa \
             -a ${RESULTS_DIR}/${sample}_depth.txt \
             -o ${RESULTS_DIR}/metabat2/${sample}_bin \
             --minContig 2500 \
             --numThreads ${THREADS}
done

# Step 3: MaxBin2 binning
echo "Running MaxBin2 binning..."
for sample in "${SAMPLES[@]}"; do
    echo "MaxBin2 binning for sample: $sample"
    
    # Create abundance file from depth
    cut -f1,3 ${RESULTS_DIR}/${sample}_depth.txt > ${RESULTS_DIR}/${sample}_abundance.txt
    
    run_MaxBin.pl -contig ${ASSEMBLY_DIR}/megahit/${sample}_contigs.fa \
                  -abund ${RESULTS_DIR}/${sample}_abundance.txt \
                  -out ${RESULTS_DIR}/maxbin2/${sample}_bin \
                  -thread ${THREADS} \
                  -min_contig_length 2500
done

# Step 4: VAMB binning (using co-assembly approach)
echo "Running VAMB binning on co-assembly..."

# Map all samples to co-assembly for VAMB
echo "Mapping all samples to co-assembly for VAMB..."
bwa index ${ASSEMBLY_DIR}/megahit/coassembly_contigs.fa

for sample in "${SAMPLES[@]}"; do
    echo "Mapping $sample to co-assembly"
    bwa mem -t ${THREADS} \
            ${ASSEMBLY_DIR}/megahit/coassembly_contigs.fa \
            ../results/01_qc/trimmed/${sample}_R1_trimmed.fastq.gz \
            ../results/01_qc/trimmed/${sample}_R2_trimmed.fastq.gz | \
    samtools view -bS -F 4 - | \
    samtools sort -o ${RESULTS_DIR}/vamb/${sample}_coassembly_sorted.bam
    
    samtools index ${RESULTS_DIR}/vamb/${sample}_coassembly_sorted.bam
done

# Create abundance matrix for VAMB
jgi_summarize_bam_contig_depths --outputDepth ${RESULTS_DIR}/vamb/coassembly_depth.txt \
                               ${RESULTS_DIR}/vamb/*_coassembly_sorted.bam

# Run VAMB
vamb --outdir ${RESULTS_DIR}/vamb/bins \
     --fasta ${ASSEMBLY_DIR}/megahit/coassembly_contigs.fa \
     --bamfiles ${RESULTS_DIR}/vamb/*_coassembly_sorted.bam \
     --mincontig 2500

# Step 5: DAS Tool for bin optimization
echo "Running DAS Tool for bin optimization..."
for sample in "${SAMPLES[@]}"; do
    echo "DAS Tool optimization for sample: $sample"
    
    # Prepare bin-to-contig files
    cd ${RESULTS_DIR}/metabat2
    ls ${sample}_bin*.fa | while read bin; do
        basename $bin .fa
        grep "^>" $bin | sed 's/>//'
    done > ${RESULTS_DIR}/das_tool/${sample}_metabat2_scaffolds2bin.txt
    
    cd ${RESULTS_DIR}/maxbin2
    ls ${sample}_bin*.fasta | while read bin; do
        basename $bin .fasta
        grep "^>" $bin | sed 's/>//'
    done > ${RESULTS_DIR}/das_tool/${sample}_maxbin2_scaffolds2bin.txt
    
    cd ${RESULTS_DIR}
    
    # Run DAS Tool
    DAS_Tool -i ${RESULTS_DIR}/das_tool/${sample}_metabat2_scaffolds2bin.txt,${RESULTS_DIR}/das_tool/${sample}_maxbin2_scaffolds2bin.txt \
             -l metabat2,maxbin2 \
             -c ${ASSEMBLY_DIR}/megahit/${sample}_contigs.fa \
             -o ${RESULTS_DIR}/das_tool/${sample}_DASToolBins \
             --threads ${THREADS} \
             --write_bins 1
done

# Step 6: Bin quality assessment with CheckM
echo "Running CheckM for bin quality assessment..."

# CheckM on MetaBAT2 bins
for sample in "${SAMPLES[@]}"; do
    if [ -d "${RESULTS_DIR}/metabat2" ] && [ "$(ls -A ${RESULTS_DIR}/metabat2/${sample}_bin*.fa 2>/dev/null)" ]; then
        echo "CheckM analysis for MetaBAT2 bins: $sample"
        checkm lineage_wf ${RESULTS_DIR}/metabat2 \
                         ${RESULTS_DIR}/checkm/${sample}_metabat2 \
                         -t ${THREADS} \
                         --file ${RESULTS_DIR}/checkm/${sample}_metabat2_quality.txt \
                         --tab_table \
                         -x fa
    fi
done

# CheckM on DAS Tool bins
for sample in "${SAMPLES[@]}"; do
    if [ -d "${RESULTS_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins" ]; then
        echo "CheckM analysis for DAS Tool bins: $sample"
        checkm lineage_wf ${RESULTS_DIR}/das_tool/${sample}_DASToolBins_DASTool_bins \
                         ${RESULTS_DIR}/checkm/${sample}_dastool \
                         -t ${THREADS} \
                         --file ${RESULTS_DIR}/checkm/${sample}_dastool_quality.txt \
                         --tab_table \
                         -x fa
    fi
done

echo "Binning completed!"
echo "Results saved in: ${RESULTS_DIR}"
echo "Check CheckM reports for bin quality:"
echo "  - MetaBAT2: ${RESULTS_DIR}/checkm/*_metabat2_quality.txt"
echo "  - DAS Tool: ${RESULTS_DIR}/checkm/*_dastool_quality.txt"
