# Metagenomic Binning Pipeline

A comprehensive pipeline for quality control, assembly, and binning of metagenomic sequencing data. This pipeline processes paired-end FASTQ files through quality control, assembly, read mapping, and multiple binning approaches with quality assessment.

## Overview

This pipeline processes metagenomic samples through the following steps:
1. **Quality Control & Primer Removal** - Remove adapters and primers using fastp
2. **Assembly** - Individual sample assembly using metaSPAdes
3. **Assembly Quality Assessment** - Evaluate assemblies with QUAST
4. **Read Mapping** - Map reads back to assemblies for coverage calculation
5. **Binning** - Multiple binning approaches (MetaBAT2, MaxBin2, VAMB) with DAS Tool optimization
6. **Quality Assessment** - Evaluate bin quality using CheckM2

## Directory Structure

```
farting_trees/
├── code/                           # All pipeline scripts
│   ├── 01_qc.sh                   # Quality control and primer removal
│   ├── 02_assembly.sh             # metaSPAdes assembly
│   ├── 03_quast.sh                # Assembly quality assessment
│   ├── 04_mapping.sh              # Read mapping for coverage
│   ├── 05_binning.sh              # Complete binning pipeline
│   └── setup_binning_envs.sh      # Environment setup script
├── data/                          # Input FASTQ files
│   ├── sample1_R1_*.fastq.gz
│   ├── sample1_R2_*.fastq.gz
│   └── ... (more samples)
└── results/                       # All pipeline outputs
    ├── 01_qc/
    ├── 02_assembly/
    ├── 03_quast/
    ├── 04_mapping/
    └── 05_binning/
```

## Prerequisites

- **Conda/Miniconda** installed and configured
- **Linux/Unix system** (tested on Linux)
- Sufficient disk space (recommend >50GB for typical datasets)
- At least 8 CPU cores recommended

## Quick Start

### 1. Setup Environments

First, set up all required conda environments:

```bash
# Setup all binning environments at once
bash code/setup_binning_envs.sh

```

This creates the following conda environments:
- `metagenome_qc`: FastP, FastQC, MultiQC
- `spades_assembly`: SPAdes assembler
- `quast`: Assembly quality assessment
- `bwa_env`: BWA and SAMtools for mapping
- `metabat218`: MetaBAT2 binning
- `maxbin2_env`: MaxBin2 binning
- `vamb_env`: VAMB binning
- `dastool_env`: DAS Tool optimization
- `checkm2_env`: CheckM2 quality assessment

### 2. Prepare Your Data

Place your paired-end FASTQ files in the `data/` directory:
```
data/
├── sample1_R1_001.fastq.gz
├── sample1_R2_001.fastq.gz
├── sample2_R1_001.fastq.gz
├── sample2_R2_001.fastq.gz
└── ...
```

### 3. Run the Complete Pipeline

```bash
# Run all steps sequentially
bash code/01_qc.sh
bash code/02_assembly.sh
bash code/03_quast.sh
bash code/04_mapping.sh
bash code/05_binning.sh all
```

## Detailed Pipeline Steps

### setup_binning_envs.sh — Environment Setup
- **Purpose:** Installs all required conda environments for the pipeline (QC, assembly, mapping, binning, etc.).
- **Outputs:** All conda environments listed in the Quick Start section.

### 01_qc.sh — Quality Control & Primer Removal
- **Purpose:** Cleans raw FASTQ files by removing adapters.
- **Main tools:** fastp (adapter/primer trimming, quality filtering), FastQC (quality reports), MultiQC (summary report)
- **Inputs:** Paired-end FASTQ files from `data/`
- **Outputs:** Cleaned FASTQ files in `results/01_qc/`, quality reports in `results/01_qc/`
- **Key steps:**
  1. Run fastp on each sample to trim adapters/primers and filter low-quality reads.
  2. Generate FastQC reports for raw and cleaned reads.
  3. Summarize QC results with MultiQC.

### 02_assembly.sh — Assembly
- **Purpose:** Assembles cleaned reads into contigs for each sample.
- **Main tools:** metaSPAdes (metagenome assembler)
- **Inputs:** Cleaned FASTQ files from `results/01_qc/`
- **Outputs:** Assembled contigs in `results/02_assembly/`
- **Key steps:**
  1. Run metaSPAdes for each sample using paired-end cleaned reads.
  2. Output contigs and assembly logs for each sample.

### 03_quast.sh — Assembly Quality Assessment
- **Purpose:** Evaluates the quality of each assembly.
- **Main tools:** QUAST
- **Inputs:** Assembled contigs from `results/02_assembly/`
- **Outputs:** QUAST reports in `results/03_quast/`
- **Key steps:**
  1. Run QUAST on each assembly to calculate N50, total length, GC content, etc.
  2. Summarize assembly statistics for all samples.

### 04_mapping.sh — Read Mapping for Coverage
- **Purpose:** Maps cleaned reads back to assemblies to estimate coverage for binning.
- **Main tools:** BWA (alignment), SAMtools (BAM processing)
- **Inputs:** Cleaned FASTQ files (`results/01_qc/`), assemblies (`results/02_assembly/`)
- **Outputs:** BAM files and coverage tables in `results/04_mapping/`
- **Key steps:**
  1. Index each assembly with BWA.
  2. Align cleaned reads to assemblies, convert SAM to BAM, sort and index BAM files.
  3. Calculate per-contig coverage for each sample.

### 05_binning.sh — Binning and Bin Optimization
- **Purpose:** Groups contigs into bins (putative genomes) using multiple algorithms and refines bins.
- **Main tools:** MetaBAT2, MaxBin2, DAS Tool, checkm2
- **Inputs:** Assemblies (`results/02_assembly/`), coverage tables (`results/04_mapping/`)
- **Outputs:** Binning results in `results/05_binning/` (one folder per tool), DAS Tool optimized bins, checkm2 report
- **Key steps:**
  1. Run MetaBAT2 and MaxBin2 on each assembly using coverage information.
  2. Collect bins from all tools.
  3. Run DAS Tool to select and optimize the best bins from all methods.
  4. Use checkm2 to assess the quality of bins.

---

For more details on each script, see the comments at the top of each `.sh` file in the `code/` directory.
