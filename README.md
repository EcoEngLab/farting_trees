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
│   ├── setup_binning_envs.sh      # Environment setup script
│   └── install_checkm2.sh         # CheckM2 installation script
├── data/                          # Input FASTQ files
│   ├── 53394_R1_*.fastq.gz
│   ├── 53394_R2_*.fastq.gz
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

# Or setup CheckM2 only (faster option)
bash code/install_checkm2.sh
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
├── 53394_R1_001.fastq.gz
├── 53394_R2_001.fastq.gz
├── 53395_R1_001.fastq.gz
├── 53395_R2_001.fastq.gz
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
