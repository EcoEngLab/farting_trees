# Farting Trees: Metagenomic Analysis Pipeline

Trees' methanogens and methanotrophs sequencing data and analysis pipeline for assembly and MAG binning.

## Project Overview

This project analyzes shotgun metagenomic data from tree samples to:
1. Assemble reads into contigs using MEGAHIT and metaSPAdes
2. Bin contigs into Metagenome-Assembled Genomes (MAGs) using MetaBAT2, MaxBin2, and VAMB
3. Characterize microbial communities, particularly methanogens and methanotrophs
4. Perform functional annotation and taxonomic classification

## Data Description

**Samples:**
- 53394_A15: Tree sample A15 (~346MB R1, ~345MB R2)
- 53395_A16: Tree sample A16 (~780MB R1, ~790MB R2)
- 53396_B6: Tree sample B6 (~517MB R1, ~525MB R2)
- 53397_B10: Tree sample B10 (~802MB R1, ~763MB R2)
- 53398_A16S: Substrate/soil sample A16S (~677MB R1, ~697MB R2)
- 53399_B10S: Substrate/soil sample B10S (~598MB R1, ~625MB R2)

Total dataset: ~6.7GB of paired-end sequencing data

## Quick Start

```bash
# 1. Setup conda environments
./scripts/run_pipeline.sh setup

# 2. Run complete pipeline
./scripts/run_pipeline.sh all

# 3. View results
firefox results/01_qc/multiqc/multiqc_report.html
```

## Pipeline Overview

### 1. **Quality Control & Preprocessing**
   - FastQC for quality assessment
   - fastp for quality trimming and adapter removal
   - MultiQC for aggregated quality reports

### 2. **Host Removal**
   - **BWA mapping** to plant/tree reference genomes
   - **Removal of host DNA** (tree/plant sequences)
   - **Statistics generation** for contamination assessment
   - Supports multiple tree species (poplar, oak, pine, spruce, etc.)

### 3. **Assembly**
   - **MEGAHIT**: Fast, memory-efficient assembly for individual samples and co-assembly
   - **metaSPAdes**: High-quality assembly (alternative/comparison)
   - **QUAST**: Assembly quality statistics
   - **BWA + Samtools**: Read mapping for coverage calculation

### 4. **Binning**
   - **MetaBAT2**: Depth-based binning using coverage information
   - **MaxBin2**: Composition-based binning using tetranucleotide frequencies
   - **VAMB**: Deep learning-based binning for complex communities
   - **DAS Tool**: Bin optimization and dereplication
   - **CheckM**: Completeness and contamination assessment

### 5. **MAG Quality Assessment**
   - Completeness and contamination metrics
   - Quality classification (high/medium/low quality)
   - Bin statistics and comparisons

### 6. **Taxonomic Classification**
   - **GTDB-Tk**: Standardized taxonomic classification
   - Phylogenetic placement of MAGs

### 7. **Functional Annotation**
   - **Prokka**: Gene prediction and basic annotation
   - **eggNOG-mapper**: Functional annotation and pathway mapping (alternative methods available)
   - Targeted analysis for methanogenic and methanotrophic pathways

### 8. **Targeted Analysis**
   - Search for key functional genes:
     - **Methanogenesis**: mcrA, mcrB, mcrG (methyl-coenzyme M reductase)
     - **Methanotrophy**: pmoA, mmoX, mxaF (methane/methanol oxidation)
   - Taxonomic filtering for known methanogens and methanotrophs
   - Pathway reconstruction and analysis

## File Structure

```
farting_trees/
├── README.md                 # This file
├── USAGE.md                 # Detailed usage instructions
├── requirements.txt         # Python dependencies
├── data/                    # Raw sequencing data
│   ├── 53394_R1_A15.fastq.gz
│   ├── 53394_R2_A15.fastq.gz
│   └── ... (other samples)
├── scripts/                 # Analysis scripts
│   ├── run_pipeline.sh      # Master pipeline script
│   ├── 01_setup_environment.sh
│   ├── 02_quality_control.sh
│   ├── 03_assembly.sh
│   ├── 04_binning.sh
│   ├── 05_annotation.sh
│   ├── 06_target_analysis.sh
│   └── advanced_analysis.py
├── config/                  # Configuration files
│   ├── analysis_config.yaml
│   └── sample_list.txt
└── results/                 # Analysis results (created during run)
    ├── 01_qc/              # Quality control
    ├── 02_assembly/        # Assembly results
    ├── 03_binning/         # Binning and MAGs
    ├── 04_annotation/      # Functional annotation
    └── 05_taxonomy/        # Targeted analysis
```

## Key Features

### Advanced Binning Strategy
- **Multi-algorithm approach**: Combines depth-based, composition-based, and ML-based binning
- **Quality optimization**: DAS Tool integration for best bin selection
- **Comprehensive assessment**: CheckM for rigorous quality control

### Targeted Microbial Analysis
- **Methanogen detection**: Focus on mcrA gene and known methanogenic taxa
- **Methanotroph detection**: Focus on pmoA/mmoX genes and methanotrophic taxa
- **Pathway analysis**: Reconstruction of methanogenic and methanotrophic pathways

### Scalable Architecture
- **Modular design**: Run individual steps or complete pipeline
- **Resource flexible**: Configurable for different computing environments
- **Quality focused**: Multiple QC checkpoints and comprehensive reporting

## Requirements

### Software Dependencies
- conda/mamba package manager
- Bioconda and conda-forge channels

### Hardware Requirements
- **Minimum**: 8 CPU cores, 16GB RAM, 100GB storage
- **Recommended**: 16+ CPU cores, 32GB+ RAM, 200GB+ storage
- **Runtime**: 6-24 hours depending on resources

### Installed via Conda
- FastQC, MultiQC, fastp, Trimmomatic
- MEGAHIT, metaSPAdes, QUAST
- MetaBAT2, MaxBin2, VAMB, DAS Tool, CheckM
- Prokka, GTDB-Tk, eggNOG-mapper
- BWA, Samtools, Bowtie2

## Expected Outputs

### Quality Reports
- MultiQC HTML report with comprehensive QC metrics
- Assembly statistics and quality assessments

### MAGs (Metagenome-Assembled Genomes)
- High-quality bins (>90% complete, <5% contamination)
- Medium-quality bins (>50% complete, <10% contamination)
- Taxonomic classification for each MAG

### Functional Analysis
- Gene predictions and annotations for all MAGs
- Functional pathway analysis
- Targeted identification of methanogens and methanotrophs

### Visualization
- Interactive plots for MAG quality assessment
- Taxonomic distribution visualizations
- Functional gene distribution analysis

## Getting Started

See [USAGE.md](USAGE.md) for detailed instructions on:
- Setting up the analysis environment
- Running the complete pipeline
- Interpreting results
- Troubleshooting common issues
- Customizing the analysis

## Citation

If you use this pipeline, please cite the relevant tools:
- MEGAHIT, metaSPAdes (assembly)
- MetaBAT2, MaxBin2, VAMB, DAS Tool (binning)
- CheckM (quality assessment)
- GTDB-Tk (taxonomic classification)
- Prokka, eggNOG-mapper (functional annotation)

## Contact

For questions about this analysis pipeline, please refer to the documentation or create an issue in the repository.
