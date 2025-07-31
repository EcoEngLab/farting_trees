# Metagenomic Analysis Pipeline - Usage Guide

## Quick Start

1. **Setup Environment**
   ```bash
   cd /home/jiayi-chen/Documents/farting_trees
   ./scripts/run_pipeline.sh setup
   ```

2. **Run Complete Pipeline**
   ```bash
   ./scripts/run_pipeline.sh all
   ```

3. **Run Individual Steps**
   ```bash
   ./scripts/run_pipeline.sh qc        # Quality control only
   ./scripts/run_pipeline.sh assembly  # Assembly only
   ./scripts/run_pipeline.sh binning   # Binning only
   ```

## Step-by-Step Instructions

### 1. Environment Setup
The pipeline requires several bioinformatics tools installed via conda:
```bash
./scripts/01_setup_environment.sh
```

This creates four conda environments:
- `metagenome_qc`: FastQC, MultiQC, fastp, Trimmomatic
- `metagenome_assembly`: MEGAHIT, metaSPAdes, QUAST, BWA, Samtools
- `metagenome_binning`: MetaBAT2, MaxBin2, VAMB, CheckM, DAS Tool
- `metagenome_annotation`: Prokka, GTDB-Tk, eggNOG-mapper

### 2. Quality Control
```bash
./scripts/02_quality_control.sh
```
- Runs FastQC on raw reads
- Performs quality trimming with fastp
- Generates MultiQC report

### 3. Assembly
```bash
./scripts/03_assembly.sh
```
- Individual sample assembly with MEGAHIT
- Co-assembly of all samples
- Alternative assembly with metaSPAdes
- Quality assessment with QUAST
- Read mapping for coverage calculation

### 4. Binning
```bash
./scripts/04_binning.sh
```
- MetaBAT2 binning (depth-based)
- MaxBin2 binning (composition-based)
- VAMB binning (deep learning)
- DAS Tool optimization
- CheckM quality assessment

### 5. Annotation
```bash
./scripts/05_annotation.sh
```
- Gene prediction with Prokka
- Taxonomic classification with GTDB-Tk
- Functional annotation with eggNOG-mapper

### 6. Targeted Analysis
```bash
./scripts/06_target_analysis.sh
```
- Search for methanogen-specific genes
- Search for methanotroph-specific genes
- Taxonomic filtering
- Generate targeted reports

## Advanced Analysis

### Python Analysis Script
```bash
# Install Python dependencies
pip install -r requirements.txt

# Run advanced analysis
python scripts/advanced_analysis.py --analysis all
```

This generates:
- Interactive quality plots
- Taxonomic overview plots
- Functional gene analysis
- Comprehensive summary report

## Expected Results Structure

```
results/
├── 01_qc/                    # Quality control results
│   ├── fastqc_raw/          # Raw read QC
│   ├── fastqc_clean/        # Trimmed read QC
│   ├── trimmed/             # Quality-trimmed reads
│   └── multiqc/             # Aggregated QC report
├── 02_assembly/             # Assembly results
│   ├── megahit/             # MEGAHIT assemblies
│   ├── spades/              # metaSPAdes assemblies
│   ├── quast/               # Assembly statistics
│   └── mapping/             # Read mapping files
├── 03_binning/              # Binning results
│   ├── metabat2/            # MetaBAT2 bins
│   ├── maxbin2/             # MaxBin2 bins
│   ├── vamb/                # VAMB bins
│   ├── das_tool/            # Optimized bins
│   └── checkm/              # Quality assessment
├── 04_annotation/           # Annotation results
│   ├── prokka/              # Gene predictions
│   ├── gtdbtk/              # Taxonomic classification
│   └── eggnog/              # Functional annotation
└── 05_taxonomy/             # Targeted analysis
    ├── methanogens/         # Methanogen analysis
    ├── methanotrophs/       # Methanotroph analysis
    └── functional_analysis/ # Functional gene analysis
```

## Key Output Files

### Quality Reports
- `results/01_qc/multiqc/multiqc_report.html` - Comprehensive QC report

### Assembly Statistics
- `results/02_assembly/quast/megahit/report.html` - MEGAHIT assembly stats
- `results/02_assembly/quast/spades/report.html` - metaSPAdes assembly stats

### MAG Quality
- `results/03_binning/checkm/*_quality.txt` - Completeness and contamination

### Taxonomic Classification
- `results/04_annotation/gtdbtk/output/gtdbtk.bac120.summary.tsv` - Taxonomy

### Target Organisms
- `results/05_taxonomy/methanogen_methanotroph_report.md` - Key findings

## Troubleshooting

### Common Issues

1. **Memory Issues**
   - Reduce thread count in scripts
   - Use MEGAHIT instead of metaSPAdes for assembly
   - Process samples individually instead of co-assembly

2. **Missing Dependencies**
   - Ensure conda environments are properly activated
   - Check conda channel priorities: `conda config --show channels`

3. **Empty Results**
   - Check input file paths and naming
   - Verify data quality with FastQC
   - Check log files for errors

### Resource Requirements

- **Minimum**: 8 CPU cores, 16GB RAM, 100GB disk space
- **Recommended**: 16+ CPU cores, 32GB+ RAM, 200GB+ disk space
- **Runtime**: 6-24 hours depending on data size and resources

## Customization

### Modifying Parameters
Edit configuration files in `config/`:
- `analysis_config.yaml` - Analysis parameters
- `sample_list.txt` - Sample names

### Adding New Samples
1. Add FASTQ files to `data/` directory
2. Update `config/sample_list.txt`
3. Modify sample lists in scripts if needed

### Different Binning Strategies
The pipeline uses three complementary binning approaches:
- **MetaBAT2**: Depth-based binning (good for abundant organisms)
- **MaxBin2**: Composition-based binning (good for low-abundance organisms)
- **VAMB**: Deep learning approach (good for complex communities)
- **DAS Tool**: Combines results for optimal bins
