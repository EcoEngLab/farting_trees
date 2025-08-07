#!/bin/bash

# Targeted analysis for methanogens and methanotrophs
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
ANNOTATION_DIR="../results/04_annotation"
BINNING_DIR="../results/03_binning"
RESULTS_DIR="../results/05_taxonomy"
THREADS=8

# Activate annotation environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenome_annotation

echo "Starting targeted analysis for methanogens and methanotrophs..."

# Create output directories
mkdir -p ${RESULTS_DIR}/{methanogens,methanotrophs,functional_analysis}

# Key functional genes for methanogens
METHANOGEN_GENES=("mcrA" "mcrB" "mcrG" "fmdA" "ftr" "mtd" "mer")

# Key functional genes for methanotrophs  
METHANOTROPH_GENES=("mmoX" "mmoY" "mmoZ" "pmoA" "pmoB" "pmoC" "mxaF" "mxaI")

# Step 1: Search for methanogen-specific genes
echo "Searching for methanogen-specific genes..."

# Create HMM database for methanogen genes (if available)
# For now, use sequence-based search with known protein sequences
cat > ${RESULTS_DIR}/methanogens/search_methanogens.py << 'EOF'
#!/usr/bin/env python3
"""
Search for methanogen-specific genes in MAG annotations
"""
import os
import glob
import pandas as pd
from Bio import SeqIO
import subprocess

def search_methanogen_genes(annotation_dir, output_dir):
    """Search for methanogen-specific genes in Prokka annotations"""
    
    methanogen_keywords = [
        'methyl-coenzyme M reductase',
        'mcrA', 'mcrB', 'mcrG',
        'formylmethanofuran dehydrogenase',
        'methenyltetrahydrofolate cyclohydrolase',
        'methane monooxygenase',
        'tetrahydrofolate dehydrogenase',
        'formate dehydrogenase'
    ]
    
    results = []
    
    # Search through all Prokka annotation files
    for tsv_file in glob.glob(f"{annotation_dir}/prokka/*/*.tsv"):
        sample_bin = os.path.basename(os.path.dirname(tsv_file))
        
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            
            for keyword in methanogen_keywords:
                matches = df[df['product'].str.contains(keyword, case=False, na=False)]
                for _, row in matches.iterrows():
                    results.append({
                        'sample_bin': sample_bin,
                        'gene_id': row['locus_tag'],
                        'product': row['product'],
                        'gene_type': 'methanogen',
                        'keyword': keyword
                    })
        except Exception as e:
            print(f"Error processing {tsv_file}: {e}")
    
    # Save results
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(f"{output_dir}/methanogen_genes.csv", index=False)
        print(f"Found {len(results)} potential methanogen genes")
    else:
        print("No methanogen genes found")

def search_methanotroph_genes(annotation_dir, output_dir):
    """Search for methanotroph-specific genes in Prokka annotations"""
    
    methanotroph_keywords = [
        'methane monooxygenase',
        'mmoX', 'mmoY', 'mmoZ',
        'particulate methane monooxygenase',
        'pmoA', 'pmoB', 'pmoC',
        'methanol dehydrogenase',
        'mxaF', 'mxaI',
        'formaldehyde dehydrogenase',
        'formate dehydrogenase'
    ]
    
    results = []
    
    # Search through all Prokka annotation files
    for tsv_file in glob.glob(f"{annotation_dir}/prokka/*/*.tsv"):
        sample_bin = os.path.basename(os.path.dirname(tsv_file))
        
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            
            for keyword in methanotroph_keywords:
                matches = df[df['product'].str.contains(keyword, case=False, na=False)]
                for _, row in matches.iterrows():
                    results.append({
                        'sample_bin': sample_bin,
                        'gene_id': row['locus_tag'],
                        'product': row['product'],
                        'gene_type': 'methanotroph',
                        'keyword': keyword
                    })
        except Exception as e:
            print(f"Error processing {tsv_file}: {e}")
    
    # Save results
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(f"{output_dir}/methanotroph_genes.csv", index=False)
        print(f"Found {len(results)} potential methanotroph genes")
    else:
        print("No methanotroph genes found")

if __name__ == "__main__":
    import sys
    annotation_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    print("Searching for methanogen genes...")
    search_methanogen_genes(annotation_dir, f"{output_dir}/methanogens")
    
    print("Searching for methanotroph genes...")
    search_methanotroph_genes(annotation_dir, f"{output_dir}/methanotrophs")
EOF

# Run the Python search script
python3 ${RESULTS_DIR}/methanogens/search_methanogens.py ${ANNOTATION_DIR} ${RESULTS_DIR}

# Step 2: Taxonomic filtering for methanogens and methanotrophs
echo "Filtering MAGs by taxonomy..."

# Create taxonomy filter script
cat > ${RESULTS_DIR}/filter_by_taxonomy.py << 'EOF'
#!/usr/bin/env python3
"""
Filter MAGs based on taxonomic classification for methanogens and methanotrophs
"""
import pandas as pd
import glob
import os

def filter_methanogens(gtdbtk_file, output_dir):
    """Filter potential methanogen MAGs"""
    
    try:
        df = pd.read_csv(gtdbtk_file, sep='\t')
        
        # Known methanogen taxa
        methanogen_taxa = [
            'Methanobrevibacter',
            'Methanococcus',
            'Methanosarcina',
            'Methanothermobacter',
            'Methanospirillum',
            'Methanopyrus',
            'Methanocaldococcus',
            'Methanobacterium',
            'Methanocorpusculum'
        ]
        
        methanogens = df[df['classification'].str.contains('|'.join(methanogen_taxa), case=False, na=False)]
        
        if not methanogens.empty:
            methanogens.to_csv(f"{output_dir}/methanogens/potential_methanogen_mags.csv", index=False)
            print(f"Found {len(methanogens)} potential methanogen MAGs")
        else:
            print("No methanogen MAGs found by taxonomy")
            
    except Exception as e:
        print(f"Error processing GTDB-Tk results: {e}")

def filter_methanotrophs(gtdbtk_file, output_dir):
    """Filter potential methanotroph MAGs"""
    
    try:
        df = pd.read_csv(gtdbtk_file, sep='\t')
        
        # Known methanotroph taxa
        methanotroph_taxa = [
            'Methylococcus',
            'Methylosinus',
            'Methylocystis',
            'Methylobacter',
            'Methylomonas',
            'Methylomicrobium',
            'Methylocaldum',
            'Methylacidiphilum'
        ]
        
        methanotrophs = df[df['classification'].str.contains('|'.join(methanotroph_taxa), case=False, na=False)]
        
        if not methanotrophs.empty:
            methanotrophs.to_csv(f"{output_dir}/methanotrophs/potential_methanotroph_mags.csv", index=False)
            print(f"Found {len(methanotrophs)} potential methanotroph MAGs")
        else:
            print("No methanotroph MAGs found by taxonomy")
            
    except Exception as e:
        print(f"Error processing GTDB-Tk results: {e}")

if __name__ == "__main__":
    import sys
    gtdbtk_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    if os.path.exists(gtdbtk_file):
        filter_methanogens(gtdbtk_file, output_dir)
        filter_methanotrophs(gtdbtk_file, output_dir)
    else:
        print(f"GTDB-Tk results file not found: {gtdbtk_file}")
EOF

# Run taxonomic filtering if GTDB-Tk results available
if [ -f "${ANNOTATION_DIR}/gtdbtk/output/gtdbtk.bac120.summary.tsv" ]; then
    python3 ${RESULTS_DIR}/filter_by_taxonomy.py \
            ${ANNOTATION_DIR}/gtdbtk/output/gtdbtk.bac120.summary.tsv \
            ${RESULTS_DIR}
fi

# Step 3: Generate comprehensive report
echo "Generating comprehensive report..."

cat > ${RESULTS_DIR}/methanogen_methanotroph_report.md << 'EOF'
# Methanogen and Methanotroph Analysis Report

## Overview
This report summarizes the identification of methanogens and methanotrophs in the tree microbiome samples.

## Methodology
1. **Functional gene search**: Searched for key genes in methanogenesis and methanotrophy pathways
2. **Taxonomic classification**: Used GTDB-Tk to identify known methanogen/methanotroph taxa
3. **Quality assessment**: Evaluated MAG completeness and contamination

## Key Functional Genes

### Methanogenesis Pathway
- **mcrA, mcrB, mcrG**: Methyl-coenzyme M reductase (terminal step)
- **fmdA**: Formylmethanofuran dehydrogenase
- **ftr**: Formyltetrahydrofolate synthetase
- **mtd**: Methylenetetrahydrofolate dehydrogenase

### Methanotrophy Pathway
- **mmoX, mmoY, mmoZ**: Soluble methane monooxygenase
- **pmoA, pmoB, pmoC**: Particulate methane monooxygenase
- **mxaF, mxaI**: Methanol dehydrogenase

## Results Summary
EOF

# Add functional gene results to report
if [ -f "${RESULTS_DIR}/methanogens/methanogen_genes.csv" ]; then
    echo "### Methanogen Genes Found" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    head -20 "${RESULTS_DIR}/methanogens/methanogen_genes.csv" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
fi

if [ -f "${RESULTS_DIR}/methanotrophs/methanotroph_genes.csv" ]; then
    echo "### Methanotroph Genes Found" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    head -20 "${RESULTS_DIR}/methanotrophs/methanotroph_genes.csv" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
fi

# Add taxonomic results to report
if [ -f "${RESULTS_DIR}/methanogens/potential_methanogen_mags.csv" ]; then
    echo "### Potential Methanogen MAGs" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    cat "${RESULTS_DIR}/methanogens/potential_methanogen_mags.csv" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
fi

if [ -f "${RESULTS_DIR}/methanotrophs/potential_methanotroph_mags.csv" ]; then
    echo "### Potential Methanotroph MAGs" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    cat "${RESULTS_DIR}/methanotrophs/potential_methanotroph_mags.csv" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
fi

echo "Targeted analysis completed!"
echo "Results saved in: ${RESULTS_DIR}"
echo "Report: ${RESULTS_DIR}/methanogen_methanotroph_report.md"
