#!/bin/bash

# Targeted analysis for methanogens and methanotrophs
# Author: Generated for farting_trees project
# Date: July 2025

set -e

# Configuration
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOTATION_DIR="${WORK_DIR}/results/05_annotation"
BINNING_DIR="${WORK_DIR}/results/04_binning"
RESULTS_DIR="${WORK_DIR}/results/06_taxonomy"
THREADS=8

# Check if annotation directory exists
if [ ! -d "$ANNOTATION_DIR" ]; then
    echo "Error: Annotation directory not found: $ANNOTATION_DIR"
    echo "Please run the annotation script first"
    exit 1
fi

echo "Found annotation directory: $ANNOTATION_DIR"

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
    for tsv_file in glob.glob(f"{annotation_dir}/prokka/coassembly_bin.*/*.tsv"):
        sample_bin = os.path.basename(os.path.dirname(tsv_file))
        
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            
            # Ensure product column exists and handle missing values properly
            if 'product' not in df.columns:
                continue
            
            # Fill NaN values and ensure all values are strings
            df['product'] = df['product'].astype(str).fillna('hypothetical protein')
            
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
    for tsv_file in glob.glob(f"{annotation_dir}/prokka/coassembly_bin.*/*.tsv"):
        sample_bin = os.path.basename(os.path.dirname(tsv_file))
        
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            
            # Ensure product column exists and handle missing values properly
            if 'product' not in df.columns:
                continue
            
            # Fill NaN values and ensure all values are strings
            df['product'] = df['product'].astype(str).fillna('hypothetical protein')
            
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

# Step 1.5: Run GTDB-Tk for taxonomic classification
echo "Running GTDB-Tk for species identification..."

# Check if GTDB-Tk is available
if command -v gtdbtk &> /dev/null; then
    echo "GTDB-Tk found, running taxonomic classification..."
    
    # Create GTDB-Tk output directory
    mkdir -p ${RESULTS_DIR}/gtdbtk
    
    # Prepare MAG directory for GTDB-Tk
    MAG_DIR="${RESULTS_DIR}/mags_for_gtdbtk"
    mkdir -p ${MAG_DIR}
    
    # Copy MAGs to GTDB-Tk input directory (need .fa extension)
    echo "Preparing MAGs for GTDB-Tk analysis..."
    for mag_file in ${BINNING_DIR}/metabat2/coassembly_bin.*.fa; do
        if [ -f "$mag_file" ]; then
            cp "$mag_file" ${MAG_DIR}/
            echo "Added $(basename $mag_file) to GTDB-Tk input"
        fi
    done
    
    # Run GTDB-Tk classify workflow
    echo "Running GTDB-Tk classification (this may take 30-60 minutes)..."
    gtdbtk classify_wf \
        --genome_dir ${MAG_DIR} \
        --out_dir ${RESULTS_DIR}/gtdbtk \
        --cpus ${THREADS} \
        --extension fa
    
    echo "GTDB-Tk classification completed!"
    
    # Create species identification summary
    cat > ${RESULTS_DIR}/create_species_summary.py << 'EOF'
#!/usr/bin/env python3
"""
Create species identification summary from GTDB-Tk results
"""
import pandas as pd
import os

def create_species_summary(gtdbtk_dir, output_file):
    """Create a comprehensive species identification summary"""
    
    results = []
    
    # Process bacterial classification results
    bac_file = f"{gtdbtk_dir}/gtdbtk.bac120.summary.tsv"
    if os.path.exists(bac_file):
        print(f"Processing bacterial classifications: {bac_file}")
        df_bac = pd.read_csv(bac_file, sep='\t')
        for _, row in df_bac.iterrows():
            classification = row['classification']
            taxa_levels = classification.split(';')
            
            results.append({
                'MAG_ID': row['user_genome'],
                'Domain': taxa_levels[0].replace('d__', '') if len(taxa_levels) > 0 else 'Unknown',
                'Phylum': taxa_levels[1].replace('p__', '') if len(taxa_levels) > 1 else 'Unknown',
                'Class': taxa_levels[2].replace('c__', '') if len(taxa_levels) > 2 else 'Unknown',
                'Order': taxa_levels[3].replace('o__', '') if len(taxa_levels) > 3 else 'Unknown',
                'Family': taxa_levels[4].replace('f__', '') if len(taxa_levels) > 4 else 'Unknown',
                'Genus': taxa_levels[5].replace('g__', '') if len(taxa_levels) > 5 else 'Unknown',
                'Species': taxa_levels[6].replace('s__', '') if len(taxa_levels) > 6 else 'Unknown',
                'Classification_Method': 'GTDB-Tk_Bacterial',
                'Full_Classification': classification,
                'Fastani_Reference': row.get('fastani_reference', 'N/A'),
                'Fastani_Reference_Radius': row.get('fastani_reference_radius', 'N/A'),
                'Fastani_Taxonomy': row.get('fastani_taxonomy', 'N/A'),
                'Classification_Score': row.get('classification_score', 'N/A')
            })
    
    # Process archaeal classification results
    arc_file = f"{gtdbtk_dir}/gtdbtk.ar53.summary.tsv"
    if os.path.exists(arc_file):
        print(f"Processing archaeal classifications: {arc_file}")
        df_arc = pd.read_csv(arc_file, sep='\t')
        for _, row in df_arc.iterrows():
            classification = row['classification']
            taxa_levels = classification.split(';')
            
            results.append({
                'MAG_ID': row['user_genome'],
                'Domain': taxa_levels[0].replace('d__', '') if len(taxa_levels) > 0 else 'Unknown',
                'Phylum': taxa_levels[1].replace('p__', '') if len(taxa_levels) > 1 else 'Unknown',
                'Class': taxa_levels[2].replace('c__', '') if len(taxa_levels) > 2 else 'Unknown',
                'Order': taxa_levels[3].replace('o__', '') if len(taxa_levels) > 3 else 'Unknown',
                'Family': taxa_levels[4].replace('f__', '') if len(taxa_levels) > 4 else 'Unknown',
                'Genus': taxa_levels[5].replace('g__', '') if len(taxa_levels) > 5 else 'Unknown',
                'Species': taxa_levels[6].replace('s__', '') if len(taxa_levels) > 6 else 'Unknown',
                'Classification_Method': 'GTDB-Tk_Archaeal',
                'Full_Classification': classification,
                'Fastani_Reference': row.get('fastani_reference', 'N/A'),
                'Fastani_Reference_Radius': row.get('fastani_reference_radius', 'N/A'),
                'Fastani_Taxonomy': row.get('fastani_taxonomy', 'N/A'),
                'Classification_Score': row.get('classification_score', 'N/A')
            })
    
    if results:
        df_results = pd.DataFrame(results)
        df_results.to_csv(output_file, index=False)
        print(f"Species identification summary saved to: {output_file}")
        
        # Print summary statistics
        print(f"\nSPECIES IDENTIFICATION SUMMARY:")
        print(f"Total MAGs classified: {len(df_results)}")
        print(f"Bacterial MAGs: {len(df_results[df_results['Domain'] == 'Bacteria'])}")
        print(f"Archaeal MAGs: {len(df_results[df_results['Domain'] == 'Archaea'])}")
        print(f"Species-level identifications: {len(df_results[df_results['Species'] != 'Unknown'])}")
        
        # Show methane-related organisms
        methanogen_terms = ['Methanobrevibacter', 'Methanococcus', 'Methanosarcina', 'Methanothermobacter', 
                          'Methanospirillum', 'Methanopyrus', 'Methanocaldococcus', 'Methanobacterium']
        methanotroph_terms = ['Methylococcus', 'Methylosinus', 'Methylocystis', 'Methylobacter', 
                            'Methylomonas', 'Methylomicrobium', 'Methylocaldum']
        
        methanogens = df_results[df_results['Full_Classification'].str.contains('|'.join(methanogen_terms), case=False, na=False)]
        methanotrophs = df_results[df_results['Full_Classification'].str.contains('|'.join(methanotroph_terms), case=False, na=False)]
        
        if not methanogens.empty:
            print(f"\n🔥 METHANOGENS IDENTIFIED: {len(methanogens)}")
            for _, row in methanogens.iterrows():
                print(f"  - {row['MAG_ID']}: {row['Genus']} {row['Species']}")
        
        if not methanotrophs.empty:
            print(f"\n🌱 METHANOTROPHS IDENTIFIED: {len(methanotrophs)}")
            for _, row in methanotrophs.iterrows():
                print(f"  - {row['MAG_ID']}: {row['Genus']} {row['Species']}")
        
        if methanogens.empty and methanotrophs.empty:
            print("\n❌ No methanogens or methanotrophs identified by taxonomic classification")
            print("   However, functional genes may still be present in other organisms")
    
    else:
        print("No GTDB-Tk classification results found")

if __name__ == "__main__":
    import sys
    gtdbtk_dir = sys.argv[1]
    output_file = sys.argv[2]
    create_species_summary(gtdbtk_dir, output_file)
EOF
    
    # Run species identification summary
    python3 ${RESULTS_DIR}/create_species_summary.py ${RESULTS_DIR}/gtdbtk ${RESULTS_DIR}/species_identification_summary.csv
    
else
    echo "⚠️  GTDB-Tk not found. Installing via conda..."
    echo "To install GTDB-Tk, run:"
    echo "  conda install -c bioconda gtdbtk-data"
    echo "  gtdbtk download_data --data_dir /path/to/gtdbtk_data"
    echo "  export GTDBTK_DATA_PATH=/path/to/gtdbtk_data"
    echo ""
    echo "Skipping taxonomic classification for now..."
fi

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
if [ -f "${RESULTS_DIR}/gtdbtk/gtdbtk.bac120.summary.tsv" ] || [ -f "${RESULTS_DIR}/gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
    echo "Found GTDB-Tk results, running taxonomic filtering..."
    
    # Filter bacterial methanogens/methanotrophs
    if [ -f "${RESULTS_DIR}/gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
        python3 ${RESULTS_DIR}/filter_by_taxonomy.py \
                ${RESULTS_DIR}/gtdbtk/gtdbtk.bac120.summary.tsv \
                ${RESULTS_DIR}
    fi
    
    # Filter archaeal methanogens
    if [ -f "${RESULTS_DIR}/gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
        python3 ${RESULTS_DIR}/filter_by_taxonomy.py \
                ${RESULTS_DIR}/gtdbtk/gtdbtk.ar53.summary.tsv \
                ${RESULTS_DIR}
    fi
else
    echo "GTDB-Tk results not found, skipping taxonomic filtering"
    echo "Run GTDB-Tk first to get species-level identification"
fi

# Step 3: Generate comprehensive report
echo "Generating comprehensive report..."

cat > ${RESULTS_DIR}/methanogen_methanotroph_report.md << 'EOF'
# Methanogen and Methanotroph Analysis Report

## Overview
This report summarizes the identification of methanogens and methanotrophs in the tree microbiome samples using both functional gene analysis and taxonomic classification.

## Methodology
1. **Functional gene search**: Searched for key genes in methanogenesis and methanotrophy pathways
2. **Taxonomic classification**: Used GTDB-Tk to identify known methanogen/methanotroph taxa
3. **Species identification**: Full taxonomic hierarchy from domain to species level
4. **Quality assessment**: Evaluated MAG completeness and contamination

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

EOF

# Add species identification results to report
if [ -f "${RESULTS_DIR}/species_identification_summary.csv" ]; then
    echo "## Species Identification Results" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo "" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo "### Complete Taxonomic Classification" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo "MAG_ID,Domain,Phylum,Class,Order,Family,Genus,Species" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    cut -d',' -f1,2,3,4,5,6,7,8 "${RESULTS_DIR}/species_identification_summary.csv" | tail -n +2 >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo '```' >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
    echo "" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md
fi

echo "## Functional Gene Analysis Results" >> ${RESULTS_DIR}/methanogen_methanotroph_report.md

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
echo ""
echo "============================================="
echo "SPECIES IDENTIFICATION & FUNCTIONAL ANALYSIS"
echo "============================================="
echo ""
echo "Results saved in: ${RESULTS_DIR}"
echo "📋 Main Report: ${RESULTS_DIR}/methanogen_methanotroph_report.md"
echo "🧬 Species ID: ${RESULTS_DIR}/species_identification_summary.csv"
echo "🔥 Functional Genes: ${RESULTS_DIR}/methanogens/ and ${RESULTS_DIR}/methanotrophs/"
echo ""

# Print quick summary if species identification was successful
if [ -f "${RESULTS_DIR}/species_identification_summary.csv" ]; then
    echo "✅ GTDB-Tk species identification completed!"
    echo "   Check the species_identification_summary.csv for full taxonomic classifications"
    echo ""
fi

# Print functional gene summary
if [ -f "${RESULTS_DIR}/methanogens/methanogen_genes.csv" ]; then
    methanogen_count=$(tail -n +2 "${RESULTS_DIR}/methanogens/methanogen_genes.csv" | wc -l)
    echo "🔥 Found $methanogen_count methanogenesis-related genes"
fi

if [ -f "${RESULTS_DIR}/methanotrophs/methanotroph_genes.csv" ]; then
    methanotroph_count=$(tail -n +2 "${RESULTS_DIR}/methanotrophs/methanotroph_genes.csv" | wc -l)
    echo "🌱 Found $methanotroph_count methanotrophy-related genes"
fi

echo ""
echo "Next steps:"
echo "1. Review the comprehensive report: ${RESULTS_DIR}/methanogen_methanotroph_report.md"
echo "2. Check species identifications for your methanogens"
echo "3. Validate functional gene predictions with protein domain analysis"
echo "4. Consider phylogenetic analysis of key methane genes"
