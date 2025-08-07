#!/usr/bin/env python3
"""
Advanced analysis and visualization for metagenomic data
Farting Trees Project
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import glob
from pathlib import Path

class MetagenomeAnalyzer:
    """
    Class for advanced metagenomic data analysis and visualization
    """
    
    def __init__(self, results_dir="../results"):
        self.results_dir = Path(results_dir)
        self.samples = ["53394_A15", "53395_A16", "53396_B6", "53397_B10", "53398_A16S", "53399_B10S"]
        
    def load_checkm_results(self):
        """Load CheckM quality results for all samples"""
        quality_data = []
        
        for sample in self.samples:
            # Try DAS Tool results first
            dastool_file = self.results_dir / f"03_binning/checkm/{sample}_dastool_quality.txt"
            metabat_file = self.results_dir / f"03_binning/checkm/{sample}_metabat2_quality.txt"
            
            if dastool_file.exists():
                df = pd.read_csv(dastool_file, sep='\t')
                df['sample'] = sample
                df['binner'] = 'DAS_Tool'
                quality_data.append(df)
            elif metabat_file.exists():
                df = pd.read_csv(metabat_file, sep='\t')
                df['sample'] = sample
                df['binner'] = 'MetaBAT2'
                quality_data.append(df)
        
        if quality_data:
            return pd.concat(quality_data, ignore_index=True)
        else:
            return pd.DataFrame()
    
    def load_gtdbtk_results(self):
        """Load GTDB-Tk taxonomic classification results"""
        gtdbtk_file = self.results_dir / "04_annotation/gtdbtk/output/gtdbtk.bac120.summary.tsv"
        
        if gtdbtk_file.exists():
            return pd.read_csv(gtdbtk_file, sep='\t')
        else:
            return pd.DataFrame()
    
    def load_functional_results(self):
        """Load functional gene analysis results"""
        methanogen_file = self.results_dir / "05_taxonomy/methanogens/methanogen_genes.csv"
        methanotroph_file = self.results_dir / "05_taxonomy/methanotrophs/methanotroph_genes.csv"
        
        functional_data = []
        
        if methanogen_file.exists():
            df = pd.read_csv(methanogen_file)
            functional_data.append(df)
        
        if methanotroph_file.exists():
            df = pd.read_csv(methanotroph_file)
            functional_data.append(df)
        
        if functional_data:
            return pd.concat(functional_data, ignore_index=True)
        else:
            return pd.DataFrame()
    
    def plot_mag_quality(self, save_path=None):
        """Create MAG quality plots"""
        quality_df = self.load_checkm_results()
        
        if quality_df.empty:
            print("No CheckM results found")
            return
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('MAG Quality Overview', 'Completeness vs Contamination', 
                          'Quality by Sample', 'Quality Distribution'),
            specs=[[{"secondary_y": True}, {}],
                   [{}, {}]]
        )
        
        # Plot 1: Quality overview
        high_quality = quality_df[(quality_df['Completeness'] >= 90) & (quality_df['Contamination'] <= 5)]
        medium_quality = quality_df[(quality_df['Completeness'] >= 50) & (quality_df['Contamination'] <= 10)]
        
        quality_counts = pd.DataFrame({
            'Sample': quality_df['sample'].unique(),
            'High Quality': [len(high_quality[high_quality['sample'] == s]) for s in quality_df['sample'].unique()],
            'Medium Quality': [len(medium_quality[medium_quality['sample'] == s]) for s in quality_df['sample'].unique()],
            'Total MAGs': [len(quality_df[quality_df['sample'] == s]) for s in quality_df['sample'].unique()]
        })
        
        fig.add_trace(
            go.Bar(x=quality_counts['Sample'], y=quality_counts['High Quality'], 
                   name='High Quality', marker_color='darkgreen'),
            row=1, col=1
        )
        fig.add_trace(
            go.Bar(x=quality_counts['Sample'], y=quality_counts['Medium Quality'], 
                   name='Medium Quality', marker_color='orange'),
            row=1, col=1
        )
        
        # Plot 2: Completeness vs Contamination scatter
        fig.add_trace(
            go.Scatter(x=quality_df['Completeness'], y=quality_df['Contamination'],
                      mode='markers', text=quality_df['sample'],
                      marker=dict(size=8, color=quality_df['sample'].astype('category').cat.codes,
                                colorscale='viridis'),
                      name='MAGs'),
            row=1, col=2
        )
        
        # Add quality thresholds
        fig.add_hline(y=5, line_dash="dash", line_color="red", row=1, col=2)
        fig.add_vline(x=90, line_dash="dash", line_color="red", row=1, col=2)
        
        # Plot 3: Quality by sample boxplot
        for sample in quality_df['sample'].unique():
            sample_data = quality_df[quality_df['sample'] == sample]
            fig.add_trace(
                go.Box(y=sample_data['Completeness'], name=sample, boxpoints='all'),
                row=2, col=1
            )
        
        # Plot 4: Quality distribution histogram
        fig.add_trace(
            go.Histogram(x=quality_df['Completeness'], name='Completeness', 
                        opacity=0.7, nbinsx=20),
            row=2, col=2
        )
        
        fig.update_layout(height=800, title_text="MAG Quality Analysis")
        fig.update_xaxes(title_text="Sample", row=1, col=1)
        fig.update_yaxes(title_text="Number of MAGs", row=1, col=1)
        fig.update_xaxes(title_text="Completeness (%)", row=1, col=2)
        fig.update_yaxes(title_text="Contamination (%)", row=1, col=2)
        
        if save_path:
            fig.write_html(save_path)
        else:
            fig.show()
    
    def plot_taxonomy_overview(self, save_path=None):
        """Create taxonomic overview plots"""
        gtdbtk_df = self.load_gtdbtk_results()
        
        if gtdbtk_df.empty:
            print("No GTDB-Tk results found")
            return
        
        # Extract taxonomic levels
        gtdbtk_df['phylum'] = gtdbtk_df['classification'].str.split(';').str[1].str.replace('p__', '')
        gtdbtk_df['class'] = gtdbtk_df['classification'].str.split(';').str[2].str.replace('c__', '')
        gtdbtk_df['order'] = gtdbtk_df['classification'].str.split(';').str[3].str.replace('o__', '')
        
        # Create treemap for phylum distribution
        phylum_counts = gtdbtk_df['phylum'].value_counts()
        
        fig = px.treemap(
            values=phylum_counts.values,
            names=phylum_counts.index,
            title="Phylum Distribution of MAGs"
        )
        
        if save_path:
            fig.write_html(save_path)
        else:
            fig.show()
    
    def plot_functional_analysis(self, save_path=None):
        """Create functional gene analysis plots"""
        functional_df = self.load_functional_results()
        
        if functional_df.empty:
            print("No functional gene results found")
            return
        
        # Count genes by type and sample
        gene_counts = functional_df.groupby(['sample_bin', 'gene_type']).size().reset_index(name='count')
        
        fig = px.bar(
            gene_counts, 
            x='sample_bin', 
            y='count', 
            color='gene_type',
            title="Functional Genes Found in MAGs",
            labels={'count': 'Number of Genes', 'sample_bin': 'MAG'}
        )
        
        fig.update_xaxes(tickangle=45)
        
        if save_path:
            fig.write_html(save_path)
        else:
            fig.show()
    
    def generate_summary_report(self):
        """Generate a comprehensive summary report"""
        print("Generating summary report...")
        
        quality_df = self.load_checkm_results()
        gtdbtk_df = self.load_gtdbtk_results()
        functional_df = self.load_functional_results()
        
        report = []
        report.append("# Farting Trees Metagenomic Analysis Report")
        report.append(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append("")
        
        # MAG quality summary
        if not quality_df.empty:
            report.append("## MAG Quality Summary")
            high_qual = quality_df[(quality_df['Completeness'] >= 90) & (quality_df['Contamination'] <= 5)]
            medium_qual = quality_df[(quality_df['Completeness'] >= 50) & (quality_df['Contamination'] <= 10)]
            
            report.append(f"- Total MAGs: {len(quality_df)}")
            report.append(f"- High-quality MAGs (>90% complete, <5% contamination): {len(high_qual)}")
            report.append(f"- Medium-quality MAGs (>50% complete, <10% contamination): {len(medium_qual)}")
            report.append("")
        
        # Taxonomic summary
        if not gtdbtk_df.empty:
            report.append("## Taxonomic Summary")
            phyla = gtdbtk_df['classification'].str.split(';').str[1].str.replace('p__', '').value_counts()
            report.append("### Top Phyla:")
            for phylum, count in phyla.head(10).items():
                report.append(f"- {phylum}: {count} MAGs")
            report.append("")
        
        # Functional genes summary
        if not functional_df.empty:
            report.append("## Functional Genes Summary")
            methanogen_genes = functional_df[functional_df['gene_type'] == 'methanogen']
            methanotroph_genes = functional_df[functional_df['gene_type'] == 'methanotroph']
            
            report.append(f"- Methanogen genes found: {len(methanogen_genes)}")
            report.append(f"- Methanotroph genes found: {len(methanotroph_genes)}")
            report.append("")
        
        # Save report
        report_path = self.results_dir / "analysis_summary_report.md"
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Summary report saved to: {report_path}")

def main():
    """Main function for running analysis"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Advanced metagenomic analysis')
    parser.add_argument('--results_dir', default='../results', help='Results directory')
    parser.add_argument('--output_dir', default='../results/plots', help='Output directory for plots')
    parser.add_argument('--analysis', choices=['quality', 'taxonomy', 'functional', 'all'], 
                       default='all', help='Type of analysis to run')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize analyzer
    analyzer = MetagenomeAnalyzer(args.results_dir)
    
    # Run analyses
    if args.analysis in ['quality', 'all']:
        analyzer.plot_mag_quality(f"{args.output_dir}/mag_quality.html")
    
    if args.analysis in ['taxonomy', 'all']:
        analyzer.plot_taxonomy_overview(f"{args.output_dir}/taxonomy_overview.html")
    
    if args.analysis in ['functional', 'all']:
        analyzer.plot_functional_analysis(f"{args.output_dir}/functional_analysis.html")
    
    if args.analysis == 'all':
        analyzer.generate_summary_report()

if __name__ == "__main__":
    main()
