#!/usr/bin/env python3
"""
Simple bin analysis script for metagenomic bins
Provides basic statistics while CheckM is unavailable
"""

import os
import sys
from pathlib import Path

def analyze_fasta(fasta_path):
    """Analyze basic statistics of a FASTA file"""
    sequences = []
    total_length = 0
    gc_content = []
    
    try:
        with open(fasta_path, 'r') as f:
            current_seq = ""
            for line in f:
                if line.startswith('>'):
                    if current_seq:
                        seq_len = len(current_seq)
                        sequences.append(seq_len)
                        total_length += seq_len
                        
                        # Calculate GC content
                        gc_count = current_seq.upper().count('G') + current_seq.upper().count('C')
                        if seq_len > 0:
                            gc_content.append((gc_count / seq_len) * 100)
                    current_seq = ""
                else:
                    current_seq += line.strip()
            
            # Don't forget the last sequence
            if current_seq:
                seq_len = len(current_seq)
                sequences.append(seq_len)
                total_length += seq_len
                gc_count = current_seq.upper().count('G') + current_seq.upper().count('C')
                if seq_len > 0:
                    gc_content.append((gc_count / seq_len) * 100)
                    
    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
        return None
    
    if not sequences:
        return None
    
    # Calculate N50
    def calculate_n50(lengths):
        lengths_sorted = sorted(lengths, reverse=True)
        total_length = sum(lengths_sorted)
        target_length = total_length * 0.5
        
        current_length = 0
        for length in lengths_sorted:
            current_length += length
            if current_length >= target_length:
                return length
        return 0
    
    return {
        'num_contigs': len(sequences),
        'total_length': total_length,
        'mean_length': sum(sequences) / len(sequences),
        'max_length': max(sequences),
        'min_length': min(sequences),
        'n50': calculate_n50(sequences),
        'mean_gc': sum(gc_content) / len(gc_content) if gc_content else 0,
        'longest_contigs': sorted(sequences, reverse=True)[:5]
    }

def find_bins(results_dir):
    """Find all bin files in the results directory"""
    bin_files = []
    results_path = Path(results_dir)
    
    # Search for different bin file patterns
    patterns = ["**/*.fa", "**/*.fasta", "**/*.fna"]
    
    for pattern in patterns:
        bin_files.extend(results_path.glob(pattern))
    
    # Filter out non-bin files
    filtered_bins = []
    skip_patterns = ['checkm', 'das_tool', 'all_metabat2_bins', 'all_maxbin2_bins', 'all_vamb_bins']
    
    for bin_file in bin_files:
        # Skip certain directories or files
        if any(skip in str(bin_file) for skip in skip_patterns):
            continue
        if 'contigs.fa' in bin_file.name or 'assembly' in bin_file.name:
            continue
        filtered_bins.append(bin_file)
    
    return filtered_bins

def main():
    results_dir = "/home/jiayi-chen/Documents/farting_trees/results/04_binning"
    
    print("=" * 80)
    print("METAGENOMIC BIN ANALYSIS REPORT")
    print("=" * 80)
    print(f"Analyzing bins in: {results_dir}")
    print()
    
    bin_files = find_bins(results_dir)
    
    if not bin_files:
        print("No bin files found!")
        print("Make sure binning has completed successfully.")
        return
    
    print(f"Found {len(bin_files)} bin files")
    print()
    
    # Analyze each bin
    all_results = []
    categories = {
        'MetaBAT2_standard': 0,
        'MetaBAT2_relaxed': 0,
        'MetaBAT2_inclusive': 0,
        'MaxBin2': 0,
        'VAMB': 0,
        'Methane_targeted': 0,
        'Other': 0
    }
    
    quality_bins = []  # Bins that might be high quality based on size
    
    for bin_file in sorted(bin_files):
        stats = analyze_fasta(bin_file)
        if stats is None:
            continue
        
        # Categorize bin
        bin_path_str = str(bin_file).lower()
        category = 'Other'
        
        if 'metabat2' in bin_path_str:
            if 'relaxed' in bin_path_str:
                category = 'MetaBAT2_relaxed'
            elif 'inclusive' in bin_path_str:
                category = 'MetaBAT2_inclusive'
            else:
                category = 'MetaBAT2_standard'
        elif 'maxbin2' in bin_path_str:
            category = 'MaxBin2'
        elif 'vamb' in bin_path_str:
            category = 'VAMB'
        elif 'methane' in bin_path_str:
            category = 'Methane_targeted'
        
        categories[category] += 1
        
        # Quality assessment based on basic metrics
        quality_score = "Unknown"
        if stats['total_length'] > 2000000 and stats['num_contigs'] < 100:
            quality_score = "Potentially High"
            quality_bins.append((bin_file.name, stats))
        elif stats['total_length'] > 1000000 and stats['num_contigs'] < 200:
            quality_score = "Potentially Medium"
        else:
            quality_score = "Low"
        
        result = {
            'bin_name': bin_file.name,
            'category': category,
            'quality_estimate': quality_score,
            **stats
        }
        all_results.append(result)
    
    # Summary by category
    print("BIN COUNT BY METHOD:")
    print("-" * 40)
    for category, count in categories.items():
        if count > 0:
            print(f"{category:20}: {count:4} bins")
    print()
    
    # Sort by total length
    all_results.sort(key=lambda x: x['total_length'], reverse=True)
    
    # Top 15 largest bins
    print("TOP 15 LARGEST BINS:")
    print("-" * 80)
    print(f"{'Rank':<4} {'Bin Name':<35} {'Size (Mb)':<10} {'Contigs':<8} {'N50 (kb)':<10} {'GC%':<6} {'Category'}")
    print("-" * 80)
    
    for i, result in enumerate(all_results[:15]):
        size_mb = result['total_length'] / 1e6
        n50_kb = result['n50'] / 1000
        print(f"{i+1:<4} {result['bin_name']:<35} {size_mb:<10.2f} {result['num_contigs']:<8} "
              f"{n50_kb:<10.1f} {result['mean_gc']:<6.1f} {result['category']}")
    
    print()
    
    # Potentially high-quality bins
    if quality_bins:
        print("POTENTIALLY HIGH-QUALITY BINS:")
        print("(>2MB size, <100 contigs - CheckM needed for completion/contamination)")
        print("-" * 80)
        for bin_name, stats in quality_bins:
            size_mb = stats['total_length'] / 1e6
            n50_kb = stats['n50'] / 1000
            print(f"{bin_name:<35} {size_mb:<8.2f} MB, {stats['num_contigs']:<4} contigs, "
                  f"N50: {n50_kb:<6.1f} kb, GC: {stats['mean_gc']:<5.1f}%")
        print()
    
    # Methane-targeted bins
    methane_results = [r for r in all_results if 'methane' in r['category'].lower()]
    if methane_results:
        print("METHANE-TARGETED BINS:")
        print("-" * 60)
        for result in methane_results:
            size_mb = result['total_length'] / 1e6
            n50_kb = result['n50'] / 1000
            print(f"{result['bin_name']:<30} {size_mb:<8.2f} MB, {result['num_contigs']:<4} contigs, "
                  f"GC: {result['mean_gc']:<5.1f}%")
        print()
    
    # Summary statistics
    total_bins = len(all_results)
    total_assembly = sum(r['total_length'] for r in all_results)
    large_bins = len([r for r in all_results if r['total_length'] > 1000000])
    
    print("SUMMARY STATISTICS:")
    print("-" * 40)
    print(f"Total bins analyzed: {total_bins}")
    print(f"Total binned assembly: {total_assembly/1e6:.1f} MB")
    print(f"Bins >1MB: {large_bins} ({large_bins/total_bins*100:.1f}%)")
    print(f"Average bin size: {total_assembly/total_bins/1e6:.2f} MB")
    print()
    
    print("NEXT STEPS:")
    print("-" * 40)
    print("1. Fix CheckM database issue for proper quality assessment")
    print("2. Focus analysis on largest bins (>1MB)")
    print("3. Run taxonomic classification on high-quality bins")
    print("4. Search for methane cycling genes in all bins")
    print("5. Use DAS Tool results for dereplicated high-quality set")
    
    # Save results as CSV manually
    output_file = Path(results_dir) / "bin_basic_analysis.csv"
    with open(output_file, 'w') as f:
        # Write header
        if all_results:
            f.write(','.join(all_results[0].keys()) + '\n')
            # Write data
            for result in all_results:
                f.write(','.join(str(v) for v in result.values()) + '\n')
    print(f"\nDetailed results saved to: {output_file}")

if __name__ == "__main__":
    main()
