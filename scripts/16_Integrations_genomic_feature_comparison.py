"""
Step 16: Statistical Comparison of Genomic Features
Analyzes GC content and contig size differences between chromosomal and 
plasmid-borne blaCTX-M-15 determinants using Welch's t-test.
"""

import pandas as pd
from scipy import stats
import os

# 1. Load data
# Assumes the summary file generated from Step 15 validation
input_file = 'results/ctxm15_chromosome_vs_plasmid.csv'
output_file = 'results/genomic_features_stats.txt'

if not os.path.exists(input_file):
    print(f"Error: {input_file} not found.")
else:
    df = pd.read_csv(input_file)

    # 2. Extract and clean columns 
    # (dropna handles the uneven sample sizes between groups)
    chrom_gc = df['Chromosome_GC_%'].dropna()
    plasmid_gc = df['Plasmid_GC_%'].dropna()

    chrom_size = df['Chromosome_Size_(Mb)'].dropna()
    plasmid_size = df['Plasmid_Size_(Mb)'].dropna()

    # 3. Perform Welch's T-test (equal_var=False)
    # This accounts for unequal variances and sample sizes
    t_stat_gc, p_val_gc = stats.ttest_ind(chrom_gc, plasmid_gc, equal_var=False)
    t_stat_size, p_val_size = stats.ttest_ind(chrom_size, plasmid_size, equal_var=False)

    # 4. Significance Helper
    def get_stars(p):
        if p < 0.001: return "***"
        elif p < 0.01: return "**"
        elif p < 0.05: return "*"
        else: return "ns"

    # 5. Formulate Output
    results = [
        "--- GC Content Analysis ---",
        f"Chromosome Mean: {chrom_gc.mean():.2f}%",
        f"Plasmid Mean:    {plasmid_gc.mean():.2f}%",
        f"P-value:         {p_val_gc:.2e} ({get_stars(p_val_gc)})",
        "\n--- Contig Size Analysis ---",
        f"Chromosome Mean: {chrom_size.mean():.2f} Mb",
        f"Plasmid Mean:    {plasmid_size.mean():.2f} Mb",
        f"P-value:         {p_val_size:.2e} ({get_stars(p_val_size)})"
    ]

    # Print to console and save to file
    with open(output_file, 'w') as f:
        for line in results:
            print(line)
            f.write(line + "\n")

    print(f"\n[Success] Statistical summary saved to {output_file}")
