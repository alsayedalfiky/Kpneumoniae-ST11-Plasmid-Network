"""
# Filter Plasmid-Borne AMR Hits
Identifies resistance genes located on contigs with 'plasmid' in their header.
This defines the subset of genomes driving plasmid-mediated resistance.
"""

import pandas as pd
import os

# Configuration
RESFINDER_RAW = "results/resfinder_filtered.tsv"
OUTPUT_LINKS = "results/physical_plasmid_links.tsv"

def main():
    if not os.path.exists(RESFINDER_RAW):
        print(f"Error: Input file {RESFINDER_RAW} not found.")
        return

    # 1. Load Raw ResFinder Data
    # Columns: 0=Genome, 1=Gene, 6=Contig Name
    df = pd.read_csv(RESFINDER_RAW, sep="\t", header=None)
    
    # 2. Filter for 'plasmid' in contig name (Column 6)
    # Case insensitive search to catch 'Plasmid', 'plasmid', 'PLASMID'
    plasmid_hits = df[df[6].str.contains('plasmid', case=False, na=False)]
    
    # 3. Export Filtering List
    # We save Genome and Gene to define our 'Plasmid-Positive' dataset
    plasmid_hits.to_csv(OUTPUT_LINKS, sep="\t", index=False, header=False)
    
    print(f"Filtering complete.")
    print(f"Found {len(plasmid_hits)} plasmid-borne AMR genes.")
    print(f"Results saved to {OUTPUT_LINKS}")

if __name__ == "__main__":
    main()
