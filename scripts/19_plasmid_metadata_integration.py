"""
Step 19: Plasmid-Metadata Integration and ST-wise Summarization
Merges Replicon profiling results with MLST and epidemiological metadata
to analyze the distribution of plasmids across different lineages.
"""

import pandas as pd
import os

# Configuration
PLASMID_FILE = "results/plasmid_summary_master.tsv"
METADATA_FILE = "results/mlst_with_metadata.tsv"
OUTPUT_MERGED = "results/plasmid_with_metadata.tsv"
OUTPUT_SUMMARY = "results/plasmid_per_ST.tsv"

# 1. Load Data
if not os.path.exists(PLASMID_FILE) or not os.path.exists(METADATA_FILE):
    print("Error: Input files missing. Ensure Step 18 is complete.")
    exit()

plasmids = pd.read_csv(PLASMID_FILE, sep="\t")
mlst_meta = pd.read_csv(METADATA_FILE, sep="\t")

# 2. Extract clean Accession ID for merging
# Standardizes 'GCF_002852675.2_ASM...' to 'GCF_002852675.2'
plasmids["Accession"] = plasmids["Genome"].str.extract(r"(GC._\d+\.\d+)")

# 3. Merge Datasets
# Left merge ensures we keep all plasmid hits and attach metadata where available
merged = plasmids.merge(mlst_meta, on="Accession", how="left")
merged.to_csv(OUTPUT_MERGED, sep="\t", index=False)
print(f"Success: Merged table saved with {len(merged)} records.")

# 4. Summarize Plasmid Distribution per Sequence Type (ST)
# Provides a frequency count of each replicon within each ST
st_summary = merged.groupby(["ST", "Replicon"]).size().reset_index(name="Count")
st_summary.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

print(f"Success: ST-wise distribution saved to {OUTPUT_SUMMARY}")
