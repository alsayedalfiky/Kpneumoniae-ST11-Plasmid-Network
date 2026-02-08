"""
Step 21: ST-wise AMR and Plasmid Burden Analysis
Aggregates isolate-level data to calculate population-level statistics.
Provides the denominators used for prevalence matrix generation.
"""

import pandas as pd

# 1. Configuration
REFINED_METADATA = "results/plasmid_with_metadata_fixed.tsv"
FILTERED_RESFINDER = "results/resfinder_filtered.tsv"
OUTPUT_FILE = "results/ST_AMR_and_Plasmid_Burden_2.tsv"

# 2. Load Data
# Using the fixed metadata (resolved STs) is crucial for accurate aggregation
meta = pd.read_csv(REFINED_METADATA, sep="\t")
resfinder = pd.read_csv(FILTERED_RESFINDER, sep="\t", header=None)
resfinder.columns = ["Genome", "Gene", "Ident", "Aln", "Cov", "PosRef", "Contig", "PosCont", "Pheno", "Acc"]

# 3. AMR Burden per ST
merged_amr = pd.merge(resfinder[["Genome", "Gene"]], meta[["Genome", "ST"]], on="Genome")
iso_amr = merged_amr.groupby(["Genome", "ST"])["Gene"].nunique().reset_index(name="Count")

st_amr = iso_amr.groupby("ST").agg(
    Isolate_Count_AMR=("Genome", "count"),
    Avg_AMR_Genes=("Count", "mean"),
    Median_AMR_Genes=("Count", "median"),
    Max_AMR_Genes=("Count", "max")
)

# 4. Plasmid Burden per ST
iso_plasmid = meta.groupby(["Genome", "ST"])["Replicon"].nunique().reset_index(name="Count")

st_plasmids = iso_plasmid.groupby("ST").agg(
    Isolate_Count_Plasmid=("Genome", "count"),
    Avg_Plasmids=("Count", "mean"),
    Max_Plasmids=("Count", "max"),
    Total_Plasmids=("Count", "sum")
)

# 5. Final Merge and Export
combined = pd.merge(st_amr, st_plasmids, on="ST").round(2)
combined.to_csv(OUTPUT_FILE, sep="\t")

print(f"Burden analysis complete. Results exported to {OUTPUT_FILE}")
