"""
Step 20: Replicon Prevalence Matrix Generation
Calculates the prevalence of plasmid replicons within each Sequence Type (ST)
and generates a pivot matrix for downstream heatmap visualization.
"""

import pandas as pd
import os

# Configuration
REPLICON_COUNTS = "results/plasmid_per_ST.tsv"
ST_BURDEN_DATA = "results/ST_AMR_and_Plasmid_Burden_2.tsv"
OUTPUT_MATRIX = "results/replicon_prevalence_matrix.tsv"

def main():
    if not os.path.exists(REPLICON_COUNTS) or not os.path.exists(ST_BURDEN_DATA):
        print("Error: Required input files for Step 20 not found.")
        return

    # 1. Load Data
    # df_replicons: ST, Replicon, Count (from Step 19)
    # df_st_info: Metadata containing total plasmid-bearing isolates per ST
    df_replicons = pd.read_csv(REPLICON_COUNTS, sep="\t")
    df_st_info = pd.read_csv(ST_BURDEN_DATA, sep="\t")

    # 2. Integrate ST-specific isolate counts
    # Rename for clarity: ST_Count is the denominator (total plasmid-bearing isolates in that ST)
    df_st_info = df_st_info[['ST', 'Isolate_Count_Plasmid']].rename(
        columns={'Isolate_Count_Plasmid': 'ST_Count'}
    )
    
    # Merge datasets on Sequence Type
    df_combined = df_replicons.merge(df_st_info, on='ST', how='left')

    # 3. Calculate Prevalence
    # Formula: (Number of isolates in ST with Replicon / Total plasmid-bearing isolates in ST) * 100
    df_combined['Prevalence_%'] = (df_combined['Count'] / df_combined['ST_Count']) * 100

    # 4. Pivot to Matrix
    # Rows = Sequence Types | Columns = Plasmid Replicons | Values = Prevalence %
    # .fillna(0) ensures that missing combinations are treated as 0% prevalence
    matrix = df_combined.pivot(
        index='ST', 
        columns='Replicon', 
        values='Prevalence_%'
    ).fillna(0)

    # 5. Save Results
    matrix.to_csv(OUTPUT_MATRIX, sep="\t")
    
    print(f"Success: Matrix generated with {matrix.shape[0]} STs and {matrix.shape[1]} replicons.")
    print(f"Matrix saved to: {OUTPUT_MATRIX}")

if __name__ == "__main__":
    main()
