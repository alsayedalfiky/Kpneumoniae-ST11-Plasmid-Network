#!/usr/bin/env python3
"""
Step 26: ST-Plasmid Co-adaptation and Enrichment Analysis
Performs statistical testing (Fisher's Exact Test) to identify significant 
associations between specific Sequence Types (STs) and Plasmid Replicons.

Input: results/plasmid_with_metadata_fixed.tsv
Method: Fisher's Exact Test with Benjamini-Hochberg (FDR) correction.
Output: results/significant_plasmid_st_associations.csv
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------
# 1. Load Data
# ---------------------------------------------------------
# Using the metadata file with resolved Kleborate STs
file_path = 'results/plasmid_with_metadata_fixed.tsv'
df = pd.read_table(file_path)

# ---------------------------------------------------------
# 2. Data Pre-processing
# ---------------------------------------------------------
# Extract simplified Replicon names (ignoring specific alleles/accessions)
# Logic: IncFII(pHN7A8)_1_... -> IncFII(pHN7A8)
df['Replicon_Clean'] = df['Replicon'].str.split('_').str[0]

# Create Presence/Absence Matrix per Isolate (Accession)
presence = df.pivot_table(index='Accession', columns='Replicon_Clean', aggfunc='size', fill_value=0)
presence = (presence > 0).astype(int)

# Map Isolates back to their Sequence Types
st_map = df[['Accession', 'ST']].drop_duplicates().set_index('Accession')
full_data = st_map.merge(presence, left_index=True, right_index=True)

# ---------------------------------------------------------
# 3. Statistical Enrichment Analysis
# ---------------------------------------------------------
sts = full_data['ST'].unique()
replicons = presence.columns
results = []

print(f"Testing associations for {len(sts)} STs and {len(replicons)} Replicons...")

for st in sts:
    # Threshold: Filter for major STs (n >= 10) to increase statistical power
    st_group = full_data[full_data['ST'] == st]
    if len(st_group) < 10: 
        continue
    
    st_mask = (full_data['ST'] == st)
    for rep in replicons:
        rep_mask = (full_data[rep] == 1)
        
        # Contingency Table for Fisher's Exact Test:
        #             |  In ST  | Not in ST
        # ------------|---------|-----------
        # Has Plasmid |    a    |     c
        # No Plasmid  |    b    |     d
        a = np.sum(st_mask & rep_mask)
        b = np.sum(st_mask & ~rep_mask)
        c = np.sum(~st_mask & rep_mask)
        d = np.sum(~st_mask & ~rep_mask)
        
        # Test only if the combination actually exists in the dataset
        if a > 0:
            odds, p = fisher_exact([[a, b], [c, d]], alternative='greater')
            results.append({
                'ST': st, 
                'Replicon': rep, 
                'ST_Total': a+b, 
                'Rep_Total': a+c, 
                'Observed': a, 
                'OddsRatio': odds, 
                'PValue': p
            })

# ---------------------------------------------------------
# 4. Multiple Testing Correction & Filtering
# ---------------------------------------------------------
res_df = pd.DataFrame(results)

# Apply False Discovery Rate (FDR) correction (Benjamini/Hochberg)
res_df['P_adj'] = multipletests(res_df['PValue'], method='fdr_bh')[1]

# Filter for statistically significant associations (P_adj < 0.05)
significant = res_df[res_df['P_adj'] < 0.05].sort_values('OddsRatio', ascending=False)

# ---------------------------------------------------------
# 5. Save Results
# ---------------------------------------------------------
significant.to_csv('results/significant_plasmid_st_associations.csv', index=False)

print(f"Analysis Complete. Found {len(significant)} significant ST-Replicon associations.")
print("\nTop 10 Highly Significant Associations:")
print(significant[['ST', 'Replicon', 'Observed', 'OddsRatio', 'P_adj']].head(10))
