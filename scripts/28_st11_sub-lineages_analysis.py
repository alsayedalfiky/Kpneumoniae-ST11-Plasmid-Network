#!/usr/bin/env python3
"""
ST11 Sub-lineage Stratification & Plasmid Association
Investigates the association of the IncFII(pHN7A8) plasmid with 
specific capsule (K-locus) and LPS (O-type) variants within ST11.

Specifically compares:
- KL64 : O2 (Emerging high-risk sub-lineage)
- KL47 : O13 (Established sub-lineage)

Method: Fisher's Exact Test for enrichment comparison.
"""

import pandas as pd
from scipy.stats import fisher_exact

# ---------------------------------------------------------
# 1. Configuration & File Paths
# ---------------------------------------------------------
MASTER_FILE = "master_results_with_metadata_focused.csv"
PLASMID_FILE = "results/plasmid_with_metadata_fixed.tsv"

# Column name mapping
ID_COL = 'strain'
ST_COL = 'Kleborate__mlst__ST'
K_COL = 'kaptive__K_locus'
O_COL = 'kaptive__O_type'
P_ACC_COL = 'Accession'
P_REP_COL = 'Replicon'

def main():
    # ---------------------------------------------------------
    # 2. Data Loading
    # ---------------------------------------------------------
    master = pd.read_csv(MASTER_FILE)
    plasmid = pd.read_table(PLASMID_FILE, sep='\t')

    # ---------------------------------------------------------
    # 3. Filter for Exact ST11 Lineage
    # ---------------------------------------------------------
    # Clean whitespace and isolate ST11 genomes
    master[ST_COL] = master[ST_COL].astype(str).str.strip()
    st11 = master[master[ST_COL].isin(["11", "ST11"])].copy()

    print(f"Total ST11 isolates identified: {len(st11)}")

    # ---------------------------------------------------------
    # 4. Stratify Sub-lineages (K/O Types)
    # ---------------------------------------------------------
    # Identify the KL47-O13 cluster
    kl47_o13 = st11[
        st11[K_COL].str.contains("KL47", case=False, na=False) & 
        st11[O_COL].str.contains("O13", case=False, na=False)
    ]

    # Identify the KL64-O2a cluster (Capsule Switch lineage)
    kl64_o2a = st11[
        st11[K_COL].str.contains("KL64", case=False, na=False) & 
        st11[O_COL].str.contains("O2", na=False)
    ]

    # ---------------------------------------------------------
    # 5. Identify Plasmid Carriers
    # ---------------------------------------------------------
    # Define the target replicon of interest
    target_rep = "IncFII(pHN7A8)"
    carrier_accessions = set(
        plasmid[plasmid[P_REP_COL].str.contains(target_rep, na=False, regex=False)][P_ACC_COL].unique()
    )

    # ---------------------------------------------------------
    # 6. Statistical Calculation Logic
    # ---------------------------------------------------------
    def calculate_stats(df, carriers):
        total = len(df)
        if total == 0: 
            return 0, 0, 0.0
        pos = df[ID_COL].isin(carriers).sum()
        neg = total - pos
        prevalence = (pos / total) * 100
        return pos, neg, prevalence

    pos64, neg64, prev64 = calculate_stats(kl64_o2a, carrier_accessions)
    pos47, neg47, prev47 = calculate_stats(kl47_o13, carrier_accessions)

    # ---------------------------------------------------------
    # 7. Output Results and Significance
    # ---------------------------------------------------------
    print(f"\n" + "="*45)
    print(f" ST11 SUB-LINEAGE ANALYSIS: {target_rep}")
    print("="*45)
    print(f"Group KL64:O2v | n = {len(kl64_o2a):<4} | Carriers: {pos64} ({prev64:.1f}%)")
    print(f"Group KL47:O13 | n = {len(kl47_o13):<4} | Carriers: {pos47} ({prev47:.1f}%)")

    if len(kl64_o2a) > 0 and len(kl47_o13) > 0:
        # Contingency table for Fisher's Exact Test
        table = [[pos64, neg64], [pos47, neg47]]
        res = fisher_exact(table)
        
        print("-" * 45)
        print(f"Fisher's Exact p-value: {res.pvalue:.4e}")
        print(f"Odds Ratio: {res.statistic:.2f}")
        print("-" * 45)

if __name__ == "__main__":
    main()
