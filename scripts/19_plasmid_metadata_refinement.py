"""
Step 19: Plasmid-Metadata Integration and ST Refinement
Consolidates replicon data with epidemiological metadata and resolves 
previously unassigned Sequence Types (STs) using Kleborate results.
"""

import pandas as pd
import os

# Configuration
PLASMID_INPUT = 'results/plasmid_with_metadata.tsv'
ST_REPORT = 'results/st_comparison_report.csv'
OUTPUT_FILE = 'results/plasmid_with_metadata_fixed.tsv'

def main():
    if not os.path.exists(PLASMID_INPUT) or not os.path.exists(ST_REPORT):
        print("Error: Required input files not found.")
        return

    # 1. Load Datasets
    # plasmid_df contains initial merge; st_report contains Kleborate 'truth'
    plasmid_df = pd.read_csv(PLASMID_INPUT, sep='\t')
    st_report = pd.read_csv(ST_REPORT)

    # 2. Extract standardized Accession IDs
    plasmid_df['clean_id'] = plasmid_df['Genome'].str.extract(r'(GC._\d+\.\d+)')

    # 3. Create mapping from Kleborate results
    # Priority is given to Kleborate for resolving 'unknown' or missing STs
    st_map = dict(zip(st_report['id'], st_report['st_kleborate']))

    def resolve_st(row):
        """Helper to fill missing STs while preserving existing valid data."""
        current_st = str(row['ST']).strip()
        kleborate_st = st_map.get(row['clean_id'])

        # Check for missing/invalid ST values
        if current_st in ['nan', '-', 'null', '0', 'Unknown', '']:
            return kleborate_st
        return current_st

    # 4. Apply refinement
    print(f"Refining ST assignments for {len(plasmid_df)} records...")
    plasmid_df['ST'] = plasmid_df.apply(resolve_st, axis=1)

    # 5. Finalize and Save
    # Drop helper column and save as a clean TSV
    plasmid_df.drop(columns=['clean_id'], inplace=True)
    plasmid_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    print(f"Success: Updated metadata saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
