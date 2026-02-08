#!/bin/bash
# Step 18: Consolidate Plasmid Results
# Aggregates 1,010 individual KMA results into a unified summary table.

# Configuration
OUTPUT="../results/plasmid_summary.tsv"

# Initialize header
echo -e "Genome\tReplicon\tIdentity\tCoverage" > "$OUTPUT"

# Merge loop
for res in ../results/plasmid_results_all/*.res; do
    genome=$(basename "$res" _plasmids.res)
    # Extract columns: Template (1), Identity (7), Coverage (8)
    awk -v id="$genome" 'NR>1 {print id "\t" $1 "\t" $7 "\t" $8}' "$res" >> "$OUTPUT"
done

echo "Successfully merged $(ls ../results/plasmid_results_all/*.res | wc -l) files."
