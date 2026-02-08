#!/bin/bash
# Step 22: Plasmid Diversity and Richness Statistics
# Calculates:
# 1. Plasmid Pool Diversity: Count of unique replicon types associated with each ST.
# 2. Average Richness: Mean number of plasmids per isolate within an ST.

INPUT_FILE="results/plasmid_with_metadata_fixed.tsv"
OUTPUT_DIVERSITY="results/plasmid_diversity_per_st.tsv"
OUTPUT_RICHNESS="results/plasmid_richness_check.tsv"

# Ensure input exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE not found."
    exit 1
fi

echo "Generating Plasmid Statistics..."

# ---------------------------------------------------------
# Part 1: Plasmid Pool Diversity (Unique Replicons per ST)
# ---------------------------------------------------------
# Logic: Extract ST & Replicon -> Unique Pairs -> Count Replicons per ST
# Note: Adjust column numbers ($8=ST, $2=Replicon) if your header changes.
echo -e "Count\tST" > "$OUTPUT_DIVERSITY"

tail -n +2 "$INPUT_FILE" | \
awk -F'\t' '{print $8 "\t" $2}' | \
sort -u | \
cut -f1 | \
sort | \
uniq -c | \
sort -nr | \
head -n 20 >> "$OUTPUT_DIVERSITY"

echo "[1/2] Diversity stats saved to $OUTPUT_DIVERSITY"

# ---------------------------------------------------------
# Part 2: Average Plasmid Richness (Validation)
# ---------------------------------------------------------
# Logic: Sum total replicons / Unique Isolates per ST
echo -e "ST\tIsolates\tTotal_Replicons\tAvg_Richness" > "$OUTPUT_RICHNESS"

tail -n +2 "$INPUT_FILE" | \
awk -F'\t' '
{
  # Count unique isolates (Accessions) per ST
  # $8 = ST, $5 = Accession (Adjust based on your file structure)
  if (!seen_iso[$8,$5]++) {
    iso_count[$8]++
  }
  # Count total replicons per ST
  rep_count[$8]++
}
END {
  for (st in iso_count) {
    # Print: ST, Isolate Count, Total Replicons, Average
    printf "%s\t%d\t%d\t%.2f\n", st, iso_count[st], rep_count[st], rep_count[st]/iso_count[st]
  }
}' | sort -t$'\t' -k2,2nr | head -n 20 >> "$OUTPUT_RICHNESS"

echo "[2/2] Richness stats saved to $OUTPUT_RICHNESS"

# ---------------------------------------------------------
# Display Top Results
# ---------------------------------------------------------
echo -e "\n--- Top 5 STs by Plasmid Diversity (Unique Types) ---"
head -n 6 "$OUTPUT_DIVERSITY" | column -t

echo -e "\n--- Top 5 STs by Isolate Count (Richness Check) ---"
head -n 6 "$OUTPUT_RICHNESS" | column -t
