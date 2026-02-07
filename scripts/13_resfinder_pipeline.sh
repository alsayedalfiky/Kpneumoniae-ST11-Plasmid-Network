#!/bin/bash
# 13_resfinder_pipeline.sh
# Automated AMR profiling and strict filtering for high-confidence gene calls

# 1. Run ResFinder 4.7.0
# Using -s Klebsiella and --acquired flags as per methods
mkdir -p results/resfinder
for f in genomes/*.fna; do
    base=$(basename "$f" .fna)
    python -m resfinder -ifa "$f" -o results/resfinder/${base} -s klebsiella --acquired
done

# 2. Consolidate Results
# Prepending Sample ID to individual tabular outputs
echo -e "Sample_ID\tGene\tIdentity\tAlignment_Length\tCoverage\tPosition\tContig\tNotes" > results/resfinder_all.tsv
for d in results/resfinder/*; do
    base=$(basename "$d")
    if [ -f "$d/ResFinder_results_tab.txt" ]; then
        awk -v id="$base" 'NR>1 {print id "\t" $0}' "$d/ResFinder_results_tab.txt" >> results/resfinder_all.tsv
    fi
done

# 3. Apply Stringent Filtering (Identity >= 95%, Coverage >= 90%)
awk -F'\t' 'NR==1 || ($3+0 >= 95 && $5+0 >= 90)' results/resfinder_all.tsv > results/resfinder_filtered.tsv

echo "Filtering complete. $(wc -l < results/resfinder_filtered.tsv) determinants retained."
