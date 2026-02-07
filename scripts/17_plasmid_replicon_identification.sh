#!/bin/bash
# Step 17: Plasmid Replicon Identification
# This script uses the KMA algorithm to query genome assemblies against 
# the PlasmidFinder database (Enterobacteriaceae scheme).

# 1. Setup Environment
# Ensure output directory exists
OUT_DIR="results/plasmidfinder"
mkdir -p "$OUT_DIR"

# Define database path (Update this path to your local database location)
DB_PATH="databases/plasmidfinder/enterobacteriaceae"

# 2. Iterative Genomic Mapping
# Processing all 1,010 genomes found in the data directory
echo "Starting PlasmidFinder analysis..."

find data/genomes/ -name "*.fna" | while read -r fna; do
    acc=$(basename "$fna" .fna)
    echo "[Processing] Genome: $acc"
    
    # Execute KMA mapping
    # -asm: Assembly mode for input contigs
    # -1t1: One-to-one mapping to prevent redundant hits
    # -cge: Output format compatible with CGE tools
    # -t 12: Utilize 12 CPU threads
    kma -i "$fna" \
        -t_db "$DB_PATH" \
        -o "${OUT_DIR}/${acc}_plasmids" \
        -asm -1t1 -cge -t 12
done

echo "Analysis complete. Results stored in $OUT_DIR."
