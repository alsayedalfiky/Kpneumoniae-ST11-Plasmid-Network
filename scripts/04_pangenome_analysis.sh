#!/bin/bash

# Define directories
INPUT_DIR="gff_for_panaroo"
OUTPUT_DIR="results/pangenome"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

echo "Starting Panaroo pangenome analysis on 1,010 isolates..."
echo "Mode: Strict (Clean and Refind)"

# Run Panaroo 
# --clean-mode strict: Removes contaminated/fragmented genes
# --remove-invalid-genes: Filters out non-standard GFF entries
# --refind-mode strict: Re-scans for genes missed by initial annotation
# -a core: Generates a core genome alignment
panaroo -i $INPUT_DIR/*.gff \
        -o $OUTPUT_DIR \
        --clean-mode strict \
        --remove-invalid-genes \
        --refind-mode strict \
        -a core \
        -t 4

echo "Pangenome analysis complete. Results saved in $OUTPUT_DIR"
