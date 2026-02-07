#!/bin/bash

# Define directories
INPUT_DIR="genomes"
OUTPUT_DIR="results"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Run MLST on all RefSeq genomes (fna) using 12 threads
# Scheme: kpneumoniae
echo "Starting MLST analysis for Klebsiella pneumoniae isolates..."

mlst --scheme kpneumoniae \
     --threads 12 \
     $(find $INPUT_DIR -name "*.fna") > $OUTPUT_DIR/mlst_results.tsv

echo "MLST analysis complete. Results saved to $OUTPUT_DIR/mlst_results.tsv"
