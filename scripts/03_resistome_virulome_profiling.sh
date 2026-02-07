#!/bin/bash

# Define directories
INPUT_DIR="fna_clean"
OUTPUT_DIR="Kleborate_Outputs"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Run Kleborate profiling
# -a: assembly files
# -o: output file prefix
# -p: species-specific profiles (K. pneumoniae species complex)
# -t: threads
echo "Running Kleborate profiling for AMR and Virulence determinants..."

kleborate -a $INPUT_DIR/*.fna \
          -o $OUTPUT_DIR/kleborate_results \
          -p kpsc \
          -t 8

echo "Profiling complete. Results saved to $OUTPUT_DIR/kleborate_results.txt"
