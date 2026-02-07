#!/bin/bash

# Define paths
ALIGNMENT="results/pangenome/core_gene_alignment.aln"
OUTPUT_PREFIX="k_pneumoniae_tree"
THREADS=6

echo "Building core-genome phylogeny for 1,010 isolates..."
echo "Model: GTR+G (General Time Reversible with Gamma distribution)"

# Run IQ-TREE
# -s: Input alignment (from Panaroo)
# -m: Substitution model (GTR+G is standard for core-genome)
# -nt: Number of threads
# -fast: Fast search mode (recommended for large 1000+ genome datasets)
# -pre: Output file prefix
iqtree -s $ALIGNMENT \
       -m GTR+G \
       -nt $THREADS \
       -fast \
       -pre $OUTPUT_PREFIX

echo "Phylogeny construction complete. Tree file saved as ${OUTPUT_PREFIX}.treefile"
