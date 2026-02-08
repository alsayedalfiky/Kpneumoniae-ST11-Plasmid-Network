#!/bin/bash
# Step 17a: Indexing the PlasmidFinder Database
# This command prepares the Enterobacteriaceae scheme for KMA mapping.

DB_DIR="/home/bioalfiky25/miniconda3/envs/amr_study/share/plasmidfinder-2.1.6/database"

cd "$DB_DIR"

# Indexing with k=16 as per standard PlasmidFinder protocol
kma_index -i enterobacteriaceae.fsa -o enterobacteriaceae -k 16

echo "Database indexing complete. Binary files generated for KMA mapping."
