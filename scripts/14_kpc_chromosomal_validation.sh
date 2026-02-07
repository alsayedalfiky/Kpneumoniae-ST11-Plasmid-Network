#!/bin/bash
# Step 14: Validation of blaKPC-2 chromosomal integration suspected cases
# This script extracts specific contigs and uses MOB-typer to verify genomic location.

# Configuration
SUSPECTS="metadata/kpc_chromosome_suspects.txt"
FNA_DIR="data/genomes"
OUT_DIR="results/validation_kpc_chromosomal"

mkdir -p "$OUT_DIR"

echo "------------------------------------------------"
echo "Starting final verification of KPC chromosomal suspects..."
echo "------------------------------------------------"

while read -r line; do
    # 1. Parse IDs (Supports 'GCF_000...' or 'GCA_000...' formats)
    full_id=$(echo "$line" | awk '{print $1}')
    acc_id=$(echo "$full_id" | grep -oP "^GC._[0-9]+\.[0-9]+")
    contig_id=$(echo "$line" | awk '{print $2}')
    
    # 2. Locate Genome File
    fna_file="${FNA_DIR}/${acc_id}.fna"
    
    # Fallback search if exact match isn't found
    if [ ! -f "$fna_file" ]; then
        fna_file=$(find "$FNA_DIR" -name "${acc_id}*" -type f | head -n 1)
    fi

    if [ -z "$fna_file" ] || [ ! -f "$fna_file" ]; then
        echo "WARNING: File not found for $acc_id"
        continue
    fi

    echo "Processing: $acc_id | Contig: $contig_id"
    
    # 3. Extract target contig using seqtk
    # Uses a unique temp name to ensure thread safety
    echo "$contig_id" > "${OUT_DIR}/${acc_id}_id.tmp"
    seqtk subseq "$fna_file" "${OUT_DIR}/${acc_id}_id.tmp" > "${OUT_DIR}/${acc_id}_contig.fasta"
    
    # 4. MOB-typer Analysis
    # Predicts if the extracted contig is chromosomal or plasmid-derived
    if [ -s "${OUT_DIR}/${acc_id}_contig.fasta" ]; then
        mob_typer --infile "${OUT_DIR}/${acc_id}_contig.fasta" --out_file "${OUT_DIR}/report_${acc_id}.txt"
        echo "SUCCESS: Validation report generated for ${acc_id}."
    fi
    
    # Internal Cleanup
    rm "${OUT_DIR}/${acc_id}_id.tmp" "${OUT_DIR}/${acc_id}_contig.fasta"
    
done < "$SUSPECTS"

echo "------------------------------------------------"
echo "Validation Complete. Results located in: $OUT_DIR"
