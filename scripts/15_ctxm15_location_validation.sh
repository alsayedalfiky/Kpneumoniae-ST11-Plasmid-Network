#!/bin/bash
# Step 15: Validation of blaCTX-M-15 Genomic Location
# This script identifies chromosomal suspects and plasmid controls, 
# then uses MOB-suite to verify their genomic context.

# --- 1. SUSPECT IDENTIFICATION ---
# Identify chromosomal suspects (94 cases) and a random subset of plasmid controls (30 cases)
# Assumes resfinder_filtered.tsv is in results/
mkdir -p metadata

# Chromosomal Suspects
awk -F'\t' '$2 == "blaCTX-M-15" && $7 ~ /[Cc]hromosome/ {print $1, $7}' results/resfinder_filtered.tsv > metadata/ctxm15_chromosomal_suspects.txt

# Plasmid Controls (Randomly sampled for validation benchmarking)
awk -F'\t' '$2 == "blaCTX-M-15" && $7 ~ /[Pp]lasmid/ {print $1, $7}' results/resfinder_filtered.tsv | shuf -n 30 > metadata/ctxm15_plasmid_controls.txt

# --- 2. CONFIGURATION ---
FNA_DIR="data/genomes"
MASTER_OUT="results/ctxm15_validation_$(date +%F)"
mkdir -p "$MASTER_OUT"

# --- 3. EXECUTION LOOP ---
for GROUP in chromosomal_suspects plasmid_controls; do
    INPUT_FILE="metadata/ctxm15_${GROUP}.txt"
    CURRENT_OUT_DIR="${MASTER_OUT}/${GROUP}"
    mkdir -p "$CURRENT_OUT_DIR"

    if [ ! -f "$INPUT_FILE" ]; then
        echo "Warning: $INPUT_FILE not found, skipping..."
        continue
    fi

    echo ">>> Processing Group: $GROUP"

    while read -r line; do
        # Parse Accession ID and Contig ID
        full_id=$(echo "$line" | awk '{print $1}')
        acc_id=$(echo "$full_id" | grep -oP "^GC._[0-9]+\.[0-9]+")
        contig_id=$(echo "$line" | awk '{print $2}')
        
        # Locate FASTA file
        fna_file="${FNA_DIR}/${acc_id}.fna"
        if [ ! -f "$fna_file" ]; then
            fna_file=$(find "$FNA_DIR" -name "${acc_id}*" -type f | head -n 1)
        fi

        if [ -z "$fna_file" ] || [ ! -f "$fna_file" ]; then
            echo "   [!] MISSING FNA: $acc_id"
            continue
        fi

        echo "   Processing: $acc_id | Contig: $contig_id"
        
        # Extract target contig
        echo "$contig_id" > "${CURRENT_OUT_DIR}/${acc_id}_id.tmp"
        seqtk subseq "$fna_file" "${CURRENT_OUT_DIR}/${acc_id}_id.tmp" > "${CURRENT_OUT_DIR}/${acc_id}_contig.fasta"
        
        # Run MOB-typer for classification
        if [ -s "${CURRENT_OUT_DIR}/${acc_id}_contig.fasta" ]; then
            mob_typer --infile "${CURRENT_OUT_DIR}/${acc_id}_contig.fasta" --out_file "${CURRENT_OUT_DIR}/report_${acc_id}.txt"
        else
            echo "   [!] EXTRACTION FAILED for $acc_id"
        fi
        
        # Clean up temp files
        rm "${CURRENT_OUT_DIR}/${acc_id}_id.tmp" "${CURRENT_OUT_DIR}/${acc_id}_contig.fasta"
        
    done < "$INPUT_FILE"
    
    echo ">>> Completed $GROUP. Reports saved in $CURRENT_OUT_DIR"
    echo "--------------------------------------------------------"
done

# --- 4. FINAL VERIFICATION SUMMARY ---
echo "Final Report Counts:"
find "$MASTER_OUT" -name "report_*.txt" | cut -d/ -f3 | uniq -c
