#Use NCBI Datasets Command Line
datasets summary genome taxon "Klebsiella pneumoniae" \
  --assembly-source refseq \
  --assembly-level complete \
  --released-after 2023-01-01 \
  --released-before 2025-01-01 \
  --exclude-atypical \
  --as-json-lines | \
jq -r '.accession' > all_accessions.txt
# Count accessions
(base) bioalfiky25@Asus-ZenBook:~/AMR_2025$ echo "Total accessions: $(wc -l all_accessions.txt)"
Total accessions: 1010 all_accessions.txt
# 4. Download genomes (in batches of 50)
split -l 50 all_accessions.txt batch_
parallel -j 4 'datasets download genome accession --inputfile {} --filename {}.zip --include genome' ::: batch_*

