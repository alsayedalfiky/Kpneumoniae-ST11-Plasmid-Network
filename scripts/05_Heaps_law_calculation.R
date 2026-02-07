#!/usr/bin/env Rscript

# Load required library
if (!require("micropan")) install.packages("micropan")
library(micropan)

# 1. Read the Panaroo output
# Panaroo gene_presence_absence.csv uses specific headers
message("Loading Panaroo gene presence/absence matrix...")
pan_data <- read.csv("gene_presence_absence.csv", 
                     check.names = FALSE, 
                     stringsAsFactors = FALSE)

# 2. Extract genome columns
# Panaroo metadata columns are usually: Gene, Non-unique Gene name, Annotation
# We extract everything from column 4 onwards
pa_matrix <- pan_data[, 4:ncol(pan_data)]

# 3. Convert to binary matrix (1 = present, 0 = absent)
message("Converting to binary matrix...")
pa_matrix_binary <- ifelse(pa_matrix != "", 1, 0)

# 4. Transpose for micropan
# Micropan requires: Rows = Genomes, Columns = Gene Clusters
pa_matrix_final <- t(pa_matrix_binary)
class(pa_matrix_final) <- "numeric"

# 5. Run Heaps' Law estimation
# Using 500 permutations to ensure a robust alpha estimate
message("Running Heaps' Law estimation (500 permutations). This may take a moment...")
heaps_results <- heaps(pa_matrix_final, n.perm = 500)

# 6. Display results
cat("\n--- Pangenome Openness (Heaps' Law) Results ---\n")
print(heaps_results)

# Note: An alpha < 1.0 indicates an 'Open' pangenome.
