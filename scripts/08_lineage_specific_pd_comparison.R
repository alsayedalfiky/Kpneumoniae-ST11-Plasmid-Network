# Lineage-Specific Phylogenetic Diversity (PD) Comparison
# Comparison of ST11 sub-lineage diversity relative to the total population

library(ape)

# 1. Load phylogeny and metadata
# Ensure the tree file and metadata CSV are in the results directory
tree <- read.tree("results/k_pneumoniae_tree.treefile")
metadata <- read.csv("results/st_comparison_report.csv")

# 2. Extract ST11 isolates
# Filtering based on Kleborate-resolved Sequence Types
st11_ids <- metadata$id[metadata$st_kleborate == 11]

# 3. Prune the phylogeny
# keep.tip generates a subtree containing only the specified lineage tips
st11_tree <- keep.tip(tree, as.character(st11_ids))

# 4. Quantify Diversity Contribution
# Faith's PD calculated as the sum of branch lengths
total_pd <- sum(tree$edge.length)    
st11_pd <- sum(st11_tree$edge.length)
div_ratio <- (st11_pd / total_pd) * 100

# 5. Output Summary Metrics
cat("\n--- Diversity Comparison Analysis ---\n")
cat("Species-wide Total PD: ", total_pd, "\n")
cat("ST11 Lineage PD:       ", st11_pd, "\n")
cat("ST11 % Contribution:   ", round(div_ratio, 4), "%\n")

# Save the pruned tree for downstream rarefaction analysis
# saveRDS(st11_tree, "results/st11_pruned_tree.rds")
