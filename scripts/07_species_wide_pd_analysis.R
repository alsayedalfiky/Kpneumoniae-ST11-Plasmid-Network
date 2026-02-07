# Species-Wide Phylogenetic Diversity (PD) Calculation
# Quantitative assessment of total branch length for the K. pneumoniae cohort

library(ape)

# 1. Load the core-genome phylogeny
# Tree file generated from IQ-TREE ML reconstruction
tree_path <- "results/k_pneumoniae_tree.treefile"
tree <- read.tree(tree_path)

# 2. Verify tree structure
print(tree)

# 3. Calculate Faith's Phylogenetic Diversity (PD)
# Defined as the sum of all branch lengths in the phylogeny
pd_value <- sum(tree$edge.length)

# 4. Output results
cat("\n--- Total Phylogenetic Diversity ---\n")
cat("Dataset: 1,010 K. pneumoniae isolates\n")
cat("Total Faith's PD:", pd_value, "\n")
