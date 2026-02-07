#!/usr/bin/env Rscript

library(ape)
library(dplyr)

# 1. Load Data
message("Loading phylogeny and metadata...")
tree <- read.tree("results/k_pneumoniae_tree.treefile")
metadata <- read.csv("results/st_comparison_report.csv")

# Ensure metadata IDs match tree tip labels
metadata <- metadata[metadata$id %in% tree$tip.label, ]

# 2. Total Species-Wide Diversity
total_pd <- sum(tree$edge.length)

# 3. Define Function for Lineage Diversity
calculate_lineage_pd <- function(target_st, label) {
  # Filter IDs (handling both Kleborate and MLST columns)
  ids <- metadata %>% 
    filter(st_kleborate == target_st | st_mlst == target_st) %>% 
    pull(id)
  
  # Prune tree to lineage
  sub_tree <- keep.tip(tree, as.character(ids))
  lineage_pd <- sum(sub_tree$edge.length)
  
  return(data.frame(
    Lineage = label,
    Count = length(ids),
    Total_PD = lineage_pd,
    Per_Isolate_PD = lineage_pd / length(ids)
  ))
}

# 4. Run Analysis for Key Groups
st11_metrics <- calculate_lineage_pd(11, "ST11")
st23_metrics <- calculate_lineage_pd(23, "ST23")

# 5. Calculate the Clonal Sweep Metric
st11_diversity_contribution <- (st11_metrics$Total_PD / total_pd) * 100

# 6. Output Results for Manuscript
cat("\n--- Evolutionary Metrics ---\n")
cat("Species-wide Faith's PD:", total_pd, "\n")
cat("ST11 Contribution to Diversity:", round(st11_diversity_contribution, 2), "%\n\n")

results_summary <- rbind(st11_metrics, st23_metrics)
print(results_summary)
