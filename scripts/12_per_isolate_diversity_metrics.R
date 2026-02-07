# Lineage-Specific Diversity Normalization
# Quantitative comparison of phylogenetic contribution per isolate (ST11 vs ST23)

library(ape)
library(dplyr)

# 1. Quantify ST11 Diversity (n=304)
# Calculation based on the previously pruned ST11 subtree
st11_exact_pd <- sum(st11_tree$edge.length)
st11_ratio <- st11_exact_pd / 304

# 2. Quantify ST23 Diversity (n=43)
# Identification of ST23 isolates using Kleborate and MLST metadata
# Intersection ensures IDs exist as tips in the phylogenetic tree
st23_ids <- intersect(
  metadata %>% filter(st_kleborate == 23 | st_mlst == 23) %>% pull(id), 
  tree$tip.label
)

# Prune the master tree to include only ST23 isolates
st23_tree <- keep.tip(tree, st23_ids)
st23_exact_pd <- sum(st23_tree$edge.length)
st23_ratio <- st23_exact_pd / 43

# 3. Output Analysis Results
# Results provide a normalized metric to compare lineages of different sample sizes
cat("\n--- Lineage-Specific Diversity Metrics ---\n")
cat("ST11 (n=304) Total PD:", round(st11_exact_pd, 4), "| Per Isolate:", round(st11_ratio, 6), "\n")
cat("ST23 (n=43)  Total PD:", round(st23_exact_pd, 4), "| Per Isolate:", round(st23_ratio, 6), "\n")

# 4. Export Summary Data
summary_stats <- data.frame(
  Lineage = c("ST11", "ST23"),
  Isolate_Count = c(304, 43),
  Total_PD = c(st11_exact_pd, st23_exact_pd),
  PD_Per_Isolate = c(st11_ratio, st23_ratio)
)

write.csv(summary_stats, "results/st_diversity_summary.csv", row.names = FALSE)
