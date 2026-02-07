# Rarefaction and Visualization of Phylogenetic Diversity (PD)
# This script performs permutation-based sampling to compare ST11 clonal expansion 
# against the species-wide diversity of K. pneumoniae.

library(picante)
library(ggplot2)
library(dplyr)

# 1. Rarefaction Function
# Iteratively samples tips and calculates PD over 100 replicates per step
run_rarefaction <- function(target_tree, group_name, steps) {
  n_reps <- 100
  valid_steps <- steps[steps <= length(target_tree$tip.label)]
  
  message(paste("Starting Rarefaction for:", group_name))
  
  pd_results <- sapply(valid_steps, function(n) {
    replicate(n_reps, {
      sampled_tips <- sample(target_tree$tip.label, n)
      sub_tree <- keep.tip(target_tree, sampled_tips)
      sum(sub_tree$edge.length)
    })
  })
  
  data.frame(
    Genomes = valid_steps,
    Mean_PD = apply(pd_results, 2, mean),
    SD = apply(pd_results, 2, sd),
    Group = group_name
  )
}

# 2. Diversity Analysis
# Requires 'tree' and 'st11_tree' to be pre-loaded in the environment
df_total <- run_rarefaction(tree, "Total Population", seq(10, 1010, by = 50))
df_st11 <- run_rarefaction(st11_tree, "ST11 Lineage", seq(10, length(st11_tree$tip.label), by = 20))

rarefaction_combined <- rbind(df_total, df_st11)

# 3. Figure Generation
p_final <- ggplot(rarefaction_combined, aes(x = Genomes, y = Mean_PD, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = Mean_PD - SD, ymax = Mean_PD + SD), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  theme_bw() +
  scale_color_manual(values = c("Total Population" = "#2c3e50", "ST11 Lineage" = "#e74c3c")) +
  scale_fill_manual(values = c("Total Population" = "#2c3e50", "ST11 Lineage" = "#e74c3c")) +
  labs(
    title = "Comparative Phylogenetic Diversity (Faith's PD)",
    subtitle = "Quantifying the Clonal Expansion of ST11 against Global K. pneumoniae Diversity",
    x = "Number of Genomes Sampled",
    y = "Total Branch Length (Evolutionary Diversity)",
    caption = "Ribbon represents Mean ± SD over 100 permutations per sampling step."
  ) +
  theme(
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

# 4. Export results
ggsave("results/Figure1_PD_Rarefaction.pdf", plot = p_final, width = 8, height = 6, dpi = 300)
write.csv(rarefaction_combined, "results/rarefaction_data_summary.csv", row.names = FALSE)
