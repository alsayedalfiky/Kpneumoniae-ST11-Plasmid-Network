#R
# Plotting just the STs on one shared axis to see the 'gap'
p_compare <- ggplot(st_only_curves, aes(x = Genomes, y = Mean_PD, color = Group)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Comparison of Evolutionary Trajectories",
    subtitle = "ST23 diversity contrasts sharply with clonal ST11/ST258 expansion",
    x = "Genomes Sampled",
    y = "Faith's PD"
  )

print(p_compare)
ggsave("Figure1_ST_Comparison_NoTotal.pdf", width = 8, height = 6)
