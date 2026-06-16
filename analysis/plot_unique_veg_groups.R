library(dplyr)
library(ggplot2)
library(tidyr)

# Load drop term data
drop_df <- read.csv("results_500m_inference/drop_term_inference_results.csv", stringsAsFactors = FALSE)
filtered_spp <- read.csv("results/filtered_native_terrestrial_species.csv", stringsAsFactors = FALSE)
conv <- read.csv("results_500m_inference/convergence_report.csv", stringsAsFactors = FALSE)
valid_species <- filtered_spp$species[filtered_spp$species %in% conv$species[conv$converged == TRUE]]

# Process data
veg_groups <- drop_df %>%
  filter(species %in% valid_species) %>%
  mutate(
    `Functional Identity` = pmax(0, unique_cwm) * 100,
    `Functional Diversity` = pmax(0, unique_fd) * 100,
    `Structure` = pmax(0, unique_structure) * 100,
    `Taxonomic` = pmax(0, unique_richness) * 100
  ) %>%
  select(species, `Functional Identity`, `Functional Diversity`, `Structure`, `Taxonomic`) %>%
  pivot_longer(cols = -species, names_to = "Veg_Group", values_to = "Unique_Dev_Perc")

# Calculate means
mean_veg <- veg_groups %>%
  group_by(Veg_Group) %>%
  summarise(Mean_Unique_Dev = mean(Unique_Dev_Perc, na.rm = TRUE)) %>%
  mutate(Plot_Group = "Unique Vegetation Variance\n(Drop-Term Analysis)")

# Order levels
mean_veg$Veg_Group <- factor(mean_veg$Veg_Group, levels = c("Functional Identity", "Functional Diversity", "Structure", "Taxonomic"))

# Plot
p <- ggplot(mean_veg, aes(x = Plot_Group, y = Mean_Unique_Dev, fill = Veg_Group)) +
  geom_bar(stat = "identity", position = position_stack(), color = "white", linewidth = 0.3, width = 0.5) +
  scale_fill_manual(values = c("Functional Identity" = "#99C24D", 
                               "Functional Diversity" = "#B5E48C", 
                               "Structure" = "#2A6F37", 
                               "Taxonomic" = "#E8D33F")) +
  theme_minimal() +
  labs(title = "Unique Explanatory Power of Vegetation Groups",
       subtitle = "Calculated via drop-term cross-validation across all species",
       x = "", y = "Average Unique Out-of-Sample Deviance (%)",
       fill = "Vegetation Group") +
  theme(legend.position = "right",
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=10, color="gray30"))

dir.create("results_500m_inference/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results_500m_inference/figures/figure_unique_veg_groups.png", p, width = 8, height = 7, dpi = 300)
cat("Successfully generated unique veg group plot.\n")
