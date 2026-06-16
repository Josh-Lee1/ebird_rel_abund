library(dplyr)
library(ggplot2)
library(tidyr)

# ── 1. Load Data ──────────────────────────────────────────────────────────────
filtered_spp <- read.csv("results/filtered_native_terrestrial_species.csv", stringsAsFactors = FALSE)
conv <- read.csv("results_500m_inference/convergence_report.csv", stringsAsFactors = FALSE)
valid_species <- filtered_spp$species[filtered_spp$species %in% conv$species[conv$converged == TRUE]]

imp_df <- read.csv("results_500m_inference/aggregated_importance.csv", stringsAsFactors = FALSE)
drop_df <- read.csv("results_500m_inference/drop_term_inference_results.csv", stringsAsFactors = FALSE)

# ── 2. Calculate Species Total Vegetation Deviance ────────────────────────────
species_weights <- drop_df %>%
  filter(species %in% valid_species) %>%
  mutate(
    Veg_Dev_Weight = (pmax(0, unique_cwm) + 
                     pmax(0, unique_structure) + 
                     pmax(0, unique_richness) + 
                     pmax(0, unique_fd)) * 100
  ) %>%
  select(species, Veg_Dev_Weight)

# ── 3. Process F-Statistics & Apply Weights ───────────────────────────────────
name_map <- c(
  "FDiv" = "Functional Divergence",
  "FEve" = "Functional Evenness",
  "FRic" = "Functional Richness",
  "gini_simpson" = "Taxonomic Diversity",
  "pielou_evenness" = "Taxonomic Evenness",
  "sp_richness" = "Taxonomic Richness",
  "mean_fhd" = "Foliage Height Diversity",
  "prop_canopy" = "Proportion Canopy",
  "cover" = "Percentage Cover",
  "LA_cwm" = "Leaf Area",
  "wooddensity_cwm" = "Stem Density",
  "height_cwm" = "Plant Height",
  "lc_built" = "Proportion Built",
  "water_occ" = "Water Availability",
  "pop_density" = "Human population",
  "elevation" = "Elevation",
  "nightlights" = "Nightlights",
  "lc_trees" = "Proportion Trees",
  "temp_annual" = "Annual Temperature",
  "precip_annual" = "Annual Precipitation"
)

imp_clean <- imp_df %>%
  filter(species %in% valid_species) %>%
  filter(category != "effort", !grepl("duration|distance|time|protocol|number", term, ignore.case = TRUE)) %>%
  group_by(species, term) %>%
  filter(timestamp == max(timestamp)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    term_raw = term %>%
      gsub("^s\\(|\\)$", "", .) %>%
      gsub("log1p\\(|log\\(", "", .) %>%
      gsub("\\)", "", .),
    var_group = case_when(
      term_raw %in% c("FDiv", "FRic", "FEve") ~ "Vegetation: Functional Diversity",
      term_raw %in% c("height_cwm", "LA_cwm", "wooddensity_cwm") ~ "Vegetation: Functional Identity",
      term_raw %in% c("cover", "prop_canopy", "mean_fhd") ~ "Vegetation: Structural Metrics",
      term_raw %in% c("sp_richness", "gini_simpson", "pielou_evenness") ~ "Vegetation: Species Diversity",
      TRUE ~ "Landscape"
    )
  ) %>%
  filter(!term_raw %in% c("LMA_cwm", "FDis", "Q")) %>%
  filter(var_group != "Landscape") %>% # REMOVE non-vegetation variables!
  mutate(term_clean = coalesce(name_map[term_raw], term_raw))

weighted_imp <- imp_clean %>%
  group_by(species) %>%
  mutate(total_f = sum(f_stat, na.rm = TRUE),
         prop_importance = f_stat / total_f) %>%
  ungroup() %>%
  left_join(species_weights, by = "species") %>%
  mutate(
    weighted_importance = prop_importance * Veg_Dev_Weight
  )

# ── 4. Calculate Community-Wide Deviance-Weighted Averages ────────────────────
avg_weighted <- weighted_imp %>%
  group_by(var_group, term_clean) %>%
  summarise(mean_weighted_imp = mean(weighted_importance, na.rm = TRUE), .groups = 'drop') %>%
  mutate(main_group = "Deviance-Weighted Importance\n(All Species)")

# ── 5. Ordering & Plotting ────────────────────────────────────────────────────
group_totals <- avg_weighted %>%
  group_by(var_group) %>%
  summarise(total = sum(mean_weighted_imp), .groups = 'drop') %>%
  arrange(desc(total))

avg_ordered <- avg_weighted %>%
  left_join(group_totals, by = "var_group") %>%
  arrange(desc(total), desc(mean_weighted_imp))

ordered_terms <- unique(avg_ordered$term_clean)

gen_shades <- function(base_color, n) {
  colorRampPalette(c("white", base_color))(n + 2)[3:(n + 2)]
}

base_colors <- list(
  "Vegetation: Functional Identity" = "#99C24D",
  "Vegetation: Structural Metrics" = "#2A6F37",
  "Vegetation: Species Diversity" = "#E8D33F",
  "Vegetation: Functional Diversity" = "#B5E48C"
)

term_colors <- c()
for (grp in group_totals$var_group) {
  terms_in_grp <- avg_ordered$term_clean[avg_ordered$var_group == grp]
  if (length(terms_in_grp) > 0) {
    shades <- rev(gen_shades(base_colors[[grp]], length(terms_in_grp)))
    names(shades) <- terms_in_grp
    term_colors <- c(term_colors, shades)
  }
}

avg_weighted$term_clean <- factor(avg_weighted$term_clean, levels = rev(ordered_terms))

p_bar <- ggplot(avg_weighted, aes(x = main_group, y = mean_weighted_imp, fill = term_clean)) +
  geom_bar(stat = "identity", position = position_stack(), color = "white", linewidth = 0.3, width = 0.5) +
  scale_fill_manual(values = term_colors) +
  theme_minimal() +
  labs(title = "Deviance-Weighted Variable Importance (Inference Models)", 
       subtitle = "Proportional F-stats weighted by each species' Absolute Unique Vegetation Deviance (% Deviance)",
       x = "", y = "Average Absolute % Deviance Explained") +
  theme(legend.position = "right", 
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=10, color="gray30"))

ggsave("results_500m_inference/figures/figure_weighted_importance_inference.png", p_bar, width = 10, height = 7, dpi=300)
cat("Plot saved to results_500m_inference/figures/figure_weighted_importance_inference.png\n")
