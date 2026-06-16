library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# 1. Load Data
df_veg <- read.csv("results_500m_inference/individual_veg_drop_results.csv", stringsAsFactors = FALSE)
df_land <- read.csv("results_500m_inference/landscape_drop_results.csv", stringsAsFactors = FALSE)

# Inner join to ensure we only use species that have both results computed
df <- df_veg %>% inner_join(df_land, by = "species")

mean_dev_full <- mean(df$dev_full, na.rm = TRUE)

# Pivot long for pure raw deviance
long_df <- df %>%
  select(species, starts_with("unique_")) %>%
  pivot_longer(cols = starts_with("unique_"), names_to = "term_raw", values_to = "unique_dev") %>%
  mutate(term_raw = gsub("unique_", "", term_raw))

# Calculate community averages for BOTH metrics
avg_unique <- long_df %>%
  group_by(term_raw) %>%
  summarise(
    raw_mean = mean(pmax(0, unique_dev) * 100, na.rm = TRUE),
    prop_mean = (mean(pmax(0, unique_dev) * 100, na.rm = TRUE) / mean_dev_full),
    .groups = 'drop'
  )

# 3. Map to pretty names and groups
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
  "landscape" = "Summed Landscape"
)

avg_unique <- avg_unique %>%
  mutate(
    var_group = case_when(
      term_raw %in% c("FDiv", "FRic", "FEve") ~ "Vegetation: Functional Diversity",
      term_raw %in% c("height_cwm", "LA_cwm", "wooddensity_cwm") ~ "Vegetation: Functional Identity",
      term_raw %in% c("cover", "prop_canopy", "mean_fhd") ~ "Vegetation: Structure",
      term_raw %in% c("sp_richness", "gini_simpson", "pielou_evenness") ~ "Vegetation: Taxonomic Diversity",
      term_raw == "landscape" ~ "Landscape Context"
    ),
    term_clean = name_map[term_raw]
  )

# Remove any NAs (in case there are unmatched terms)
avg_unique <- avg_unique %>% filter(!is.na(var_group))

# Ordering & Plotting (using raw_mean for sorting)
group_totals <- avg_unique %>%
  group_by(var_group) %>%
  summarise(total = sum(raw_mean), .groups = 'drop') %>%
  arrange(desc(total))

avg_ordered <- avg_unique %>%
  left_join(group_totals, by = "var_group") %>%
  arrange(desc(total), desc(raw_mean))

ordered_terms <- unique(avg_ordered$term_clean)

gen_shades <- function(base_color, n) {
  if (n == 1) return(base_color)
  colorRampPalette(c("white", base_color))(n + 2)[3:(n + 2)]
}

base_colors <- list(
  "Landscape Context" = "#E9C46A",
  "Vegetation: Functional Identity" = "#3D5A80",
  "Vegetation: Structure" = "#E07A5F",
  "Vegetation: Taxonomic Diversity" = "#81B29A",
  "Vegetation: Functional Diversity" = "#98C1D9"
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

avg_ordered$term_clean <- factor(avg_ordered$term_clean, levels = rev(ordered_terms))

# Create the two plot versions
p_raw <- ggplot(avg_ordered %>% mutate(main_group = "Pure Unique Deviance\n(All Species)"), aes(x = main_group, y = raw_mean, fill = term_clean)) +
  geom_bar(stat = "identity", position = position_stack(), color = "white", linewidth = 0.3, width = 0.5) +
  scale_fill_manual(values = term_colors) +
  theme_minimal() +
  labs(x = "", y = "Average Unique Out-of-Sample Deviance (%)", fill = NULL) +
  theme(legend.position = "right", 
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

p_prop <- ggplot(avg_ordered %>% mutate(main_group = "Proportional Unique Deviance\n(All Species)"), aes(x = main_group, y = prop_mean, fill = term_clean)) +
  geom_bar(stat = "identity", position = position_stack(), color = "white", linewidth = 0.3, width = 0.5) +
  scale_fill_manual(values = term_colors) +
  theme_minimal() +
  labs(x = "", y = "Unique Deviance Explained (% of Total Model Power)", fill = NULL) +
  theme(legend.position = "right", 
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

dir.create("results_500m_inference/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results_500m_inference/figures/figure_unique_veg_and_landscape_raw.png", p_raw, width = 10, height = 7, dpi=300)
ggsave("results_500m_inference/figures/figure_unique_veg_and_landscape_proportional.png", p_prop, width = 10, height = 7, dpi=300)

cat("Successfully generated both raw and proportional unique veg + landscape plots.\n")
