setwd("c:/Users/90957140/OneDrive - Western Sydney University/R/bird_FD_ebirdabund")

library(dplyr)
library(ggplot2)
library(tidyr)

# Load data
imp <- read.csv("results_500m_inference/aggregated_importance.csv", stringsAsFactors = FALSE)
conv <- read.csv("results_500m_inference/convergence_report.csv", stringsAsFactors = FALSE)
valid_species <- conv$species[conv$converged == TRUE]

imp <- imp %>% filter(species %in% valid_species)

# Define term categories explicitly
func_id_terms <- c("s(LA_cwm)", "s(wooddensity_cwm)", "s(height_cwm)")
func_div_terms <- c("s(FDiv)", "s(FEve)", "s(FRic)")
struct_terms <- c("s(mean_fhd)", "s(prop_canopy)", "s(cover)")
taxa_terms <- c("s(gini_simpson)", "s(pielou_evenness)", "s(sp_richness)")

effort_terms <- c("s(log1p(effort_distance_km))", "s(log1p(duration_minutes))", "s(day_of_year)", 
                  "s(time_observations_started)", "s(observer_expertise)", "s(number_observers)",
                  "s(nightlights)", "s(pop_density)")
env_terms <- c("s(elevation)", "s(lc_water)", "s(water_occ)", "s(lc_built)", "s(lc_trees)", "s(temp_annual)", "s(precip_annual)")

imp <- imp %>%
  mutate(MacroCategory = case_when(
    term %in% func_id_terms ~ "Functional Identity",
    term %in% func_div_terms ~ "Functional Diversity",
    term %in% struct_terms ~ "Structure",
    term %in% taxa_terms ~ "Taxonomic",
    term %in% env_terms ~ "Landscape & Climate",
    TRUE ~ "Effort" # Everything else is Effort
  ))

# Aggregate F-statistics
agg_imp <- imp %>%
  group_by(species, MacroCategory) %>%
  summarise(Total_F = sum(f_stat, na.rm = TRUE), .groups = "drop") %>%
  group_by(species) %>%
  mutate(Prop_F = Total_F / sum(Total_F, na.rm = TRUE)) %>%
  ungroup()

# Get mean proportion for sorting (sorting by Functional + Structure + Taxa to keep 'total vegetation' order)
sort_order <- agg_imp %>%
  filter(MacroCategory %in% c("Functional Identity", "Functional Diversity", "Structure", "Taxonomic")) %>%
  group_by(species) %>%
  summarise(Total_Veg_Prop = sum(Prop_F, na.rm = TRUE)) %>%
  arrange(desc(Total_Veg_Prop)) %>%
  pull(species)

agg_imp$species <- factor(agg_imp$species, levels = sort_order)
agg_imp$MacroCategory <- factor(agg_imp$MacroCategory, levels = c("Functional Identity", "Functional Diversity", "Structure", "Taxonomic", "Landscape & Climate", "Effort"))

# Plot
p <- ggplot(agg_imp, aes(x = species, y = Prop_F, fill = MacroCategory)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("Functional Identity" = "#3D5A80", "Functional Diversity" = "#98C1D9", "Structure" = "#E07A5F", "Taxonomic" = "#81B29A", "Landscape & Climate" = "#F2CC8F", "Effort" = "#BDBDBD")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Species (Ordered by Vegetation Importance)",
       y = "Proportion of Total F-statistic",
       fill = "Category")

dir.create("results_500m_inference/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results_500m_inference/figures/figure_1_macro_categories_split.png", p, width = 10, height = 6, dpi=300)

suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))
p_plotly <- ggplotly(p)
saveWidget(p_plotly, "results_500m_inference/figures/figure_1_macro_categories_split_interactive.html", selfcontained = FALSE)
