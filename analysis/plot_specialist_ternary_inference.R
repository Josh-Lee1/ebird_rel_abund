setwd("c:/Users/90957140/OneDrive - Western Sydney University/R/bird_FD_ebirdabund")

library(dplyr)
library(ggplot2)
library(ggtern)

# Load data
filtered_spp <- read.csv("results/filtered_native_terrestrial_species.csv", stringsAsFactors = FALSE)
conv <- read.csv("results_500m_inference/convergence_report.csv", stringsAsFactors = FALSE)
valid_species <- filtered_spp$species[filtered_spp$species %in% conv$species[conv$converged == TRUE]]

drop_df <- read.csv("results_500m_inference/individual_veg_drop_results.csv", stringsAsFactors = FALSE)
traits <- read.csv("data/all_bird_traits.csv", stringsAsFactors = FALSE)

traits_clean <- traits %>% 
  select(COMMON.NAME, Trophic.Niche, Habitat) %>%
  rename(species = COMMON.NAME, Diet = Trophic.Niche) %>%
  distinct(species, .keep_all = TRUE)

# Find vegetation-reliant species (> 0 unique deviance)
reliant_spp <- drop_df %>%
  filter(species %in% valid_species) %>%
  mutate(Total_Veg_Dev_Raw = pmax(0, unique_FDiv) + pmax(0, unique_FEve) + pmax(0, unique_FRic) + pmax(0, unique_LA_cwm) + pmax(0, unique_wooddensity_cwm) + pmax(0, unique_height_cwm) + pmax(0, unique_cover) + pmax(0, unique_mean_fhd) + pmax(0, unique_prop_canopy) + pmax(0, unique_sp_richness) + pmax(0, unique_pielou_evenness) + pmax(0, unique_gini_simpson)) %>%
  filter(Total_Veg_Dev_Raw > 0) %>%
  pull(species)

# Prepare data for ternary plot using ONLY reliant species
bubble_df <- drop_df %>% 
  filter(species %in% reliant_spp) %>%
  left_join(traits_clean, by = "species") %>%
  filter(!is.na(Diet) & Diet != "")

tern_df <- bubble_df %>%
  mutate(
    Func_Dev = pmax(0, unique_FDiv) + pmax(0, unique_FEve) + pmax(0, unique_FRic) + pmax(0, unique_LA_cwm) + pmax(0, unique_wooddensity_cwm) + pmax(0, unique_height_cwm),
    Struct_Dev = pmax(0, unique_cover) + pmax(0, unique_mean_fhd) + pmax(0, unique_prop_canopy),
    Rich_Dev = pmax(0, unique_sp_richness) + pmax(0, unique_pielou_evenness) + pmax(0, unique_gini_simpson)
  ) %>%
  mutate(
    Total_Sum = Func_Dev + Struct_Dev + Rich_Dev,
    Func_Prop = Func_Dev / Total_Sum,
    Struct_Prop = Struct_Dev / Total_Sum,
    Rich_Prop = Rich_Dev / Total_Sum
  ) %>%
  filter(!is.na(Total_Sum) & Total_Sum > 0.0001)

cat(sprintf("Plotting %d species on the ternary plot (out of %d reliant species)\n", nrow(tern_df), length(reliant_spp)))

p3 <- ggtern(data = tern_df, aes(x = Struct_Prop, y = Func_Prop, z = Rich_Prop, color = Diet)) +
  geom_point(alpha = 0.5, size = 2.5) +
  theme_rgbw() +
  scale_color_brewer(palette = "Set2") +
  theme(
    tern.axis.title.L = element_text(color = "#E07A5F", face = "bold"),
    tern.axis.line.L = element_line(color = "#E07A5F", linewidth = 1),
    tern.axis.text.L = element_text(color = "#E07A5F"),
    tern.axis.ticks.major.L = element_line(color = "#E07A5F"),
    tern.axis.arrow.L = element_line(color = "#E07A5F", linewidth = 1),
    tern.axis.arrow.text.L = element_text(color = "#E07A5F"),
    tern.panel.grid.major.L = element_line(color = "#E07A5F", linetype = "dashed", linewidth = 0.4),
    tern.panel.grid.minor.L = element_blank(),
    
    tern.axis.title.T = element_text(color = "#3D5A80", face = "bold"),
    tern.axis.line.T = element_line(color = "#3D5A80", linewidth = 1),
    tern.axis.text.T = element_text(color = "#3D5A80"),
    tern.axis.ticks.major.T = element_line(color = "#3D5A80"),
    tern.axis.arrow.T = element_line(color = "#3D5A80", linewidth = 1),
    tern.axis.arrow.text.T = element_text(color = "#3D5A80"),
    tern.panel.grid.major.T = element_line(color = "#3D5A80", linetype = "dashed", linewidth = 0.4),
    tern.panel.grid.minor.T = element_blank(),
    
    tern.axis.title.R = element_text(color = "#81B29A", face = "bold"),
    tern.axis.line.R = element_line(color = "#81B29A", linewidth = 1),
    tern.axis.text.R = element_text(color = "#81B29A"),
    tern.axis.ticks.major.R = element_line(color = "#81B29A"),
    tern.axis.arrow.R = element_line(color = "#81B29A", linewidth = 1),
    tern.axis.arrow.text.R = element_text(color = "#81B29A"),
    tern.panel.grid.major.R = element_line(color = "#81B29A", linetype = "dashed", linewidth = 0.4),
    tern.panel.grid.minor.R = element_blank(),
    
    legend.key = element_blank()
  ) +
  labs(x = "Structure", y = "Functional Traits", z = "Taxonomic Diversity",
       color = "Trophic Niche")

dir.create("results_500m_inference/figures/dropterm_deepdive", showWarnings = FALSE, recursive = TRUE)
ggsave("results_500m_inference/figures/dropterm_deepdive/specialist_ternary_inference.png", p3, width = 9, height = 7, dpi=600)
