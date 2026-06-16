setwd("c:/Users/90957140/OneDrive - Western Sydney University/R/bird_FD_ebirdabund")

library(dplyr)
library(ggplot2)
library(patchwork)
library(broom)

# Load data
drop_df <- read.csv("results_500m_inference/individual_veg_drop_results.csv", stringsAsFactors = FALSE)
traits <- read.csv("data/all_bird_traits.csv", stringsAsFactors = FALSE)

# Extract and clean traits
traits_clean <- traits %>%
  select(
    species = COMMON.NAME,
    Diet = Trophic.Niche,
    Lifestyle = Primary.Lifestyle,
    BodyMass = Mass,
    RangeSize = X208_Area_of_suitable_climate_inside_range_14,
    HWI = Hand.Wing.Index,
    Habitat = Habitat
  ) %>%
  distinct(species, .keep_all = TRUE)

# Merge and compute Functional Importance
analysis_df <- drop_df %>%
  mutate(
    Total_Veg_Dev_Raw = pmax(0, unique_FDiv) + pmax(0, unique_FEve) + pmax(0, unique_FRic) + pmax(0, unique_LA_cwm) + pmax(0, unique_wooddensity_cwm) + pmax(0, unique_height_cwm) + pmax(0, unique_cover) + pmax(0, unique_mean_fhd) + pmax(0, unique_prop_canopy) + pmax(0, unique_sp_richness) + pmax(0, unique_pielou_evenness) + pmax(0, unique_gini_simpson),
    Func_Dev_Raw = pmax(0, unique_FDiv) + pmax(0, unique_FEve) + pmax(0, unique_FRic) + pmax(0, unique_LA_cwm) + pmax(0, unique_wooddensity_cwm) + pmax(0, unique_height_cwm)
  ) %>%
  filter(Total_Veg_Dev_Raw > 0) %>% # Subset to reliant species
  mutate(Func_Dev = (Func_Dev_Raw / Total_Veg_Dev_Raw) * 100) %>% # Percentage of VEGETATION variation from functional traits
  left_join(traits_clean, by = "species") %>%
  filter(!is.na(Func_Dev))

# Run simple linear models to get p-values for the plots
# Diet
model_diet <- lm(Func_Dev ~ Diet, data = analysis_df)
pval_diet <- glance(model_diet)$p.value

# Lifestyle
model_lifestyle <- lm(Func_Dev ~ Lifestyle, data = analysis_df)
pval_lifestyle <- glance(model_lifestyle)$p.value

# Habitat
model_habitat <- lm(Func_Dev ~ Habitat, data = analysis_df %>% filter(!is.na(Habitat) & Habitat != ""))
pval_habitat <- glance(model_habitat)$p.value

# Body Mass
model_mass <- lm(Func_Dev ~ log10(BodyMass), data = analysis_df %>% filter(!is.na(BodyMass) & BodyMass > 0))
pval_mass <- glance(model_mass)$p.value

# Hand Wing Index
model_hwi <- lm(Func_Dev ~ HWI, data = analysis_df %>% filter(!is.na(HWI)))
pval_hwi <- glance(model_hwi)$p.value

# Plot A: Diet
p_diet <- ggplot(analysis_df %>% filter(!is.na(Diet) & Diet != ""), aes(x = reorder(Diet, Func_Dev, FUN = median), y = Func_Dev)) +
  geom_boxplot(fill = "#3D5A80", alpha = 0.7) +
  theme_minimal() +
  coord_flip() +
  labs(title = "A) Functional Importance by Trophic Niche",
       x = "Trophic Niche",
       y = "% of Veg Variance from Func Traits")

# Plot B: Habitat
p_hab <- ggplot(analysis_df %>% filter(!is.na(Habitat) & Habitat != ""), aes(x = reorder(Habitat, Func_Dev, FUN = median), y = Func_Dev)) +
  geom_boxplot(fill = "#3D5A80", alpha = 0.7) +
  theme_minimal() +
  coord_flip() +
  labs(title = "B) Functional Importance by Habitat",
       x = "Primary Habitat",
       y = "% of Veg Variance from Func Traits")

# Plot C: Body Mass
p_mass <- ggplot(analysis_df %>% filter(!is.na(BodyMass) & BodyMass > 0), aes(x = BodyMass, y = Func_Dev)) +
  geom_point(color = "#3D5A80", alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "C) Functional Importance vs Body Mass",
       x = "Body Mass (g, log scale)",
       y = "% of Veg Variance from Func Traits")

# Plot D: Hand Wing Index
p_hwi <- ggplot(analysis_df %>% filter(!is.na(HWI)), aes(x = HWI, y = Func_Dev)) +
  geom_point(color = "#3D5A80", alpha = 0.6) +
  theme_minimal() +
  labs(title = "D) Functional Importance vs Hand Wing Index",
       x = "Hand Wing Index (Dispersal Ability)",
       y = "% of Veg Variance from Func Traits")

# Combine plots using patchwork
combined_plot <- (p_diet | p_hab) / (p_mass | p_hwi)

dir.create("results_500m_inference/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results_500m_inference/figures/figure_4_trait_analysis.png", combined_plot, width = 14, height = 10, dpi = 300)

cat(sprintf("\n--- P-VALUES FOR RESULTS REPORT ---\n"))
cat(sprintf("Diet ANOVA p-value: %.3f\n", pval_diet))
cat(sprintf("Lifestyle ANOVA p-value: %.3f\n", pval_lifestyle))
cat(sprintf("Habitat ANOVA p-value: %.3f\n", pval_habitat))
cat(sprintf("Body Mass LM p-value: %.3f\n", pval_mass))
cat(sprintf("Hand Wing Index LM p-value: %.3f\n", pval_hwi))

cat("\nSuccessfully generated Figure 4. Check results_500m_inference/figures/figure_4_trait_analysis.png\n")
