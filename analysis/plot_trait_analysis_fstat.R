setwd("c:/Users/90957140/OneDrive - Western Sydney University/R/bird_FD_ebirdabund")

library(dplyr)
library(ggplot2)
library(patchwork)
library(broom)

# Load data
df <- read.csv("results_500m_inference/aggregated_importance.csv", stringsAsFactors = FALSE)
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

# Define functional terms
func_terms <- c("s(LA_cwm)", "s(wooddensity_cwm)", "s(height_cwm)", "s(FDiv)", "s(FEve)", "s(FRic)")

# Compute F-Statistic Proportions
analysis_df <- df %>%
  filter(is_vegetation == TRUE) %>%
  group_by(species) %>%
  summarise(
    Total_Veg_F = sum(f_stat, na.rm = TRUE),
    Func_F_Raw = sum(f_stat[term %in% func_terms], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Total_Veg_F > 0) %>% # Subset to reliant species
  mutate(Func_Dev = (Func_F_Raw / Total_Veg_F) * 100) %>% # Percentage of VEGETATION F-stat from functional traits
  left_join(traits_clean, by = "species") %>%
  filter(!is.na(Func_Dev))

# Run simple linear models to get p-values for the plots
# Diet
model_diet <- lm(Func_Dev ~ Diet, data = analysis_df)
pval_diet <- glance(model_diet)$p.value

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
  geom_boxplot(fill = "#006E90", alpha = 0.7) +
  theme_minimal() +
  coord_flip() +
  labs(title = "A) Functional Importance by Diet",
       x = "Trophic Niche",
       y = "% of Veg F-statistic from Func Traits")

# Plot B: Habitat
p_hab <- ggplot(analysis_df %>% filter(!is.na(Habitat) & Habitat != ""), aes(x = reorder(Habitat, Func_Dev, FUN = median), y = Func_Dev)) +
  geom_boxplot(fill = "#99C24D", alpha = 0.7) +
  theme_minimal() +
  coord_flip() +
  labs(title = "B) Functional Importance by Habitat",
       x = "Primary Habitat",
       y = "% of Veg F-statistic from Func Traits")

# Plot C: Body Mass
p_mass <- ggplot(analysis_df %>% filter(!is.na(BodyMass) & BodyMass > 0), aes(x = BodyMass, y = Func_Dev)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "C) Functional Importance vs Body Mass",
       x = "Body Mass (g, log scale)",
       y = "% of Veg F-statistic from Func Traits")

# Plot D: Hand Wing Index
p_hwi <- ggplot(analysis_df %>% filter(!is.na(HWI)), aes(x = HWI, y = Func_Dev)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "D) Functional Importance vs Hand Wing Index",
       x = "Hand Wing Index (Dispersal Ability)",
       y = "% of Veg F-statistic from Func Traits")

# Combine plots using patchwork
combined_plot <- (p_diet | p_hab) / (p_mass | p_hwi)

dir.create("results_500m_inference/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results_500m_inference/figures/figure_4_trait_analysis_fstat.png", combined_plot, width = 14, height = 10, dpi = 300)

cat(sprintf("\n--- P-VALUES FOR F-STAT RESULTS ---\n"))
cat(sprintf("Diet ANOVA p-value: %.3f\n", pval_diet))
cat(sprintf("Habitat ANOVA p-value: %.3f\n", pval_habitat))
cat(sprintf("Body Mass LM p-value: %.3f\n", pval_mass))
cat(sprintf("Hand Wing Index LM p-value: %.3f\n", pval_hwi))

cat("\nSuccessfully generated Figure 4. Check results_500m_inference/figures/figure_4_trait_analysis_fstat.png\n")
