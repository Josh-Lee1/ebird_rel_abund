setwd("c:/Users/90957140/OneDrive - Western Sydney University/R/bird_FD_ebirdabund")

library(dplyr)
library(car)   # for VIF

# ── 1. Load benchmark species data ────────────────────────────────────────────
cat("Loading benchmark species data...\n")
bench <- read.csv("data/cleaned_benchmark_data.csv", stringsAsFactors = FALSE)

# Remove records with negative or zero cover (data artefacts)
bench <- bench %>% filter(CoverScore > 0)

# ── 2. Calculate diversity metrics per plot ────────────────────────────────────
# Pielou's Evenness = Shannon H / log(S)
# Gini-Simpson = 1 - sum(p_i^2)
# where p_i = relative cover of species i

cat("Calculating diversity metrics per plot...\n")

diversity_metrics <- bench %>%
  group_by(CensusKey) %>%
  summarise(
    sp_richness_bench = n_distinct(aligned_name[CoverScore > 0]),
    
    # Proportional cover for each species
    .groups = "drop"
  )

# Do it in one step for speed
diversity_metrics <- bench %>%
  group_by(CensusKey) %>%
  summarise(
    # Raw counts
    sp_richness_bench = n_distinct(aligned_name),
    total_cover       = sum(CoverScore, na.rm = TRUE),
    
    # Shannon H and Pielou's Evenness
    shannon_H = {
      p <- CoverScore / sum(CoverScore)
      p <- p[p > 0]
      -sum(p * log(p))
    },
    pielou_evenness = {
      p <- CoverScore / sum(CoverScore)
      p <- p[p > 0]
      H <- -sum(p * log(p))
      S <- length(p)
      if (S > 1) H / log(S) else NA_real_
    },
    
    # Gini-Simpson
    gini_simpson = {
      p <- CoverScore / sum(CoverScore)
      1 - sum(p^2)
    },
    
    .groups = "drop"
  )

cat("Diversity metrics calculated for", nrow(diversity_metrics), "plots.\n")
cat("Sample:\n")
print(head(diversity_metrics, 5))
cat("\nSummary stats:\n")
print(summary(diversity_metrics[, c("sp_richness_bench","shannon_H","pielou_evenness","gini_simpson")]))

# ── 3. Load GEDI structural/FD data and join ─────────────────────────────────
cat("\nLoading GEDI metrics...\n")
gedi <- read.csv("data/bench_GEDI_struct_metrics.csv", stringsAsFactors = FALSE)
cat("GEDI sites:", nrow(gedi), "\n")
cat("Matching key in GEDI: 'site'\n")

gedi_enriched <- gedi %>%
  left_join(diversity_metrics, by = c("site" = "CensusKey"))

n_matched <- sum(!is.na(gedi_enriched$shannon_H))
cat(sprintf("Joined: %d / %d GEDI sites matched to diversity metrics.\n", n_matched, nrow(gedi)))

# Save enriched plot file
write.csv(gedi_enriched, "data/bench_GEDI_struct_metrics_enriched.csv", row.names = FALSE)
cat("Saved: data/bench_GEDI_struct_metrics_enriched.csv\n")

# ── 4. VIF Analysis ───────────────────────────────────────────────────────────
cat("\n=== VIF Analysis ===\n")

# Variables to test (all potential plot predictors)
vif_vars <- c(
  # FD (5)
  "FRic", "FEve", "FDis", "FDiv", "Q",
  # Structural (3)
  "cover", "mean_fhd", "prop_canopy",
  # Diversity (3)
  "sp_richness", "pielou_evenness", "gini_simpson"
)

# Load CWM traits for VIF too
cwm <- read.csv("data/cleaned_veg_w_fd.csv", stringsAsFactors = FALSE)
cat("CWM columns:", paste(names(cwm), collapse=", "), "\n")

# We need to join CWM to GEDI enriched
# We need to join CWM to GEDI enriched - only bring in the _cwm columns to avoid
# column name collisions (CWM file also has FRic, FEve etc.)
cwm_key <- if ("CensusKey" %in% names(cwm)) "CensusKey" else if ("site" %in% names(cwm)) "site" else names(cwm)[1]
cat("CWM join key:", cwm_key, "\n")

# Only keep the _cwm trait columns plus the key
cwm_cols_only <- grep("_cwm$", names(cwm), value = TRUE)
cat("CWM trait columns found:", paste(cwm_cols_only, collapse=", "), "\n")
cwm_slim <- cwm[, c(cwm_key, cwm_cols_only)]
names(cwm_slim)[1] <- "site"  # normalise key name

gedi_full <- gedi_enriched %>%
  left_join(cwm_slim, by = "site")

cwm_cols <- cwm_cols_only
vif_vars_full <- c(vif_vars, cwm_cols)
vif_vars_available <- vif_vars_full[vif_vars_full %in% names(gedi_full)]

# Build a complete-cases matrix for VIF
vif_df <- gedi_full[, vif_vars_available, drop = FALSE]
vif_df <- vif_df[complete.cases(vif_df), ]
cat(sprintf("\nVIF analysis on %d complete-case plots.\n", nrow(vif_df)))

# VIF via a dummy OLS regression: use a numeric index as response,
# regress all predictors against each other
vif_df$dummy_y <- seq_len(nrow(vif_df))
vif_formula <- as.formula(paste("dummy_y ~", paste(vif_vars_available, collapse=" + ")))
vif_model <- lm(vif_formula, data = vif_df)
vif_scores <- car::vif(vif_model)

vif_result <- data.frame(
  variable = names(vif_scores),
  VIF = round(as.numeric(vif_scores), 2)
) %>% arrange(desc(VIF))

cat("\nVIF Scores (all variables):\n")
print(vif_result)

# ── 5. Identify variables to drop ────────────────────────────────────────────
cat("\n=== Selection for 3-per-group ===\n")

# FD group: keep 3 with lowest VIF
fd_vif <- vif_result %>% filter(variable %in% c("FRic","FEve","FDis","FDiv","Q"))
fd_keep <- head(fd_vif %>% arrange(VIF), 3)$variable
fd_drop <- setdiff(c("FRic","FEve","FDis","FDiv","Q"), fd_keep)
cat("FD group VIFs:\n"); print(fd_vif)
cat("FD metrics to KEEP:", paste(fd_keep, collapse=", "), "\n")
cat("FD metrics to DROP:", paste(fd_drop, collapse=", "), "\n")

# CWM group: keep 3 with lowest VIF
cwm_vif <- vif_result %>% filter(variable %in% cwm_cols)
cwm_keep <- head(cwm_vif %>% arrange(VIF), 3)$variable
cwm_drop <- setdiff(cwm_cols, cwm_keep)
cat("\nCWM group VIFs:\n"); print(cwm_vif)
cat("CWM traits to KEEP:", paste(cwm_keep, collapse=", "), "\n")
cat("CWM traits to DROP:", paste(cwm_drop, collapse=", "), "\n")

# Diversity group: sp_richness, pielou_evenness, gini_simpson
div_vif <- vif_result %>% filter(variable %in% c("sp_richness","pielou_evenness","gini_simpson"))
cat("\nDiversity group VIFs:\n"); print(div_vif)

# ── 6. Save summary ───────────────────────────────────────────────────────────
write.csv(vif_result, "results_500m/vif_scores.csv", row.names = FALSE)
cat("\nSaved: results_500m/vif_scores.csv\n")

# Print final variable set for config
final_vars <- c(fd_keep, cwm_keep, "sp_richness", "pielou_evenness", "gini_simpson",
                "cover", "mean_fhd", "prop_canopy")
cat("\n=== FINAL VARIABLE SET (copy into config_500m.yaml) ===\n")
cat(paste0("    - ", final_vars, collapse="\n"), "\n")
