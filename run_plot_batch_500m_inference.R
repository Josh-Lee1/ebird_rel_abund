#!/usr/bin/env Rscript
# =============================================================================
# Inference Multi-species batch runner
# =============================================================================
# Filters to terrestrial species only, uses 500m config, and runs in parallel.
#
# Usage:
#   Rscript run_plot_batch_500m_inference.R
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(mgcv)
  library(dplyr)
  library(yaml)
  library(geodata)
  library(parallel)
})

Sys.setenv(EBIRD_CONFIG = "config_500m_inference.yaml")
Sys.setenv(EBIRD_RESULTS = "results_500m_inference")

CONFIG_FILE <- Sys.getenv("EBIRD_CONFIG")
cfg <- yaml::read_yaml(CONFIG_FILE)
plot_cfg <- cfg$plot_data

RESULTS_DIR <- Sys.getenv("EBIRD_RESULTS")
CONV_FILE   <- file.path(RESULTS_DIR, "convergence_report.csv")
IMP_FILE    <- file.path(RESULTS_DIR, "aggregated_importance.csv")
SMOOTH_DIR  <- file.path(RESULTS_DIR, "smooths")
MODEL_DIR   <- file.path(RESULTS_DIR, "models")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SMOOTH_DIR,  showWarnings = FALSE)
dir.create(MODEL_DIR,   showWarnings = FALSE)

# ── Species list (filtered for terrestrial only) ─────────────────────────────
species_df <- read.csv(cfg$species_list_file, stringsAsFactors = FALSE)
threshold  <- cfg$reporting_rate_threshold
species_df <- species_df[species_df$reporting_rate >= threshold, ]

# Filter to only native terrestrial species to remove waterbirds
terrestrial <- read.csv("results/filtered_native_terrestrial_species.csv", stringsAsFactors = FALSE)
species_list <- intersect(species_df$common_name, terrestrial$species)

message(sprintf("Species to process: %d (terrestrial only)", length(species_list)))

# ── Initialise output CSVs ───────────────────────────────────────────────────
conv_cols <- c("species", "converged", "n_checklists", "n_positive",
               "pct_positive", "deviance_explained", "n_smooths",
               "max_edf", "status", "error_message", "timestamp")
imp_cols  <- c("species", "term", "edf", "f_stat", "p_value",
               "category", "is_vegetation", "rank", "timestamp")

if (!file.exists(CONV_FILE)) {
  write.csv(data.frame(matrix(nrow = 0, ncol = length(conv_cols), dimnames = list(NULL, conv_cols))), CONV_FILE, row.names = FALSE)
}
if (!file.exists(IMP_FILE)) {
  write.csv(data.frame(matrix(nrow = 0, ncol = length(imp_cols), dimnames = list(NULL, imp_cols))), IMP_FILE, row.names = FALSE)
}

safe_name <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# Resume: skip already done
done <- read.csv(CONV_FILE, stringsAsFactors = FALSE)
species_list <- setdiff(species_list, done$species)
if (length(species_list) == 0) {
  message("All species completed.")
  quit(save = "no")
}

# ── Shared Data Load (One Time) ──────────────────────────────────────────────
aus <- geodata::gadm(country = cfg$study_polygon$country, level = 1, path = cfg$covariate_cache)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
polygon <- sf::st_transform(sf::st_buffer(sf::st_transform(region_sf, cfg$study_polygon$proj_processing), cfg$study_polygon$buffer_metres), cfg$study_polygon$proj_output)

plots_sf <- load_plot_data(plot_cfg$file, plot_cfg$cwm_file, plot_cfg$variables)
buffered_plots <- buffer_plots(plots_sf, buffer_m = plot_cfg$buffer_m)

ebd_path <- file.path(cfg$raw_data, cfg$ebd_files[[1]])
smp_path <- file.path(cfg$raw_data, cfg$sampling_files[[1]])

cov <- prepare_covariates(polygon, cache_dir = cfg$covariate_cache, nightlights_path = "data/covariates/nightlights/World_Atlas_2015.tif")
prepare_sampling_master(smp_path, polygon, cov, cfg$zerofill_cache)

# ── Parallel Processing Function ─────────────────────────────────────────────
process_species <- function(sp) {
  t_sp <- proc.time()[["elapsed"]]
  sp_safe <- safe_name(sp)
  
  res <- tryCatch({
    df <- load_ebird_with_plots(
      polygon = polygon, ebird_zip = ebd_path, sampling_txt = smp_path,
      species = sp, cache_dir = cfg$zerofill_cache,
      plots_buffered = buffered_plots, plot_vars = plot_cfg$variables
    )
    df_samp <- subsample_hex(df)
    n_total <- nrow(df_samp)
    n_positive <- sum(df_samp$observation_count > 0)
    
    if (n_positive < 50L) stop(sprintf("Insufficient detections: %d", n_positive))
    
    # Fit GAM (using our new model.R inference settings)
    mod <- fit_gam(df_samp, extra_smooth_cols = plot_cfg$variables)
    
    gam_sum <- suppressWarnings(summary(mod))
    dev_expl <- gam_sum$dev.expl
    max_edf <- max(gam_sum$s.table[, "edf"], na.rm = TRUE)
    
    saveRDS(mod, file.path(MODEL_DIR, paste0(sp_safe, ".rds")))
    
    # Save smooth plots
    n_smooths <- length(mod$smooth)
    ncol_plot <- 4
    nrow_plot <- ceiling(n_smooths / ncol_plot)
    png(file.path(SMOOTH_DIR, paste0(sp_safe, ".png")), width = 1200, height = max(400, nrow_plot * 300), res = 150)
    par(mfrow = c(nrow_plot, ncol_plot), mar = c(4, 4, 2, 1))
    for (j in seq_len(n_smooths)) plot(mod, select = j, scheme = 1, shade = TRUE, shade.col = "lightblue", all.terms = FALSE)
    dev.off()
    
    imp <- extract_variable_importance(mod, veg_vars = plot_cfg$variables)
    imp$rank <- seq_len(nrow(imp))
    
    conv_row <- data.frame(
      species = sp, converged = mod$converged, n_checklists = n_total, n_positive = n_positive,
      pct_positive = round(100 * n_positive / n_total, 1), deviance_explained = round(dev_expl, 4),
      n_smooths = n_smooths, max_edf = round(max_edf, 2), status = "ok", error_message = NA_character_,
      timestamp = format(Sys.time()), stringsAsFactors = FALSE
    )
    
    imp_out <- data.frame(
      species = sp, term = imp$term, edf = imp$edf, f_stat = imp$f_stat, p_value = imp$p_value,
      category = imp$category, is_vegetation = imp$is_vegetation, rank = imp$rank,
      timestamp = format(Sys.time()), stringsAsFactors = FALSE
    )
    
    list(success = TRUE, sp = sp, conv = conv_row, imp = imp_out, elapsed = proc.time()[["elapsed"]] - t_sp)
  }, error = function(e) {
    conv_row <- data.frame(
      species = sp, converged = FALSE, n_checklists = NA_integer_, n_positive = NA_integer_,
      pct_positive = NA_real_, deviance_explained = NA_real_, n_smooths = NA_integer_,
      max_edf = NA_real_, status = "failed", error_message = conditionMessage(e),
      timestamp = format(Sys.time()), stringsAsFactors = FALSE
    )
    list(success = FALSE, sp = sp, conv = conv_row, error = conditionMessage(e))
  })
  return(res)
}

# ── Run Parallel ─────────────────────────────────────────────────────────────
message(sprintf("Starting parallel processing across %d cores...", cfg$n_cores))
cl <- makeCluster(cfg$n_cores)

# Export necessary variables and libraries to the cluster
clusterEvalQ(cl, {
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(mgcv)
  library(dplyr)
})
clusterExport(cl, c("polygon", "ebd_path", "smp_path", "cfg", "plot_cfg", "buffered_plots", "MODEL_DIR", "SMOOTH_DIR", "safe_name", "process_species"))

results <- parLapply(cl, species_list, process_species)
stopCluster(cl)

# ── Aggregate and Save Results ───────────────────────────────────────────────
message("Writing results to disk...")
n_ok <- 0; n_fail <- 0
for (res in results) {
  write.table(res$conv, CONV_FILE, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE, quote = TRUE)
  if (res$success) {
    write.table(res$imp, IMP_FILE, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE, quote = TRUE)
    n_ok <- n_ok + 1
    message(sprintf("  ✓ %s", res$sp))
  } else {
    n_fail <- n_fail + 1
    message(sprintf("  ✗ %s — %s", res$sp, res$error))
  }
}

message(sprintf("\nDone. %d ok | %d failed", n_ok, n_fail))
