#!/usr/bin/env Rscript
Sys.setenv(EBIRD_CONFIG  = "config_500m_inference.yaml")
Sys.setenv(EBIRD_RESULTS = "results_500m_inference")

suppressPackageStartupMessages({
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(mgcv)
  library(dplyr)
  library(yaml)
  library(geodata)
  library(parallel)
})

CONFIG_FILE <- Sys.getenv("EBIRD_CONFIG")
cfg <- yaml::read_yaml(CONFIG_FILE)
plot_cfg <- cfg$plot_data

# Load filtered species list and convergence report
filtered_spp <- read.csv("results/filtered_native_terrestrial_species.csv", stringsAsFactors = FALSE)$species
conv <- read.csv("results_500m_inference/convergence_report.csv", stringsAsFactors = FALSE)
run_species <- filtered_spp[filtered_spp %in% conv$species[conv$converged == TRUE]]

OUT_FILE <- "results_500m_inference/landscape_drop_results.csv"
if (file.exists(OUT_FILE)) {
  done_df <- read.csv(OUT_FILE, stringsAsFactors = FALSE)
  run_species <- setdiff(run_species, done_df$species)
  cat(sprintf("Resuming: %d species remaining (%d already done)\n", length(run_species), nrow(done_df)))
} else {
  out_cols <- c("species", "dev_no_land", "unique_landscape")
  write.table(data.frame(matrix(ncol = length(out_cols), nrow = 0, dimnames = list(NULL, out_cols))), 
              OUT_FILE, sep = ",", row.names = FALSE, quote = TRUE)
}

if (length(run_species) == 0) {
  cat("All species completed. Merging into drop_term_inference_results.csv...\n")
  old_df <- read.csv("results_500m_inference/drop_term_inference_results.csv", stringsAsFactors = FALSE)
  new_df <- read.csv(OUT_FILE, stringsAsFactors = FALSE)
  if (!"unique_landscape" %in% names(old_df)) {
    merged <- left_join(old_df, new_df, by = "species")
    write.csv(merged, "results_500m_inference/drop_term_inference_results.csv", row.names = FALSE)
    cat("Merged unique_landscape into drop_term_inference_results.csv\n")
  }
  quit(save = "no")
}

cat("Running Landscape Drop-Term CV Analysis for", length(run_species), "species...\n")

aus       <- geodata::gadm(country = cfg$study_polygon$country, level = 1, path = cfg$covariate_cache)
region_sf <- sf::st_as_sf(aus[aus$NAME_1 == cfg$study_polygon$region, ])
polygon   <- sf::st_transform(
  sf::st_buffer(sf::st_transform(region_sf, cfg$study_polygon$proj_processing), cfg$study_polygon$buffer_metres),
  cfg$study_polygon$proj_output
)

plots_sf <- load_plot_data(plot_csv=plot_cfg$file, cwm_csv=plot_cfg$cwm_file, plot_vars=plot_cfg$variables)
buffered_plots <- buffer_plots(plots_sf, buffer_m = plot_cfg$buffer_m)

ebd_path <- file.path(cfg$raw_data, cfg$ebd_files[[1]])
smp_path <- file.path(cfg$raw_data, cfg$sampling_files[[1]])

cov <- prepare_covariates(polygon, cache_dir = cfg$covariate_cache, nightlights_path = "data/covariates/nightlights/World_Atlas_2015.tif")
prepare_sampling_master(smp_path, polygon, cov, cfg$zerofill_cache)

fd_vars    <- c("FRic", "FEve", "FDiv")
cwm_vars   <- c("wooddensity_cwm", "LA_cwm", "height_cwm")
rich_vars  <- c("sp_richness", "pielou_evenness", "gini_simpson")
struct_vars <- c("cover", "mean_fhd", "prop_canopy")

process_species_cv <- function(sp) {
  t_sp <- proc.time()[["elapsed"]]
  
  df <- tryCatch({
    load_ebird_with_plots(
      polygon = polygon, ebird_zip = ebd_path, sampling_txt = smp_path,
      species = sp, cache_dir = cfg$zerofill_cache, plots_buffered = buffered_plots, plot_vars = plot_cfg$variables
    )
  }, error = function(e) return(list(success = FALSE, sp = sp, error = conditionMessage(e))))
  
  df_samp <- subsample_hex(df)
  n_pos <- sum(df_samp$observation_count > 0)
  if (n_pos < 50) return(list(success = FALSE, sp = sp, error = "insufficient detections"))
  
  hab_cols <- grep("^(lc_|elevation|precip_|temp_|pop_|water_|clay|tree_height|nightlights|palsar_hv)", names(df_samp), value = TRUE)
  hab_cols <- setdiff(hab_cols, c("lc_shrubs", "palsar_hv", "clay", "lc_grassland", "lc_cropland", "lc_built", "pop_density", "tree_height"))
  hab_cols <- hab_cols[vapply(hab_cols, function(col) length(unique(stats::na.omit(df_samp[[col]]))) >= 4L, logical(1))]
  
  valid_smooths <- function(vars) {
    vars[vars %in% names(df_samp) & vapply(vars, function(col) length(unique(stats::na.omit(df_samp[[col]]))) >= 4L, logical(1))]
  }
  
  all_veg_valid <- c(valid_smooths(fd_vars), valid_smooths(cwm_vars), valid_smooths(rich_vars), valid_smooths(struct_vars))
  
  form_full      <- build_gam_formula(df_samp, hab_cols, extra_smooth_cols = all_veg_valid)
  form_no_land   <- build_gam_formula(df_samp, hab_cols = character(0), extra_smooth_cols = all_veg_valid)
  
  set.seed(42)
  my_folds <- sample(rep_len(seq_len(5), nrow(df_samp)))
  
  tryCatch({
    cv_full      <- evaluate_model_cv(df_samp, formula = form_full, fold_ids = my_folds, gamma = 1.0)
    cv_no_land   <- evaluate_model_cv(df_samp, formula = form_no_land, fold_ids = my_folds, gamma = 1.0)
    
    res <- data.frame(
      species = sp,
      dev_no_land = cv_no_land$summary["holdout_dev_expl"]
    )
    res$unique_landscape <- cv_full$summary["holdout_dev_expl"] - res$dev_no_land
    
    list(success = TRUE, sp = sp, res = res, elapsed = proc.time()[["elapsed"]] - t_sp)
  }, error = function(e) {
    list(success = FALSE, sp = sp, error = conditionMessage(e))
  })
}

message(sprintf("Starting parallel processing across %d cores...", cfg$n_cores))
cl <- makeCluster(cfg$n_cores)
clusterEvalQ(cl, {
  devtools::load_all("ebirdabund", quiet = TRUE)
  library(sf)
  library(mgcv)
  library(dplyr)
})
clusterExport(cl, c(
  "polygon", "ebd_path", "smp_path", "cfg", "plot_cfg", "buffered_plots", 
  "fd_vars", "cwm_vars", "rich_vars", "struct_vars", "process_species_cv"
))

results <- parLapply(cl, run_species, process_species_cv)
stopCluster(cl)

message("Writing results to disk...")
n_ok <- 0; n_fail <- 0
for (res in results) {
  if (res$success) {
    write.table(res$res, OUT_FILE, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE, quote = TRUE)
    n_ok <- n_ok + 1
    message(sprintf("  ✓ %s", res$sp))
  } else {
    n_fail <- n_fail + 1
    message(sprintf("  ✗ %s — %s", res$sp, res$error))
  }
}

message(sprintf("\nDone. %d ok | %d failed", n_ok, n_fail))

# Merge into main dataset
old_df <- read.csv("results_500m_inference/drop_term_inference_results.csv", stringsAsFactors = FALSE)
new_df <- read.csv(OUT_FILE, stringsAsFactors = FALSE)
if (!"unique_landscape" %in% names(old_df)) {
  merged <- left_join(old_df, new_df, by = "species")
  write.csv(merged, "results_500m_inference/drop_term_inference_results.csv", row.names = FALSE)
  cat("Merged unique_landscape into drop_term_inference_results.csv\n")
}
