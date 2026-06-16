# Floristic plot ↔ eBird checklist spatial intersection utilities.
#
# Loads vegetation/functional-diversity plot data, buffers each plot to
# a circular polygon, spatially joins eBird checklists that fall within
# those buffers, and provides simple diagnostics on coverage.

# ---------------------------------------------------------------------------
# Load plot data
# ---------------------------------------------------------------------------
# Read the primary plot CSV (bench_GEDI_struct_metrics.csv) and optionally
# join CWM trait columns from a second CSV (cleaned_veg_w_fd.csv).
#
# plot_csv  : path to the primary plot CSV (must contain site, lat, lon)
# cwm_csv   : optional path to CWM CSV (joined by site = CensusKey)
# plot_vars : optional character vector of variable names to keep
#
# Returns an sf POINT object (CRS 4326) with one row per plot.
#' @export
load_plot_data <- function(plot_csv, cwm_csv = NULL, plot_vars = NULL) {
  if (!file.exists(plot_csv)) stop("Plot CSV not found: ", plot_csv)

  plots <- read.csv(plot_csv, stringsAsFactors = FALSE)

  required <- c("site", "lat", "lon")
  missing  <- setdiff(required, names(plots))
  if (length(missing) > 0L) {
    stop("Plot CSV missing required columns: ", paste(missing, collapse = ", "))
  }


  # Join CWM traits if a second file is provided
  if (!is.null(cwm_csv)) {
    if (!file.exists(cwm_csv)) stop("CWM CSV not found: ", cwm_csv)
    cwm <- read.csv(cwm_csv, stringsAsFactors = FALSE)
    if (!"CensusKey" %in% names(cwm)) {
      stop("CWM CSV must contain a 'CensusKey' column for joining.")
    }
    # Keep only CensusKey and requested plot_vars from CWM that are NOT already in plots
    if (!is.null(plot_vars)) {
      novel_vars <- setdiff(intersect(names(cwm), plot_vars), names(plots))
      cwm_keep <- c("CensusKey", novel_vars)
      cwm <- cwm[, cwm_keep, drop = FALSE]
    }
    
    # Join on site = CensusKey; keep all plot rows (left join)
    plots <- merge(plots, cwm, by.x = "site", by.y = "CensusKey", all.x = TRUE)
  }

  # Subset to requested variables (always keep site, lat, lon)
  if (!is.null(plot_vars)) {
    keep <- unique(c("site", "lat", "lon", plot_vars))
    avail <- intersect(keep, names(plots))
    miss  <- setdiff(keep, names(plots))
    if (length(miss) > 0L) {
      warning("Requested plot variables not found (ignored): ",
              paste(miss, collapse = ", "))
    }
    plots <- plots[, avail, drop = FALSE]
  }

  plots_sf <- sf::st_as_sf(plots, coords = c("lon", "lat"),
                            crs = 4326, remove = FALSE)

  message(sprintf("Loaded %d floristic plots from %s",
                  nrow(plots_sf), basename(plot_csv)))
  plots_sf
}

# ---------------------------------------------------------------------------
# Buffer plots
# ---------------------------------------------------------------------------
# Create circular buffer polygons around each plot point.
#
# plots_sf  : sf POINT object returned by load_plot_data()
# buffer_m  : buffer radius in metres
# proj_crs  : projected CRS for buffering (default EPSG:3577 GDA94/Albers)
#
# Returns an sf POLYGON object (CRS 4326) with the same attributes as the
# input plus plot_lon / plot_lat columns recording original point coordinates.
#' @export
buffer_plots <- function(plots_sf, buffer_m, proj_crs = 3577) {
  if (!inherits(plots_sf, "sf")) stop("`plots_sf` must be an sf object.")

  # Preserve the original point coordinates for later distance calculations
  plots_sf$plot_lon <- plots_sf$lon
  plots_sf$plot_lat <- plots_sf$lat

  plots_proj <- sf::st_transform(plots_sf, proj_crs)
  buffered   <- sf::st_buffer(plots_proj, dist = buffer_m)
  buffered   <- sf::st_transform(buffered, 4326)

  message(sprintf("Buffered %d plots with radius %d m", nrow(buffered), buffer_m))
  buffered
}

# ---------------------------------------------------------------------------
# Intersect checklists with buffered plots
# ---------------------------------------------------------------------------
# Spatially join eBird checklists to buffered plot polygons.
# When a checklist falls inside multiple overlapping buffers the nearest
# plot (by Euclidean distance to the original plot point) is kept.
#
# checklists_df  : data.frame with latitude and longitude columns
#                  (e.g. from sampling_master.rds)
# plots_buffered : sf POLYGON object returned by buffer_plots()
# plot_vars      : character vector of plot variable names to attach
# method         : "nearest" (default) — keep closest plot when ambiguous
#
# Returns a plain data.frame with all original checklist columns plus
# plot_id (= site) and the requested plot_vars columns.
#' @export
intersect_checklists_with_plots <- function(checklists_df, plots_buffered,
                                            plot_vars, method = "nearest") {
  if (!all(c("latitude", "longitude") %in% names(checklists_df))) {
    stop("`checklists_df` must contain 'latitude' and 'longitude' columns.")
  }

  checklists_sf <- sf::st_as_sf(
    checklists_df,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

  # Spatial join: keep checklists that fall within a buffered plot
  joined <- sf::st_join(checklists_sf, plots_buffered, join = sf::st_within)

  # Drop checklists that matched no plot (inner join)
  joined <- joined[!is.na(joined$site), ]

  if (nrow(joined) == 0L) {
    message("No checklists intersected any plot buffer.")
    cols <- c(names(checklists_df), "plot_id", plot_vars)
    return(data.frame(matrix(nrow = 0, ncol = length(cols),
                             dimnames = list(NULL, cols))))
  }

  # Resolve duplicates: when a checklist lands in multiple buffers, keep

  # the one closest to the original plot point.
  joined_df <- sf::st_drop_geometry(joined)

  if (method == "nearest") {
    # Identify checklist ID column (from sampling master)
    id_col <- if ("checklist_id" %in% names(joined_df)) {
      "checklist_id"
    } else {
      stop("Cannot find 'checklist_id' column in checklists_df.")
    }

    # Compute distance from each checklist to its matched plot point
    joined_df$.dist <- sqrt(
      (joined_df$longitude - joined_df$plot_lon)^2 +
      (joined_df$latitude  - joined_df$plot_lat)^2
    )

    # Keep closest plot per checklist
    joined_df <- joined_df[order(joined_df[[id_col]], joined_df$.dist), ]
    joined_df <- joined_df[!duplicated(joined_df[[id_col]]), ]
    joined_df$.dist <- NULL
  }

  # Rename site -> plot_id; select final columns
  joined_df$plot_id <- joined_df$site

  keep_cols <- c(names(checklists_df), "plot_id",
                 intersect(plot_vars, names(joined_df)))
  out <- joined_df[, keep_cols, drop = FALSE]
  rownames(out) <- NULL

  message(sprintf("Matched %d checklists to %d unique plots",
                  nrow(out), length(unique(out$plot_id))))
  out
}

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
# Compute and print summary statistics for the plot–checklist intersection.
#
# intersection_df : data.frame returned by intersect_checklists_with_plots()
# n_total_plots   : total number of plots in the dataset (for coverage calc)
#
# Returns an invisible list of summary values.
#' @export
plot_intersection_diagnostics <- function(intersection_df, n_total_plots) {
  n_matched_plots     <- length(unique(intersection_df$plot_id))
  n_matched_checklists <- nrow(intersection_df)
  per_plot            <- table(intersection_df$plot_id)
  n_empty_plots       <- n_total_plots - n_matched_plots

  message("=== Plot-checklist intersection diagnostics ===")
  message(sprintf("  Plots with >= 1 checklist : %d / %d (%.1f%%)",
                  n_matched_plots, n_total_plots,
                  100 * n_matched_plots / max(n_total_plots, 1)))
  message(sprintf("  Total matched checklists  : %d", n_matched_checklists))
  message(sprintf("  Checklists per plot — min: %d  median: %.1f  mean: %.1f  max: %d",
                  min(per_plot), stats::median(per_plot),
                  mean(per_plot), max(per_plot)))
  message(sprintf("  Plots with 0 checklists   : %d", n_empty_plots))

  invisible(list(
    n_matched_plots      = n_matched_plots,
    n_matched_checklists = n_matched_checklists,
    checklists_per_plot  = as.integer(per_plot),
    min_per_plot         = as.integer(min(per_plot)),
    median_per_plot      = stats::median(per_plot),
    mean_per_plot        = mean(per_plot),
    max_per_plot         = as.integer(max(per_plot)),
    n_empty_plots        = n_empty_plots
  ))
}
