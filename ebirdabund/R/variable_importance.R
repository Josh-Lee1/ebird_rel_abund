# Extract and classify variable importance from a fitted mgcv GAM.
#
# model    : a fitted mgcv::gam (or bam) object.
# veg_vars : optional character vector of vegetation variable names
#            (e.g. CWM trait columns from floristic plots). Used to tag
#            terms with category = "vegetation".
#
# Returns a data.frame sorted by F-statistic descending, with columns:
#   term, edf, f_stat, p_value, category
#   (plus is_vegetation logical if veg_vars is supplied)
#' @export
extract_variable_importance <- function(model, veg_vars = NULL) {
  gam_sum  <- suppressWarnings(summary(model))
  gam_anov <- suppressWarnings(anova(model))

  # ── Smooth terms ──────────────────────────────────────────────────────────
  rows <- data.frame(
    term    = character(0),
    edf     = numeric(0),
    f_stat  = numeric(0),
    p_value = numeric(0),
    stringsAsFactors = FALSE
  )

  if (!is.null(gam_sum$s.table) && nrow(gam_sum$s.table) > 0) {
    stbl <- as.data.frame(gam_sum$s.table)
    rows <- data.frame(
      term    = rownames(stbl),
      edf     = stbl$edf,
      f_stat  = stbl[, "F"],
      p_value = stbl[, "p-value"],
      stringsAsFactors = FALSE
    )
  }

  # ── Parametric terms (e.g. protocol_type) ─────────────────────────────────
  if (!is.null(gam_anov$pTerms.table) && nrow(gam_anov$pTerms.table) > 0) {
    ptbl <- as.data.frame(gam_anov$pTerms.table)
    prows <- data.frame(
      term    = rownames(ptbl),
      edf     = ptbl$df,
      f_stat  = ptbl$F,
      p_value = ptbl[, "p-value"],
      stringsAsFactors = FALSE
    )
    rows <- rbind(rows, prows)
  }

  # ── Classify each term ────────────────────────────────────────────────────
  # Strip the s(...) wrapper to get the bare variable name for matching.
  # mgcv terms may include parameters like s(FEve,k=3) — strip those too.
  bare_name <- gsub("^s\\(([^,)]+).*\\)$", "\\1", rows$term)
  # Also strip any transform wrappers like log(), log1p()
  bare_name <- gsub("^log1p\\((.+)\\)$", "\\1", bare_name)
  bare_name <- gsub("^log\\((.+)\\)$", "\\1", bare_name)

  effort_pattern <- paste0(
    "^(day_of_year|time_obs|duration|effort_distance|",
    "number_observers|protocol|observer_expertise)"
  )
  landscape_pattern <- paste0(
    "^(lc_|elevation|precip_|temp_|pop_|water_|clay|",
    "tree_height|nightlights|palsar)"
  )

  classify_term <- function(bare) {
    if (grepl(effort_pattern, bare))    return("effort")
    if (grepl(landscape_pattern, bare)) return("landscape")
    if (!is.null(veg_vars) && bare %in% unlist(veg_vars)) return("vegetation")
    "other"
  }

  rows$category <- vapply(bare_name, classify_term, character(1),
                           USE.NAMES = FALSE)

  # ── Optional vegetation flag ──────────────────────────────────────────────
  if (!is.null(veg_vars)) {
    rows$is_vegetation <- bare_name %in% unlist(veg_vars)
  }

  # Sort by F-statistic descending

  rows <- rows[order(-rows$f_stat), ]
  rownames(rows) <- NULL
  rows
}

# Save a variable-importance data.frame to CSV with species and timestamp.
#
# importance_df : data.frame from extract_variable_importance().
# species_name  : character, species common name.
# output_path   : file path for the output CSV.
#' @export
save_importance_csv <- function(importance_df, species_name, output_path) {
  importance_df$species   <- species_name
  importance_df$timestamp <- Sys.time()

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(importance_df, output_path, row.names = FALSE)
  message(sprintf(
    "Variable importance saved: %s (%d terms for '%s').",
    output_path, nrow(importance_df), species_name
  ))
  invisible(importance_df)
}
