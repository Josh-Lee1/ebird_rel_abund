library(mgcv)
library(dplyr)
library(tidyr)

# Load the filtered 218 species list
filtered_spp <- read.csv("results/filtered_native_terrestrial_species.csv", stringsAsFactors = FALSE)
# Note: The user said there are 218 species in the 500m model. 
# We should cross-reference this list to find the actual files.
# The filename format is `common_name` in snake_case.
safe_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9_]+", "_", gsub("'", "", x)))
}

# The true list of 500m converged species is in the convergence report.
# We will just pull the 500m convergence report.
conv_500 <- read.csv("results_500m/convergence_report.csv", stringsAsFactors = FALSE)
conv_500 <- conv_500 %>% 
  filter(status == "ok", converged == TRUE) %>%
  arrange(desc(timestamp)) %>%
  group_by(species) %>% slice(1) %>% ungroup()

# Only keep terrestrial species
conv_500 <- conv_500 %>% filter(species %in% filtered_spp$species)

cat("Found", nrow(conv_500), "converged terrestrial species in 500m model.\n")

# Target a subset of up to 10 for concurvity check to be fast.
test_species <- head(conv_500$species, 10)
cat("Running concurvity on", length(test_species), "species...\n")

concurvity_results <- list()

for (sp in test_species) {
  file_name <- paste0(safe_name(sp), ".rds")
  file_path <- file.path("results_500m", "models", file_name)
  
  if (file.exists(file_path)) {
    tryCatch({
      m <- readRDS(file_path)
      # Calculate worst-case concurvity pairwise
      conc <- concurvity(m, full = FALSE)$worst
      
      # Convert to long format
      conc_df <- as.data.frame(conc)
      conc_df$term1 <- rownames(conc_df)
      rownames(conc_df) <- NULL
      
      long_conc <- conc_df %>%
        pivot_longer(cols = -term1, names_to = "term2", values_to = "concurvity_score") %>%
        filter(term1 != term2) %>%
        filter(!grepl("para", term1) & !grepl("para", term2)) # ignore parametric part (intercept)
      
      long_conc$species <- sp
      concurvity_results[[sp]] <- long_conc
    }, error = function(e) {
      cat("Error for", sp, ":", e$message, "\n")
    })
  } else {
    cat("File not found:", file_path, "\n")
  }
}

all_conc <- bind_rows(concurvity_results)

# Filter to interactions between FD/CWM terms and Structural terms
fd_cwm_patterns <- c("FDiv", "FRic", "FEve", "FDis", "Q", "_cwm")
struct_patterns <- c("tree_height", "prop_canopy", "mean_fhd", "cover")

# Function to check if a term matches any pattern
matches_any <- function(term, patterns) {
  any(sapply(patterns, function(p) grepl(p, term)))
}

if(nrow(all_conc) > 0) {
  all_conc <- all_conc %>%
    rowwise() %>%
    mutate(
      is_fd1 = matches_any(term1, fd_cwm_patterns),
      is_struct1 = matches_any(term1, struct_patterns),
      is_fd2 = matches_any(term2, fd_cwm_patterns),
      is_struct2 = matches_any(term2, struct_patterns)
    ) %>%
    ungroup() %>%
    filter((is_fd1 & is_struct2) | (is_struct1 & is_fd2)) %>%
    select(species, term1, term2, concurvity_score)
  
  cat("\n=== Concurvity between Functional and Structural Terms (Worst-case) ===\n")
  summary_conc <- all_conc %>%
    group_by(term1, term2) %>%
    summarise(
      mean_concurvity = mean(concurvity_score, na.rm=TRUE),
      max_concurvity = max(concurvity_score, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(max_concurvity))
    
  print(head(summary_conc, 15))
  
  write.csv(all_conc, "results_500m/concurvity_summary.csv", row.names = FALSE)
  cat("\nSaved full pairwise scores to results_500m/concurvity_summary.csv\n")
} else {
  cat("\nNo concurvity scores extracted.\n")
}
