# Vegetation Functional Diversity and Bird Abundance in New South Wales

**How do vegetation functional diversity, functional identity, structural complexity and taxonomic diversity explain bird relative abundance across NSW, Australia?**

This repository contains the data, R code, and results for a study that models the relative abundance of 213 native terrestrial bird species across New South Wales (NSW), Australia, and uses drop-term cross-validation to partition the unique explanatory power of four vegetation variable groups: functional diversity (FD), functional identity (CWM traits), vegetation structure, and plant taxonomic diversity.

---

## Study System

Bird abundance models are fit to citizen-science eBird observations across NSW using Generalised Additive Models (GAMs) following the [eBird Best Practices](https://doi.org/10.5281/zenodo.3620739) framework. Each model integrates two data streams:

1. **eBird checklist observations** — zero-filled detection histories from effort-filtered complete checklists in NSW, ACT, VIC, QLD, and SA (training buffer).
2. **Floristic plot covariates** — vegetation diversity, structure, and functional trait metrics derived from the NSW Vegetation Information System benchmark census plots, enriched with GEDI LiDAR structural metrics.

For each species, 12 vegetation variables are included as GAM smooth terms alongside landscape covariates (climate, land cover, elevation, nightlights). Drop-term 5-fold cross-validation then quantifies the unique out-of-sample deviance explained by each vegetation variable and variable group.

---

## Repository Structure

```
bird_FD_ebirdabund/
│
├── ebirdabund/                   # R package: core modelling functions
│   ├── R/                        # batch.R, covariates.R, model.R, load_ebird.R, etc.
│   ├── DESCRIPTION
│   └── NAMESPACE
│
├── data/                         # Input data (see Data Requirements below)
│   ├── bench_GEDI_struct_metrics.csv          # Raw GEDI structural metrics per plot
│   ├── bench_GEDI_struct_metrics_enriched.csv # + taxonomic diversity (computed by pipeline)
│   ├── cleaned_veg_w_fd.csv                   # CWM traits and FD indices per plot
│   ├── all_bird_traits.csv                    # AVONET / EltonTraits bird trait data
│   ├── Australian_Bird_Data_Version_1.xlsx    # Australian bird ecological traits
│   ├── nsw_ebird_taxonomy.csv                 # eBird common ↔ scientific name lookup
│   ├── botw_name_aliases.csv                  # eBird ↔ BirdLife name alias map
│   └── covariates/                            # Landscape raster layers (auto-cached)
│
├── analysis/                     # Pipeline analysis scripts (see Workflow below)
│   ├── compute_diversity_and_vif.R
│   ├── run_drop_term_inference.R
│   ├── run_individual_veg_drop_terms.R
│   ├── run_landscape_drop.R
│   ├── plot_macro_category_importance.R
│   ├── plot_specialist_ternary_inference.R
│   ├── plot_unique_individual_veg.R
│   ├── plot_trait_analysis.R
│   └── check_concurvity.R
│
├── results_500m_inference/       # All model outputs (500 m plot buffer, inference run)
│   ├── convergence_report.csv         # Per-species convergence and deviance metrics
│   ├── aggregated_importance.csv      # Per-species GAM term F-statistics
│   ├── drop_term_inference_results.csv# Drop-term CV: 4 vegetation groups × species
│   ├── individual_veg_drop_results.csv# Drop-term CV: 12 individual vars × species
│   ├── landscape_drop_results.csv     # Drop-term CV: landscape context × species
│   ├── collinearity_vif.csv           # VIF scores for vegetation predictors
│   ├── collinearity_correlation_plot.png
│   ├── smooths/                       # Per-species GAM smooth diagnostic plots
│   └── figures/                       # Publication figures (see below)
│
├── results/
│   └── filtered_native_terrestrial_species.csv  # Final species list (native terrestrial)
│
├── nsw_species_list.R            # Generates nsw_species_list.csv from eBird EBD
├── nsw_species_list.csv          # Species list with reporting rates
├── run_plot_batch_500m_inference.R  # Main GAM fitting batch runner
├── config_500m_inference.yaml    # Configuration for the inference run
└── botw_name_aliases.csv         # eBird–BirdLife taxonomy alias map
```

---

## Publication Figures

The four final publication figures are produced by the scripts listed, reading from `results_500m_inference/`:

| Figure | Script | Output file |
|--------|--------|-------------|
| **Fig. 1** — Macro-category variable importance stacked bar (F-statistics) | `analysis/plot_macro_category_importance.R` | `results_500m_inference/figures/figure_1_macro_categories_split.png` |
| **Fig. 2** — Unique deviance by individual vegetation variable | `analysis/plot_unique_individual_veg.R` | `results_500m_inference/figures/figure_unique_individual_veg_raw.png` |
| **Fig. 3** — Ternary plot: structure vs functional traits vs taxonomic diversity by trophic niche | `analysis/plot_specialist_ternary_inference.R` | `results_500m_inference/figures/dropterm_deepdive/specialist_ternary_inference.png` |
| **Fig. 4** — Functional importance vs species traits (diet, habitat, body mass, HWI) | `analysis/plot_trait_analysis.R` | `results_500m_inference/figures/figure_4_trait_analysis.png` |

---

## Data Requirements

### ✅ Included in this repository
| File | Description | Size |
|------|-------------|------|
| `data/bench_GEDI_struct_metrics.csv` | GEDI LiDAR structural metrics per floristic plot | ~5 MB |
| `data/bench_GEDI_struct_metrics_enriched.csv` | Above + taxonomic diversity metrics (pipeline output) | ~5.2 MB |
| `data/cleaned_veg_w_fd.csv` | Community-weighted mean (CWM) traits and FD indices per plot | ~13 MB |
| `data/all_bird_traits.csv` | AVONET / EltonTraits bird functional traits | ~620 KB |
| `data/nsw_ebird_taxonomy.csv` | eBird taxonomy lookup table | ~23 KB |
| `nsw_species_list.csv` | NSW species list with reporting rates | ~36 KB |
| `results/filtered_native_terrestrial_species.csv` | Final modelled species list | ~24 KB |
| `results_500m_inference/*.csv` | All model result summary tables | <1 MB each |
| `results_500m_inference/figures/` | All publication figures | ~2 MB total |

### ⬇️ Must be downloaded separately

#### eBird observation data
eBird EBD files must be requested from the [eBird Data Access portal](https://ebird.org/data/download) (free, requires registration). The pipeline uses the following files (placed in `data/`):
- `ebd_AU-NSW_unv_smp_relFeb-2026/ebd_AU-NSW_unv_smp_relFeb-2026.txt` + `_sampling.txt`
- `ebd_AU-ACT_unv_smp_relMar-2026/` (ACT buffer region)
- `ebd_AU-VIC_unv_smp_relMar-2026/` (VIC buffer region)
- `ebd_AU-QLD_unv_smp_relMar-2026/` (QLD buffer region)
- `ebd_AU-SA_unv_smp_relMar-2026/` (SA buffer region)

These files are not redistributable under the eBird data use policy.

### 📧 Available on request from the lead author

| File | Reason not included |
|------|---------------------|
| `data/cleaned_benchmark_data.csv` (~520 MB) | Floristic plot raw census data — too large for GitHub; contact lead author |
| `data/botw_species/BOTW_2025.gpkg` | BirdLife range polygons — licence-restricted; contact lead author |
| `results_500m_inference/models/*.rds` | Fitted GAM objects, 213 files, ~850 MB total — see below |

### 🔄 Fitted model files
The 213 fitted GAM `.rds` files (~4 MB each, ~850 MB total) cannot be hosted on GitHub due to file-size limits. They can be:
- **Reproduced** by running the pipeline (Step 3 below) — requires ~24–48 h of compute on 6 cores.
- **Archived** on a data repository (e.g. Zenodo, OSF) — contact the lead author if you need pre-fitted models without re-running.

> **Note:** RDS files are already gzip-compressed by R internally; additional ZIP compression yields negligible size reduction.

---

## Analysis Workflow

Run the following steps in order. All scripts assume the working directory is the repository root.

### Step 0 — Install the `ebirdabund` package

```r
# From local source (recommended)
devtools::load_all("ebirdabund")

# Or install from GitHub if published:
# pak::pak("username/ebirdabund")
```

### Step 1 — Generate the NSW species list

Reads the eBird EBD sampling file and computes reporting rates across effort-filtered complete checklists. Outputs `nsw_species_list.csv`.

```r
source("nsw_species_list.R")
```

### Step 2 — Compute vegetation diversity metrics and check collinearity

Calculates Pielou's evenness, Gini-Simpson diversity, and species richness per floristic plot; joins to GEDI structural metrics; runs VIF analysis to select the final 12-variable set. Outputs `data/bench_GEDI_struct_metrics_enriched.csv` and `results_500m_inference/collinearity_vif.csv`.

```r
source("analysis/compute_diversity_and_vif.R")
```

### Step 3 — Fit GAMs for all species (inference run)

Fits negative-binomial GAMs for all native terrestrial species using the 500 m plot buffer configuration. Runs in parallel (6 cores by default; set in `config_500m_inference.yaml`). Saves fitted models to `results_500m_inference/models/`, GAM smooth plots to `results_500m_inference/smooths/`, and aggregated F-statistics to `results_500m_inference/aggregated_importance.csv`.

```r
source("run_plot_batch_500m_inference.R")
```

> Requires eBird EBD data and floristic plot data. Allow 24–48 h on a 6-core machine.

### Step 4 — Drop-term cross-validation (vegetation groups)

5-fold cross-validation comparing the full model against models with each vegetation group dropped. Quantifies the unique out-of-sample deviance explained by functional diversity, functional identity, structure, and taxonomic diversity.

```r
source("analysis/run_drop_term_inference.R")
```

### Step 5 — Drop-term cross-validation (individual variables)

As above but dropping each of the 12 vegetation variables individually (13 models per species). Outputs `results_500m_inference/individual_veg_drop_results.csv`.

```r
source("analysis/run_individual_veg_drop_terms.R")
```

### Step 6 — Drop-term cross-validation (landscape context)

Quantifies the unique contribution of landscape covariates (climate, land cover, elevation) by comparison to a vegetation-only model. Merges result into `results_500m_inference/drop_term_inference_results.csv`.

```r
source("analysis/run_landscape_drop.R")
```

### Step 7 — Generate publication figures

Each script can be run independently once Steps 3–6 are complete:

```r
source("analysis/plot_macro_category_importance.R")      # Fig. 1
source("analysis/plot_unique_individual_veg.R")           # Fig. 2
source("analysis/plot_specialist_ternary_inference.R")    # Fig. 3
source("analysis/plot_trait_analysis.R")                  # Fig. 4
```

---

## Model Description

Each species is modelled with a **negative-binomial GAM** (via `mgcv::bam`) with:

- **Effort covariates** (as smooths): duration, distance, time of day, day of year, observer expertise, number of observers
- **Landscape covariates** (as smooths): elevation, land cover classes, annual temperature, annual precipitation, human population density, water occurrence, nightlights
- **Vegetation covariates** (as smooths, 500 m buffer around each plot):
  - *Functional diversity*: FRic (Functional Richness), FEve (Functional Evenness), FDiv (Functional Divergence)
  - *Functional identity*: Height CWM, Leaf Area CWM, Wood Density CWM
  - *Vegetation structure*: Percentage cover, proportion canopy, mean foliage height diversity (FHD)
  - *Taxonomic diversity*: Species richness, Pielou's evenness, Gini-Simpson diversity

Spatial subsampling uses `subsample_hex()` before model fitting to reduce spatial autocorrelation. Final predictions use standardised effort (median covariate values).

Variable importance is assessed via:
1. **F-statistics** from the GAM summary (`s.table`) — proportional contribution of each term.
2. **Drop-term 5-fold cross-validation** — unique out-of-sample deviance explained by each variable/group.

---

## Key Output Files

| File | Description |
|------|-------------|
| `results_500m_inference/convergence_report.csv` | Species-level model convergence, deviance explained, N checklists, N detections |
| `results_500m_inference/aggregated_importance.csv` | Per-term F-statistics and EDF for all species |
| `results_500m_inference/drop_term_inference_results.csv` | Unique deviance by 4 vegetation groups + landscape |
| `results_500m_inference/individual_veg_drop_results.csv` | Unique deviance for each of the 12 vegetation variables |
| `results_500m_inference/landscape_drop_results.csv` | Unique deviance attributable to landscape context |
| `results_500m_inference/collinearity_vif.csv` | VIF scores used to select the final variable set |

---

## R Dependencies

```r
# Core pipeline
install.packages(c(
  "mgcv",      # GAM fitting
  "sf",        # Spatial operations
  "terra",     # Raster handling
  "data.table",# Fast CSV reading
  "yaml",      # Configuration files
  "geodata",   # GADM boundary download
  "parallel",  # Parallel processing
  "devtools"   # Loading the ebirdabund package
))

# Figures
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "ggtern",    # Ternary plots (Fig. 3)
  "broom",
  "plotly",
  "htmlwidgets"
))

# Optional (for eBird range masking)
install.packages("ebirdst")
```

Set your eBird Status & Trends API key before using `ebirdst`:

```r
Sys.setenv(EBIRDST_KEY = "your_key_here")
```

---

## Citation

If you use this code or data, please cite:

> *[Author et al. (in prep)] Vegetation functional diversity and bird abundance across New South Wales, Australia.*

The `ebirdabund` modelling framework follows:

> Strimas-Mackey, M., Hochachka, W.M., Ruiz-Gutierrez, V., Robinson, O.J., Miller, E.T., Auer, T., Kelling, S., Fink, D., Johnston, A. (2023). *Best Practices for Using eBird Data*. Version 2.0. Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.5281/zenodo.3620739

---

## Contact

For access to restricted data files (`cleaned_benchmark_data.csv`, `botw_species/BOTW_2025.gpkg`, pre-fitted model `.rds` files), please contact the lead author.
