# Data Archive: Vegetation Functional Diversity and Bird Abundance in NSW, Australia

This archive contains supplementary data files for:

> *[Author et al. (in prep)] Vegetation functional diversity explains bird relative abundance across New South Wales, Australia.*

Source code and analysis scripts are available at:
**https://github.com/Josh-Lee1/ebird_rel_abund**

---

## Files in This Archive

| File | Size | Description |
|------|------|-------------|
| `models_500m_inference.zip` | ~805 MB | 213 fitted GAM objects (`.rds` format), one per native terrestrial bird species modelled for NSW |
| `cleaned_benchmark_data.csv` | ~518 MB | NSW Vegetation Information System floristic plot census data: species-level cover scores per census plot, used to compute functional diversity, taxonomic diversity, and community-weighted mean trait variables |
| `nsw_abundance_stack_3km.tif` | ~145 MB | Multi-band GeoTIFF (WGS84, 3 km grid): relative abundance predictions for all modelled species across NSW. Band names correspond to eBird common names. |
| `nsw_abundance_se_stack_3km.tif` | ~136 MB | Multi-band GeoTIFF: standard error of relative abundance predictions (same species order as abundance stack) |

---

## Model Files (`models_500m_inference.zip`)

Each `.rds` file is a fitted `mgcv::bam()` object for one bird species, named by eBird common name in snake_case (e.g. `superb_fairywren.rds`). Models were fit using the 500 m plot buffer inference configuration (`config_500m_inference.yaml`).

**To load a model in R:**
```r
mod <- readRDS("superb_fairywren.rds")
summary(mod)
mgcv::plot.gam(mod, pages = 1)
```

**Model structure:**
- Family: Negative binomial (`nb()`)
- Covariates: effort terms + landscape covariates + 12 vegetation variables (see paper)
- Package: `mgcv` version ≥ 1.9

**Species list:** See `results_500m_inference/convergence_report.csv` in the GitHub repository for the full list of 213 modelled species, their convergence status, and deviance explained.

---

## Floristic Plot Data (`cleaned_benchmark_data.csv`)

NSW Vegetation Information System (VIS) benchmark plot census records. Each row is one species observation in one census plot.

| Column | Description |
|--------|-------------|
| `CensusKey` | Unique plot identifier (join key for `bench_GEDI_struct_metrics_enriched.csv`) |
| `aligned_name` | Aligned plant species name |
| `CoverScore` | Braun-Blanquet cover score for this species in this plot |
| *(additional columns)* | Site metadata (coordinates, date, vegetation community) |

**Note:** This dataset contains floristic records for vegetation plots across NSW. It is derived from the NSW Vegetation Information System. Please contact the lead author for further information on data access terms.

---

## Abundance Raster Stacks (`nsw_abundance_stack_3km.tif`, `nsw_abundance_se_stack_3km.tif`)

- **Projection:** WGS84 (EPSG:4326)
- **Resolution:** ~3 km × 3 km
- **Extent:** NSW + 100 km buffer
- **Bands:** One band per modelled species (213 bands total)
- **Values:** Relative abundance (expected count under standardised effort conditions)

**To read in R:**
```r
library(terra)
abd <- terra::rast("nsw_abundance_stack_3km.tif")
names(abd)        # species names
plot(abd[["superb_fairywren"]])
```

---

## Reproducing the Analysis

To reproduce the models and rasters from raw inputs:

1. Clone the code repository: `git clone https://github.com/Josh-Lee1/ebird_rel_abund`
2. Download eBird EBD data for NSW/ACT/VIC/QLD/SA from https://ebird.org/data/download
3. Obtain `cleaned_benchmark_data.csv` (this archive) and place in `data/`
4. Follow the step-by-step workflow in the repository README

---

## Citation

If you use these data or models, please cite:

> *[Author et al. (in prep)] Vegetation functional diversity explains bird relative abundance across New South Wales, Australia.*

And the underlying modelling framework:

> Strimas-Mackey, M., et al. (2023). *Best Practices for Using eBird Data*. Version 2.0. Cornell Lab of Ornithology. https://doi.org/10.5281/zenodo.3620739

---

## Contact

**Lead author:** [Name]  
**Institution:** Western Sydney University  
**Email:** [email]
