# ==============================================================================
# Environmental-variable selection: Spearman correlation + VIF
# ==============================================================================
# Supplementary analysis script for:
# Sánchez-González, L. A., Prieto-Torres, D. A., & Espinosa-Chávez, O. J. (2026).
# Ecological Niche Differences Underlie the Assembly of Endemic Birds in
# Acrobiont Forests of Northern Mesoamerica. Journal of Biogeography, 53, e70141.
# https://doi.org/10.1111/jbi.70141
#
# Purpose
# -------
# This script reproduces the first environmental-variable selection approach
# described in the article: (1) remove strongly correlated predictors using a
# Spearman correlation threshold of |r| >= 0.8, and (2) remove remaining
# multicollinear predictors using a Variance Inflation Factor (VIF) threshold
# of 10.
#
# The published study evaluated a second variable-selection approach based on a
# PCA-derived set of six variables. That second approach is not reconstructed
# here because the present script corresponds specifically to the Spearman + VIF
# workflow. If the original PCA-selection script is archived, it should be
# provided separately rather than inferred from the publication.
#
# Reproducibility notes
# ---------------------
# 1. All paths are relative to the repository root; no setwd() calls are used.
# 2. Run this script from the repository root, or define NICHE_PROJECT_ROOT.
# 3. Occurrence coordinates must be in decimal degrees and should use WGS84,
#    consistent with the article.
# 4. Present-day climatic predictors should correspond to the 19 WorldClim 1.4
#    bioclimatic variables at 30 arc-seconds (~1 km), as used in the study.
# 5. `usdm::vifcor()` is used with method = "spearman" and th = 0.8. It removes
#    one variable from each highly correlated pair through the package's
#    stepwise procedure. `usdm::vifstep()` is then applied with th = 10.
#
# Suggested repository structure
# ------------------------------
# repository/
# |-- scripts/
# |   `-- 01_select_variables_spearman_vif.R
# |-- data/
# |   |-- occurrences/
# |   `-- climate/present/
# `-- results/
#     `-- variable_selection/
# ==============================================================================

# -------------------------------
# 0. User configuration
# -------------------------------

project_root <- Sys.getenv("NICHE_PROJECT_ROOT", unset = ".")

# Change this value when processing a different species.
species_name <- "Atlapetes_pileatus"

# Occurrence files may be CSV or DBF and must contain longitude and latitude
# columns. Edit the column names below only if your archived data use different
# names.
occurrence_file <- file.path(
  project_root,
  "data",
  "occurrences",
  paste0(species_name, ".csv")
)

longitude_column <- "lon"
latitude_column <- "lat"

# Directory containing the 19 WorldClim 1.4 bioclimatic layers (.asc).
climate_dir <- file.path(project_root, "data", "climate", "present")

results_dir <- file.path(
  project_root,
  "results",
  "variable_selection",
  species_name
)

# Thresholds reported in the article.
correlation_threshold <- 0.8
vif_threshold <- 10

# -------------------------------
# 1. Package requirements
# -------------------------------

required_packages <- c(
  "corrplot",
  "foreign",
  "raster",
  "usdm"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 2. Helper function
# -------------------------------

read_occurrences <- function(path, lon_col, lat_col) {
  if (!file.exists(path)) {
    stop("Occurrence file not found: ", path)
  }

  extension <- tolower(tools::file_ext(path))

  dat <- switch(
    extension,
    "csv" = utils::read.csv(path, stringsAsFactors = FALSE),
    "dbf" = foreign::read.dbf(path, as.is = TRUE),
    stop("Unsupported occurrence file format: .", extension)
  )

  missing_columns <- setdiff(c(lon_col, lat_col), names(dat))
  if (length(missing_columns) > 0) {
    stop(
      "Occurrence file is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  lon <- suppressWarnings(as.numeric(as.character(dat[[lon_col]])))
  lat <- suppressWarnings(as.numeric(as.character(dat[[lat_col]])))

  out <- data.frame(
    species = species_name,
    lon = lon,
    lat = lat,
    stringsAsFactors = FALSE
  )

  out <- out[
    is.finite(out$lon) & is.finite(out$lat) &
      out$lon >= -180 & out$lon <= 180 &
      out$lat >= -90 & out$lat <= 90,
    ,
    drop = FALSE
  ]

  if (nrow(out) == 0) {
    stop("No valid geographic coordinates were found in: ", path)
  }

  out
}

# -------------------------------
# 3. Read occurrence data
# -------------------------------

occurrences <- read_occurrences(
  occurrence_file,
  longitude_column,
  latitude_column
)

# -------------------------------
# 4. Read present-day climatic variables
# -------------------------------

climate_files <- list.files(
  climate_dir,
  pattern = "\\.asc$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(climate_files) == 0) {
  stop("No .asc climatic layers were found in: ", climate_dir)
}

if (length(climate_files) != 19) {
  warning(
    "The article used 19 WorldClim 1.4 bioclimatic variables, but ",
    length(climate_files),
    " .asc files were found. Verify the climate directory before interpreting results."
  )
}

climate_present <- raster::stack(climate_files)

# -------------------------------
# 5. Extract climatic values at occurrence localities
# -------------------------------

environment_values <- raster::extract(
  climate_present,
  occurrences[, c("lon", "lat")]
)

environment_values <- as.data.frame(environment_values, check.names = FALSE)

occurrence_environment <- data.frame(
  occurrences,
  environment_values,
  check.names = FALSE
)

# Retain only records with complete climatic information.
occurrence_environment <- stats::na.omit(occurrence_environment)

if (nrow(occurrence_environment) < 3) {
  stop("Too few occurrence localities contain complete climatic information.")
}

predictor_names <- names(climate_present)
environment_complete <- occurrence_environment[
  ,
  predictor_names,
  drop = FALSE
]

# Save occurrence records that contain complete climatic information.
utils::write.csv(
  occurrence_environment,
  file.path(results_dir, "occurrences_with_climate.csv"),
  row.names = FALSE
)

utils::write.csv(
  occurrence_environment[, c("species", "lon", "lat")],
  file.path(results_dir, "occurrences_complete.csv"),
  row.names = FALSE
)

# -------------------------------
# 6. Spearman correlation analysis (|r| < 0.8)
# -------------------------------

spearman_matrix <- stats::cor(
  environment_complete,
  method = "spearman",
  use = "complete.obs"
)

utils::write.csv(
  spearman_matrix,
  file.path(results_dir, "spearman_correlation_matrix_all_variables.csv"),
  row.names = TRUE
)

# Save a graphical representation of the complete Spearman correlation matrix.
grDevices::png(
  filename = file.path(results_dir, "spearman_correlation_matrix.png"),
  width = 2600,
  height = 2400,
  res = 220
)

corrplot::corrplot(
  spearman_matrix,
  method = "color",
  type = "lower",
  order = "AOE",
  diag = FALSE,
  addCoef.col = "black",
  number.cex = 0.55,
  tl.col = "black",
  tl.srt = 45,
  mar = c(0, 0, 2, 0),
  title = paste0(
    gsub("_", " ", species_name),
    " - Spearman correlation"
  )
)

grDevices::dev.off()

# Identify and remove variables involved in strong pairwise correlations.
# `vifcor()` uses the selected correlation method to find pairs at or above the
# threshold and iteratively removes one variable from each pair.
correlation_selection <- usdm::vifcor(
  environment_complete,
  th = correlation_threshold,
  method = "spearman"
)

environment_after_correlation <- usdm::exclude(
  environment_complete,
  correlation_selection
)

# -------------------------------
# 7. VIF filtering (VIF < 10)
# -------------------------------

vif_before_filtering <- usdm::vif(environment_after_correlation)

vif_selection <- usdm::vifstep(
  environment_after_correlation,
  th = vif_threshold,
  method = "spearman"
)

selected_environment <- usdm::exclude(
  environment_after_correlation,
  vif_selection
)

selected_variables <- names(selected_environment)

if (length(selected_variables) < 2) {
  stop("Fewer than two environmental variables remain after filtering.")
}

# -------------------------------
# 8. Verify the final variable set
# -------------------------------

final_spearman_matrix <- stats::cor(
  selected_environment,
  method = "spearman",
  use = "complete.obs"
)

final_vif <- usdm::vif(selected_environment)

# Maximum absolute off-diagonal Spearman correlation.
cor_check <- abs(final_spearman_matrix)
diag(cor_check) <- NA_real_
max_abs_spearman <- max(cor_check, na.rm = TRUE)
max_final_vif <- max(final_vif$VIF, na.rm = TRUE)

if (max_abs_spearman >= correlation_threshold) {
  warning(
    "At least one retained variable pair has |Spearman r| >= ",
    correlation_threshold,
    ". Inspect the saved correlation matrix."
  )
}

if (max_final_vif >= vif_threshold) {
  warning(
    "At least one retained variable has VIF >= ",
    vif_threshold,
    ". Inspect the saved VIF table."
  )
}

# -------------------------------
# 9. Save variable-selection outputs
# -------------------------------

utils::write.csv(
  data.frame(
    species = species_name,
    variable = selected_variables,
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "selected_variables.csv"),
  row.names = FALSE
)

utils::write.csv(
  final_spearman_matrix,
  file.path(results_dir, "spearman_correlation_matrix_selected_variables.csv"),
  row.names = TRUE
)

utils::write.csv(
  vif_before_filtering,
  file.path(results_dir, "vif_after_spearman_filter.csv"),
  row.names = FALSE
)

utils::write.csv(
  final_vif,
  file.path(results_dir, "vif_selected_variables.csv"),
  row.names = FALSE
)

capture.output(
  cat("Species:", gsub("_", " ", species_name), "\n"),
  cat("Occurrence records with complete climate data:", nrow(occurrence_environment), "\n"),
  cat("Initial number of climatic variables:", ncol(environment_complete), "\n"),
  cat("Spearman threshold: |r| <", correlation_threshold, "\n"),
  cat("VIF threshold: VIF <", vif_threshold, "\n\n"),
  cat("Spearman-correlation filtering:\n"),
  print(correlation_selection),
  cat("\nVIF filtering after Spearman filtering:\n"),
  print(vif_selection),
  cat("\nFinal retained variables:\n"),
  cat(paste(selected_variables, collapse = ", "), "\n"),
  cat("\nMaximum final |Spearman r|:", round(max_abs_spearman, 4), "\n"),
  cat("Maximum final VIF:", round(max_final_vif, 4), "\n"),
  file = file.path(results_dir, "variable_selection_summary.txt")
)

capture.output(
  utils::sessionInfo(),
  file = file.path(results_dir, "sessionInfo.txt")
)

message(
  "Variable selection completed for ",
  species_name,
  ". Retained variables: ",
  paste(selected_variables, collapse = ", "),
  ". Results saved to: ",
  normalizePath(results_dir, mustWork = FALSE)
)
