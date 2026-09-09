# ==============================================================================
# Ecological niche overlap, equivalency, and similarity tests
# ==============================================================================
# Supplementary analysis script for:
# Sánchez-González, L. A., Prieto-Torres, D. A., & Espinosa-Chávez, O. J. (2026).
# Ecological Niche Differences Underlie the Assembly of Endemic Birds in
# Acrobiont Forests of Northern Mesoamerica. Journal of Biogeography, 53, e70141.
# https://doi.org/10.1111/jbi.70141
#
# Purpose
# -------
# This script performs a pairwise comparison of the environmental niches of two
# species using the PCA-env framework implemented in the R package `ecospat`.
# It calculates Schoener's D niche overlap and runs niche equivalency and
# background similarity tests following Broennimann et al. (2012).
#
# Reproducibility notes
# ---------------------
# 1. All paths are relative to the repository root. Do not use setwd().
# 2. Run this script from the repository root, or define the environment variable
#    NICHE_PROJECT_ROOT with the path to the repository.
# 3. The input occurrence files must contain columns named `lon` and `lat` in
#    decimal degrees (WGS84 in the original analyses).
# 4. The script accepts occurrence files in either CSV or DBF format.
# 5. Environmental predictors are screened using the thresholds reported in the
#    article: Spearman |r| < 0.8 followed by VIF < 10. The same retained set is
#    then used for both species in each pairwise PCA-env comparison.
# 6. The original workflow used 0.008333 degrees (~30 arc-seconds) as the spatial
#    sampling resolution, R = 300 for the PCA-env grid, and 1,000 randomizations
#    for equivalency and similarity tests.
# 7. This script implements the PCA-env analyses with `ecospat` only. The
#    `humboldt` package mentioned in the article is intentionally not included.
# 8. Randomization-based tests are stochastic. If the original analysis used a
#    fixed random seed, enter it below. Otherwise leave `random_seed <- NULL`.
#
# Suggested repository structure
# ------------------------------
# repository/
# |-- scripts/
# |   `-- niche_overlap_equivalency_similarity.R
# |-- data/
# |   |-- occurrences/
# |   |-- climate/present/
# |   `-- calibration_areas/
# `-- results/
#
# IMPORTANT
# ---------
# This script is written as a transparent pairwise template. Change only the
# configuration block below when running a different species comparison.
# ============================================================================== 

# -------------------------------
# 0. User configuration
# -------------------------------

project_root <- Sys.getenv("NICHE_PROJECT_ROOT", unset = ".")

# Species names are used for labels and output file names.
species_1_name <- "Atlapetes_pileatus"
species_2_name <- "Selasphorus_heloisa"

# Input files. Change these names if the repository uses different file names.
species_1_file <- file.path(
  project_root, "data", "occurrences", paste0(species_1_name, ".dbf")
)
species_2_file <- file.path(
  project_root, "data", "occurrences", paste0(species_2_name, ".dbf")
)

climate_dir <- file.path(project_root, "data", "climate", "present")
calibration_dir <- file.path(project_root, "data", "calibration_areas")

species_1_M_file <- file.path(
  calibration_dir, paste0(species_1_name, ".shp")
)
species_2_M_file <- file.path(
  calibration_dir, paste0(species_2_name, ".shp")
)

comparison_id <- paste0(species_1_name, "_vs_", species_2_name)
results_dir <- file.path(
  project_root, "results", "niche_comparisons", comparison_id
)

# Analysis parameters reported in the article / used in the original workflow.
correlation_threshold <- 0.8
vif_threshold <- 10
sampling_resolution <- 0.008333
pca_grid_resolution <- 300
n_randomizations <- 1000
random_seed <- NULL

# -------------------------------
# 1. Package requirements
# -------------------------------

required_packages <- c(
  "ade4",
  "ecospat",
  "foreign",
  "raster",
  "sf",
  "sp",
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

# Create the output directory if it does not exist.
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Use a fixed seed only when one is known from the original analysis.
if (!is.null(random_seed)) {
  set.seed(random_seed)
}

# -------------------------------
# 2. Helper functions
# -------------------------------

read_occurrences <- function(path) {
  if (!file.exists(path)) {
    stop("Occurrence file not found: ", path)
  }

  extension <- tolower(tools::file_ext(path))

  dat <- switch(
    extension,
    "dbf" = foreign::read.dbf(path, as.is = TRUE),
    "csv" = utils::read.csv(path, stringsAsFactors = FALSE),
    stop("Unsupported occurrence file format: .", extension)
  )

  required_columns <- c("lon", "lat")
  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0) {
    stop(
      "Occurrence file ", basename(path),
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat$lon <- suppressWarnings(as.numeric(as.character(dat$lon)))
  dat$lat <- suppressWarnings(as.numeric(as.character(dat$lat)))

  dat <- dat[
    is.finite(dat$lon) & is.finite(dat$lat),
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) {
    stop("No valid longitude/latitude records were found in: ", path)
  }

  dat
}

read_calibration_area <- function(path, target_crs) {
  if (!file.exists(path)) {
    stop("Calibration-area shapefile not found: ", path)
  }

  area_sf <- sf::st_read(path, quiet = TRUE)

  if (is.na(sf::st_crs(area_sf))) {
    stop("The calibration-area shapefile has no defined CRS: ", path)
  }

  target_sf_crs <- sf::st_crs(target_crs)

  if (!is.na(target_sf_crs) && sf::st_crs(area_sf) != target_sf_crs) {
    area_sf <- sf::st_transform(area_sf, target_sf_crs)
  }

  # `raster::crop()` and `raster::mask()` operate reliably with Spatial objects.
  methods::as(area_sf, "Spatial")
}

species_label <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

# -------------------------------
# 3. Read occurrence data
# -------------------------------

species_1 <- read_occurrences(species_1_file)
species_2 <- read_occurrences(species_2_file)

# Keep coordinates explicit. Raster extraction expects x = longitude and
# y = latitude when using geographic coordinates.
all_occurrences <- rbind(
  species_1[, c("lon", "lat")],
  species_2[, c("lon", "lat")]
)

# -------------------------------
# 4. Read present-day climatic layers
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

climate_present <- raster::stack(climate_files)

# Extract climatic values at all occurrence localities. These values are used
# to diagnose multicollinearity before the PCA-env niche comparison.
occurrence_environment <- raster::extract(
  climate_present,
  all_occurrences[, c("lon", "lat")]
)
occurrence_environment <- as.data.frame(occurrence_environment)
occurrence_environment <- stats::na.omit(occurrence_environment)

if (nrow(occurrence_environment) < 2) {
  stop("Too few occurrence localities contain complete climatic information.")
}

# -------------------------------
# 5. Environmental-variable screening
# -------------------------------

# The article reports two criteria for the first variable-selection approach:
# Spearman |r| < 0.8 and VIF < 10. Here they are applied sequentially to the
# climatic values extracted from the pooled occurrence localities of the two
# species being compared. This produces one common predictor set for the
# pairwise PCA-env analysis.

# Save the complete Spearman correlation matrix for inspection.
correlation_matrix <- stats::cor(
  occurrence_environment,
  method = "spearman",
  use = "complete.obs"
)

utils::write.csv(
  correlation_matrix,
  file.path(results_dir, "spearman_correlation_matrix_all_variables.csv"),
  row.names = TRUE
)

# Step 1: remove variables involved in strong pairwise correlations.
# `vifcor()` supports Spearman correlations directly and applies the requested
# correlation threshold through a stepwise procedure.
correlation_selection <- usdm::vifcor(
  occurrence_environment,
  th = correlation_threshold,
  method = "spearman"
)

environment_after_correlation <- usdm::exclude(
  occurrence_environment,
  correlation_selection
)

# Step 2: remove remaining variables with VIF >= 10.
vif_initial <- usdm::vif(environment_after_correlation)

vif_selection <- usdm::vifstep(
  environment_after_correlation,
  th = vif_threshold,
  method = "spearman"
)

selected_environment <- usdm::exclude(
  environment_after_correlation,
  vif_selection
)

env_names <- names(selected_environment)

if (length(env_names) < 2) {
  stop("Fewer than two environmental variables remain after filtering.")
}

# Verify that the final set meets both published thresholds.
final_correlation_matrix <- stats::cor(
  selected_environment,
  method = "spearman",
  use = "complete.obs"
)
final_vif <- usdm::vif(selected_environment)

utils::write.csv(
  final_correlation_matrix,
  file.path(results_dir, "spearman_correlation_matrix_selected_variables.csv"),
  row.names = TRUE
)

utils::write.csv(
  final_vif,
  file.path(results_dir, "vif_selected_variables.csv"),
  row.names = FALSE
)

utils::write.csv(
  data.frame(variable = env_names),
  file.path(results_dir, "selected_environmental_variables.csv"),
  row.names = FALSE
)

capture.output(
  cat("Spearman-correlation filtering (threshold =", correlation_threshold, "):\n"),
  print(correlation_selection),
  cat("\nVIF values after Spearman filtering:\n"),
  print(vif_initial),
  cat("\nVIF filtering (threshold =", vif_threshold, "):\n"),
  print(vif_selection),
  cat("\nFinal retained variables:\n"),
  cat(paste(env_names, collapse = ", "), "\n"),
  file = file.path(results_dir, "environmental_variable_selection.txt")
)

# -------------------------------
# 6. Restrict climatic layers to each species' accessible area (M)
# -------------------------------

raster_crs <- raster::projection(climate_present)

M_species_1 <- read_calibration_area(species_1_M_file, raster_crs)
M_species_2 <- read_calibration_area(species_2_M_file, raster_crs)

# Species 1
cropped_env_1 <- raster::crop(
  climate_present,
  raster::extent(M_species_1)
)
masked_env_1 <- raster::mask(cropped_env_1, M_species_1)

# Species 2
cropped_env_2 <- raster::crop(
  climate_present,
  raster::extent(M_species_2)
)
masked_env_2 <- raster::mask(cropped_env_2, M_species_2)

# Apply the same final predictor set to both accessible areas (M). Subsetting by
# name is used here because the final set reflects exclusions from both the
# Spearman-correlation and VIF steps.
missing_env_1 <- setdiff(env_names, names(masked_env_1))
missing_env_2 <- setdiff(env_names, names(masked_env_2))

if (length(missing_env_1) > 0 || length(missing_env_2) > 0) {
  stop(
    "One or more retained environmental variables are missing from the climate stack. ",
    "Species 1 missing: ", paste(missing_env_1, collapse = ", "),
    "; Species 2 missing: ", paste(missing_env_2, collapse = ", ")
  )
}

selected_env_1 <- masked_env_1[[env_names]]
selected_env_2 <- masked_env_2[[env_names]]

# -------------------------------
# 7. Build environmental background data for each species
# -------------------------------

# Generate one point for each valid raster cell in M, using the first retained
# variable as a spatial template, and then extract all retained variables.
background_points_1 <- raster::rasterToPoints(
  selected_env_1[[1]],
  spatial = TRUE
)
background_points_2 <- raster::rasterToPoints(
  selected_env_2[[1]],
  spatial = TRUE
)

background_values_1 <- raster::extract(selected_env_1, background_points_1)
background_values_2 <- raster::extract(selected_env_2, background_points_2)

background_1 <- data.frame(
  sp::coordinates(background_points_1),
  background_values_1,
  check.names = FALSE
)
background_2 <- data.frame(
  sp::coordinates(background_points_2),
  background_values_2,
  check.names = FALSE
)

names(background_1)[1:2] <- c("lon", "lat")
names(background_2)[1:2] <- c("lon", "lat")

background_1 <- stats::na.omit(background_1)
background_2 <- stats::na.omit(background_2)

# -------------------------------
# 8. Match occurrences to environmental cells
# -------------------------------

occ_species_1 <- species_1[, c("lon", "lat")]
occ_species_2 <- species_2[, c("lon", "lat")]

# `ecospat.sample.envar()` samples environmental values at occurrence locations
# and removes duplicate records at the selected spatial resolution.
occ_env_1 <- stats::na.exclude(
  ecospat::ecospat.sample.envar(
    dfsp = occ_species_1,
    colspxy = 1:2,
    colspkept = 1:2,
    dfvar = background_1,
    colvarxy = 1:2,
    colvar = "all",
    resolution = sampling_resolution
  )
)

occ_env_2 <- stats::na.exclude(
  ecospat::ecospat.sample.envar(
    dfsp = occ_species_2,
    colspxy = 1:2,
    colspkept = 1:2,
    dfvar = background_2,
    colvarxy = 1:2,
    colvar = "all",
    resolution = sampling_resolution
  )
)

# Ensure that occurrence and background data contain the same columns and order.
expected_columns <- c("lon", "lat", env_names)

occ_env_1 <- occ_env_1[, expected_columns, drop = FALSE]
occ_env_2 <- occ_env_2[, expected_columns, drop = FALSE]
background_1 <- background_1[, expected_columns, drop = FALSE]
background_2 <- background_2[, expected_columns, drop = FALSE]

# Add an indicator distinguishing occurrence records (1) from background cells (0).
occ_env_1$species_occ <- 1
background_1$species_occ <- 0
occ_env_2$species_occ <- 1
background_2$species_occ <- 0

data_species_1 <- rbind(occ_env_1, background_1)
data_species_2 <- rbind(occ_env_2, background_2)

# -------------------------------
# 9. PCA-env
# -------------------------------

# PCA is calculated using the combined environmental space of both species.
pca_input <- rbind(
  data_species_1[, env_names, drop = FALSE],
  data_species_2[, env_names, drop = FALSE]
)

pca_env <- ade4::dudi.pca(
  pca_input,
  center = TRUE,
  scale = TRUE,
  scannf = FALSE,
  nf = 2
)

# Save PCA variable contributions.
grDevices::png(
  filename = file.path(results_dir, "pca_variable_contributions.png"),
  width = 2000,
  height = 1600,
  res = 200
)
ecospat::ecospat.plot.contrib(
  contrib = pca_env$co,
  eigen = pca_env$eig
)
grDevices::dev.off()

# PCA scores for the complete environmental space.
scores_global <- pca_env$li

# Project occurrences and species-specific backgrounds onto the PCA axes.
scores_occ_1 <- ade4::suprow(
  pca_env,
  data_species_1[
    data_species_1$species_occ == 1,
    env_names,
    drop = FALSE
  ]
)$li

scores_occ_2 <- ade4::suprow(
  pca_env,
  data_species_2[
    data_species_2$species_occ == 1,
    env_names,
    drop = FALSE
  ]
)$li

scores_background_1 <- ade4::suprow(
  pca_env,
  data_species_1[
    data_species_1$species_occ == 0,
    env_names,
    drop = FALSE
  ]
)$li

scores_background_2 <- ade4::suprow(
  pca_env,
  data_species_2[
    data_species_2$species_occ == 0,
    env_names,
    drop = FALSE
  ]
)$li

# -------------------------------
# 10. Occupancy-density grids in environmental space
# -------------------------------

grid_species_1 <- ecospat::ecospat.grid.clim.dyn(
  glob = scores_global,
  glob1 = scores_background_1,
  sp = scores_occ_1,
  R = pca_grid_resolution,
  th.sp = 0
)

grid_species_2 <- ecospat::ecospat.grid.clim.dyn(
  glob = scores_global,
  glob1 = scores_background_2,
  sp = scores_occ_2,
  R = pca_grid_resolution,
  th.sp = 0
)

# -------------------------------
# 11. Schoener's D niche overlap and niche dynamics
# -------------------------------

overlap <- ecospat::ecospat.niche.overlap(
  grid_species_1,
  grid_species_2,
  cor = TRUE
)

schoener_D <- unname(overlap$D)
schoener_D_rounded <- round(schoener_D, 3)

niche_dynamics <- ecospat::ecospat.niche.dyn.index(
  grid_species_1,
  grid_species_2,
  intersection = NA
)

# Save a summary table with the observed overlap value.
utils::write.csv(
  data.frame(
    species_1 = species_label(species_1_name),
    species_2 = species_label(species_2_name),
    schoener_D = schoener_D
  ),
  file.path(results_dir, "niche_overlap_summary.csv"),
  row.names = FALSE
)

# Save the niche-dynamics plot.
grDevices::png(
  filename = file.path(results_dir, "niche_overlap_pca_env.png"),
  width = 1800,
  height = 1600,
  res = 200
)
ecospat::ecospat.plot.niche.dyn(
  grid_species_1,
  grid_species_2,
  quant = 0.25,
  colZ1 = "Green",
  colZ2 = "Red",
  title = paste0(
    "Niche overlap: ",
    species_label(species_1_name),
    " vs. ",
    species_label(species_2_name),
    " | Schoener's D = ",
    schoener_D_rounded
  ),
  name.axis1 = "PC1",
  name.axis2 = "PC2"
)
grDevices::dev.off()

# -------------------------------
# 12. Niche equivalency test
# -------------------------------

equivalency_test <- ecospat::ecospat.niche.equivalency.test(
  grid_species_1,
  grid_species_2,
  rep = n_randomizations
)

# -------------------------------
# 13. Background similarity tests
# -------------------------------

# Similarity tests are directional; therefore, both directions are evaluated.
similarity_1_to_2 <- ecospat::ecospat.niche.similarity.test(
  grid_species_1,
  grid_species_2,
  rep = n_randomizations,
  rand.type = 2
)

similarity_2_to_1 <- ecospat::ecospat.niche.similarity.test(
  grid_species_2,
  grid_species_1,
  rep = n_randomizations,
  rand.type = 2
)

# -------------------------------
# 14. Save equivalency and similarity plots
# -------------------------------

grDevices::png(
  filename = file.path(results_dir, "equivalency_similarity_tests.png"),
  width = 3000,
  height = 1100,
  res = 200
)
graphics::par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

ecospat::ecospat.plot.overlap.test(
  equivalency_test,
  "D",
  "Equivalency"
)

ecospat::ecospat.plot.overlap.test(
  similarity_1_to_2,
  "D",
  paste0(
    "Similarity: ",
    species_label(species_1_name),
    " -> ",
    species_label(species_2_name)
  )
)

ecospat::ecospat.plot.overlap.test(
  similarity_2_to_1,
  "D",
  paste0(
    "Similarity: ",
    species_label(species_2_name),
    " -> ",
    species_label(species_1_name)
  )
)

grDevices::dev.off()

# -------------------------------
# 15. Save analysis objects and session information
# -------------------------------

saveRDS(
  list(
    species_1 = species_1_name,
    species_2 = species_2_name,
    environmental_variables = env_names,
    pca = pca_env,
    grid_species_1 = grid_species_1,
    grid_species_2 = grid_species_2,
    overlap = overlap,
    niche_dynamics = niche_dynamics,
    equivalency_test = equivalency_test,
    similarity_1_to_2 = similarity_1_to_2,
    similarity_2_to_1 = similarity_2_to_1
  ),
  file.path(results_dir, "niche_comparison_objects.rds")
)

capture.output(
  utils::sessionInfo(),
  file = file.path(results_dir, "sessionInfo.txt")
)

message(
  "Analysis completed: ", comparison_id,
  " | Schoener's D = ", schoener_D_rounded,
  " | Results saved to: ", normalizePath(results_dir, mustWork = FALSE)
)
