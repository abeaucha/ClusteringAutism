#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# process_mouse_images.R
# Author: Jacob Ellegood, Antoine Beauchamp
# Created: April 29th, 2025
# 
# Adapted from Clustering_Pipeline.R


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "Clustering_Functions_AB.R"))


# Main -----------------------------------------------------------------------

# Parameter set ID
params_id <- "001"

# Image resolution 
resolution_um <- 200
resolution_mm <- resolution_um/1000

# Cluster parameters
cluster_nk_max <- 2
cluster_metric <- "correlation"
cluster_K <- 10
cluster_sigma <- 0.5
cluster_t <- 20

# Centroid parameters
centroid_method <- "mean"

# Pipeline directory
pipeline_dir <- "data/mouse/derivatives/"

# Export parameter set to metadata
metadata <- file.path(pipeline_dir, "metadata.csv")
df_metadata <- tibble(dataset = "MICe",
                      resolution = resolution_mm,
                      cluster_resolution = resolution_mm,
                      cluster_nk_max = nk_max,
                      cluster_metric = cluster_metric,
                      cluster_K = cluster_K,
                      cluster_sigma = cluster_sigma,
                      cluster_t = cluster_t, 
                      centroid_method = centroid_method,
                      id = params_id)
write_csv(x = df_metadata, file = metadata)

# Update pipeline directory with parameter set ID
pipeline_dir <- file.path(pipeline_dir, params_id)

# Registration directory
registration_dir <- "data/mouse/registration/"

# List of scanbase files
scanbase_files <- list(scans = "scanbase_40um - Scans_31Jan22.csv",
                       studies = "scanbase_40um - Studies_Feb2023.csv",
                       genotypes = "scanbase_40um - Genotypes_Feb2023.csv")

# Prepend path to scanbase files
scanbase_files <- map(.x = scanbase_files, 
                      .f = function(x) {
                        file.path(registration_dir, "resources", x)
                      })

# List of paths
paths <- list(
  jacobians = file.path(registration_dir, "jacobians"),
  effect_sizes = file.path(pipeline_dir, "effect_sizes"),
  clusters = file.path(pipeline_dir, "clusters", paste0("resolution_", resolution_mm)),
  centroids = file.path(pipeline_dir, "centroids", paste0("resolution_", resolution_mm))
)

# Create output directories
for (path in paths[-1]) {dir.create(path, recursive = TRUE)}

# Generate list of models
model_list <- MakeModelList(
  scanbase_scans_file = scanbase_files[["scans"]],
  scanbase_studies_file = scanbase_files[["studies"]],
  scanbase_genotypes_file = scanbase_files[["genotypes"]],
  scanbase_focus = "Autism",
  scanbase_sample_size_threshold = 6,
  base_directory = pipeline_dir
)


# Effect sizes ----------------------------------------------------------------

message("Computing absolute effect size images...")

effect_size_data_matrix_abs <- MakeEffectSizeMaps(
  jdtype = "Absolute", model_list = model_list, 
  resolution = as.character(resolution_um),
  dir_determinants = paths[["jacobians"]],
  output_phenotype_dir = paths[["effect_sizes"]],
  base_directory = pipeline_dir, boot = "N"
)

message("Computing relative effect size images...")

effect_size_data_matrix_rel <- MakeEffectSizeMaps(
  jdtype = "Relative", model_list = model_list,
  resolution = as.character(resolution_um),
  dir_determinants = paths[["jacobians"]],
  output_phenotype_dir = paths[["effect_sizes"]],
  base_directory = pipeline_dir, boot = "N"
)

# Replace NaN with 0
effect_size_data_matrix_abs[is.nan(effect_size_data_matrix_abs)] <- 0
effect_size_data_matrix_rel[is.nan(effect_size_data_matrix_rel)] <- 0


# Clustering ------------------------------------------------------------------

message("Generating clusters...")

W <- SNFCombine(Data1 = effect_size_data_matrix_abs,
                Data2 = effect_size_data_matrix_rel,
                K = cluster_K, alpha = cluster_sigma, 
                T = cluster_t, distfunc = "cor",
                output_dir = paths[["clusters"]])


Clusters <- CreateClusters(num_clusters = as.character(cluster_nk_max),
                           cluster_method = "spectral",
                           output_dir = paths[["clusters"]],
                           NameFile = file.path(registration_dir, "Names_Paper.csv"))


# Centroids -------------------------------------------------------------------

message("Generating cluster centroids...")
jacobians <- c("absolute", "relative")
for (jtype in jacobians) {
  CreateClusterAnatomyMaps(num_clusters = as.character(cluster_nk_max),
                           cluster_method = "spectral", 
                           average_kind = centroid_method,
                           volume_type = jtype, 
                           resolution = as.character(resolution_um),
                           dir_determinants = paths[["jacobians"]],
                           output_dir = paths[["centroids"]])
}

