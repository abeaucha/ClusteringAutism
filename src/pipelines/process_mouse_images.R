#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# process_mouse_images.R
# Author: Jacob Ellegood, Antoine Beauchamp
# Created: April 29th, 2025
# 
# Adapted from Clustering_Pipeline.R


# Packages -------------------------------------------------------------------

# library(tidyverse)
# library(magrittr)
# library(RMINC)
# library(glue)
# library(MRIcrotome)
# library(grid)
# library(SNFtool)
# library(tmod)
# library(Exact)


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "Clustering_Functions_AB.R"))


# Main -----------------------------------------------------------------------

# Parameter set ID
params_id <- "100"

# Registration directory
registration_dir <- "data/mouse/registration/"

# Pipeline directory
pipeline_dir <- "data/mouse/derivatives/"
pipeline_dir <- file.path(pipeline_dir, params_id)

# Image resolution 
resolution_um <- 200
resolution_mm <- resolution_um/1000

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


# Effect sizes

effect_size_data_matrix_abs <- MakeEffectSizeMaps(
  jdtype = "Absolute", model_list = model_list, 
  resolution = as.character(resolution_um),
  dir_determinants = paths[["jacobians"]],
  output_phenotype_dir = paths[["effect_sizes"]],
  base_directory = pipeline_dir, boot = "N"
)

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


# Clusters


W <- SNFCombine(Data1 = effect_size_data_matrix_abs,
                Data2 = effect_size_data_matrix_rel,
                K = 10, alpha = 0.5, T = 20, distfunc = "cor",
                output_dir = paths[["clusters"]])


Clusters <- CreateClusters(num_clusters="10", cluster_method="spectral",
                           output_dir = paths[["clusters"]],
                           NameFile = file.path(registration_dir, "Names_Paper.csv"))


# Centroids

jacobians <- c("absolute", "relative")
for (jtype in jacobians) {
  CreateClusterAnatomyMaps(num_clusters = "10", cluster_method = "spectral",
                           average_kind = "mean", volume_type = jtype, 
                           resolution = as.character(resolution_um), 
                           dir_determinants = paths[["jacobians"]],
                           output_dir = paths[["centroids"]])
}



# Organization

# file.rename(file.path(paths[["clusters"]], "Clusters.csv"),
#             file.path(paths[["clusters"]], "clusters.csv"))
# 
# file.rename(file.path(paths[["clusters"]], "WMatrix.RData"),
#             file.path(paths[["clusters"]], "affinity.RData"))
# 
# 
# 
# jtype <- "absolute"
# centroid_dir_j <- file.path(paths[["centroids"]], jtype)
# if (!(dir.exists(centroid_dir_j))) {
#   dir.create(centroid_dir_j, recursive = TRUE)
# }

