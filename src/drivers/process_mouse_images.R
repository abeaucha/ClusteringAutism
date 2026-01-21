#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# process_mouse_images.R
# Author: Jacob Ellegood, Antoine Beauchamp
# Created: April 29th, 2025
#
#

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--pipeline-dir", type = "character"),
  make_option("--registration-dir", type = "character"),
  make_option("--mask", type = "character"),
  make_option("--cluster-resolution", type = "numeric"),
  make_option("--cluster-nk-max", type = "numeric"),
  make_option("--cluster-metric", type = "character"),
  make_option("--cluster-K", type = "numeric"),
  make_option("--cluster-sigma", type = "numeric"),
  make_option("--cluster-t", type = "numeric"),
  make_option("--cluster-file", type = "character"),
  make_option("--cluster-affinity-file", type = "character"),
  make_option("--centroid-method", type = "character")
)


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "processing.R"))


# Main -----------------------------------------------------------------------

args <- parse_args(OptionParser(option_list = option_list))
pipeline_dir <- args[["pipeline-dir"]]
registration_dir <- args[["registration-dir"]]
mask <- args[["mask"]]
cluster_resolution <- args[["cluster-resolution"]]
cluster_nk_max <- args[["cluster-nk-max"]]
cluster_metric <- args[["cluster-metric"]]
cluster_K <- args[["cluster-K"]]
cluster_sigma <- args[["cluster-sigma"]]
cluster_t <- args[["cluster-t"]]
cluster_file <- args[["cluster-file"]]
affinity_file <- args[["cluster-affinity-file"]]
centroid_method <- args[["centroid-method"]]

resolution_mm <- cluster_resolution
resolution_um <- cluster_resolution * 1000

# List of scanbase files
scanbase_files <- list(
  scans = "scanbase_40um - Scans_31Jan22.csv",
  studies = "scanbase_40um - Studies_Feb2023.csv",
  genotypes = "scanbase_40um - Genotypes_Feb2023.csv"
)

# Prepend path to scanbase files
scanbase_files <- map(.x = scanbase_files, .f = function(x) {
  file.path(registration_dir, "resources", x)
})

# List of paths
paths <- list(
  jacobians = file.path(registration_dir, "jacobians"),
  effect_sizes = file.path(pipeline_dir, "effect_sizes"),
  clusters = file.path(
    pipeline_dir,
    "clusters",
    paste0("resolution_", resolution_mm)
  ),
  centroids = file.path(
    pipeline_dir,
    "centroids",
    paste0("resolution_", resolution_mm)
  )
)

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

jacobians <- c("absolute", "relative")
list_es <- vector(mode = "list", length = length(jacobians))
names(list_es) <- jacobians
for (jtype in jacobians) {
  # Compute effect sizes
  message(paste("Computing", jtype, "effect size images..."))
  es_matrix <- MakeEffectSizeMaps(
    jdtype = str_to_title(jtype),
    model_list = model_list,
    resolution = as.character(resolution_um),
    dir_determinants = paths[["jacobians"]],
    output_phenotype_dir = paths[["effect_sizes"]],
    base_directory = pipeline_dir,
    boot = "N"
  )

  # Replace NaN with 0
  es_matrix[is.nan(es_matrix)] <- 0

  # Add matrix to list
  list_es[[jtype]] <- es_matrix
}


# Clustering ------------------------------------------------------------------

message("Generating clusters...")

# Run similarity network fusion
outfile <- file.path(paths[["clusters"]], affinity_file)
W <- similarity_network(
  x = list_es,
  metric = cluster_metric,
  K = cluster_K,
  sigma = cluster_sigma,
  t = cluster_t,
  outfile = outfile
)

# Identify clusters from fused affinity matrix
outfile <- file.path(paths[["clusters"]], cluster_file)
clusters <- create_clusters(W = W, nk = cluster_nk_max)
clusters[["ID"]] <- paste0(clusters[["ID"]], ".mnc")
write.csv(x = clusters, file = outfile, row.names = FALSE)


# Centroids -------------------------------------------------------------------

for (jtype in jacobians) {
  message(paste("Generating", jtype, "centroid images..."))
  clusters_jtype <- clusters %>%
    as_tibble() %>%
    mutate(
      ID = file.path(paths[["effect_sizes"]], resolution_um, ID),
      ID = str_replace(
        ID,
        ".mnc",
        paste0("_ES_", str_to_title(jtype), "_", resolution_um, ".mnc")
      )
    ) %>%
    column_to_rownames("ID")

  outdir <- file.path(paths[["centroids"]], jtype)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  for (j in 1:ncol(clusters_jtype)) {
    compute_cluster_centroids(
      i = j,
      clusters = clusters_jtype,
      mask = mask,
      method = centroid_method,
      outdir = outdir
    )
  }
}
