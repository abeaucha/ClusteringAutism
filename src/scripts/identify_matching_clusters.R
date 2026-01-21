#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# identify_matching_clusters.R
# Author: Antoine Beauchamp
# Created: June 10th, 2025
#
# Identify matching clusters between data sets

# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "analysis.R"))


# Main -----------------------------------------------------------------------

# Parameter set ID
params_id <- 375

# Pipeline directory
pipeline_dir <- "data/cross_species/v3/"

# Import cluster similarity
similarity <- import_similarity(
  param_id = params_id,
  pipeline_dir = pipeline_dir
)

# Import cluster similarity permutations
permutations <- import_similarity_permutations(
  param_id = params_id,
  pipeline_dir = pipeline_dir
)

# Compute significance of correlations
df_sim_pvals <- compute_similarity_significance(
  similarity = similarity,
  permutations = permutations
)

# Identify matching clusters
df_matching <- df_sim_pvals %>%
  filter(pval < 0.05)

# Export matching clusters
outfile <- "matching_clusters.csv"
outfile <- file.path(pipeline_dir, params_id, outfile)
write_csv(x = df_matching, file = outfile)
