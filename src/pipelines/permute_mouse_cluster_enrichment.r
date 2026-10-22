#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# permute_mouse_cluster_enrichment.R
# Author: Antoine Beauchamp
# Created: October 20th, 2026
#
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option(
    "--pipeline-dir",
    type = "character",
    default = "data/mouse/derivatives/"
  ),
  make_option(
    "--params-id", 
  type = "character"
  ),
  make_option(
    "--background-file",
    type = "character",
    default = "data/enrichment/sagittal_gene_table_normalized_filtered.csv"
  ),
  make_option(
    "--modules-file",
    type = "character",
    default = "data/enrichment/Human_Reactome_June_01_2025_symbol.gmt"
  ),
  make_option(
    "--output-dir", 
    type = "character"
    ),
  make_option(
    "--models", 
    type = "character", 
    default = "model_names.csv"
    ),
  make_option(
    "--gene-score", 
    type = "numeric",
    default = 950
    ),
  make_option(
    "--n-permutations",
    type = "numeric", 
    default = 1000L
    ),
  make_option(
    "--n-workers",
    type = "numeric",
    default = 4L
    )
)

# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")

# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "enrichment.R"))


# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

pipeline_dir <- args[["pipeline-dir"]]

# Get pipeline parameters
metadata <- file.path(pipeline_dir, "metadata.csv")
params <- fetch_params_metadata(
  metadata = metadata,
  params = list(id = args[["params-id"]])
)

# Complete pipeline directory path with parameter set id
pipeline_dir <- file.path(pipeline_dir, args[["params-id"]])

# Get max number of clusters
nk_max <- params[["cluster_nk_max"]]

# Get cluster resolution
cluster_resolution <- params[["cluster_resolution"]]
cluster_resolution <- sprintf("%.1f", cluster_resolution)

# Define output directory
if (is.null(args[["output-dir"]])) {
  output_dir <- file.path(pipeline_dir, "enrichment", "StringDB_12.0_Bader_2025")
} else {
  output_dir <- args[["output-dir"]]
}
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

# Cluster directory path
clusters_dir <- file.path(pipeline_dir, "clusters", paste0("resolution_", cluster_resolution))

# Import clusters
clusters_file <- file.path(clusters_dir, "clusters.csv")
clusters <- read_csv(file = clusters_file, show_col_types = FALSE)

# Import model names dictionary
models_file <- file.path(pipeline_dir, args[["models"]])
models <- get_model_genes(models_file)

# Join with model gene names (and filter for single genes)
clusters <- clusters |>
    rename(file = ID) |>
    right_join(models, by = "file") |>
    select(-file) |>
    rename(model = ID)

# Import Bader Reactome pathway modules
modules <- importMsigDBGMT(args[["modules-file"]])

# Construct background set of genes
message("Building background gene set ...")
background_set <- read_csv(file = args[["background-file"]], show_col_types = FALSE) |>
  pull(msg.genes.acronym) |>
  as.character() |> 
  unique() |> 
  fetch_gene_homologs(species = "mouse") |> 
  filter(!is.na(GeneHuman)) |> 
  pull(GeneHuman)

# Compute gene interaction network for all genes in clusters
message("Building gene interaction network ...")
gene_network <- build_gene_interaction_network(
  genes = unique(clusters$gene),
  score = args[["gene-score"]]
)

# Build data frame of human homologs for all target genes
gene_network_hom <- fetch_gene_homologs(
  genes = unique(c(gene_network[["gene_A"]], 
                   gene_network[["gene_B"]])), 
  species = "mouse"
)

# Get human gene homologues for gene network
gene_network <- gene_network |> 
  left_join(select(gene_network_hom, GeneInput, GeneHuman),
            by = c("gene_A" = "GeneInput"), 
            relationship = "many-to-many") |> 
  select(-gene_A) |> 
  rename(gene_A = GeneHuman) |> 
  left_join(select(gene_network_hom, GeneInput, GeneHuman),
            by = c("gene_B" = "GeneInput"),
            relationship = "many-to-many") |> 
  select(-gene_B) |> 
  rename(gene_B = GeneHuman)


message("Running permutations for cluster enrichment analysis ...")
ti <- Sys.time()
permutation_status <- run_cluster_enrichment_permutations(
  clusters = clusters,
  gene_network = gene_network,
  background_set = background_set,
  modules = modules,
  string_score = args[["gene-score"]],
  nk_max = nk_max,
  output_dir = output_dir,
  B = args[["n-permutations"]], 
  n_workers = args[["n-workers"]]
)
tf <- Sys.time()
print(difftime(tf, ti, units = "hours"))

if (any(!permutation_status$success)) {
  warning("Some permutations failed; inspect permutation_status for details.")
}
