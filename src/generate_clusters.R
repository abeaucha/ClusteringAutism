#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# generate_clusters.R
# Author: Antoine Beauchamp, Jacob Ellegood
# Created: December 24th, 2023
#
# Identify clusters of patients based on effect size matrices.
#
# Description
# -----------
# This script identifies clusters of human participants using matrices of
# absolute and relative voxel-wise effect sizes. Absolute and relative
# participant affinity matrices are calculated using the effect size matrices.
# The affinity matrices are then combined into a single participant affinity
# matrix using similarity network fusion. Clusters are identified based on this
# fused affinity matrix using spectral clustering.


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--file1",
              type = "character",
              help = paste("Path to file (.csv) containing first effect size",
                           "matrix.")),
  make_option("--file2",
              type = "character",
              help = paste("Path to file (.csv) containing second effect size",
                           "matrix.")),
  make_option("--rownames",
              type = "character",
              help = "Column in the input files containing row names."),
  make_option("--nk-max",
              type = "numeric",
              default = 10,
              help = paste("Maximum number of clusters to identify. Solutions",
                           "will be obtained for nk = 2 to nk = --nk-max",
                           "[default %default]")),
  make_option("--metric",
              type = "character",
              default = "correlation",
              help = paste("Similarity metric used to compute the SNF affinity",
                           "matrices. [default %default]")),
  make_option("--K",
              type = "numeric",
              default = 10,
              help = paste("Number of nearest-neighbours used to compute",
                           "the SNF affinity matrices. [default %default]")),
  make_option("--sigma",
              type = "numeric",
              default = 0.5,
              help = paste("Variance for the local model in the SNF",
                           "affinity matrices. [default %default]")),
  make_option("--t",
              type = "numeric",
              default = 20,
              help = paste("Number of iterations for the diffusion",
                           "process in SNF. [default %default]")),
  make_option("--cluster-file",
              type = "character",
              default = "clusters.csv",
              help = paste("Path to the file (.csv) in which to save the",
                           "cluster assignments. [default %default]")),
  make_option("--affinity-file",
              type = "character",
              help = paste("Path to file (.csv) in which to save the SNF",
                           "affinity matrix. If NULL, the affinity matrix is",
                           "not saved. [default %default]")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbosity option. [default %default]"))  
)


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "processing.R"))


# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
file1 <- args[["file1"]]
file2 <- args[["file2"]]
row_names <- args[["rownames"]]
SNF_K <- args[["K"]]
SNF_sigma <- args[["sigma"]]
SNF_t <- args[["t"]]
SNF_metric <- args[["metric"]]
nk_max <- args[["nk-max"]]
cluster_file <- args[["cluster-file"]]
affinity_file <- args[["affinity-file"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# Check required arguments
args_req <- c("file1", "file2", "rownames")
for (arg in args_req) {
  if (is.null(args[[arg]])) {
    arg <- paste0("--", arg)
    stop("Argument ", arg, " must be specified.")
  }
}

# Check that input files are CSV
for (file in c("file1", "file2")) {
  ext <- tools::file_ext(args[[file]])
  if (ext != "csv") {
    stop("Input file must be a CSV.")
  }
}

# At least 2 clusters must be specified
if (nk_max < 2) {
  stop("Argument --nk-max must be greater than 1")
}

# Create outdir if needed
outdir <- dirname(cluster_file)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Import matrices
if (verbose) {message("Importing data...")}
x1 <- as_tibble(data.table::fread(file1, header = TRUE))
x2 <- as_tibble(data.table::fread(file2, header = TRUE))

if (nrow(x1) != nrow(x2)) {
  stop("Input matrices have different numbers of rows.")
}

if (nrow(x1) < SNF_K + 1) {
  stop("Argument K must be less than the number of input matrix rows.")
}

if (!is.null(row_names)) {
  if (!(row_names %in% colnames(x1))) {
    stop(paste0("Row names column '", row_names,
                "' not found in input file ", file1))
  } else if (!(row_names %in% colnames(x2))) {
    stop(paste0("Row names column '", row_names,
                "' not found in input file ", file2))
  } else {
    x1 <- column_to_rownames(x1, var = row_names)
    x2 <- column_to_rownames(x2, var = row_names)
    if (sum(rownames(x1) != rownames(x2)) != 0) {
      stop("Input matrix row names don't match.")
    }
  }
}

colnames(x1) <- NULL
colnames(x2) <- NULL

numeric_test_1 <- any(!purrr::map_lgl(.x = x1, .f = is.numeric))
if (numeric_test_1) {
  stop(paste0("Input data contains non-numeric columns: ", file1, ".\n",
              "A single column containing row names can be specified using ",
              "the rownames argument."))
}

numeric_test_2 <- any(!purrr::map_lgl(.x = x2, .f = is.numeric))
if (numeric_test_2) {
  stop(paste0("Input data contains non-numeric columns: ", file2, ".\n",
              "A single column containing row names can be specified using ",
              "the rownames argument."))
}

x1 <- as.matrix(x1)
x2 <- as.matrix(x2)

# Run similarity network fusion
if (verbose) {message("Running similarity network fusion...")}
W <- similarity_network(x1 = x1, x2 = x2,
                        K = SNF_K, sigma = SNF_sigma,
                        t = SNF_t, metric = SNF_metric,
                        outfile = affinity_file)

# Identify clusters
if (verbose) {message("Identifying clusters...")}
clusters <- create_clusters(W = W, nk = nk_max,
                            outfile = cluster_file)
