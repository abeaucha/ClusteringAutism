# ----------------------------------------------------------------------------
# cluster_human_data.R
# Author: Antoine Beauchamp, Jacob Ellegood
# Created: May 17th, 2022
#
# Description
# -----------
# 


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SNFtool))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--file1',
              type = 'character',
              help = "Path to .csv file containing first effect size matrix."),
  make_option('--file2',
              type = 'character',
              help = "Path to .csv file containing second effect size matrix."),
  make_option('--rownames',
              type = 'character',
              help = "Name of column in .csv file containing row names."),
  make_option('--nk-max',
              type = 'numeric',
              default = 10,
              help = paste("Maximum number of clusters to identify.", 
                           "[default %default]")),
  make_option('--metric',
              type = 'character',
              default = 'correlation',
              help = paste("Distance metric used to compute the SNF affinity",
                           "matrices. [default %default]")),
  make_option('--K',
              type = 'numeric',
              default = 10,
              help = paste("Number of nearest-neighbours used to compute",
                           "the SNF affinity matrices. [default %default]")),
  make_option('--sigma',
              type = 'numeric',
              default = 0.5,
              help = paste("Variance for the local model in the SNF",
                           "affinity matrices. [default %default]")),
  make_option('--t',
              type = 'numeric',
              default = 20,
              help = paste("Number of iterations for the diffusion",
                           "process in SNF. [default %default]")),
  make_option('--cluster-file',
              type = 'character',
              help = paste("Path to .csv file in which to save the cluster", 
                           "assignments.")),
  make_option('--affinity-file',
              type = 'character',
              help = paste("Path to .csv file in which to save the SNF",
                           "affinity matrix [default %default]")),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))  
)


# Functions ------------------------------------------------------------------

#' Run similarity network fusion (SNF)
#'
#' @param x1 (matrix) Input matrix
#' @param x2 (matrix) Input matrix
#' @param metric (character scalar) Distance metric used to compute the 
#' affinity matrices.
#' @param K (numeric scalar) Number of nearest-neighbours used to 
#' compute the SNF affinity matrices. 
#' @param sigma (numeric scalar) Variance for the local model in the 
#' SNF affinity matrices.
#' @param t (numeric scalar) Number of iterations for the diffusion
#' process in SNF.
#' @param outfile (character scalar) Path to file in which to save
#' affinity matrix.
#'
#' @return (matrix) SNF affinity matrix.
SNF_combine <- function(x1, x2, metric = "correlation", K = 10, 
                        sigma = 0.5, t = 20, outfile = NULL){
  
  if (metric == "correlation") { 
    d1 <- (1-cor(t(x1)))
    d2 <- (1-cor(t(x2)))
  } else if (metric == "euclidean") {
    d1 <- (dist2(as.matrix(x1), as.matrix(x1)))^(1/2)
    d2 <- (dist2(as.matrix(x2), as.matrix(x2)))^(1/2)
  } else {
    stop(paste("Argument metric must be one of {correlation, euclidean}:", metric))
  }
  
  W1 <- affinityMatrix(d1, K = K, sigma = sigma)
  W2 <- affinityMatrix(d2, K = K, sigma = sigma)
  
  W <- SNF(list(W1, W2), K = K, t = t)
  
  if (!is.null(outfile)){
    save(W, file = outfile)
  }
  
  return(W)
  
}

#' Create clusters from SNF affinity matrix
#'
#' @param W (matrix) SNF affinity matrix.
#' @param nk (numeric scalar) Maximum number of clusters to use in 
#' clustering.
#' @param outfile (character scalar) Optional path to .csv file in 
#' which to save cluster assignments.
#'
#' @return (data.frame) Cluster assignments.
create_clusters <- function(W, nk = 10, outfile = NULL) {
  
  for(k in 2:nk) {
    group <- spectralClustering(affinity = W, K = k)
    group_name <- paste0('nk', k)
    assign(group_name, group)
    if (k == 2) {
      if (is.null(rownames(W))) {
        ids <- as.character(1:nrow(W))
      } else {
        ids <- rownames(W)
      }
      all_clusters <- data.frame(ids, group, stringsAsFactors = F)
      colnames(all_clusters) <- c('ID', group_name)
    } else {
      group <- data.frame(group)
      colnames(group) <- group_name
      all_clusters <- cbind(all_clusters, group)
    }
  }
  
  if (!is.null(outfile)) {
    write.csv(x = all_clusters, file = outfile, row.names = FALSE)
  }
  
  return(all_clusters)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
file1 <- args[['file1']]
file2 <- args[['file2']]
row_names <- args[['rownames']]
SNF_K <- args[['K']]
SNF_sigma <- args[['sigma']]
SNF_t <- args[['t']]
SNF_metric <- args[['metric']]
nk_max <- args[['nk-max']]
cluster_file <- args[['cluster-file']]
affinity_file <- args[['affinity-file']]
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

if (nk_max < 2) {
  stop("Argument nk-max must be greater than 1")
}

#Create outdir if needed
outdir <- dirname(cluster_file)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

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
    stop(paste0("Row names column '", row_names, "' not found in input file ", file1))
  } else if (!(row_names %in% colnames(x2))) {
    stop(paste0("Row names column '", row_names, "' not found in input file ", file2))
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
              "A single column containing row names can be specified using the rownames argument."))
}

numeric_test_2 <- any(!purrr::map_lgl(.x = x2, .f = is.numeric))
if (numeric_test_2) {
  stop(paste0("Input data contains non-numeric columns: ", file2, ".\n",
              "A single column containing row names can be specified using the rownames argument."))
}

x1 <- as.matrix(x1)
x2 <- as.matrix(x2)

if (verbose) {message("Running similarity network fusion...")}

W <- SNF_combine(x1 = x1, x2 = x2,
                 K = SNF_K, sigma = SNF_sigma,
                 t = SNF_t, metric = SNF_metric,
                 outfile = affinity_file)

if (verbose) {message("Assigning clusters...")}

clusters <- create_clusters(W = W, nk = nk_max,
                            outfile = cluster_file)
