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
  make_option('--nclusters',
              type = 'numeric',
              default = 10,
              help = paste("Maximum number of clusters to identify.", 
                           "[default %default]")),
  make_option('--metric',
              type = 'character',
              default = 'cor',
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
  make_option('--outfile',
              type = 'character',
              help = paste("Path to .csv file in which to save the cluster", 
                           "assignments.")),
  make_option('--wfile',
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
SNF_combine <- function(x1, x2, metric = 'cor', K = 10, 
                        sigma = 0.5, t = 20, outfile = NULL){
  
  if(metric == 'cor'){
    d1 <- (1-cor(t(x1)))
    d2 <- (1-cor(t(x2)))
  } else {
    d1 <- (dist2(as.matrix(x1), as.matrix(x1)))^(1/2)
    d2 <- (dist2(as.matrix(x2), as.matrix(x2)))^(1/2)
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
    group_name <- paste0('Group', k)
    assign(group_name, group)
    if (k == 2) {
      all_clusters <- data.frame(rownames(W), group, stringsAsFactors = F)
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
outfile <- args[['outfile']]
wfile <- args[['wfile']]
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

if (verbose) {message("Importing data...")}

x1 <- data.table::fread(file1, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames(row_names) %>% 
  as.matrix() 
colnames(x1) <- NULL
rownames(x1) <- basename(rownames(x1))

x2 <- data.table::fread(file2, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames(row_names) %>% 
  as.matrix() 
colnames(x2) <- NULL
rownames(x2) <- basename(rownames(x2))

if (verbose) {message("Running similarity network fusion...")}

W <- SNF_combine(x1 = x1, x2 = x2,
                 K = args[['K']], alpha = args[['alpha']],
                 t = args[['t']], metric = args[['metric']],
                 outfile = wfile)

if (verbose) {message("Assigning clusters...")}

clusters <- create_clusters(W = W, nk = args[['nclusters']], 
                            outfile = outfile)
