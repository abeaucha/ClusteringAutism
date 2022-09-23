# ----------------------------------------------------------------------------
# human_cluster_signatures.R
# Author: Antoine Beauchamp
# Created: May 18th, 2022
#
# Create gene expression signatures for human imaging clusters.

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(tcltk))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--cluster-dir',
              type = 'character',
              help = "Path to directory containing cluster mask images."),
  make_option('--expr-dir',
              type = 'character',
              help = paste("Path to directory containing gene expression",
                           "data sets.")),
  make_option('--coordinates',
              type = 'character',
              help = ""),
  make_option('--template',
              type = 'character',
              help = "Path to MINC file containing human imaging template."),
  make_option('--sign',
              type = 'character',
              help = ""),
  make_option('--outfile',
              type = 'character',
              default = 'cluster_signatures.csv',
              help = paste("Path to the .csv file in which to export cluster",
                           "signatures. [default %default]")),
  make_option('--parallel',
              type = 'character',
              default = 'false',
              help = "Option to run in parallel. [default %default]"),
  make_option('--nproc',
              type = 'numeric',
              help = paste("Number of processors to use in parallel.",
                           "Ignored if --parallel is false.")),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))
)


# Functions ------------------------------------------------------------------


#' Convert a set of world coordinates to voxel coordinates
#'
#' @param coords (character scalar)
#' @param template (character scalar) 
#'
#' @return (matrix)
world_to_voxel <- function(coords, template) {
  
  coords <- suppressMessages(read_csv(coords)) %>% 
    select(label, x, y, z) %>% 
    column_to_rownames('label') %>% 
    as.matrix()
  
  coords_voxel <- mincConvertWorldMatrix(world_matrix = t(coords),
                                         file = template,
                                         nearest_voxel = TRUE)
  coords_voxel <- t(coords_voxel)
  colnames(coords_voxel) <- c('x', 'y', 'z')
  rownames(coords_voxel) <- rownames(coords)
  
  return(coords_voxel)
  
}

#'  Create an expression signature for a cluster
#'
#' @param infile (data.frame) A data.frame row with two columns 
#' containing 1. the path to the cluster mask image, and 2. the path to 
#' the gene expression .csv file.
#' @param sample_coordinates (matrix) AHBA microarray sample 
#' voxel coordinates.
#' @param sign (character scalar) One of {'positive', 'negative'} 
#' indicating whether to use only negative or positive mask values. 
#' All values are used if NULL.  
#'
#' @return (data.frame) A data.frame row containing the expression 
#' signature.
create_cluster_signature <- function(infiles, sample_coordinates, sign = NULL) {
  
  cluster <- mincGetVolume(infiles[[1, 1]])
  cluster <- floor(cluster)
  if (!is.null(sign)) {
    if (sign == 'positive') {
      cluster[cluster < 0] <- 0
    } else if (sign == 'negative') {
      cluster[cluster > 0] <- 0
    } else {
      stop()
    }
  }
  cluster <- abs(cluster)
  cluster <- mincArray(cluster)
  
  cluster_voxels <- which(cluster == 1, arr.ind = TRUE)
  
  cluster_coordinates <- cluster_voxels %>% 
    as_tibble() %>% 
    unite(coords, 1:3, sep = '-') %>% 
    pull(coords)
  
  sample_coordinates <- sample_coordinates[,3:1] %>% 
    as_tibble() %>% 
    unite(coords, 1:3, sep = '-') %>% 
    pull(coords)
  
  samples_in_cluster <- sample_coordinates %in% cluster_coordinates
  
  expression <- data.table::fread(infiles[[1, 2]], header = TRUE) %>% 
    as_tibble()
  
  signature <- expression[samples_in_cluster,] %>% 
    summarise_all(.funs = mean)
  
  return(signature)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
cluster_dir <- args[['cluster-dir']]
expr_dir <- args[['expr-dir']]
coords <- args[['coordinates']]
template <- args[['template']]
outfile <- args[['outfile']]
sign <- args[['sign']]
inparallel <- ifelse(args[['parallel']] == 'true', TRUE, FALSE)
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

#Create outdir if needed
outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

if (verbose) {message("Getting voxel coordinates for microarray samples...")}

sample_coordinates <- world_to_voxel(coords = coords,
                                     template = template)

if (verbose) {message("Identifying input files...")}

#Get cluster and expression input files
expr_files <- Sys.glob(file.path(expr_dir, '*.csv'))
cluster_files <- list.files(cluster_dir, full.names = TRUE)
infiles <- expand_grid(cluster_file = cluster_files, 
                       expr_file = expr_files)

if (verbose) {message("Computing cluster signatures...")}

#Compute cluster signatures
pb <- txtProgressBar(max = nrow(infiles), style = 3)
progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
if (inparallel) {
  nproc <- args[['nproc']]
  cl <- makeSOCKcluster(nproc)
  registerDoSNOW(cl)
  opts <- list(progress=progress)
  signatures <- foreach(i = 1:nrow(infiles),
                        .packages = c('tidyverse', 'RMINC'),
                        .combine = 'bind_rows', .options.snow=opts) %dopar% {
                          create_cluster_signature(infile = infiles[i,],
                                                   sample_coordinates = sample_coordinates,
                                                   sign = sign)
                        }
  close(pb)
  stopCluster(cl)
} else {
  signatures <- foreach(i = 1:nrow(infiles),
                        .packages = c('tidyverse', 'RMINC'),
                        .combine = 'bind_rows') %do% {
                          progress(n = i)
                          create_cluster_signature(infile = infiles[i,],
                                                   sample_coordinates = sample_coordinates,
                                                   sign = sign)
                        }
  close(pb)
}

#Include file info
signatures <- bind_cols(infiles, signatures)
if (is.null(sign)) {
  signatures$sign <- NA
} else {
  signatures$sign <- sign
}

#Write to file
data.table::fwrite(signatures,
                   file = outfile)
