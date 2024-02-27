#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# create_cluster_centroids.R
# Author: Antoine Beauchamp
# Created: May 18th, 2022
#
# Create cluster centroid images.
#
# Description
# -----------
# This script creates a cluster centroid image for each cluster in the input
# file. The centroid images are computed by aggregating the voxel-wise values
# for all images in a cluster.


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--imgdir",
              type = "character",
              help = paste("Path to the directory containing the individual",
                           "images (.mnc).")),
  make_option("--cluster-file",
              type = "character",
              help = "Path to the file (.csv) containing cluster assignments."),
  make_option("--mask",
              type = "character",
              help = "Path to a mask image (.mnc)."),
  make_option("--outdir",
              type = "character",
              help = "Path to the output directory."),
  make_option("--method",
              type = "character",
              default = "mean",
              help = paste("One of {mean, median} specifying how to",
                           "compute the centroids.",
                           "[default %default]")),
  make_option("--execution",
              type = "character",
              default = "local",
              help = paste("[default %default]")),
  make_option("--nproc",
              type = "numeric",
              default = 1,
              help = paste("Number of processors to use in parallel.",
                           "Executed serially if 1.",
                           "[default %default]")),
  make_option("--slurm-njobs",
              type = "numeric",
              help = "Number of jobs to deploy on Slurm."),
  make_option("--slurm-mem",
              type = "character",
              help = paste("Memory per CPU core")),
  make_option("--slurm-time",
              type = "numeric",
              help = paste("Walltime in minutes")),
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

# TODO remove lines when script works
# REMOVE THESE LINES WHEN FINISHED
args[["imgdir"]] <- "data/test/human/derivatives/v2/547/effect_sizes/resolution_3.0/absolute/"
args[["cluster-file"]] <- "data/test/human/derivatives/v2/547/clusters/resolution_3.0/clusters.csv"
args[["mask"]] <- "data/human/registration/v2/reference_files/mask_3.0mm.mnc"
args[["outdir"]] <- "data/test/human/derivatives/v2/547/centroids/resolution_3.0/absolute/"
args[["method"]] <- "mean"
args[["execution"]] <- "local"
args[["nproc"]] <- 8

imgdir <- args[["imgdir"]]
clusterfile <- args[["cluster-file"]]
mask <- args[["mask"]]
outdir <- args[["outdir"]]
method <- args[["method"]]
execution <- args[["execution"]]
nproc <- args[["nproc"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# Check required arguments
args_req <- c("cluster-file", "imgdir", "mask", "outdir")
for (arg in args_req) {
  if (is.null(args[[arg]])) {
    arg <- paste0("--", arg)
    stop("Argument ", arg, " must be specified.")
  }
}

# Check that cluster file is a CSV
if (tools::file_ext(clusterfile) != "csv") {
  stop(paste("Clusters file must be a CSV:",
             clusterfile))
}

# Check that imgdir contains images
imgfiles <- list.files(imgdir, full.names = TRUE, pattern = "*.mnc")
if (length(imgfiles) == 0) {
  stop(paste("No MINC files were found in the image directory:",
             imgdir))
}

# Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Check execution option
if (execution == "local") {
  resources <- list() 
} else if (execution == "slurm") {
  njobs <- args[["slurm-njobs"]]
  resources <- list(memory = args[["slurm-mem"]],
                    walltime = args[["slurm-time"]]*60)
} else {
  stop()
}

# Import cluster information
if (verbose) {message("Importing cluster information...")}
df_clusters <- as_tibble(data.table::fread(clusterfile, header = TRUE)) %>% 
  rename(file = ID) %>% 
  mutate(file = file.path(imgdir, file))


# Create centroid images for all clusters
if (verbose) {message("Creating centroid images...")}
sink(nullfile(), type = "output")
cols_iter <- colnames(df_clusters)[colnames(df_clusters) != "file"]
for (col in cols_iter) {
  clusters <- df_clusters[, c("file", col)]
  colnames(clusters) <- c("file", "k")
  centroids <- compute_cluster_centroids(clusters = clusters,
                                         mask = mask,
                                         outdir = outdir,
                                         method = method,
                                         execution = execution,
                                         nproc = nproc,
                                         njobs = njobs,
                                         resources = resources)
}
sink(NULL)
