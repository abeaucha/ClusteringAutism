#!Rscript
# ----------------------------------------------------------------------------
# compute_cluster_centroid.R
# Author: Antoine Beauchamp
# Created: May 18th, 2022
#
# Compute the centroid image for all clusters in a cluster solution. 
#
# Description
# -----------


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(RMINC))


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
  make_option("--nk",
              type = "numeric",
              default = 2,
              help = paste("Cluster solution for which to compute",
                           "the centroids. [default %default]")),
  make_option("--method",
              type = "character",
              default = "mean",
              help = paste("One of {mean, median} specifying how to",
                           "compute the centroids.",
                           "[default %default]")),
  make_option("--nproc",
              type = "numeric",
              default = 1,
              help = paste("Number of processors to use in parallel.",
                           "Executed serially if 1.",
                           "[default %default]")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbosity option. [default %default]"))
) 


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
imgdir <- args[["imgdir"]]
clusterfile <- args[["cluster-file"]]
mask <- args[["mask"]]
outdir <- args[["outdir"]]
nk <- args[["nk"]]
method <- args[["method"]]
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
  stop("Clusters file must be a CSV.")
}

# Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Import cluster information
if (verbose) {message("Importing cluster information...")}
df_clusters <- data.table::fread(clusterfile, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames("ID")

# Create centroids
if (verbose) {message("Creating centroid images...")}
sink(nullfile(), type = "output")
imgfiles <- list.files(imgdir, full.names = TRUE, pattern = "*.mnc")

# Get the right cluster solution
# Iterate over the clusters in the solution
# Generate the centroid image for each cluster

# for (j in 1:ncol(df_clusters)) {
#   
#   krange <- sort(unique(df_clusters[,j]))
#   for (k in krange) {
#     
#     if (verbose) {
#       message(paste("nk =", max(krange), ":", "k =", k))
#     }
#     
#     rows_k <- df_clusters[,j] == k
#     id_k <- rownames(df_clusters)[rows_k]
#     files_k <- file.path(imgdir, id_k)
#     file_exists <- files_k %in% imgfiles
#     if (sum(!file_exists) != 0) {
#       warning(paste(sum(!(file_exists)), "files not found in directory", 
#                     imgdir, "\nThese files will be ommitted from",
#                     "the calculation."))
#     }
#     
#     # Centroid function
#     if (method == "mean") {
#       centroid <- mean
#     } else if (method == "median") {
#       centroid <- median
#     } else {
#       stop("Argument --method must be one of {mean, median}.")
#     }
#     
#     # Create centroid image
#     cluster_map <- mcMincApply(filenames = files_k,
#                                fun = centroid,
#                                mask = mask,
#                                cores = nproc)
#     
#     # Export image
#     outfile <- paste0("cluster_map_nk_", max(krange), "_k_", k, ".mnc")
#     outfile <- file.path(outdir, outfile)
#     mincWriteVolume(cluster_map,
#                     output.filename = outfile,
#                     clobber = TRUE)
#   }
# }
sink(NULL)
