# ----------------------------------------------------------------------------
# create_cluster_maps.R
# Author: Antoine Beauchamp
# Created: May 18th, 2022
#
# Create representative voxel-wise maps for clustered images.
#
# Description
# -----------
# This script creates a representative voxel-wise map for each cluster. 
# The representative cluster maps are computed by aggregating the voxel-wise
# values for all images in a cluster. 


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--cluster-file",
              type = "character",
              help = "Path to .csv file containing cluster assignment data."),
  make_option("--imgdir",
              type = "character",
              help = "Path to directory containing images to use."),
  make_option("--mask",
              type = "character",
              help = "Path to the mask file to use."),
  make_option("--outdir",
              type = "character",
              help = paste("Path to output directory.")),
  make_option("--method",
              type = "character",
              default = "mean",
              help = paste("Method used to create the representative cluster",
                           "maps. [default %default]")),
  make_option("--nproc",
              type = "numeric",
              default = 2,
              help = paste("Number of processors to use in parallel.")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbosity option. [default %default]"))
) 


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
clusterfile <- args[["cluster-file"]]
imgdir <- args[["imgdir"]]
mask <- args[["mask"]]
outdir <- args[["outdir"]]
method <- args[["method"]]
nproc <- args[["nproc"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

#Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Import cluster information
if (verbose) {message("Importing cluster information...")}
df_clusters <- data.table::fread(clusterfile, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames("ID")

#Create cluster maps
if (verbose) {message("Creating centroid images...")}
sink(nullfile(), type = "output")
imgfiles <- list.files(imgdir, full.names = TRUE, pattern = "*.mnc")
for (j in 1:ncol(df_clusters)) {
  
  krange <- sort(unique(df_clusters[,j]))
  for (k in krange) {
    
    if (verbose) {
      message(paste("nk =", max(krange), ":", "k =", k))
    }
    
    rows_k <- df_clusters[,j] == k
    id_k <- rownames(df_clusters)[rows_k]
    files_k <- file.path(imgdir, id_k)
    file_exists <- files_k %in% imgfiles
    if (sum(!file_exists) != 0) {
      warning(paste(sum(!(file_exists)), "files not found in directory", 
                    imgdir, "\nThese files will be ommitted from",
                    "the calculation."))
    }
    
    #Centroid function
    if (method == "mean") {
      centroid <- mean
    } else if (method == "median") {
      centroid <- median
    } else {
      stop("Argument --method must be one of {mean, median}.")
    }
    
    #Create centroid image
    cluster_map <- mcMincApply(filenames = files_k,
                               fun = centroid,
                               mask = mask,
                               cores = nproc)
    
    #Export image
    outfile <- paste0("cluster_map_nk_", max(krange), "_k_", k, ".mnc")
    outfile <- file.path(outdir, outfile)
    mincWriteVolume(cluster_map,
                    output.filename = outfile,
                    clobber = TRUE)
  }
}
sink(NULL)
