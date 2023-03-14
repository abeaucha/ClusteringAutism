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
  make_option("--method",
              type = "character",
              default = "mean",
              help = paste("Method used to create the representative cluster",
                           "maps. [default %default]")),
  make_option("--outdir",
              type = "character",
              help = paste("Path to output directory.")),
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
method <- args[["method"]]
outdir <- args[["outdir"]]
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
if (verbose) {message("Creating cluster maps...")}

sink(file = "tmp.log", type = "output")
imgfiles <- list.files(imgdir, full.names = T)
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
    
    if (method == "mean") {
      cluster_map <- mincMean(filenames = files_k,
                              mask = mask)
      cluster_map <- cluster_map[,1]
    } else if (method == "median") {
      cluster_map <- mincApply(filenames = files_k,
                               function.string = quote(median(x)),
                               mask = mask)
    } else {
      stop("Argument method must be one of {mean, median}.")
    }
    
    class(cluster_map) <- class(mincGetVolume(files_k[1]))
    attributes(cluster_map) <- attributes(mincGetVolume(files_k[1]))
    
    outfile <- paste0("cluster_map_nk_", max(krange), "_k_", k, ".mnc")
    outfile <- file.path(outdir, outfile)
    
    mincWriteVolume(cluster_map,
                    output.filename = outfile,
                    clobber = TRUE)
  }
}
sink(file = NULL)
system("rm tmp.log")