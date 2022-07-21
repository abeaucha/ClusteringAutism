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
  make_option('--clusterfile',
              type = 'character',
              help = "Path to .csv file containing cluster assignment data."),
  make_option('--imgdir',
              type = 'character',
              help = "Path to directory containing images to use."),
  make_option('--method',
              type = 'character',
              default = 'mean',
              help = paste("Method used to create the representative cluster",
                           "maps. [default %default]")),
  make_option('--jacobians',
              type = 'character',
              help = paste("Optional flag to indicate type of jacobians",
                           "used in outfile.")),
  make_option('--outdir',
              type = 'character',
              help = paste("Path to output directory.")),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))
) 


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
clusterfile <- args[['clusterfile']]
imgdir <- args[['imgdir']]
method <- args[['method']]
outdir <- args[['outdir']]
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

#Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Import cluster information
if (verbose) {message("Importing cluster information...")}

df_clusters <- data.table::fread(clusterfile, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames('ID')

#Create cluster maps
if (verbose) {message("Creating cluster maps...")}

sink(file = 'tmp.log', type = 'output')
imgfiles <- list.files(imgdir, full.names = T)
for (j in 1:ncol(df_clusters)) {
  
  krange <- sort(unique(df_clusters[,j]))
  
  for (k in krange) {
    
    if (verbose) {
      message(paste('nk =', max(krange), ':', 'k =', k))
    }
    
    rows_k <- df_clusters[,j] == k
    id_k <- rownames(df_clusters)[rows_k]
    
    files_k <- character(length(id_k))
    for (f in 1:length(id_k)) {
      files_k[f] <- str_subset(imgfiles, id_k[f])
      
    }
    
    cluster_map <- mincSummary(filenames = files_k,
                               method = method,
                               grouping = NULL)
    
    if (method == 'mean') {
      cluster_map <- cluster_map[,1]
    } else if (method == 'median') {
      cluster_map <- cluster_map
    } else {
      stop()
    }
    
    class(cluster_map) <- class(mincGetVolume(files_k[1]))
    attributes(cluster_map) <- attributes(mincGetVolume(files_k[1]))
    
    res <- max(minc.separation.sizes(files_k[1]))
    
    outfile <- paste0('Group_', k,
                      '_Clusternum_', max(krange), '_ES')
    
    if (!is.null(args[['jacobians']])) {
      outfile <- paste0(outfile, '_', args[['jacobians']])
    }
    
    outfile <- paste0(outfile, '_', res,
                      '_', method, '.mnc')
    
    outfile <- file.path(outdir, outfile)
    
    mincWriteVolume(cluster_map,
                    output.filename = outfile,
                    clobber = TRUE)
  }
}
sink(file = NULL)
system('rm tmp.log')