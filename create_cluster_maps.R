# ----------------------------------------------------------------------------
# create_cluster_maps.R
# Author: Antoine Beauchamp
# Created: May 18th, 2022
#
# Description
# -----------
#

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RMINC))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--clusterfile',
              type = 'character',
              help = "Path to CSV file containing cluster data."),
  make_option('--imgdir',
              type = 'character',
              help = ""),
  make_option('--method',
              type = 'character',
              default = 'mean',
              help = "[default %default]"),
  make_option('--outdir',
              type = 'character',
              help = ""),
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

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

if (verbose) {message("Importing cluster information...")}

df_clusters <- as.data.frame(data.table::fread(clusterfile, header = TRUE))
rownames(df_clusters) <- df_clusters[,'ID']
df_clusters <- df_clusters[,colnames(df_clusters) != 'ID']

if (verbose) {message("Creating cluster maps...")}

sink(file = 'tmp.log', type = 'output')
for (j in 1:ncol(df_clusters)) {
  
  krange <- sort(unique(df_clusters[,j]))
  
  for (k in krange) {
    
    if (verbose) {
      message(paste('nk =', max(krange), ':', 'k =', k))
    }
    
    rows_k <- df_clusters[,j] == k
    files_k <- rownames(df_clusters)[rows_k]
    files_k <- file.path(imgdir, files_k)
  
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
    
    outfile <- paste0('Group_', k,
                      '_Clusternum_', max(krange), 
                      '_ES_', method, '.mnc')
    outfile <- file.path(outdir, outfile)
    
    mincWriteVolume(cluster_map,
                    output.filename = outfile,
                    clobber = TRUE)
  }
}
sink(file = NULL)
system('rm tmp.log')