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
suppressPackageStartupMessages(library(doParallel))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--clusterdir',
              type = 'character',
              help = "Path to directory containing cluster mask images."),
  make_option('--exprdir',
              type = 'character',
              help = paste("Path to directory containing gene expression",
                           "data sets.")),
  make_option('--metadata',
              type = 'character',
              help = "Path to .csv file containing AHBA sample metadata."),
  make_option('--template',
              type = 'character',
              help = "Path to MINC file containing human imaging template."),
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
                           "Ignored if --parallel is false."))
)


# Functions ------------------------------------------------------------------

#' Get microarray sample coordinates
#'
#' @param metadata (character scalar) Path to the .csv file containing AHBA
#' sample metadata.
#' @param template (character scalar) Path to the .mnc file containing human
#' imaging template.
#'
#' @return (character vector) AHBA microarray sample voxel coordinates in the 
#' form "x-y-z".
get_sample_coordinates <- function(metadata, template) {
  
  metadata <- suppressMessages(read_csv(metadata))
  donors <- unique(metadata[['Donor']])
  for (i in 1:length(donors)) {
    
    url <- paste0("https://raw.githubusercontent.com/gdevenyi/",
                  "AllenHumanGeneMNI/master/transformed-points/recombine/", 
                  donors[i], 
                  "_SampleAnnot.csv")
    
    df_tmp <- suppressMessages(read_csv(url)) %>% 
      mutate(Donor = donors[i]) %>% 
      unite(SampleID, 
            structure_id, slab_num, well_id, 
            sep = '-', remove = FALSE) %>% 
      column_to_rownames('SampleID') %>% 
      select(contains('mni_nlin')) %>% 
      as.matrix()
    
    if (i == 1) {
      sample_coordinates_world <- df_tmp
    } else {
      sample_coordinates_world <- rbind(sample_coordinates_world,
                                        df_tmp)
    }
  }
  
  ind_match_metadata <- match(rownames(sample_coordinates_world),
                              metadata[['SampleID']])
  sample_coordinates_world <- sample_coordinates_world[ind_match_metadata,]
  
  sample_coordinates_voxel <- mincConvertWorldMatrix(world_matrix = t(sample_coordinates_world),
                                                     file = template,
                                                     nearest_voxel = TRUE)
  sample_coordinates_voxel <- t(sample_coordinates_voxel)
  colnames(sample_coordinates_voxel) <- c('x', 'y', 'z')
  
  sample_coordinates <- sample_coordinates_voxel %>% 
    as_tibble() %>% 
    unite(coords, x, y, z, sep = '-') %>% 
    pull(coords)
  
  return(sample_coordinates)
  
}


#'  Create an expression signature for a cluster
#'
#' @param infile (data.frame) A data.frame row with two columns 
#' containing 1. the path to the cluster mask image, and 2. the path to 
#' the gene expression .csv file.
#' @param sample_coordinates (character vector) AHBA microarray sample 
#' voxel coordinates in the form "x-y-z".
#'
#' @return (data.frame) A data.frame row containing the expression 
#' signature.
create_cluster_signature <- function(infile, sample_coordinates) {
  
  cluster <- mincArray(mincGetVolume(infile[,1][[1]]))
  
  cluster_voxels <- which(cluster == 1, arr.ind = TRUE)
  
  cluster_coordinates <- cluster_voxels %>% 
    as_tibble() %>% 
    unite(coords, 1:3, sep = '-') %>% 
    pull(coords)
  
  expression <- data.table::fread(infile[, 2][[1]], header = TRUE) %>% 
    as_tibble()
  
  ind <- sample_coordinates %in% cluster_coordinates
  
  signature <- expression[ind,] %>% 
    summarise_all(.funs = mean)
  
  return(signature)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
cluster_dir <- args[['clusterdir']]
expr_dir <- args[['exprdir']]
metadata <- args[['metadata']]
template <- args[['template']]
inparallel <- ifelse(args[['parallel']] == 'true', TRUE, FALSE)

#Get AHBA microarray sample coordinates
sample_coordinates <- get_sample_coordinates(metadata = metadata, 
                                             template = template)

#Get cluster and expression input files
expr_files <- Sys.glob(file.path(expr_dir, '*.csv'))

cluster_files <- list.files(cluster_dir, full.names = TRUE)
infiles <- expand_grid(clusterfile = cluster_files, 
                       exprfile = expr_files)

#Option to run in parallel
if (inparallel) {
  nproc <- args[['nproc']]
  cl <- makeCluster(nproc)
  registerDoParallel(cl)
  signatures <- foreach(i = 1:nrow(infiles),
                        .packages = c('tidyverse', 'RMINC'),
                        .combine = 'bind_rows') 
  %dopar% {create_cluster_signature(infile = infiles[i,],
                                    sample_coordinate = sample_coordinates)
  }
  stopCluster(cl)
} else {
  signatures <- foreach(i = 1:nrow(infiles),
                        .packages = c('tidyverse', 'RMINC'),
                        .combine = 'bind_rows') 
  %do% {create_cluster_signature(infile = infiles[i,],
                                 sample_coordinate = sample_coordinates)
    }
}

#Write to file
data.table::fwrite(signatures,
                   file = args[['outfile']])
