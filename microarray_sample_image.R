# ----------------------------------------------------------------------------
# template.R
# Author: Antoine Beauchamp
# Created:
#
# Brief description
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--metadata',
              type = 'character',
              help = "Path to .csv file containing AHBA sample metadata."),
  make_option('--template',
              type = 'character',
              help = "Path to MINC file containing human imaging template."),
  make_option('--outdir',
              type = 'character',
              help = "Path to directory in which to save output files."),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))
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
  
  return(sample_coordinates_voxel)
  
}

#' Title
#'
#' @param coords 
#' @param template 
#'
#' @return
coords_to_img <- function(coords, template) {
  
  template <- template %>% 
    mincGetVolume() %>% 
    mincArray
  
  img <- array(data = 0, dim = dim(template))
  for (s in 1:nrow(coords)){
    i <- coords[s,3]
    j <- coords[s,2]
    k <- coords[s,1]
    img[i,j,k] <- as.integer(s)
  }
  attributes(img) <- attributes(template)
  
  return(img)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
metadata <- args[['metadata']]
template <- args[['template']]
outdir <- args[['outdir']]
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

if (verbose) {message("Downloading AHBA microarray sample coordinates...")}

# metadata <- 'data/human/SampleInformation_pipeline_v1.csv'
# template <- 'data/human/atlas/mni_icbm152_t1_tal_nlin_sym_09c.mnc'

get_sample_coordinates(metadata = metadata,
                       template = template)

sample_coordinates <- get_sample_coordinates(metadata = metadata,
                                             template = template)

if (verbose) {message("Mapping AHBA samples to image space...")}

sample_labels <- coords_to_img(coords = sample_coordinates,
                               template = template)

defs <- tibble(sample_id = rownames(sample_coordinates),
               label = 1:length(sample_id))

if (verbose) {message("Writing sample image to file...")}

outfile_img <- file.path(outdir, 'AHBA_microarray_samples.mnc')
mincWriteVolume(buffer = sample_labels,
                output.filename = outfile_img,
                clobber = TRUE)

if (verbose) {message("Writing sample definitions to file...")}

outfile_defs <- file.path(outdir, 'AHBA_microarray_samples_defs.csv')
write_csv(x = defs, file = outfile_defs)

quit()

