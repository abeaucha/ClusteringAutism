# ----------------------------------------------------------------------------
# coordinates_to_minc.R
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
  make_option('--coordinates',
              type = 'character',
              help = "Path to .csv file containing AHBA microarray coordinates."),
  make_option('--template',
              type = 'character',
              help = "Path to MINC file containing human imaging template."),
  make_option('--outfile',
              type = 'character',
              help = "Path to .mnc file in which to save output image."),
  make_option('--type',
              type = 'character',
              default = 'labels',
              help = "[default %default]"),
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


#' Convert voxel coordinates to MINC image
#'
#' @param coords (matrix)
#' @param template (character scalar)
#' @param type (character scalar)
#'
#' @return
coords_to_img <- function(coords, template, type = 'labels') {
  
  template <- template %>% 
    mincGetVolume() %>% 
    mincArray()
  
  img <- array(data = 0, dim = dim(template))
  for (s in 1:nrow(coords)){
    i <- coords[s,3]
    j <- coords[s,2]
    k <- coords[s,1]
    if (type == 'mask') {
      img[i,j,k] <- 1
    } else if (type == 'labels') {
      img[i,j,k] <- as.integer(rownames(coords)[s])
    } else {
      stop()
    }
  }
  attributes(img) <- attributes(template)
  
  return(img)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
coords <- args[['coordinates']]
template <- args[['template']]
outfile <- args[['outfile']]
type <- args[['type']]
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

if (verbose) {message("Converting world coordinates to voxel...")}

voxel_coordinates <- world_to_voxel(coords = coords,
                                    template = template)

if (verbose) {message("Creating image from coordinates...")}

coords_img <- coords_to_img(coords = voxel_coordinates,
                            template = template,
                            type = type)

if (verbose) {message("Writing image to file...")}

mincWriteVolume(buffer = coords_img,
                output.filename = outfile,
                clobber = TRUE)
