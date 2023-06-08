# ----------------------------------------------------------------------------
# fetch_microarray_coordinates.R
# Author: Antoine Beauchamp
# Created: June 23rd, 2022
#
# Brief description
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--metadata",
              type = "character",
              help = "Path to .csv file containing AHBA sample metadata."),
  make_option("--outfile",
              type = "character",
              help = "Path to the .csv file in which to save coordinates."),
  make_option("--labels",
              type = "character",
              default = "true",
              help = "Option to add labels."),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbosity option. [default %default]"))
) 


# Functions ------------------------------------------------------------------

#' Get microarray sample coordinates
#'
#' @param metadata (character scalar) Path to the .csv file containing 
#' AHBA sample metadata.
#'
#' @return (tibble) 
fetch_microarray_coordinates <- function(metadata) {
  
  metadata <- suppressMessages(read_csv(metadata))
  donors <- unique(metadata[["Donor"]])
  for (i in 1:length(donors)) {
    
    url <- paste0("https://raw.githubusercontent.com/gdevenyi/",
                  "AllenHumanGeneMNI/master/transformed-points/recombine/", 
                  donors[i], 
                  "_SampleAnnot.csv")
    
    coords_tmp <- suppressMessages(read_csv(url)) %>% 
      unite(label, 
            structure_id, slab_num, well_id, 
            sep = "-", remove = TRUE) %>%
      mutate(t = 0) %>% 
      select(x = mni_nlin_x,
             y = mni_nlin_y,
             z = mni_nlin_z,
             t,
             label)
      
    if (i == 1) {
      coords <- coords_tmp
    } else {
      coords <- rbind(coords,
                      coords_tmp)
    }
  }
  
  return(coords)
  
}

#' Make labels for microarray samples
#'
#' @param coords (tibble)
#'
#' @return (list)
make_labels <- function(coords) {
  
  coords <- coords %>%
    mutate(sample_id = label,
           label = 1:nrow(.))
  
  defs <- coords %>% 
    select(sample_id, label)
  
  coords <- coords %>% 
    select(-sample_id)
  
  return(list(coords = coords,
              defs = defs))
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
metadata <- args[["metadata"]]
outfile <- args[["outfile"]]
labels <- ifelse(args[["labels"]] == "true", TRUE, FALSE)
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

if (verbose) {message("Downloading AHBA microarray sample coordinates...")}

coords <- fetch_microarray_coordinates(metadata = metadata)

if (labels) {
  label_list <- make_labels(coords = coords)
  coords <- label_list[["coords"]]
  defs <- label_list[["defs"]]
  write_csv(x = defs, file = str_replace(outfile, ".csv", "_defs.csv"))
} else {
  coords <- coords %>% 
    mutate(label = 0)
}

write_csv(x = coords, file = outfile)
