# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Functions ------------------------------------------------------------------

#' Title
#'
#' @param metadata 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
fetch_params_metadata <- function(metadata, ...) {
  
  kwargs <- list(...)
  
  if (!file.exists(metadata)) {
    stop(paste("Metadata file not found:", metadata))
  }
  
  df_metadata <- read_csv(metadata, show_col_types = FALSE)
  
  if (length(kwargs) == 0) {
    message("No parameters provided. Returning all metadata.")
    return(df_metadata)
  } else {
    df_params = as_tibble(kwargs)
  }
  
  df_match <- inner_join(df_metadata, df_params, by = colnames(df_params))
  nmatch <- nrow(df_match)
  if (nmatch == 0) {
    return(NULL)
  } else {
    return(df_match)
  }
}


#' Resample a MINC image
#'
#' @param infile (character scalar) Path to the image file to resample.
#' @param isostep (numeric scalar) Isotropic dimension of voxels in 
#' the resampled image (mm).
#' @param outdir (character scalar) Path to the directory in which 
#' to save the resampled image. If None, resampled image will be 
#' written to the directory of the input file.
#' @param suffix (character scalar) Suffix to append to output file name 
#' before the file extension.
#'
#' @return (character scalar) Path to the resampled image.
resample_image <- function(infile, isostep, outdir = NULL, suffix = NULL){
  return(outfile)
}

