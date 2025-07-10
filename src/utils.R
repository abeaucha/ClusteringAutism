# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Functions ------------------------------------------------------------------

#' Fetch pipeline metadata for a set of pipeline parameters
#'
#' @param metadata (character scalar) Path to the metadata file (.csv)
#' @param params (list) List of named parameters
#'
#' @return (tbl, data.frame) Data frame containing parameter set metadata
fetch_params_metadata <- function(metadata, params = NULL) {
  
  # Check existence of metadata file
  if (!file.exists(metadata)) {
    stop(paste("Metadata file not found:", metadata))
  }
  
  # Import metadata
  df_metadata <- read_csv(metadata, show_col_types = FALSE) %>% 
    mutate(id = as.character(id))
  
  # If no parameters provided, return all metadata
  if (is.null(params) | length(params) == 0) {
    message("No parameter list provided. Returning all metadata.")
    return(df_metadata)
  } else {
    df_params = as_tibble(params)
  }
  
  df_match <- inner_join(df_metadata, df_params, by = colnames(df_params))
  return(df_match)
}


#' Fetch ID for a set of pipeline parameters
#'
#' @param metadata (character scalar) Path to the metadata file (.csv)
#' @param params (list) List of named parameters
#'
#' @return (character vector or NULL) Pipeline parameter IDs
fetch_params_id <- function(metadata, params = NULL) {
  df_metadata <- fetch_params_metadata(metadata = metadata, params = params)
  if (nrow(df_metadata) > 0) {
    return(df_metadata[["id"]])
  } else {
    return(NULL)
  }
}


#' Set the unique ID for a set of parameters
#'
#' @param metadata (character scalar) Path to the metadata file (.csv)
#' @param params (list) List of named parameters
#' @param params_id (character scalar or NULL) Optional ID to be specified
#'
#' @returns (character scalar)
set_params_id <- function(metadata, params, params_id = NULL) {
  
  df_params <- as_tibble(params)
  if (file.exists(metadata)) {
    df_metadata <- read_csv(metadata, show_col_types = FALSE) %>% 
      mutate(id = as.character(id))
    df_match <- inner_join(df_metadata, df_params, by = colnames(df_params))
    nmatch <- nrow(df_match)
    if (nmatch == 1) {
      message(paste("Parameters already identified:", df_match[["id"]]))
      df_params[["id"]] <- df_match[["id"]]
    } else if (nmatch == 0) {
      df_params[["id"]] <- ifelse(is.null(params_id), sample(100:999, 1), params_id)
      df_params[["id"]] <- as.character(df_params[["id"]])
      df_metadata <- bind_rows(df_metadata, df_params)
      write_csv(df_metadata, file = metadata)
    } else {
      stop("Multiple IDs have been identified for the parameters")
    }
  } else {
    df_params[["id"]] <- ifelse(is.null(params_id), sample(100:999, 1), params_id)
    df_params[["id"]] <- as.character(df_params[["id"]])
    write_csv(df_params, file = metadata)
  }
  
  return(df_params[["id"]])
  
}


#' Create unique directory for a set of pipeline parameters
#'
#' @param params (list)
#' @param outdir (character scalar) 
#' @param metadata (character scalar) Metadata file base name.
#' @param params_id (character scalar or NULL) 
#'
#' @returns (character scalar) Path to the new directory
mkdir_from_params <- function(params, outdir, metadata = "metadata.csv", 
                              params_id = NULL) {
  
  if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
  
  metadata <- file.path(outdir, metadata)
  
  # Check existence of metadata file
  params_id_check <- try(fetch_params_id(metadata = metadata, params = params), 
                         silent = TRUE)
  if (class(params_id_check) == "try-error") {
    params_id_check <- NULL
  }
  
  if (is.null(params_id_check)) {
    # If file does not exist or file exists and params not ID'd
    params_id <- set_params_id(metadata = metadata, 
                               params = params,
                               params_id = params_id)
  } else {
    if (length(params_id_check) > 1) {
      stop(paste("Multiple entries found for the parameters:",
                 str_flatten(params_id_check, collapse = " ")))
    } else {
      # File exists and params identified
      message(paste("Parameters already identified:", params_id_check))
      params_id <- params_id_check
    }
  }
  
  outdir <- file.path(outdir, params_id)
  if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
  
  return(outdir)
  
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
  
  # Append suffix if specified
  if (!is.null(suffix)) {
    outfile <- paste0(tools::file_path_sans_ext(infile), 
                      suffix, ".", tools::file_ext(infile))
  } else {
    outfile <- infile
  }
  
  # Create output directory if needed
  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) {
      dir.create(outdir, 
                 showWarnings = FALSE, 
                 recursive = TRUE)
    }
    outfile <- basename(outfile)
    outfile <- file.path(outdir, outfile)
  } else {
    if (is.null(suffix)) {
      stop("Parameters outdir and suffix cannot both be NULL.")
    }
  }
  
  cmd_autocrop <- paste("autocrop", "-quiet", "-clobber",
                        "-isostep", isostep, infile, outfile)
  
  system(cmd_autocrop)
  
  return(outfile)
}
