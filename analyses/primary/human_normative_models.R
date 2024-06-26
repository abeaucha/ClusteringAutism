#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(MRIcrotome))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))

fit_predict_mean <- function(y, demographics, newdata, batch = NULL, df = 3) {
  
  # Residualize using batch variable if specified
  if (!is.null(batch)) {
    batch <- demographics %>%
      select(all_of(batch)) %>%
      unite(col = batch) %>%
      pull(batch)
    y <- residuals(lm(y ~ batch))
    names(y) <- NULL
  }
  
  ind_fit <- demographics[["DX"]] == "Control"
  
  # Training data frame
  df_fit <- demographics[ind_fit, c("Age", "Sex")]
  df_fit[["y"]] <- y[ind_fit]
  
  model_fit <- lm(y ~ Sex + ns(Age, df = df), data = df_fit)
  
  y_pred <- predict(model_fit, newdata = newdata)
  
  return(y_pred)
  
}


# Main -----------------------------------------------------------------------

# Output directory
output_dir <- "outputs/human_normative_models/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Pipeline directory
version <- "v3"
pipeline_dir <- "../../data/human/derivatives/"
pipeline_dir <- file.path(pipeline_dir, version)

# Registration directory
registration_dir <- "../../data/human/registration/v3/"

# Parameter set IDs
param_ids <- c("547", "664")

# Get pipeline parameters
metadata <- file.path(pipeline_dir, "metadata.csv")
# params <- fetch_params_metadata(metadata,
# id = param_ids)
params <- read_csv(metadata) %>% filter(id %in% param_ids)

params <- params %>% 
  select(id, dataset, resolution, es_df, es_batch)

params <- params %>% 
  bind_rows(mutate(params, es_batch = NA))

for (i in 1:nrow(params)) {
  
  print(paste("i", "=", i))
  
  # Parameter set ID
  param_id <- params[[i, "id"]]
  
  dataset <- params[[i, "dataset"]]
  
  # Image resolution
  resolution <- params[[i, "resolution"]]
  resolution <- sprintf("%.1f", resolution)
  
  # Effect size batch variables
  batch <- params[[i, "es_batch"]]  
  if (is.na(batch)) {batch <- NULL}
  
  # Spline dof
  df <- params[[i, "es_df"]]
  
  # Directory to Jacobian images
  jacobians_dir <- file.path(registration_dir, "jacobians_resampled/resolution_3.0/")
  
  # Image mask
  mask <- paste0("mask_", resolution, "mm.mnc")
  mask <- file.path(registration_dir, "reference_files", mask)
  
  # Path to demographics
  demographics <- file.path(pipeline_dir, param_id, "demographics.csv")
  # demographics <- read_csv(demographics, show_col_types = FALSE) %>%
  demographics <- read_csv(demographics) %>% 
    filter(!is.na(DX),
           !is.na(Age),
           !is.na(Sex))
  
  if (!is.null(batch)) {
    batch <- str_split(batch, pattern = "-")[[1]]
    batch_check <- batch %in% colnames(demographics)
    if (!all(batch_check)) {
      stop("Batch columns not found in demographics:\n",
           str_flatten(batch, collapse = "\n"))
    }
  }
  
  # Jacobians
  jacobians <- c("absolute", "relative")
  sex <- c("Male", "Female")
  ages <- seq(5, 30, by = 1)
  for (j in 1:length(jacobians)) {
    for (s in 1:length(sex)) {
      
      print(paste(jacobians[j], sex[s], sep = ","))
      
      # Image files
      imgdir <- file.path(jacobians_dir, jacobians[j])
      imgfiles <- list.files(imgdir, full.names = TRUE)
      
      imgs_in_demographics <- basename(imgfiles) %in% demographics[["file"]]
      imgfiles <- imgfiles[imgs_in_demographics]
      row_match <- match(basename(imgfiles), demographics[["file"]])
      demographics <- demographics[row_match,]
      
      df_pred <- tibble(Age = ages, Sex = sex[s])
      
      voxels <- qMincApply(filenames = imgfiles,
                           fun = fit_predict_mean,
                           newdata = df_pred,
                           demographics = demographics,
                           batch = batch,
                           df = df,
                           mask = mask,
                           batches = 50,
                           source = file.path(SRCPATH, "processing.R"),
                           registry_name = "registry_normative_growth",
                           cleanup = FALSE,
                           return_raw = TRUE,
                           resources = list(memory = "8G",
                                            walltime = 30*60))
      voxels <- simplify_masked(voxels[["vals"]])
      voxels <- t(voxels)
      colnames(voxels) <- paste0("V", 1:ncol(voxels))

      df_voxels <- bind_cols(df_pred, as_tibble(voxels))
      
      outfile <- paste("normative_models", dataset, jacobians[j], str_to_lower(sex[s]), sep = "_")
      if (!is.null(batch)) {outfile <- paste(outfile, "wbatch", sep = "_")}
      outfile <- paste0(outfile, ".csv")
      outfile = file.path(output_dir, outfile)
      write_csv(x = df_voxels, path = outfile)
      
    }
  }
}