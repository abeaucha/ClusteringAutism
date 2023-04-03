# ----------------------------------------------------------------------------
# normative_growth_model.R
# Author: Antoine Beauchamp
# Created: February 16th, 2023
#
# Brief description
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(doSNOW))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--demographics",
              type = "character",
              help = paste("Help message")),
  make_option("--imgdir",
              type = "character",
              help = paste("Help message")),
  make_option("--mask",
              type = "character"),
  make_option("--key",
              type = "character",
              default = "file",
              help = paste("Help message")),
  make_option("--df",
              type = "numeric",
              help = paste("Help message")),
  make_option("--combat",
              type = "character",
              default = "true",
              help = paste("Help message")),
  make_option("--combat-batch",
              type = "character",
              help = paste("Help message")),
  make_option("--outdir",
              type = "character",
              help = paste("Help message")),
  make_option("--matrix-file",
              type = "character",
              default = "effect_sizes.csv",
              help = paste("Help")),
  make_option("--parallel",
              type = "character",
              default = "false",
              help = "Option to run in parallel. [default %default]"),
  make_option("--nproc",
              type = "numeric",
              help = paste("Number of processors to use in parallel.",
                           "Ignored if --parallel is false.")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = "Verbosity [default %default]")
) 


# Functions ------------------------------------------------------------------

# Processing functions
source("src/processing.R")


#' Function description
#'
#' @param x (class) Parameter description.
#'
#' @return (class) Return description.
fit_predict_model <- function(y, demographics, df) {
  
  if (length(y) != nrow(demographics)) {
    stop()
  }
  
  ind_fit <- demographics[["DX"]] == "Control"
  ind_pred <- !ind_fit
  
  df_fit <- demographics[ind_fit, c("Age", "Sex")] %>% 
    mutate(y = y[ind_fit])
  
  df_pred <- demographics[ind_pred, c("Age", "Sex")] %>% 
    mutate(y = y[ind_pred])
  
  model_fit <- lm(y ~ Sex + ns(Age, df = df), data = df_fit)
  model_pred <- predict(model_fit, 
                        newdata = df_pred, 
                        interval = "prediction",
                        level = pnorm(q = 1) - pnorm(q = -1))
  
  df_pred <- df_pred %>% 
    mutate(y_pred = model_pred[,"fit"],
           y_lwr = model_pred[,"lwr"],
           y_upr = model_pred[,"upr"],
           y_sd = y_pred - y_lwr)
  
  return(df_pred)
}


zscore <- function(x){
  cols_check <- c("y", "y_pred", "y_sd")
  if (any(!(cols_check %in% colnames(x)))){
    stop()
  }
  x <- mutate(x, z = (y - y_pred)/y_sd)
  return(x)  
}


compute_normative_zscore <- function(y, demographics, df) {
  
  y_pred <- fit_predict_model(y = y, 
                              demographics = demographics,
                              df = df)
  z <- pull(zscore(y_pred), "z")
  return(z)
  
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

imgdir <- args[["imgdir"]]
demographics <- args[["demographics"]]
mask <- args[["mask"]]
key <- args[["key"]]
df <- args[["df"]]
combat <- ifelse(args[["combat"]] == "true", TRUE, FALSE)
combat_batch <- args[["combat-batch"]]
outdir <- args[["outdir"]]
matrix_file <- args[["matrix-file"]]
inparallel <- ifelse(args[["parallel"]] == "true", TRUE, FALSE)
nproc <- args[["nproc"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc_wfile.csv"
# imgdir <- "data/human/derivatives/POND_SickKids/jacobians/resolution_3.0/absolute/"
# mask <- "data/human/registration/reference_files/mask_3.0mm.mnc"
# key <- "file"
# df <- 3
# combat <- TRUE
# combat_batch <- "Site-Scanner"
# outdir <- "tmp/"
# outfile <- "effect_sizes.csv"
# verbose <- TRUE
# inparallel = TRUE
# nproc <- 4

if (inparallel) {
  if (is.null(nproc)) {
    stop("Specify the number of processors to use in parallel.")
  }
}

#Import demographics data
if (verbose) {message("Importing data...")}
demographics <- as_tibble(data.table::fread(demographics, header = TRUE))

#Import images as tibble
imgfiles <- list.files(imgdir, full.names = TRUE)
voxels <- build_voxel_matrix(imgfiles = imgfiles,
                             mask = mask,
                             file_col = TRUE,
                             sort = TRUE,
                             inparallel = inparallel,
                             nproc = nproc)
voxels[["file"]] <- basename(voxels[["file"]])
if (key != "file") {
  voxels[[key]] <- voxels[["file"]]
  voxels <- select(voxels, -file)
}

#Check existence of key column in demographics
if (!(key %in% colnames(demographics))) {
  stop(paste("demographics data is missing key column:", key))
}

#Check existence of key column in voxels
if (!(key %in% colnames(voxels))) {
  stop(paste("voxels data is missing key column:", key))
}

#Remove entries with missing diagnosis, age, or sex
demographics <- demographics %>% 
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex),
         !is.na(Site),
         !is.na(Scanner)) 

#Align voxels and demographics rows
if (verbose) {message("Aligning voxel and demographics data...")}
voxels <- semi_join(voxels, demographics, by = key)
row_match <- match(voxels[[key]], demographics[[key]])
demographics <- demographics[row_match,]
voxels <- select(voxels, -all_of(key))

#Apply ComBat normalization if specified
if (combat) {
  if (verbose) {message("Running ComBat normalization...")}
  if (is.null(combat_batch)) {
    stop("Argument --combat-batch must be specified for ComBat normalization.")
  }
  combat_batch <- str_split(combat_batch, pattern = "-")[[1]]
  batches_exist <- all(combat_batch %in% colnames(demographics))
  if (!batches_exist) {
    stop(paste("ComBat batch variables not found in --demographics:", 
               str_flatten(combat_batch, collapse = ", ")))
  }
  for (i in 1:length(combat_batch)) {
    message(paste("Normalizing using batch variable:", combat_batch[i]))
    batch <- demographics[[combat_batch[i]]]
    voxels <- t(as.matrix(voxels))
    voxels <- sva::ComBat(dat = voxels, batch = batch)
    voxels <- as_tibble(t(voxels))
  }
}

#Compute normative z-scores voxelwise
if (verbose) {message("Applying normative growth modelling...")}
pb <- txtProgressBar(max = ncol(voxels), style = 3)
progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
if (inparallel) {
  cl <- makeSOCKcluster(nproc)
  registerDoSNOW(cl)
  opts <- list(progress=progress)
  zscores <- foreach(j = 1:ncol(voxels), 
                     .packages = c("tidyverse", "splines"), 
                     .options.snow = opts) %dopar% {
                       compute_normative_zscore(y = voxels[[j]], 
                                                demographics = demographics, 
                                                df = df)
                     }
  close(pb)
  stopCluster(cl)
} else {
  zscores <- foreach(j = 1:ncol(voxels), 
                     .packages = c("tidyverse", "splines")) %do% {
                       progress(n = j)
                       compute_normative_zscore(y = voxels[[j]],
                                                demographics = demographics,
                                                df = df)
                     }
  close(pb)
}

#Store zscore data in matrix
voxels_norm <- matrix(data = 0, nrow = length(zscores[[1]]), ncol = length(zscores))
for (j in 1:length(zscores)) {
  voxels_norm[,j] <- zscores[[j]]
}
colnames(voxels_norm) <- 1:length(zscores)

#Create output directory if necessary
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Export z-scores to images
if (verbose) {message("Exporting normalized images...")}
outfiles <- demographics[demographics[["DX"]] != "Control", key][[1]]
outfiles <- file.path(outdir, outfiles)
outfiles <- matrix_to_images(x = voxels_norm, outfiles = outfiles, mask = mask)

#Convert z-scores to data frame
df_voxels_norm <- as_tibble(voxels_norm)
df_voxels_norm[[key]] <- basename(outfiles)

#Write effect sizes to file
if (verbose) {message("Exporting normalized voxel matrix...")}
matrix_file <- basename(matrix_file)
matrix_file <- file.path(outdir, matrix_file)
data.table::fwrite(df_voxels_norm, matrix_file)

