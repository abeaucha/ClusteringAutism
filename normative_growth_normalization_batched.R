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
suppressPackageStartupMessages(library(RMINC))


# Command line arguments -----------------------------------------------------

option_list <- list(
#   make_option("--demographics",
#               type = "character",
#               help = paste("Help message")),
#   make_option("--imgdir",
#               type = "character",
#               help = paste("Help message")),
#   make_option("--mask",
#               type = "character"),
#   make_option("--key",
#               type = "character",
#               default = "file",
#               help = paste("Help message")),
#   make_option("--df",
#               type = "numeric",
#               help = paste("Help message")),
#   make_option("--outdir",
#               type = "character",
#               help = paste("Help message")),
#   make_option("--matrix-file",
#               type = "character",
#               default = "effect_sizes.csv",
#               help = paste("Help")),
  make_option("--nproc",
              type = "numeric",
              help = paste("Number of processors to use in parallel.",
                           "Ignored if --parallel is false."))
#   make_option("--verbose",
#               type = "character",
#               default = "true",
#               help = "Verbosity [default %default]")
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
# 
# imgdir <- args[["imgdir"]]
# demographics <- args[["demographics"]]
# mask <- args[["mask"]]
# key <- args[["key"]]
# df <- args[["df"]]
# outdir <- args[["outdir"]]
# matrix_file <- args[["matrix-file"]]
nproc <- args[["nproc"]]
# verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# demographics <- "data/tmp_16/POND_SickKids/DBM_input_demo_passedqc_wfile.csv"
demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc_wfile.csv"
# imgdir <- "data/tmp_16/POND_SickKids/jacobians/resolution_0.8/absolute/"
imgdir <- "data/human/derivatives/POND_SickKids/jacobians/resolution_3.0/absolute/"
# mask <- "data/human/registration/reference_files/mask_0.8mm.mnc"
mask <- "data/human/registration/reference_files/mask_3.0mm.mnc"
key <- "file"
df <- 3
outdir <- "tmp/"
matrix_file <- "effect_sizes.csv"
verbose <- TRUE
nproc <- 8

if (is.null(nproc)) {
  stop("Specify the number of processors to use in parallel.")
}

#Import demographics data
if (verbose) {message("Importing demographics information...")}
demographics <- as_tibble(data.table::fread(demographics, header = TRUE))

#Check existence of key column in demographics
if (!(key %in% colnames(demographics))) {
  stop(paste("demographics data is missing key column:", key))
}

#Remove entries with missing diagnosis, age, or sex
demographics <- demographics %>%
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex),
         !is.na(Site),
         !is.na(Scanner))

#Image files
imgfiles <- list.files(imgdir, full.names = TRUE)

#Match image files to demographics
if (verbose) {message("Matching image files to demographics...")}
imgs_in_demographics <- basename(imgfiles) %in% demographics[[key]]
imgfiles <- imgfiles[imgs_in_demographics]
row_match <- match(basename(imgfiles), demographics[[key]])
demographics <- demographics[row_match,]



mask_vol <- mincGetVolume(mask)
mask_ind <- which(mask_vol == 1)
batch_ind <- mask_ind[1:1000]
mask_batch <- numeric(length(mask_vol))
mask_batch[batch_ind] <- 1
attributes(mask_batch) <- attributes(mask_vol)
mask_batch_file <- "tmp/mask_batch.mnc"
mincWriteVolume(buffer = mask_batch,
                output.filename = mask_batch_file,
                like.filename = mask,
                clobber = TRUE)

ind_mask_tmp <- sample(x = which(maskvol == 1), size = 1000, replace = FALSE)
mask_tmp <- numeric(length(maskvol))
mask_tmp[ind_mask_tmp] <- 1
attributes(mask_tmp) <- attributes(maskvol)
outfile <- "masktmp.mnc"
mincWriteVolume(buffer = mask_tmp, output.filename = outfile, like.filename = mask, clobber = TRUE)

#Run normative growth modelling
if (verbose) {message("Evaluating normative growth models...")}
ti <- Sys.time()
voxels <- mcMincApply(filenames = imgfiles, 
                      fun = compute_normative_zscore,
                      demographics = demographics,
                      df = df,
                      mask = mask,
                      cores = nproc, 
                      return_raw = TRUE)
tf <- Sys.time()
tdiff <- tf-ti

ti <- Sys.time()
message("Converting voxel list to matrix")
voxels <- simplify_masked(voxels[["vals"]])
voxels <- asplit(voxels, MARGIN = 2)
tf <- Sys.time()
tdiff <- tf-ti

#Export images
if (verbose) {message("Exporting normalized images...")}
outfiles <- demographics[demographics[["DX"]] != "Control", key][[1]]
outfiles <- file.path(outdir, outfiles)
out <- parallel::mcmapply(vector_to_image, voxels, outfiles, 
                          MoreArgs = list(mask = mask),
                          SIMPLIFY = TRUE, mc.cores = nproc)


#THis is a much faster way to import (I think?)
#It also import the data in the orientation desired by ComBat
# imgdir <- "data/human/registration/jacobians_resampled/resolution_0.8/absolute/"
# imgfiles <- list.files(imgdir, full.names = TRUE)
# mask <- "data/human/registration/reference_files/mask_0.8mm.mnc"
# tmp <- pMincApply(filenames = imgfiles[1:100],
#                   fun = function(x) {return(x)},
#                   mask = mask,
#                   cores = nproc,
#                   local = TRUE)

