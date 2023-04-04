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

demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc_wfile.csv"
# imgdir <- "data/human/registration/jacobians_resampled/resolution_0.8/absolute/"
imgdir <- "data/human/registration/jacobians_resampled/resolution_3.0/absolute/"
# mask <- "data/human/registration/reference_files/mask_0.8mm.mnc"
mask <- "data/human/registration/reference_files/mask_3.0mm.mnc"
key <- "file"
df <- 3
outdir <- "tmp_batch/"
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

#Batch size
batch_size <- 5000

#Mask voxels
mask_vol <- mincGetVolume(mask)
mask_ind <- which(mask_vol == 1)
mask_size <- length(mask_ind)

#Number of batches
nbatches <- ceiling(mask_size/batch_size)

#Create batches
batches <- rep(1:nbatches, each = batch_size)
batches <- batches[1:mask_size]
batches <- factor(batches)
batches <- split(x = mask_ind, f = batches)

#For loop here
for (b in 1:nbatches) {
  
  message("On batch ", b, " of ", nbatches)
  
  #Create batch mask
  batch_ind <- batches[[b]]
  batch_mask <- numeric(length(mask_vol))
  batch_mask[batch_ind] <- 1
  attributes(batch_mask) <- attributes(mask_vol)
  batch_mask_file <- "batch_mask.mnc"
  batch_mask_file <- file.path(outdir, batch_mask_file)
  sink(nullfile())
  mincWriteVolume(buffer = batch_mask,
                  output.filename = batch_mask_file,
                  like.filename = mask,
                  clobber = TRUE)
  sink(NULL)
  
  #Run normative growth modelling
  batch_voxels <- mcMincApply(filenames = imgfiles, 
                              fun = compute_normative_zscore,
                              demographics = demographics,
                              df = df,
                              mask = batch_mask_file,
                              cores = nproc, 
                              return_raw = TRUE)
  batch_voxels <- simplify_masked(batch_voxels[["vals"]])
  
  #Export batched voxels to csv
  batch_voxels <- t(batch_voxels)
  colnames(batch_voxels) <- 1:ncol(batch_voxels)
  batch_voxels <- as_tibble(batch_voxels)
  batch_file <- paste0("batch_", b, ".csv")
  batch_file <- file.path(outdir, batch_file)
  data.table::fwrite(x = batch_voxels, file = batch_file, col.names = FALSE)
  
}

outfiles <- demographics[demographics[["DX"]] != "Control", key][[1]]
outfiles <- file.path(outdir, outfiles)
for (i in 1:length(outfiles)) {
  
  img <- numeric(length = length(mask_vol))
  
  for (b in 1:nbatches) {
    batch_file <- paste0("batch_", b, ".csv")
    batch_file <- file.path(outdir, batch_file)
    img_batch <- data.table::fread(batch_file, header = FALSE, skip = i-1, nrows = 1)
    img_batch <- as.matrix(img_batch)
    img_batch <- img_batch[1,]
    names(img_batch) <- NULL
    img[batches[[b]]] <- img_batch
  }
  
  attributes(img) <- attributes(mask_vol)
  outfile <- outfiles[i]
  sink(nullfile())
  mincWriteVolume(buffer = img, 
                  output.filename = outfile, 
                  like.filename = file.path(imgdir, basename(outfile)), 
                  clobber = TRUE)
  sink(NULL)
  
  file1 <- file.path("tmp", basename(outfile))
  file2 <- file.path("tmp_batch", basename(outfile))

  img1 <- mincGetVolume(file1)
  img2 <- mincGetVolume(file2)

  img1 <- img1[mask_vol > 0.5]
  img2 <- img2[mask_vol > 0.5]

  n <- 9
  ind_match <- round(img1, n) == round(img2, n)
  print(sum(!ind_match))

}


