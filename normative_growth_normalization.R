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
suppressPackageStartupMessages(library(tcltk))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--demographics",
              type = "character",
              help = paste("Help message")),
  make_option("--voxels",
              type = "character",
              help = paste("Help message")),
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
  make_option("--outfile",
              type = "character",
              help = paste("Help message")),
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
              default = "false",
              help = "Verbosity [default %default]")
) 


# Functions ------------------------------------------------------------------

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
  model_pred <- predict(model_fit, df_pred, se = TRUE)
  
  df_pred <- df_pred %>% 
    mutate(y_pred = model_pred[["fit"]],
           y_se = model_pred[["se.fit"]])
  
  return(df_pred)
}

zscore <- function(x){
  cols_check <- c("y", "y_pred", "y_se")
  if (any(!(cols_check %in% colnames(x)))){
    stop()
  }
  x <- mutate(x, z = (y - y_pred)/y_se)
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
demographics <- args[["demographics"]]
voxels <- args[["voxels"]]
key <- args[["key"]]
df <- args[["df"]]
combat <- ifelse(args[["combat"]] == "true", TRUE, FALSE)
combat_batch <- args[["combat-batch"]]
outfile <- args[["outfile"]]
inparallel <- ifelse(args[["parallel"]] == "true", TRUE, FALSE)
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc.csv"
# voxels <- "data/human/derivatives/POND_SickKids/jacobians_3mm/absolute/jacobians.csv"
# outfile <- "data/human/derivatives/POND_SickKids/effect_sizes/ES_normative.csv"
# key <- "file"
# df <- 5
# combat <- TRUE
# combat_batch <- "Site-Scanner"
# inparallel <- FALSE

#Import data
if (verbose) {message("Importing data...")}
demographics <- as_tibble(data.table::fread(demographics, header = TRUE))
demographics <- demographics %>% rename(file = File) #REMOVE THIS LINE ONCE FIXED
voxels <- as_tibble(data.table::fread(voxels, header = TRUE))

#Check existence of key column
if (!(key %in% colnames(demographics))) {
  stop(paste("demographics data is missing key column:", key))
}

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
    stop(paste("ComBat batch variables not found in --demographics:", str_flatten(combat_batch, collapse = ", ")))
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
  nproc <- args[["nproc"]]
  cl <- makeSOCKcluster(nproc)
  registerDoSNOW(cl)
  opts <- list(progress=progress)
  zscores <- foreach(j = 1:ncol(voxels), 
                     .packages = c("tidyverse", "splines"),
                     .combine = "cbind",
                     .options.snow = opts) %dopar% {
                       compute_normative_zscore(y = voxels[[j]],
                                                demographics = demographics,
                                                df = df)
                     } 
  close(pb)
  stopCluster(cl)
} else {
  zscores <- foreach(j = 1:ncol(voxels), 
                     .packages = c("tidyverse", "splines"),
                     .combine = "cbind") %do% {
                       progress(n = j)
                       compute_normative_zscore(y = voxels[[j]],
                                                demographics = demographics,
                                                df = df)
                     }
  close(pb)
}


#Convert z-scores to data frame
colnames(zscores) <- 1:ncol(voxels)
zscores <- as_tibble(zscores)
zscores[[key]] <- demographics %>%
  filter(DX != "Control") %>%
  pull(file)

#Create output directory if necessary
outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Write effect sizes to file
if (verbose) {message("Writing out normalized data...")}
data.table::fwrite(zscores, outfile)
