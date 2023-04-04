# ----------------------------------------------------------------------------
# combat_normalization.R
# Author: Antoine Beauchamp
# Created: February 15th, 2023
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
  make_option("--imgdir",
              type = "character",
              help = paste("Help message")),
  make_option("--demographics",
              type = "character",
              help = paste("Help message")),
  make_option("--mask",
              type = "character"),
  make_option("--outdir",
              type = "character",
              help = paste("Help message")),
  make_option("--key",
              type = "character",
              default = "file",
              help = paste("Help message")),
  make_option("--batch",
              type = "character",
              help = paste("Help message")),
  make_option("--nproc",
              type = "numeric",
              help = "Number of processors to use in parallel."),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = "Verbosity [default %default]")
) 


# Functions ------------------------------------------------------------------

source("src/processing.R")


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
imgdir <- args[["imgdir"]]
demographics <- args[["demographics"]]
mask <- args[["mask"]]
key <- args[["key"]]
batch <- args[["batch"]]
outdir <- args[["outdir"]]
nproc <- args[["nproc"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc_wfile.csv"
# imgdir <- "data/human/registration/jacobians_resampled/resolution_3.0/absolute/"
# mask <- "data/human/registration/reference_files/mask_3.0mm.mnc"
# outdir <- "data/human/test/"
# key <- "file"
# batch <- "Site-Scanner"
# nproc <- 8
# verbose <- TRUE

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

#Check existence of key column
if (!(key %in% colnames(demographics))) {
  stop(paste("demographics data is missing key column:", key))
}

#Check existence of batch columns
if (!is.null(batch)) {
  batch <- str_split(batch, pattern = "-")[[1]]
  batch_check <- batch %in% colnames(demographics)
  if (!all(batch_check)) {
    stop("Batch columns not found in demographics:\n", str_flatten(batch, collapse = "\n"))
  }
}

#Image files
imgfiles <- list.files(imgdir, full.names = TRUE)

#Match image files to demographics
if (verbose) {message("Matching image files to demographics...")}
imgs_in_demographics <- basename(imgfiles) %in% demographics[[key]]
imgfiles <- imgfiles[imgs_in_demographics]
row_match <- match(basename(imgfiles), demographics[[key]])
demographics <- demographics[row_match,]

#Import images
if (verbose) {message("Importing images...")}
voxels <- import_images(imgfiles = imgfiles, 
                        mask = mask, 
                        output_format = "matrix", 
                        margin = 2,
                        inparallel = TRUE, 
                        nproc = nproc)
gc()

#Pull batch information
batch <- demographics %>% 
  select(all_of(batch)) %>% 
  unite(col = batch) %>% 
  pull(batch)

#Run ComBat normalization
voxels <- sva::ComBat(dat = voxels, batch = batch)

#Export images
if (verbose) {message("Exporting ComBat normalized images...")}
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
outfiles <- file.path(outdir, basename(imgfiles))
matrix_to_images(x = voxels, outfiles = outfiles, mask = mask,
                 margin = 2, nproc = nproc)

