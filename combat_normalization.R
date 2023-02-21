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


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--arg',
              type = 'character',
              help = paste("Help message")),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))
) 


# Functions ------------------------------------------------------------------

#' Function description
#'
#' @param x (class) Parameter description.
#'
#' @return (class) Return description.
template <- function(x){
  return()
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
demographics <- args[["demographics"]]
infile <- args[["infile"]]
outfile <- args[["outfile"]]

# if (verbose) {message("")}

demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc.csv"
voxels <- "data/human/derivatives/POND_SickKids/jacobians_3mm/absolute/jacobians.csv"
outfile <- "data/human/derivatives/POND_SickKids/combat/jacobians_normalized.csv"
key <- "file"

outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Import data
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
dat <- semi_join(dat, demographics, by = "file")

ind_match <- match(dat[["file"]], demographics[["file"]])
demographics <- demographics[ind_match,]

batch <- "Site"
batch <- demographics[[batch]]

dat <- dat %>% 
  column_to_rownames(var = "file") %>% 
  as.matrix() %>% 
  t()

dat_norm <- sva::ComBat(dat = dat,
                        batch = batch)

dat_norm <- t(dat_norm)

df_dat_norm <- as_tibble(dat_norm, rownames = "file")

data.table::fwrite(x = df_dat_norm, file = outfile)

