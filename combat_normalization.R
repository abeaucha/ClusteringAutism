# ----------------------------------------------------------------------------
# template.R
# Author: Antoine Beauchamp
# Created:
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
# args <- parse_args(OptionParser(option_list = option_list))
# arg <- args[['arg']]
# verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

# if (verbose) {message("")}

demographics <- "data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc.csv"
infile <- "data/human/derivatives/POND_SickKids/jacobians_3mm/absolute/jacobians.csv"
outfile <- "data/human/derivatives/POND_SickKids/combat/jacobians_normalized.csv"

outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Import demographics data
demographics <- as_tibble(data.table::fread(demographics, header = TRUE))

ind_file_col <- str_detect(colnames(demographics), "(f|F)ile")
if (any(ind_file_col)) {
  colnames(demographics)[ind_file_col] <- str_to_lower(colnames(demographics)[ind_file_col])
} else {
  stop()
}

#Remove entries with missing diagnosis, age, or sex
demographics <- demographics %>% 
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex),
         !is.na(Site),
         !is.na(Scanner))


dat <- as_tibble(data.table::fread(infile, header = TRUE))
ind_file_col <- str_detect(colnames(dat), "(f|F)ile")
if (any(ind_file_col)) {
  colnames(dat)[ind_file_col] <- str_to_lower(colnames(dat)[ind_file_col])
} else {
  stop()
}

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

outfile <- file.path(outdir, outfile)
data.table::fwrite(x = df_dat_norm, file = outfile)

