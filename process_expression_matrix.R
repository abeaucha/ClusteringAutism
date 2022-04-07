# ----------------------------------------------------------------------------
# process_expression_matrix.R
# Antoine Beauchamp
# Created: August 25th, 2021
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--infile",
              type = "character",
              help = paste("Path to CSV file containing labelled expression",
                           "matrix.")),
  make_option("--transpose",
              type = "character",
              default = "false",
              help = paste("Transpose expression matrix.",
                           "[default %default]")),
  make_option("--scale",
              type = "character",
              default = "true",
              help = paste("Scale expression matrix.",
                           "[default %default]")),
  make_option("--aggregate",
              type = "character",
              default = "false",
              help = paste("Aggregate expression data under a set of",
                           "atlas labels. [default %default]")),
  make_option("--outdir",
              default = "data/",
              type = "character",
              help = paste("Output directory.",
                           "[default %default]")),
  make_option("--verbose",
              default = "true",
              type = "character",
              help = "[default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))

if (!(args[["transpose"]] %in% c("true", "false"))) {
  stop('Argument --transpose must be one of [true, false]')
}

if (!(args[["scale"]] %in% c("true", "false"))) {
  stop('Argument --scale must be one of [true, false]')
}

if (!(args[["aggregate"]] %in% c("true", "false"))) {
  stop('Argument --aggregate must be one of [true, false]')
}


# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset('--file=') %>% 
  str_remove('--file=') %>% 
  dirname()

path_processing_tools <- file.path(working_dir,
                                   script_dir,
                                   'functions',
                                   'processing_tools.R')
source(path_processing_tools)


# Main -----------------------------------------------------------------------

verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)
transpose <- ifelse(args[['transpose']] == 'true', TRUE, FALSE)
normalize <- ifelse(args[['scale']] == 'true', TRUE, FALSE)
aggregate <- ifelse(args[['aggregate']] == 'true', TRUE, FALSE)
infile <- args[['infile']]

if (!transpose & !normalize & !aggregate){
  stop("One of {--transpose, --scale, --aggregate} must be true.")
}


if (verbose) {message("Processing data from file: ", infile)}

if(verbose){message("Importing...")}

#Import data
dfExpression <- suppressMessages(data.table::fread(infile,
                                                   header = TRUE)) %>%
  as_tibble()

if (('Gene' %in% colnames(dfExpression)) & !transpose) {
  stop("Detected column named 'Gene' in matrix. Consider transposing.")
}

if (transpose) {
  
  if(verbose){message("Transposing...")}
  
  dfExpression <- dfExpression %>% 
    column_to_rownames('Gene') %>% 
    as.matrix() %>% 
    t() %>% 
    as_tibble() 
}

#Extract genes list from data
genes <- colnames(dfExpression)[!str_detect(colnames(dfExpression), 'Region')]

#Do any columns contain region labels?
containsLabels <- any(str_detect(colnames(dfExpression), 'Region'))

if (normalize) {
  
  if(verbose){message("Normalizing...")}
  
  #Extract labels from data frame
  if (containsLabels) {
    dfLabels <- dfExpression %>% select(contains('Region'))
  }
  
  #Normalize data
  dfExpression <- dfExpression %>% 
    select(all_of(genes)) %>% 
    as.matrix() %>% 
    scaler(axis = 'rows') %>% 
    scaler(scale = FALSE, axis = 'columns') %>% 
    as_tibble()
  
  if (containsLabels) {
    dfExpression <- dfExpression %>% 
      bind_cols(dfLabels)
  }
  
  outfile <- infile %>% 
    basename() %>% 
    str_replace('.csv', '_scaled.csv')
  
}

if (aggregate) {
  
  if(verbose){message("Aggregating...")}
  
  # #Aggregate mouse expression data under label set
  dfExpression <- dfExpression %>% 
    group_by(Region) %>% 
    summarise_all(mean) %>% 
    ungroup()
    
  outfile <- infile %>% 
    basename() %>% 
    str_replace('voxel', 'ROI')
  
  if (normalize) {
    outfile <- str_replace(outfile, ".csv", "_scaled.csv")
  }
  
}

outfile <- file.path(args[['outdir']], outfile)

if (verbose) {message(paste("Writing to file:", outfile, "..."))}

data.table::fwrite(dfExpression,
                   file = outfile)
