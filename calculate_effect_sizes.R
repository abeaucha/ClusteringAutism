# ----------------------------------------------------------------------------
# calculate_effect_sizes.R
# Author: Antoine Beauchamp
# Created: May 9th, 2022
#
# Description
# -----------
#

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--demographics',
              type = 'character',
              help = "Path to CSV file containing demographics data."),
  make_option('--imgdir',
              type = 'character',
              help = ""),
  make_option('--outdir',
              type = 'character',
              help = "")
) 




#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

demofile <- args[['demographics']]
imgdir <- args[['imgdir']]
outdir <- args[['outdir']]

demofile <- 'data/human/registration/DBM_input_demo_passedqc.csv'
demographics <- data.table::fread(demofile, header = TRUE) %>% 
  as_tibble()

demographics <- demographics %>% 
  filter(!is.na(DX))

controls <- demographics %>% 
  filter(DX == 'Control')

participants <- demographics %>% 
  filter(DX != 'Control')

i <- 1
participants[i, 'Extract_ID'][[1]]




