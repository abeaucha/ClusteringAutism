# ----------------------------------------------------------------------------
# permute_cluster_labels.R
# Author: Antoine Beauchamp
# Created:February 1st, 2023
#
# Permute the cluster labels for a given cluster number.


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--infile",
              type = "character",
              help = paste("Path to .csv file containing cluster",
                           "labels.")),
  make_option("--outdir",
              type = "character",
              help = paste("Path to directory in which to save",
                           "permutated cluster labels.")),
  make_option("--nk",
              type = "integer",
              help = paste("Number of clusters to use for",
                           "permutation.")),
  make_option("--np",
              type = "integer",
              help = paste("Number of permutations."))
) 


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
arg <- args[["arg"]]
infile <- args[["infile"]]
outdir <- args[["outdir"]]
nk <- args[["nk"]]
np <- args[["np"]]

#Create output directory if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

#Import cluster labels
df_clusters <- read_csv(infile, show_col_types = FALSE)

#Select desired number of clusters
kcol <- paste0("Group", nk)
cols <- c("ID", kcol)

#Permutation
for (p in 1:np) {
  set.seed(p)
  
  #Permute cluster labels
  df_permute <- df_clusters[,cols]
  df_permute[[kcol]] <- sample(x = df_clusters[[kcol]],
                               size = nrow(df_clusters),
                               replace = FALSE)
  
  #Write out permuted labels
  outfile <- basename(infile) %>% 
    str_remove(".csv") %>% 
    str_c("nk", nk, "permutation", p, sep = "_") %>% 
    str_c(".csv")
  outfile <- file.path(outdir, outfile)  
  write_csv(x = df_permute,
            file = outfile)
}
