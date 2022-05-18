# ----------------------------------------------------------------------------
# clustering_temp.R
# Author: Antoine Beauchamp, Jacob Ellegood
# Created: May 17th, 2022
#
# Description
# -----------
#

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SNFtool))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--file1',
              type = 'character',
              help = ""),
  make_option('--file2',
              type = 'character',
              help = ""),
  make_option('--rownames',
              type = 'character',
              help = ""),
  make_option('--nclusters',
              type = 'numeric',
              default = 10,
              help = ""),
  make_option('--K',
              type = 'numeric',
              default = 10,
              help = ""),
  make_option('--t',
              type = 'numeric',
              default = 20,
              help = ""),
  make_option('--alpha',
              type = 'numeric',
              default = 0.5,
              help = ""),
  make_option('--metric',
              type = 'character',
              default = 'cor',
              help = ""),
  make_option('--outfile',
              type = 'character',
              help = ""),
  make_option('--wsave',
              type = 'character',
              default = 'false',
              help = ""),
  make_option('--wfile',
              type = 'character',
              help = ""),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))  
)


# Functions ------------------------------------------------------------------

SNF_combine <-function(x1, x2, K=10, alpha=0.5, t=20, metric='cor', 
                       save_output = FALSE, outfile = 'W_matrix.RData'){
  
  if(metric == 'cor'){
    d1 <- (1-cor(t(x1)))
    d2 <- (1-cor(t(x2)))
  } else {
    d1 <- (dist2(as.matrix(x1), as.matrix(x1)))^(1/2)
    d2 <- (dist2(as.matrix(x2), as.matrix(x2)))^(1/2)
  }
  
  W1 <- affinityMatrix(d1, K = K, sigma = alpha)
  W2 <- affinityMatrix(d2, K = K, sigma = alpha)
  
  W <- SNF(list(W1, W2), K = K, t = t)
  
  if (save_output){
    save(W, file = outfile)
  }
  return(W)
}

create_clusters <- function(W, nk = 10, method = 'spectral', 
                            save_output = FALSE, outfile = 'clusters.csv') {
  if(method == 'spectral') {
    for(k in 2:nk) {
      group <- spectralClustering(affinity = W, K = k)
      group_name <- paste0('Group', k)
      assign(group_name, group)
      if (k == 2) {
        all_clusters <- data.frame(rownames(W), group, stringsAsFactors = F)
        colnames(all_clusters) <- c('ID', group_name)
      } else {
        group <- data.frame(group)
        colnames(group) <- group_name
        all_clusters <- cbind(all_clusters, group)
      }
    }
  }
  
  if (save_output) {
    write.csv(x = all_clusters, file = outfile, row.names = FALSE)
  }
  
  return(all_clusters)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
file1 <- args[['file1']]
file2 <- args[['file2']]
row_names <- args[['rownames']]
outfile <- args[['outfile']]
save_w <- ifelse(args[['wsave']] == 'true', TRUE, FALSE)
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

if (verbose) {message("Importing data...")}

x1 <- data.table::fread(file1, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames(row_names) %>% 
  as.matrix() 
colnames(x1) <- NULL
rownames(x1) <- basename(rownames(x1))

x2 <- data.table::fread(file2, header = TRUE) %>% 
  as_tibble() %>% 
  column_to_rownames(row_names) %>% 
  as.matrix() 
colnames(x2) <- NULL
rownames(x2) <- basename(rownames(x2))

if (verbose) {message("Running similarity network fusion...")}

W <- SNF_combine(x1 = x1, x2 = x2,
                 K = args[['K']], alpha = args[['alpha']],
                 t = args[['t']], metric = args[['metric']],
                 save_output = save_w, outfile = args[['wfile']])

if (verbose) {message("Assigning clusters...")}

clusters <- create_clusters(W = W, nk = args[['nclusters']],
                            method = 'spectral', save_output = TRUE,
                            outfile = outfile)
