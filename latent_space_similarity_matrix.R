# ----------------------------------------------------------------------------
# latent_space_similarity_matrix.R
# Author: Antoine Beauchamp
# Created: June 6th, 2022
#
# Create a similarity matrix from latent space signatures.
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--mouse',
              type = 'character',
              help = paste("Path to .csv file containing mouse latent space",
                           "signatures.")),
  make_option('--human',
              type = 'character',
              help = paste("Path to .csv file containing human latent space",
                           "signatures.")),
  make_option('--metric',
              type = 'character',
              default = 'correlation',
              help = paste("Metric used to compute the similarity matrix.",
                           "[default %default]")),
  make_option('--outfile',
              type = 'character',
              help = paste("Path to .csv file containing similarity matrix.")),
  make_option('--save-intermediate',
              type = 'character',
              default = 'false',
              help = paste("Option to save intermediate latent space",
                           "similarity matrices. File names will be ",
                           "adapted from --outfile. [default %default]")),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))
) 


# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset("--file=") %>% 
  str_remove("--file=") %>% 
  dirname()

path_func <- file.path(working_dir,
                       script_dir,
                       "functions",
                       "buildSimilarityMatrix.R")
source(path_func)


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
mousefile <- args[['mouse']]
humanfile <- args[['human']]
metric <- args[['metric']]
outfile <- args[['outfile']]
save_intermediate <- ifelse(args[['save-intermediate']] == 'true', TRUE, FALSE)
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

#Metric error catch
metric_choices <- c('correlation', 'cosine', 'euclidean')
if (!(metric %in% metric_choices)) {
  stop(paste("Argument --metric must be one of",
             "[correlation, cosine, euclidean]"))
}

#Create outdir if needed
outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

if (verbose) {message("Preparing data...")}

df_mouse <- data.table::fread(mousefile, header = TRUE) %>% 
  as_tibble()

df_human <- data.table::fread(humanfile, header = TRUE) %>% 
  as_tibble()

df_mouse <- df_mouse %>%
  mutate(latent_space = exprfile %>% 
           str_extract('transform_[0-9]+') %>% 
           str_extract('[0-9]+') %>% 
           as.numeric(),
         nk = clusterfile %>% 
           str_extract('Clusternum_[0-9]+') %>% 
           str_extract('[0-9]+') %>% 
           as.numeric(),
         k = clusterfile %>% 
           str_extract('Group_[0-9]+') %>% 
           str_extract('[0-9]+') %>% 
           as.numeric()) %>% 
  arrange(latent_space, nk, k) %>% 
  unite(cluster_id, nk, k, sep = '-')

df_human <- df_human %>%
  mutate(latent_space = exprfile %>% 
           str_extract('transform_[0-9]+') %>% 
           str_extract('[0-9]+') %>% 
           as.numeric(),
         nk = clusterfile %>% 
           str_extract('Clusternum_[0-9]+') %>% 
           str_extract('[0-9]+') %>% 
           as.numeric(),
         k = clusterfile %>% 
           str_extract('Group_[0-9]+') %>% 
           str_extract('[0-9]+') %>% 
           as.numeric()) %>% 
  arrange(latent_space, nk, k) %>% 
  unite(cluster_id, nk, k, sep = '-')

column_names <- unique(df_mouse[,'cluster_id'][[1]])
row_names <- unique(df_human[,'cluster_id'][[1]])

mouse_latent_space_ids <- df_mouse[,'latent_space'][[1]] %>%
  unique() %>% 
  sort()

human_latent_space_ids <- df_human[,'latent_space'][[1]] %>%
  unique() %>% 
  sort()

if(length(mouse_latent_space_ids) != length(human_latent_space_ids)) {
  stop("Mouse and human expression files contain a different number of latent spaces.")
} else {
  if (sum(mouse_latent_space_ids != human_latent_space_ids) != 0) {
    stop("Mouse and human expression files contain different latent spaces.")
  } else {
    latent_space_ids <- intersect(mouse_latent_space_ids,
                                  human_latent_space_ids)
  }
}

if (verbose) {message("Computing similarity matrices...")}

list_sim <- vector(mode = 'list', length = length(latent_space_ids))
for (i in 1:length(list_sim)){
  
  #Prepare mouse data
  mat_mouse <- df_mouse %>% 
    filter(latent_space == latent_space_ids[i]) %>% 
    select(-clusterfile,
           -exprfile,
           -latent_space) %>% 
    column_to_rownames('cluster_id') %>% 
    as.matrix() %>% 
    t() 
  
  #Prepare human data  
  mat_human <- df_human %>% 
    filter(latent_space == latent_space_ids[i]) %>% 
    select(-clusterfile,
           -exprfile,
           -latent_space) %>% 
    column_to_rownames('cluster_id') %>% 
    as.matrix() %>% 
    t()
  
  #Similarity matrix
  mat_sim <- buildSimilarityMatrix(x1 = mat_human,
                                   x2 = mat_mouse,
                                   method = metric)
    
  #Add to list
  list_sim[[i]] <- mat_sim
  
  #Save to file if desired
  if (save_intermediate) {
    outfile_tmp <- outfile %>% 
      str_replace('.csv', 
                  str_c('_', latent_space_ids[i], '.csv'))
    mat_sim %>% 
      as_tibble(rownames = 'cluster_id') %>% 
      data.table::fwrite(file = outfile_tmp)
  }
  
  if (i != 1) {  
    ncol_current <- ncol(list_sim[[i]])
    ncol_prev <- ncol(list_sim[[i-1]])
    if (ncol_current != ncol_prev) {
      stop(paste("Mismatch between number of mouse clusters across latent",
                 "spaces. Interrupted at iteration", i))
    }
    
    nrow_current <- nrow(list_sim[[i]])
    nrow_prev <- nrow(list_sim[[i-1]])
    if (nrow_current != nrow_prev) {
      stop(paste("Mismatch between number of human clusters across latent",
                 "spaces. Interrupted at iteration", i))
    }
  }
}

array_sim <- array(unlist(list_sim),
                   dim = c(nrow_current,
                           ncol_current,
                           length(list_sim)))

avg_sim <- rowMeans(array_sim, dims = 2)
rownames(avg_sim) <- row_names
colnames(avg_sim) <- column_names

avg_sim %>% 
  as_tibble(rownames = 'cluster_id') %>% 
  data.table::fwrite(file = outfile)
