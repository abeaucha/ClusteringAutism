# ----------------------------------------------------------------------------
# mouse_cluster_enrichment.R
# Author: Antoine Beauchamp
# Created: 
#
# 
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

# option_list <- list(
#   make_option("--arg",
#               type = "character",
#               help = paste("Help message")),
#   make_option("--verbose",
#               type = "character",
#               default = "true",
#               help = paste("Verbosity option. [default %default]"))
# ) 


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "enrichment.R"))


# Main -----------------------------------------------------------------------

#Parse command line args
# args <- parse_args(OptionParser(option_list = option_list))
# arg <- args[["arg"]]
# verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# if (verbose) {message("")}

# Command line args? 
pipeline_dir <- file.path(PROJECTPATH, "data/mouse/derivatives/v3/")
enrichment_dir <- file.path(PROJECTPATH, "data/enrichment/")
params_id <- 107
models <- "model_names.csv"
nk_max <- 10
gene_scores <- 950
stringdb_versions <- "12.0"
bader_versions <- 2023

# Get pipeline parameters
metadata <- file.path(pipeline_dir, "metadata.csv")
params <- fetch_params_metadata(metadata = metadata, id = params_id)

# Update pipeline directory with parameter set ID
pipeline_dir <- file.path(pipeline_dir, params_id)

# Get cluster resolution
cluster_resolution <- params[["cluster_resolution"]]
cluster_resolution <- sprintf("%.1f", cluster_resolution)

# Update cluster directory
cluster_dir <- file.path(pipeline_dir, "clusters", 
                         paste0("resolution_", cluster_resolution))

# Model names dictionary
models <- file.path(pipeline_dir, "model_names.csv")
models <- get_model_genes(models)

# Import clusters
clusters <- file.path(cluster_dir, "clusters.csv")
clusters <- read_csv(clusters, show_col_types = FALSE) 

# Join with model gene names (and filter for single genes)
clusters <- clusters %>%
  rename(file = ID) %>% 
  right_join(models, by = "file") %>% 
  select(-file) %>% 
  rename(model = ID)

# Background set for enrichment analysis
background_set <- file.path(enrichment_dir, "sagittal_gene_table_normalized_filtered.csv")
background_set <- read_csv(background_set, show_col_types = FALSE) %>% 
  pull(msg.genes.acronym) %>% 
  as.character()

params <- expand_grid(stringdb = stringdb_versions,
                      bader = bader_versions,
                      score = gene_scores)
for (i in 1:nrow(params)) {
  
  bader_version <- params[[i, "bader"]]
  stringdb_version <- params[[i, "stringdb"]]
  gene_score <- params[[i, "score"]]
  
  message(
    paste(paste("StringDB version:", stringdb_version),
          paste("Bader version:", bader_version),
          paste("Gene score:", gene_score),
          sep = "\n")
  )
  
  output_dirs <- c("NeighbourhoodEnrichment", "NeighbourhoodInfo")
  output_dirs <- file.path(pipeline_dir, "enrichment", 
                           paste("StringDB", stringdb_version, "Bader", 
                                 bader_version, sep = "_"), 
                           gene_score, output_dirs)
  for (path in output_dirs) {
    if (!dir.exists(path)) {dir.create(path, recursive = TRUE)}
  }
  
  if (bader_version == 2020) {
    bader_modules <- file.path(enrichment_dir, "Human_Reactome_March_01_2020_symbol.gmt")
  } else if (bader_version == 2023) {
    bader_modules <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  
  } else {
    stop()
  }
  
  
  for (nk in 2:nk_max) {
    
    # Identify unique genes in each cluster
    cluster_genes_nk <- clusters %>%
      select(cluster = paste0("nk", nk), gene) %>% 
      unite(dup_select, "gene", "cluster", remove = FALSE) %>% 
      filter(!duplicated(dup_select)) %>% 
      select(cluster, gene)
    
    for (k in 1:nk) {
      
      message(paste0("On cluster ", nk, "-", k))
      
      cluster_genes_k <- cluster_genes_nk %>% 
        filter(cluster == k) %>% 
        pull(gene)
      
      message("Building gene neighbourhood...")
      neighbourhood_k <- get_gene_neighbourhood(genes = cluster_genes_k, 
                                               score = gene_score, 
                                               stringdb_version = stringdb_version)
      
      outfile <- paste("cluster_neighbourhood", nk, k, gene_score, sep = "_")
      outfile <- paste0(outfile, ".csv")
      outfile <- file.path(output_dirs[2], outfile)
      write_csv(neighbourhood_k, file = outfile)
      
      target_set <- unique(c(neighbourhood_k[["gene_A"]], neighbourhood_k[["gene_B"]]))
      
      message("Calculating pathway enrichment...")
      enrichment <- get_neighbourhood_enrichment(target = target_set,
                                                 background = background_set, 
                                                 modules = bader_modules)
      
      enrichment[["nk"]] <- nk
      enrichment[["k"]] <- k
      
      outfile <- paste("cluster_pathway_enrichment", nk, k, gene_score, sep = "_")
      outfile <- paste0(outfile, ".csv")
      outfile <- file.path(output_dirs[1], outfile)
      write_csv(enrichment, file = outfile)
      
    }
  }
}