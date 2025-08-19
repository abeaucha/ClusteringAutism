#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# human_cluster_enrichment.R
# Authors: Antoine Beauchamp
# Created: November 15th, 2023


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--pipeline-dir",
              type = "character",
              default = "data/human/derivatives/"),
  make_option("--enrichment-dir",
              type = "character",
              default = "data/enrichment/"),
  make_option("--params-id",
              type = "character"),
  make_option("--variants",
              type = "character",
              default = "cluster_variants.csv"),
  make_option("--gene-score",
              type = "numeric",
              default = 950),
  make_option("--stringdb-version",
              type = "character",
              default = "12.0"),
  make_option("--bader-version",
              type = "numeric",
              default = 2025)
)


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "enrichment.R"))


# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
pipeline_dir <- args[["pipeline-dir"]]
enrichment_dir <- args[["enrichment-dir"]]
params_id <- args[["params-id"]]
variants <- args[["variants"]]
gene_score <- args[["gene-score"]]
stringdb_version <- args[["stringdb-version"]]
bader_version <- args[["bader-version"]]

if (is.null(params_id)) {
  stop("Specify the parameter set ID for the pipeline outputs.")
}

# Get pipeline parameters
metadata <- file.path(pipeline_dir, "metadata.csv")
params <- fetch_params_metadata(metadata = metadata, 
                                params = list(id = params_id))

# Path to pipeline version directory
pipeline_dir <- file.path(pipeline_dir, params_id)
if (!dir.exists(pipeline_dir)) {stop("Pipeline outputs not found.")}

# Get max number of clusters
nk_max <- params[["cluster_nk_max"]]

# Import variants
variants <- file.path(pipeline_dir, variants)
df_variants <- read_csv(variants, show_col_types = FALSE)

# Remove entries without genes
df_variants <- df_variants %>% 
  filter(!is.na(gene))

# Identify genes for those patients with single gene mutations
# Version 2
# df_variants[["gene"]] <- c(NA, NA, 
#                            "Kcnq2", "Gria2",
#                            "Rai1", "Fmr1",
#                            "Tnpo3", "Stxbp1",
#                            NA, "Rbm10", 
#                            NA, NA, 
#                            "Fgf14", "Tsc1",
#                            "Meis2", "Rfx3",
#                            "Pten", "Rtl9",
#                            "Nlgn1")

# Version 3
# df_variants[["gene"]] <- c("Pten", "Rtl9", 
#                            "Nlgn1", NA,
#                            "Rbm10", NA,
#                            NA, "Fgf14",
#                            "Tsc1", "Meis2",
#                            NA, "Kcnq2",
#                            "Gria2", "Rfx3",
#                            "Rai1", "Fmr1",
#                            "Tnpo3", "Stxbp1",
#                            NA) 

# Import background gene set
background_set <- "sagittal_gene_table_normalized_filtered.csv"
background_set <- file.path(enrichment_dir, background_set)
background_set <- read_csv(background_set, show_col_types = FALSE) %>% 
  pull(msg.genes.acronym) %>% 
  as.character()

message(
  paste(paste("StringDB version:", stringdb_version),
        paste("Bader version:", bader_version),
        paste("Gene score:", gene_score),
        sep = "\n")
)
  
# Create output directory
output_dirs <- c("NeighbourhoodEnrichment", "NeighbourhoodInfo")
output_dirs <- file.path(pipeline_dir, "enrichment", 
                        paste("StringDB", stringdb_version, "Bader", 
                              bader_version, sep = "_"), 
                        gene_score, output_dirs)
for (path in output_dirs) {
  if (!dir.exists(path)) {dir.create(path, recursive = TRUE)}
}

# Set path to Bader modules file
if (bader_version == 2020) {
  bader_modules <- file.path(enrichment_dir, "Human_Reactome_March_01_2020_symbol.gmt")
} else if (bader_version == 2023) {
  bader_modules <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  
} else if (bader_version == 2025) {
  bader_modules <- file.path(enrichment_dir, "Human_Reactome_June_01_2025_symbol.gmt")  
} else {
  stop(paste("Bader modules version", bader_version, "not found."))
}

# Iterate over cluster solutions
for (nk in 2:nk_max) {
  
  # Identify unique genes in each cluster
  cluster_genes_nk <- df_variants %>% 
    select(cluster = paste0("nk", nk), gene) %>% 
    unite(dup_select, "gene", "cluster", remove = FALSE) %>% 
    filter(!duplicated(dup_select)) %>% 
    select(cluster, gene)
  
  # Iterate over clusters in solution
  for (k in 1:nk) {
    
    message(paste0("On cluster ", nk, "-", k))
    
    cluster_genes_k <- cluster_genes_nk %>% 
      filter(cluster == k) %>% 
      pull(gene)
    
    if (length(cluster_genes_k) == 0) {
      message("\tNo genes in cluster.")
      next
    }
    
    message("\tBuilding gene neighbourhood...")
    neighbourhood_k <- get_gene_neighbourhood(genes = cluster_genes_k,
                                              score = gene_score,
                                              stringdb_version = stringdb_version)
    
    outfile <- paste("cluster_neighbourhood", nk, k, gene_score, sep = "_")
    outfile <- paste0(outfile, ".csv")
    outfile <- file.path(output_dirs[2], outfile)
    write_csv(neighbourhood_k, file = outfile)
    
    # Construct target set for cluster
    target_set <- unique(c(neighbourhood_k[["gene_A"]], 
                           neighbourhood_k[["gene_B"]]))
    
    # Evaluate pathway enrichment for cluster
    message("\tCalculating pathway enrichment...")
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