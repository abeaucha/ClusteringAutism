#!/usr/bin/env Rscript
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

# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")

# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "enrichment.R"))



# string_score <- 950
# nk_max = 2
# permute = FALSE
# seed = NULL
# nk <- 2; k <- 1

compute_cluster_enrichment <- function(
    clusters,
    gene_network,
    background_set,
    modules,
    string_score = 950,
    nk_max = 10,
    output_dir = "data/enrichment/StringDB_12.0_Bader_2025",
    permute = FALSE,
    seed = NULL
) {
  
    # Add STRING score to output dir
    output_dir <- file.path(output_dir, string_score)

    # If permutation, reshuffle cluster assignments
    if (permute) {
        output_dir <- file.path(output_dir, "permutations")

        if (is.null(seed)) {
            seed <- sample(1:1E4, size = 1)
        }

        output_dir <- file.path(output_dir, paste0("seed_", seed))

        set.seed(seed)
        idx_sample <- sample(
            seq_len(nrow(clusters)),
            size = nrow(clusters),
            replace = FALSE
        )

        clusters$model <- clusters$model[idx_sample]
        clusters$gene <- clusters$gene[idx_sample]
    }

    # Construct output directories
    output_dirs <- c(
        file.path(output_dir, "NeighbourhoodInfo"),
        file.path(output_dir, "NeighbourhoodEnrichment")
    )
    for (path in output_dirs) {
        if (!dir.exists(path)) {
            dir.create(path, recursive = TRUE)
        }
    }

    # Iterate over cluster solutions
    for (nk in 2:nk_max) {
        # Identify unique genes in each cluster
        cluster_genes_nk <- clusters |>
            select(cluster = paste0("nk", nk), gene) |>
            unite(dup_select, "gene", "cluster", remove = FALSE) |>
            filter(!duplicated(dup_select)) |>
            select(cluster, gene)

        # Iterate over clusters in solution
        for (k in 1:nk) {
            message(paste0("On cluster ", nk, "-", k))

            # Pull genes for cluster
            cluster_genes_k <- cluster_genes_nk %>%
                filter(cluster == k) %>%
                pull(gene)

            # Extract gene neighbourhood for cluster
            neighbourhood_k <- gene_network %>%
                filter(gene_in %in% cluster_genes_k)

            # Export gene neighbourhood
            if (!permute) {
                outfile <- paste(
                    "cluster_neighbourhood",
                    nk,
                    k,
                    string_score,
                    sep = "_"
                )
                outfile <- paste0(outfile, ".csv")
                outfile <- file.path(output_dir, "NeighbourhoodInfo", outfile)
                write_csv(neighbourhood_k, file = outfile)
            }

            # Construct target set for cluster
            target_set <- unique(c(
                neighbourhood_k[["gene_A"]],
                neighbourhood_k[["gene_B"]]
            ))

            # Evaluate pathway enrichment for cluster
            enrichment <- gene_set_enrichment(
              target = target_set,
              background = background_set,
              modules = modules
            )
            
            enrichment[["nk"]] <- nk
            enrichment[["k"]] <- k

            outfile <- paste(
                "cluster_pathway_enrichment",
                nk,
                k,
                string_score,
                sep = "_"
            )
            outfile <- paste0(outfile, ".csv")
            outfile <- file.path(output_dir, "NeighbourhoodEnrichment", outfile)
            write_csv(enrichment, file = outfile)
        }
    }
}


# Main -----------------------------------------------------------------------

pipeline_dir <- "data/mouse/derivatives/v3/"
params_id <- "107"
nk_max <- 10
background_file <- "data/enrichment/sagittal_gene_table_normalized_filtered.csv"
modules_file <- "data/enrichment/Human_Reactome_June_01_2025_symbol.gmt"
string_score <- 950
permute <- FALSE

# Complete pipeline directory path with parameter set id
pipeline_dir <- file.path(pipeline_dir, params_id)

# Cluster directory path
clusters_dir <- file.path(pipeline_dir, "clusters", "resolution_0.2")

# Import clusters
clusters_file <- file.path(clusters_dir, "clusters.csv")
clusters <- read_csv(file = clusters_file, show_col_types = FALSE)

# Import model names dictionary
models_file <- file.path(pipeline_dir, "model_names.csv")
models <- get_model_genes(models_file)

# Join with model gene names (and filter for single genes)
clusters <- clusters %>%
    rename(file = ID) %>%
    right_join(models, by = "file") %>%
    select(-file) %>%
    rename(model = ID)

# Import Bader Reactome pathway modules
modules <- importMsigDBGMT(modules_file)

# Construct background set of genes
message("Building background gene set ...")
background_set <- read_csv(file = background_file, show_col_types = FALSE) |>
  pull(msg.genes.acronym) |>
  as.character() %>% 
  unique() |> 
  fetch_gene_homologs(species = "mouse") |> 
  filter(!is.na(GeneHuman)) |> 
  pull(GeneHuman)

# Compute gene interaction network for all genes in clusters
message("Building gene interaction network ...")
gene_network <- build_gene_interaction_network(
  genes = unique(clusters$gene),
  score = string_score
)

# Extract a vector of all possible target genes
# target_genes_all <- unique(c(gene_network[["gene_A"]], gene_network[["gene_B"]]))

# Build data frame of human homologs for all target genes
gene_network_hom <- fetch_gene_homologs(
  genes = unique(c(gene_network[["gene_A"]], 
                   gene_network[["gene_B"]])), 
  species = "mouse"
)

# Get human gene homologues for gene network
gene_network <- gene_network |> 
  left_join(select(gene_network_hom, GeneInput, GeneHuman),
            by = c("gene_A" = "GeneInput"), 
            relationship = "many-to-many") |> 
  select(-gene_A) |> 
  rename(gene_A = GeneHuman) |> 
  left_join(select(gene_network_hom, GeneInput, GeneHuman),
            by = c("gene_B" = "GeneInput"),
            relationship = "many-to-many") |> 
  select(-gene_B) |> 
  rename(gene_B = GeneHuman)


# Define output directory
# output_dir <- file.path(pipeline_dir, "enrichment", "StringDB_12.0_Bader_2025")

output_dir <- file.path("tmp", "enrichment", "StringDB_12.0_Bader_2025")

# for permutation sample in number of samples
# Run cluster enrichment

permute <- FALSE
seed <- NULL

compute_cluster_enrichment(
    clusters = clusters,
    gene_network = gene_network,
    background_set = background_set,
    modules = modules,
    string_score = string_score,
    nk_max = nk_max,
    output_dir = output_dir,
    permute = TRUE,
    seed = seed
)

# Questions
# Should background_set be the file or the vector?
# Should modules be the file or the list?
# 
# df_true <- read_csv("tmp/enrichment/StringDB_12.0_Bader_2025/950/NeighbourhoodEnrichment/cluster_pathway_enrichment_6_1_950.csv")
# df_permute <- read_csv("tmp/enrichment/StringDB_12.0_Bader_2025/950/permutations/seed_5555/NeighbourhoodEnrichment/cluster_pathway_enrichment_6_1_950.csv")
# 
#  
# colnames(df_true)
# 
# true_sig <- df_true |> 
#   filter(adj.P.Val < 0.01) |> 
#   pull(Title)
# 
# permute_sig <- df_permute |> 
#   filter(adj.P.Val < 0.01) |> 
#   pull(Title)
# 
# length(true_sig)
# length(permute_sig)
# 
# length(intersect(true_sig, permute_sig))
