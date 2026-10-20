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


compute_cluster_enrichment <- function(
    clusters,
    neighbourhoods,
    background_set,
    modules,
    string_score = 950,
    nk_max = 10,
    output_dir = "data/enrichment/StringDB_12.0_Bader_2025",
    permute = FALSE,
    seed = NULL
) {
    output_dir <- file.path(output_dir, string_score)

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

            cluster_genes_k <- cluster_genes_nk %>%
                filter(cluster == k) %>%
                pull(gene)

            # Construct gene neighbourhood for cluster
            message("\tBuilding gene neighbourhood...")
            # neighbourhood_k <- get_gene_neighbourhood(
            #     genes = cluster_genes_k,
            #     score = string_score
            # )

            # Extract gene neighbourhood for cluster
            neighbourhood_k <- neighbourhoods %>%
                filter(gene_in %in% cluster_genes_k)

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
            message("\tCalculating pathway enrichment...")
            enrichment <- get_neighbourhood_enrichment(
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
nk_max <- 2
background_file <- "data/enrichment/sagittal_gene_table_normalized_filtered.csv"
modules_file <- "data/enrichment/Human_Reactome_June_01_2025_symbol.gmt"
string_score <- 950

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

# Import background set of genes
background_genes_mouse <- read_csv(file = background_file, show_col_types = FALSE) |>
    pull(msg.genes.acronym) |>
    as.character()

background_genes <- get_homologs(genes = background_genes_mouse,
                                 species = "mouse")[["human_genes"]]

# Import Bader Reactome pathway modules
modules <- importMsigDBGMT(modules_file)

# Compute gene interaction network for all genes in clusters
neighbourhoods <- get_gene_neighbourhood(
    genes = unique(clusters$gene),
    score = string_score
)

neighbourhoods_all_genes <- unique(c(
  neighbourhoods[["gene_A"]],
  neighbourhoods[["gene_B"]]
))

genes = unique(c(neighbourhoods[["gene_A"]], 
                 neighbourhoods[["gene_B"]]))
species <- "mouse"

list_neighbourhood_homologs <- get_homologs(genes = unique(c(neighbourhoods[["gene_A"]], 
                                                             neighbourhoods[["gene_B"]])),
                                            species = "mouse")

list_neighbourhood_homologs$input_genes %>% length()
list_neighbourhood_homologs$mouse_genes %>% length()
list_neighbourhood_homologs$human_genes %>% length()

# Define output directory
# output_dir <- file.path(pipeline_dir, "enrichment", "StringDB_12.0_Bader_2025")

output_dir <- file.path("tmp", "enrichment", "StringDB_12.0_Bader_2025")

# for permutation sample in number of samples
# Run cluster enrichment

permute <- FALSE
seed <- NULL

compute_cluster_enrichment(
    clusters = clusters,
    neighbourhoods = neighbourhoods,
    background_set = background_genes,
    modules = modules,
    string_score = string_score,
    nk_max = nk_max,
    output_dir = output_dir,
    permute = permute,
    seed = seed
)

# Questions
# Should background_set be the file or the vector?
# Should modules be the file or the list?
