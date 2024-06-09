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
# suppressPackageStartupMessages(library(magrittr))
# suppressPackageStartupMessages(library(glue))
# suppressPackageStartupMessages(library(tmod))


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


#' Filter mouse models
#'
#' @param clusters (character scalar)
#' @param models (character scalar)
#'
#' @return
filter_models <- function(clusters, models) {
  
  # Genes to rename
  genes_rename <- tibble(Gene = c("Andr", "Caspr2", "Dat", "Mor", "Nl1", "Nl3",
                                  "Nrxn1a", "Sert", "Pcdh", "Chd7;En1Cre",
                                  "Ube3a.2", "FusDelta14", "Nr1a", "Snf2L", "Snf2H"),
                         GeneNew = c("Ar", "Cntnap2", "Slc6a3", "Oprm1", "Nlgn1",
                                     "Nlgn3", "Nrxn1", "Slc6a4", "Pcdhga3", "Chd7",
                                     "Ube3a", "Fus", "Nmdar1", "Smarca1", "Smarca5"))
  
  # Genes to exclude
  genes_exclude <- c("15q11-13", "16p11.2", "22q11.2", "XO", "Btbr", "Balbc", 
                     "MAR", "15q25", "TCDD", "VPA", "BtbrTT")
  
  # Import models 
  models <- read_csv(models, show_col_types = FALSE)
  
  # Import clusters
  clusters <- read_csv(clusters, show_col_types = FALSE) %>% 
    rename(file = ID) %>% 
    left_join(models, by = "file") %>% 
    select(-file) %>% 
    rename(Model = ID)
  
  # Create genes columns
  clusters <- clusters %>% 
    mutate(Gene = Model %>%
             str_split("\\(") %>% 
             map_chr(.f = function(x){x[[1]]}))
  
  # Clean up gene names
  clusters <- clusters %>% 
    left_join(genes_rename, by = "Gene") %>% 
    mutate(GeneNew = ifelse(is.na(GeneNew), Gene, GeneNew),
           GeneNew = ifelse(Model == "itsn1(+/+);itsn2(-/-)", "itsn2", GeneNew),
           GeneNew = ifelse(Model == "Snf2H(+/+);Snf2L(-/-);emxcre", "Snf2l", GeneNew),
           GeneNew = ifelse(Model == "Gsk3(a)", "Gsk3A", GeneNew),
           GeneNew = ifelse(Model == "Gsk3(B)", "Gsk3B", GeneNew)) %>% 
    filter(!(GeneNew %in% genes_exclude)) %>% 
    select(-Gene) %>% 
    rename(Gene = GeneNew)  
  
  return(clusters)
  
}


# Main -----------------------------------------------------------------------

#Parse command line args
# args <- parse_args(OptionParser(option_list = option_list))
# arg <- args[["arg"]]
# verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# if (verbose) {message("")}

# Directory containing mouse cluster files
# Cluster_Dir <- "/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters_Paper/"
cluster_dir <- file.path(PROJECTPATH, "data/mouse/derivatives/v3/107/clusters/resolution_0.2/")

# Command line args? 
pipeline_dir <- file.path(PROJECTPATH, "data/mouse/derivatives/v3/")
params_id <- 107
models <- "model_names.csv"

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

# Enrichment output directory
enrichment_dir <- file.path(PROJECTPATH, "data/mouse/enrichment/")
if (!dir.exists(enrichment_dir)) {
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
}

# Model names dictionary
models <- file.path(pipeline_dir, models)

# Clusters
clusters <- file.path(cluster_dir, "clusters.csv")

# Filter models in clusters
clusters <- filter_models(clusters = clusters, models = models)


# Parameters for enrichment pipeline

# Max number of clusters
nk_max <- 10

# Gene score for StringDB
gene_scores <- seq(400, 950, by = 50)

# Version of StringDB
# stringdb_versions <- c("11.5", "12.0")
stringdb_versions <- "12.0"

# Version of Bader pathways
# bader_versions <- c(2020, 2023)
bader_versions <- 2023

# Background set for enrichment analysis
background_set <- file.path(enrichment_dir, "sagittal_gene_table_normalized_filtered.csv")

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
  
  output_dir <- paste("StringDB", stringdb_version, "Bader", 
                      bader_version, sep = "_")
  output_dir <- file.path(enrichment_dir, output_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  if (bader_version == 2020) {
    bader_list <- file.path(enrichment_dir, "Human_Reactome_March_01_2020_symbol.gmt")
  } else if (bader_version == 2023) {
    bader_list <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  
  } else {
    stop()
  }
  
  

  # Function new version of GetGeneNeighbourhood
    
  GeneScore <- gene_score
  total_clusters <- nk_max
  
  # StringDB version URL
  if (stringdb_version == "11.5") {
    stringdb_url <- "https://version-11-5.string-db.org/"
  } else if (stringdb_version == "12.0") {
    stringdb_url <- "https://string-db.org/"
  } else {
    stop()
  }
  
  # for (nk in 2:nk_max) {}
  nk <- 2
  
  all_df <- clusters %>% 
    select(Model, gene = Gene, cluster = paste0("nk", nk))
  
  # Remove duplicates where multiple genes are in the same cluster
  all_df <- all_df %>% 
    unite(dup_select, "gene", "cluster", remove = FALSE) %>% 
    mutate(is_dup = duplicated(dup_select)) %>% 
    filter(!is_dup) %>% 
    select(Model, gene, cluster)
  
  all_df <- all_df %>% 
    mutate(queryIndex = 1:nrow(.)-1)

  # StringDB API query to get gene IDs
  api_query_string_ids <- paste0(stringdb_url, "api/tsv/get_string_ids?identifiers=", 
                                 paste(all_df[["gene"]], collapse = "%0D"), 
                                 "&species=10090&limit=1")
  
  # Fetch gene IDs
  string_ids <- read_tsv(api_query_string_ids, show_col_types = FALSE) %>% 
    select(queryIndex, stringID = stringId, preferredName) 
  
  # Join String IDs to clusters
  all_df <- all_df %>% 
    left_join(string_ids, by = "queryIndex")
  
  # Remove organism ID
  all_df <- all_df %>% 
    mutate(fixed_id = stringID %>% 
             strsplit("\\.") %>% 
             map_chr(`[`, 2))
  
  # StringDB API query for gene neighbourhoods
  api_query_interactions <- paste0(stringdb_url, "api/tsv/interaction_partners?identifiers=", 
                                   paste(all_df[["stringID"]], collapse = "%0D"),
                                   "&species=10090&required_score=", gene_score)
  
  # Fetch gene neighbourhoods
  df_interactions_k <- read_tsv(api_query_interactions, show_col_types = FALSE) %>% 
    rename(stringID_A = stringId_A, stringID_B = stringId_B) %>% 
    left_join(all_df %>% 
              select(cluster, stringID),
              by = c("stringID_A" = "stringID"))  
  

  
    base_set <- unique(c(gx_data$preferredName_A, gx_data$preferredName_B))
  
  GetGeneNeighbourhood(GeneScore = gene_score, 
                       total_clusters = nk_max, 
                       stringdb_version = stringdb_version,
                       output_dir = file.path(output_dir, "NeighbourhoodInfo"))
  
  GetNeighbourhoodEnrichment(GeneScore = gene_score,
                             total_clusters = nk_max,
                             Bader_List = bader_list,
                             background_set = background_set,
                             Neighbour_dir = file.path(output_dir, "NeighbourhoodInfo"),
                             output_dir = file.path(output_dir, "NeighbourhoodEnrichment"))
  
}



