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
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(MRIcrotome))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggplotify))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--gene-space",
              type = "character",
              default = "average_latent_space",
              help = paste("[default %default]")),
  make_option("--jacobians",
              type = "character",
              default = "absolute",
              help = paste("[default %default]")),
  make_option("--threshold",
              type = "numeric",
              default = 0.5,
              help = paste("[default %default]")),
  make_option("--metric",
              type = "character",
              default = "correlation",
              help = paste("[default %default]")),
  make_option("--latent-space-id",
              type = "numeric",
              default = 100,
              help = paste("[default %default]")),
  make_option("--clustering",
              type = "character",
              default = "false",
              help = paste("[default %default]")),
  make_option("--kcut",
              type = "numeric",
              default = 8,
              help = paste("[default %default]")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbosity option. [default %default]"))
) 


# Functions ------------------------------------------------------------------

source("../functions/buildSimilarityMatrix.R")
source("../functions/tree_tools.R")

#' Process cluster signatures
#'
#' @param signatures (character scalar) Path to the .csv file containing cluster expression signatures
#' @param expr (character scalar) Name of the .csv containing the expression data using to computed the cluster signatures
#'
#' @return (matrix) Expression signatures for each cluster
process_signatures <- function(signatures, expr){
  mat_expr <- suppressMessages(read_csv(signatures)) %>% 
    filter(str_detect(exprfile, expr)) %>% 
    mutate(nk = clusterfile %>% 
             str_extract('Clusternum_[0-9]+') %>% 
             str_extract('[0-9]+') %>% 
             as.integer(),
           k = clusterfile %>% 
             str_extract('Group_[0-9]+') %>% 
             str_extract('[0-9]+') %>% 
             as.integer()) %>% 
    arrange(nk, k) %>% 
    select(-exprfile, -clusterfile) %>% 
    unite(col = 'cluster_id', nk, k, sep = '-') %>% 
    column_to_rownames('cluster_id') %>% 
    as.matrix() %>% 
    t() 
  return(mat_expr)
}

generate_heatmap <- function(x, palette) {
  
  df_annotations <- tibble(cluster_ids = colnames(x)) %>% 
    separate(cluster_ids, into = c("nk", "k"), remove = FALSE) %>% 
    mutate(nk = factor(nk, levels = 1:10),
           k = factor(k, levels = 1:10)) %>% 
    column_to_rownames(var = "cluster_ids")
  
  df_colour_palette <- tibble(i = factor(1:10),
                              colour = palette)
  
  df_annotation_colours <- df_annotations %>% 
    left_join(df_colour_palette, by = c("nk" = "i")) %>% 
    rename(nk_col = colour) %>% 
    left_join(df_colour_palette, by = c("k" = "i")) %>% 
    rename(k_col = colour) 
  
  nk_colours <- df_annotation_colours$nk_col
  names(nk_colours) <- df_annotation_colours$nk
  
  k_colours <- df_annotation_colours$k_col
  names(k_colours) <- df_annotation_colours$k
  
  annotation_colours <- list(nk = nk_colours,
                             k = k_colours)
  
  p <- pheatmap(mat = x, 
                cluster_cols = F, cluster_rows = F, 
                annotation_row = df_annotations,
                annotation_col = df_annotations,
                annotation_colors = annotation_colours, 
                silent = TRUE, 
                na_col = "black")
  
  return(p)
}


#' Remove empty matrix cells
#'
#' @param x (matrix) A matrix with cells containing NA values
#'
#' @return (matrix) A matrix without any NA values
remove_empty_cells <- function(x) {
  
  df_empty <- as_tibble(which(is.na(x), arr.ind = TRUE))
  
  empty_rows <- df_empty %>% 
    group_by(row) %>% 
    count() %>% 
    filter(n == ncol(x)) %>% 
    pull(row)
  
  empty_cols <- df_empty %>% 
    group_by(col) %>% 
    count() %>% 
    filter(n == nrow(x)) %>% 
    pull(col)
  
  rows_not_empty <- !(1:nrow(x) %in% empty_rows)
  cols_not_empty <- !(1:ncol(x) %in% empty_cols)
  
  x <- x[rows_not_empty, cols_not_empty]
  
  return(x)
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
gene_space <- args[["gene-space"]]
jacobians <- args[["jacobians"]]
threshold <- args[["threshold"]]
metric <- args[["metric"]]
clustering <- ifelse(args[['clustering']] == 'true', TRUE, FALSE)
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

if (jacobians == "absolute") {
  jacobians <- "abs"
} else if (jacobians == "relative") {
  jacobians <- "rel"
} else {
  stop()
}

if (gene_space == "input_space") {
  
  mouse_file <- str_c("../data/mouse/cluster_signatures/input_space/mouse_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_inputspace.csv")
  human_file <- str_c("../data/human/cluster_signatures/input_space/human_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_inputspace.csv")
  
  mouse_expr_file <- "MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv"
  human_expr_file <- "data/human/expression/input_space//HumanExpressionMatrix_samples_pipeline_abagen_homologs_scaled.csv"
  
  mat_mouse <- process_signatures(signatures = mouse_file, expr = mouse_expr_file)
  mat_human <- process_signatures(signatures = human_file, expr = human_expr_file)
  
  mat_sim <- buildSimilarityMatrix(x1 = mat_human,
                                   x2 = mat_mouse,
                                   method = metric)
  
} else if (gene_space == "latent_space") {
  
  mouse_file <- str_c("../data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_latentspace100.csv")
  human_file <- str_c("../data/human/cluster_signatures/latent_space_100/human_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_latentspace100.csv")
  
  latent_space_id <- args[["latent-space-id"]]
  
  mouse_expr_file <- str_c("MLP_labels67_layers3_units200_L20.0_mousetransform_", latent_space_id, ".csv")
  human_expr_file <- str_c("MLP_labels67_layers3_units200_L20.0_humantransform_", latent_space_id, ".csv")
  
  mat_mouse <- process_signatures(signatures = mouse_file, expr = mouse_expr_file)
  mat_human <- process_signatures(signatures = human_file, expr = human_expr_file)
  
  mat_sim <- buildSimilarityMatrix(x1 = mat_human,
                                   x2 = mat_mouse,
                                   method = metric)
  
} else if (gene_space == "average_latent_space") {
  
  sim_file <- str_c("../data/similarity/latent_space_100/similarity_hm_", jacobians, "_mean_threshold", threshold, "_latentspace100.csv")
  mat_sim <- read_csv(file = sim_file,
                      show_col_types = FALSE) %>% 
    column_to_rownames('cluster_id') %>% 
    as.matrix()
  
} else {
  stop()
}

if (clustering) {
  
  kcut <- args[["kcut"]]
  
  mat_sim_complete <- remove_empty_cells(x = mat_sim)
  
  similarity_heatmap_clustered <- pheatmap(mat = mat_sim_complete, 
                                           cutree_rows = kcut, 
                                           cutree_cols = kcut)
  
  human_labels <- sort(cutree(similarity_heatmap_clustered$tree_row, k = kcut))
  mouse_labels <- sort(cutree(similarity_heatmap_clustered$tree_col, k = kcut))
  
  annotation_row <- tibble(meta_k = human_labels,
                           cluster_id = names(human_labels)) %>% 
    column_to_rownames("cluster_id") %>% 
    mutate(meta_k = factor(meta_k))
  
  annotation_col <- tibble(meta_k = mouse_labels,
                           cluster_id = names(mouse_labels)) %>% 
    column_to_rownames("cluster_id") %>% 
    mutate(meta_k = factor(meta_k))
  
  similarity_heatmap <- as.ggplot(pheatmap(mat_sim_complete,
                                           cutree_rows = kcut, cutree_cols = kcut,
                                           annotation_row = annotation_row,
                                           annotation_col = annotation_col, 
                                           annotation_names_row = F,
                                           annotation_names_col = F, 
                                           silent = T))
  
} else {
  
  palette <- mako(n = 10, direction = -1, begin = 0.3)
  similarity_heatmap <- as.ggplot(generate_heatmap(x = mat_sim, 
                                                   palette = palette)) 
}

padding <- 0.2

similarity_heatmap <- similarity_heatmap + 
  labs(title = "Cluster similarity",
       subtitle = str_c("Gene space: ", gene_space, "\n",
                        "Threshold: ", threshold, "\n",
                        "Jacobians: ", jacobians, "\n",
                        "Metric: ", metric, "\n",
                        "Rows: Human", "\n",
                        "Columns: Mouse", "\n")) +
  theme(plot.margin = margin(t = padding, 
                             r = padding, 
                             b = padding, 
                             l = padding, 
                             unit = "in"))

outfile <- str_c("ClusterSimilarity", gene_space, sep = "_")

if (gene_space == "latent_space") {
  outfile <- str_c(outfile, latent_space_id, sep = "_")  
}

outfile <- str_c(outfile, jacobians, threshold, metric, sep = "_")

if (clustering) {
  outfile <- str_c(outfile, "clustered", sep = "_")
}

outfile <- str_c(outfile, ".pdf")

pdf(file = outfile,
    width = unit(10, "inch"),
    height = unit(10, "inch"))

similarity_heatmap

dev.off()
