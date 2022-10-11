#' Import and process cluster signatures
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


import_similarity_matrix <- function(gene_space = "average_latent_space", 
                                     jacobians = "abs", 
                                     threshold_method = "intensity", 
                                     threshold = 0.5, mask_type = "symmetric", 
                                     latent_space_id = 1) {
  if (gene_space == "input_space") {
    
    # mouse_file <- str_c("../../data/mouse/cluster_signatures/input_space/mouse_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_inputspace.csv")
    # human_file <- str_c("../../data/human/cluster_signatures/input_space/human_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_inputspace.csv")
    # 
    # mouse_expr_file <- "MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv"
    # human_expr_file <- "data/human/expression/input_space//HumanExpressionMatrix_samples_pipeline_abagen_homologs_scaled.csv"
    # 
    # mat_mouse <- process_signatures(signatures = mouse_file, expr = mouse_expr_file)
    # mat_human <- process_signatures(signatures = human_file, expr = human_expr_file)
    # 
    # mat_sim <- buildSimilarityMatrix(x1 = mat_human,
    #                                  x2 = mat_mouse,
    #                                  method = metric)
    
  } else if (gene_space == "latent_space") {
    
    # mouse_file <- str_c("../../data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_latentspace100.csv")
    # human_file <- str_c("../../data/human/cluster_signatures/latent_space_100/human_cluster_signatures_", jacobians, "_mean_threshold", threshold, "_latentspace100.csv")
    # 
    # latent_space_id <- args[["latent-space-id"]]
    # 
    # mouse_expr_file <- str_c("MLP_labels67_layers3_units200_L20.0_mousetransform_", latent_space_id, ".csv")
    # human_expr_file <- str_c("MLP_labels67_layers3_units200_L20.0_humantransform_", latent_space_id, ".csv")
    # 
    # mat_mouse <- process_signatures(signatures = mouse_file, expr = mouse_expr_file)
    # mat_human <- process_signatures(signatures = human_file, expr = human_expr_file)
    # 
    # mat_sim <- buildSimilarityMatrix(x1 = mat_human,
    #                                  x2 = mat_mouse,
    #                                  method = metric)
    
  } else if (gene_space == "average_latent_space") {
    
    sim_file <- str_c("../../data/similarity_matrix/latent_space_100/similarity_hm_", jacobians, 
                      "_mean_threshold_", threshold_method, "_", threshold, "_", mask_type, "_latentspace100.csv")
    mat_sim <- read_csv(file = sim_file,
                        show_col_types = FALSE) %>% 
      column_to_rownames("cluster_id") %>% 
      as.matrix()
    
  } else {
    stop()
  }
  
  return(mat_sim)
  
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

#' Generate a set of meta-clusters
#'
#' @param x (matrix) The similarity matrix to cluster
#' @param kcut (numeric scalar) The number of clusters at which to cut
#' the dendrogram
#' @param axis (character scalar) One of "columns" or "rows" indicating
#' which axis to cluster
#'
#' @return (data.frame)
generate_meta_clusters <- function(x, kcut, axis = "columns") {
  
  x <- remove_empty_cells(x = x)
  clustered_heatmap <- pheatmap(mat = x, silent = TRUE)
  
  if (axis == "columns") {
    field <- "tree_col"
  } else if (axis == "rows") {
    field <- "tree_row"
  }
  meta_cluster_labels <- sort(cutree(clustered_heatmap[[field]], k = kcut))
  
  df_meta_clusters <- tibble(meta_k = meta_cluster_labels,
                             cluster_id = names(meta_cluster_labels)) %>% 
    separate(col = "cluster_id", 
             into = c("nk", "k"),
             sep = "-", 
             remove = FALSE) %>% 
    mutate(nk = as.integer(nk),
           k = as.integer(k)) %>% 
    arrange(meta_k, nk, k)
  
  return(df_meta_clusters)
  
}


#' Create reduced human atlas definitions
#'
#' @param tree (data.tree) Tree whose leaves will be new labels
#' @param defs (data.frame) Definitions containing labels for all human microarray samples
#' @param simplify (logical scalar) Option to simplify returned data frame
#'
#' @return (data.frame) Reduced atlas label definitions
reduce_human_defs <- function(tree, defs, simplify = FALSE) {
  
  defs_reduced <- tree$Get(attribute = "samples", 
                           filterFun = isLeaf, 
                           simplify = FALSE) %>% 
    map_dfr(.f = function(x){tibble(sample_id = x)},
            .id = "name") %>% 
    inner_join(defs,
               by = "sample_id") %>% 
    rename(sample_label = label) %>% 
    group_by(name) %>% 
    mutate(label = min(sample_label)) %>% 
    ungroup()
  
  if(simplify) {
    defs_reduced <- defs_reduced %>% 
      select(name, label) %>% 
      distinct()
  }
  
  return(defs_reduced)
  
}

#' Create reduced human atlas labels
#'
#' @param tree 
#' @param labels 
#' @param defs 
#'
#' @return (mincSingleDim) Reduced labels
reduce_human_labels <- function(tree, labels, defs) {
  
  defs_reduced <- reduce_human_defs(tree = tree,
                                    defs = defs,
                                    simplify = FALSE)
  
  ind_match <- match(x = labels, 
                     table = defs_reduced[["sample_label"]])
  labels_reduced <- defs_reduced[ind_match,"label"][[1]]
  labels_reduced[is.na(labels_reduced)] <- 0
  attributes(labels_reduced) <- attributes(labels)
  
  return(labels_reduced)
  
}

#' Calculate cluster region fractions
#'
#' @param cluster (mincSingleDim)
#' @param labels (mincSingleDim)
#' @param defs (data.frame)
#' @param mask_type (character scalar)
#'
#' @return (data.frame)
calc_cluster_region_fractions <- function(cluster, labels, defs, mask_type = "symmetric") {
  # if (mask_type == "positive") {
  #   cluster[cluster < 0] <- 0
  # } else if (mask_type == "negative") {
  #   cluster[cluster > 0] <- 0
  #   cluster <- abs(cluster)
  # } else if (mask_type == "symmetric") {
  #   cluster <- abs(cluster)
  # }
  voxels_cluster <- cluster == 1
  voxels_nonzero <- labels != 0
  labels_cluster <- labels[voxels_cluster & voxels_nonzero]
  if (length(labels_cluster) == 0){
    stop(paste("There are no voxels in the cluster with mask_type:",
               mask_type))
  } else {
    out <- table(labels_cluster) %>% 
      as_tibble() %>% 
      rename(label = labels_cluster,
             nvoxels = n) %>% 
      mutate(label = as.integer(label)) %>% 
      right_join(defs, by = "label") %>%
      mutate(nvoxels = ifelse(is.na(nvoxels), 0, nvoxels),
             nvoxels_cluster = sum(voxels_cluster),
             nvoxels_expr = length(labels_cluster),
             fvoxels_cluster = nvoxels/nvoxels_cluster,
             fvoxels_expr = nvoxels/nvoxels_expr) %>% 
      select(name,
             label,
             nvoxels,
             nvoxels_cluster,
             nvoxels_expr,
             fvoxels_cluster,
             fvoxels_expr)
  }
  return(out)
}


#' Import a set of cluster masks
#'
#' @param files (character vector)
#' @param df (data.frame)
#'
#' @return (list)
import_cluster_masks <- function(files, df) {
  masks <- vector(mode = "list", length = nrow(df))
  names(masks) <- df[["cluster_id"]]
  for (i in 1:length(masks)) {
    nk <- df[[i,"nk"]]
    k <- df[[i, "k"]]
    infile <- files %>% 
      str_subset(str_c("Clusternum_", nk)) %>%
      str_subset(str_c("Group_", k))
    masks[[i]] <- mincGetVolume(infile)
    masks[[i]] <- floor(masks[[i]])
  }
  return(masks)
}



prepare_spider_data <- function(mouse, human, homologues) {
  ind_match <- match(human[["name"]], table = homologues[["human"]])
  human[["name"]] <- homologues[["mouse"]][ind_match]
  human <- human %>% 
    mutate(name = factor(name, levels = homologues[["mouse"]]),
           species = "human")
  mouse <- mouse %>% 
    mutate(name = factor(name, levels = homologues[["mouse"]]),
           species = "mouse")
  out <- bind_rows(human, mouse)
  out <- out %>% 
    unite(col = group_id, 
          species, cluster_id, 
          sep = "-", remove = FALSE)
  return(out)
}

plot_heatmap <- function(x, clustering = FALSE, cutree_rows = 1, cutree_cols = 1, 
                         gene_space, jacobians, mask_type, threshold_method,
                         threshold, metric, header = TRUE, padding = 0.2, 
                         save = FALSE, outdir = "plots/", 
                         fig_width = 10, fig_height = 10) {
  
  palette <- mako(n = 10, direction = -1, begin = 0.3)
  
  if (clustering) {
    
    annotation_row <- generate_meta_clusters(x = x, 
                                             kcut = cutree_rows, 
                                             axis = "rows") %>% 
      select(cluster_id, meta_k) %>% 
      column_to_rownames("cluster_id") %>% 
      mutate(meta_k = factor(meta_k))
    
    annotation_col <- generate_meta_clusters(x = x, 
                                             kcut = cutree_cols, 
                                             axis = "columns") %>% 
      select(cluster_id, meta_k = meta_k) %>% 
      column_to_rownames("cluster_id") %>% 
      mutate(meta_k = factor(meta_k))
    
    if (cutree_rows != cutree_cols) {
      annotation_row <- annotation_row %>% 
        rename(meta_k_row = meta_k)
      annotation_col <- annotation_col %>% 
        rename(meta_k_col = meta_k)
    }
    
    x <- remove_empty_cells(x = x)
    p_heatmap <- as.ggplot(pheatmap(mat = x,
                                    cutree_rows = cutree_rows, 
                                    cutree_cols = cutree_cols,
                                    annotation_row = annotation_row,
                                    annotation_col = annotation_col, 
                                    annotation_names_row = F,
                                    annotation_names_col = F, 
                                    silent = TRUE))
    
  } else {
    
    df_annotations <- tibble(cluster_ids = colnames(x)) %>% 
      separate(cluster_ids, into = c('nk', 'k'), remove = FALSE) %>% 
      mutate(nk = factor(nk, levels = 1:10),
             k = factor(k, levels = 1:10)) %>% 
      column_to_rownames(var = 'cluster_ids')
    
    df_colour_palette <- tibble(i = factor(1:10),
                                colour = palette)
    
    df_annotation_colours <- df_annotations %>% 
      left_join(df_colour_palette, by = c('nk' = 'i')) %>% 
      rename(nk_col = colour) %>% 
      left_join(df_colour_palette, by = c('k' = 'i')) %>% 
      rename(k_col = colour) 
    
    nk_colours <- df_annotation_colours$nk_col
    names(nk_colours) <- df_annotation_colours$nk
    
    k_colours <- df_annotation_colours$k_col
    names(k_colours) <- df_annotation_colours$k
    
    annotation_colours <- list(nk = nk_colours,
                               k = k_colours)
    
    p_heatmap <- as.ggplot(pheatmap(mat = x, 
                                    cluster_cols = F, cluster_rows = F, 
                                    annotation_row = df_annotations,
                                    annotation_col = df_annotations,
                                    annotation_colors = annotation_colours, 
                                    silent = TRUE, 
                                    na_col = 'black'))
    
  }
  
  if (header) {
    p_heatmap <- p_heatmap +
      labs(title = "Cluster similarity",
           subtitle = str_c("Gene space: ", gene_space, "\n",
                            "Jacobians: ", jacobians, "\n",
                            "Mask type: ", mask_type, "\n",
                            "Thresholding method: ", threshold_method, "\n",
                            "Threshold: ", threshold, "\n",
                            "Metric: ", metric, "\n",
                            "Rows: Human", "\n",
                            "Columns: Mouse", "\n"))
      
  }
    
  p_heatmap <- p_heatmap +
    theme(plot.margin = margin(t = padding, 
                               r = padding, 
                               b = padding, 
                               l = padding, 
                               unit = "in"))
  
  if (save) {
    outfile <- str_c("ClusterSimilarity", 
                     gene_space, 
                     "jacobians", jacobians,
                     "threshold", threshold_method,
                     threshold,
                     mask_type,
                     sep = "_")
    if (clustering) {
      outfile <- str_c(outfile,
                       "clustered",
                       "rowcut", cutree_rows,
                       "colcut", cutree_cols,
                       sep = "_")
    }
    outfile <- str_c(outfile, ".pdf")
    outfile <- file.path(outdir, outfile)
    pdf(file = outfile,
        width = unit(fig_width, "inch"),
        height = unit(fig_height, "inch"))
    print(p_heatmap)
    dev.off()
  }
  
  return(p_heatmap)
  
}

plot_cluster_maps <- function(x, species, kcut, meta_k, gene_space, jacobians,
                              mask_type, threshold_method, threshold, 
                              save = FALSE, outdir = "plots/",
                              fig_width = 10, fig_height = 10) {
  
  if (species == "mouse") {
    axis <- "columns"
    anat_file <- "../../data/mouse/atlas/DSURQE_CCFv3_average_200um.mnc"
    cluster_dir <- str_c("../../data/mouse/clustering/cluster_maps/", 
                         jacobians, "/resolution_200/mean/")
    cluster_mask_dir <- str_c("../../data/mouse/clustering/cluster_masks/", 
                              jacobians, "/resolution_200/mean/threshold_", 
                              threshold_method, "/symmetric/threshold_", threshold)
    slice_begin <- 10
    slice_end <- 60
    anat_low <- 700
    anat_high <- 1400 
  } else if (species == "human") {
    axis <- "rows"
    anat_file <- "../../data/human/registration/reference_files/model_1.0mm.mnc"
    cluster_dir <- str_c("../../data/human/clustering/cluster_maps/", jacobians,
                         "/resolution_1.0/mean/")
    cluster_mask_dir <- str_c("../../data/human/clustering/cluster_masks/", jacobians,
                              "/resolution_1.0/mean/threshold_", threshold_method,
                              "/symmetric/threshold_", threshold)
    slice_begin <- 50
    slice_end <- 200
    anat_low <- 3
    anat_high <- 7
  } else {
    stop("species must be one of {mouse, human}")
  }
  
  df_clusters <- generate_meta_clusters(x = x, kcut = kcut, axis = axis)
  df_clusters <- df_clusters[df_clusters[["meta_k"]] == meta_k,]
  
  anat_vol <- mincArray(mincGetVolume(anat_file))
  
  #Cluster maps and masks
  cluster_map_files <- list.files(cluster_dir, full.names = TRUE)
  cluster_mask_files <- list.files(cluster_mask_dir, full.names = TRUE)
  
  #Import cluster maps in meta-cluster
  meta_cluster_files <- character(nrow(df_clusters))
  cluster_maps <- vector(mode = "list", length = nrow(df_clusters))
  names(cluster_maps) <- df_clusters$cluster_id
  for (i in 1:length(cluster_maps)) {
    nk <- df_clusters[[i, "nk"]]
    k <- df_clusters[[i, "k"]]
    map_infile <- cluster_map_files %>%
      str_subset(str_c("Clusternum_", nk)) %>%
      str_subset(str_c("Group_", k))
    mask_infile <- cluster_mask_files %>%
      str_subset(str_c("Clusternum_", nk)) %>%
      str_subset(str_c("Group_", k))
    cluster_maps[[i]] <- mincGetVolume(map_infile)
    cluster_mask <- floor(mincGetVolume(mask_infile))
    if (mask_type == "positive") {
      cluster_mask[cluster_mask < 0] <- 0
    } else if (mask_type == "negative") {
      cluster_mask[cluster_mask > 0] <- 0
      cluster_mask <- abs(cluster_mask)
    } else if (mask_type == "symmetric"){
      cluster_mask <- abs(cluster_mask)
    } else {
      stop("mask_type must be one of {'symmetric', 'positive', 'negative'}")
    }
    cluster_maps[[i]][cluster_mask == 0] <- 0
  }
  
  for (i in 1:length(cluster_maps)) {
    overlay_low_test <- min(abs(cluster_maps[[i]][cluster_maps[[i]] != 0]))
    overlay_high_test <- max(abs(cluster_maps[[i]][cluster_maps[[i]] != 0]))
    if (i == 1) {
      overlay_low <- overlay_low_test
      overlay_high <- overlay_high_test
    } else {
      if (overlay_low_test < overlay_low) {
        overlay_low <- overlay_low_test
      }
      if (overlay_high_test > overlay_high) {
        overlay_high <- overlay_high_test
      }
    }
  }
  
  overlay_low <- round(overlay_low, 1)
  overlay_high <- round(0.6*overlay_high, 1)
  
  #Base slice series anatomy
  ss_cluster_maps <- sliceSeries(nrow = 8, 
                                 ncol = 1, 
                                 begin = slice_begin, 
                                 end = slice_end) %>% 
    anatomy(anat_vol, 
            low = anat_low, 
            high = anat_high)
  
  #Add cluster maps
  for (i in 1:length(cluster_maps)) {
    if(i == 1) {
      ss_cluster_maps <- ss_cluster_maps %>% 
        overlay(mincArray(cluster_maps[[i]]), 
                low = overlay_low, 
                high = overlay_high, 
                symmetric = TRUE)
    } else {
      ss_cluster_maps <- ss_cluster_maps %>% 
        sliceSeries() %>% anatomy() %>% 
        overlay(mincArray(cluster_maps[[i]]), 
                low = overlay_low, 
                high = overlay_high, 
                symmetric = TRUE)
    }
    ss_cluster_maps <- ss_cluster_maps %>% 
      addtitle(names(cluster_maps)[i])
  }
  
  #Add legend
  ss_cluster_maps <- ss_cluster_maps %>% 
    legend("Effect size")
  
  if (save) {
    outfile <- str_c("ClusterMaps",
                     species,
                     gene_space,
                     "jacobians", jacobians,
                     "threshold", threshold_method,
                     threshold, mask_type,
                     "kcut", kcut,
                     "meta_k", meta_k,
                     sep = "_")
    outfile <- str_c(outfile, ".pdf")
    outfile <- file.path(outdir, outfile)
    pdf(file = outfile,
        width = unit(fig_width, "inch"),
        height = unit(fig_height, "inch"))
    draw(ss_cluster_maps)   
    dev.off()
  }
  
  return(draw(ss_cluster_maps))
  
}


plot_spider <- function(x, mouse_kcut, human_kcut, mouse_meta_k, human_meta_k, spokes, gene_space,
                        jacobians, mask_type, threshold_method, threshold, 
                        save = FALSE, outdir = "plots/", 
                        fig_width = 10, fig_height = 10) {
  
  if (spokes == "neuro-coarse") {
    infile <- "../../data/MouseHumanMatches_H10M09.csv"
  } else if (spokes == "neuro-mid") {
    infile <- "../../data/MouseHumanMatches_H88M67.csv"
  } else {
    stop("spokes must be one of {'neuro-coarse', 'neuro-mid'}")
  }
  df_neuro_pairs <- read_csv(infile, show_col_types = FALSE)
  
  df_human_meta_cluster <- generate_meta_clusters(x = x,
                                                  kcut = human_kcut,
                                                  axis = "rows") %>% 
    filter(meta_k == human_meta_k)
  
  
  df_mouse_meta_cluster <- generate_meta_clusters(x = x,
                                                  kcut = mouse_kcut,
                                                  axis = "columns") %>% 
    filter(meta_k == mouse_meta_k)
  
  #Human labels for microarray sample
  human_label_file <- "../../data/human/expression/AHBA_microarray_labels_studyspace_1.0mm.mnc"
  human_labels <- mincGetVolume(human_label_file)
  
  #Human definitions for microarray sample labels
  human_defs_file <- "../../data/human/expression/AHBA_microarray_coordinates_mni_defs.csv"
  human_defs <- read_csv(human_defs_file, show_col_types = FALSE)
  ind_atlas <- human_defs$label %in% human_labels
  human_defs <- human_defs[ind_atlas,]
  
  #Import human hierarchy
  human_tree_file <- "../../data/human/expression/HumanExpressionTree.RData"
  load(human_tree_file)
  tree_human <- Clone(treeHumanExpr)
  rm(treeHumanExpr)
  
  #Prune human tree to coarse levels
  tree_human_reduced <- Clone(tree_human)
  pruneAnatTree(tree_human_reduced, 
                nodes = df_neuro_pairs[["Human"]], 
                method = "BelowNode")
  
  human_labels_reduced <- reduce_human_labels(tree = tree_human_reduced,
                                              labels = human_labels,
                                              defs = human_defs)
  
  human_defs_reduced <- reduce_human_defs(tree = tree_human_reduced,
                                          defs = human_defs,
                                          simplify = TRUE)
  
  human_defs_reduced <- human_defs_reduced %>% 
    filter(name %in% df_neuro_pairs[["Human"]])
  
  human_labels_reduced[!(human_labels_reduced %in% human_defs_reduced[["label"]])] <- 0
  
  human_cluster_mask_dir <- str_c("../../data/human/clustering/cluster_masks/", 
                                  jacobians, "/resolution_1.0/mean/threshold_", 
                                  threshold_method, "/symmetric", 
                                  "/threshold_", threshold)
  human_cluster_mask_files <- list.files(human_cluster_mask_dir, 
                                         full.names = TRUE)
  human_masks <- import_cluster_masks(files = human_cluster_mask_files,
                                      df = df_human_meta_cluster)
  
  df_human_cluster_fractions <- map_dfr(.x = human_masks,
                                        .f = calc_cluster_region_fractions,
                                        .id = "cluster_id",
                                        labels = human_labels_reduced,
                                        defs = human_defs_reduced,
                                        mask_type = mask_type)
  
  mouse_label_file <- "../../data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc"
  mouse_labels <- mincGetVolume(mouse_label_file)
  mouse_labels_vol <- mincArray(mouse_labels)
  mouse_defs_file <- "../../data/mouse/atlas/DSURQE_40micron_R_mapping_long.csv"
  mouse_defs <- read_csv(mouse_defs_file, show_col_types = FALSE) %>% 
    select(name = Structure, label = Label)
  
  mouse_tree_file <- "../../data/mouse/expression/MouseExpressionTree_DSURQE.RData"
  load(mouse_tree_file)
  tree_mouse <- Clone(treeMouseExpr)
  rm(treeMouseExpr)
  
  tree_mouse_reduced <- Clone(tree_mouse)
  pruneAnatTree(tree_mouse_reduced, 
                nodes = df_neuro_pairs[["Mouse"]], 
                method = "BelowNode")
  
  mouse_labels_reduced <- hanatToAtlas(tree_mouse_reduced, mouse_labels_vol)
  mouse_defs_reduced <- hanatToAtlasDefs(tree_mouse_reduced)
  mouse_defs_reduced <- mouse_defs_reduced %>% 
    rename(name = Structure,
           label = Label)
  
  mouse_defs_reduced <- mouse_defs_reduced %>% 
    filter(name %in% df_neuro_pairs[["Mouse"]])
  
  mouse_labels_reduced[!(mouse_labels_reduced %in% mouse_defs_reduced[["label"]])] <- 0
  
  mouse_cluster_mask_dir <- str_c("../../data/mouse/clustering/cluster_masks/", 
                                  jacobians, "/resolution_200/mean/threshold_", 
                                  threshold_method, "/symmetric", 
                                  "/threshold_", threshold)
  mouse_cluster_mask_files <- list.files(mouse_cluster_mask_dir, 
                                         full.names = TRUE)
  
  mouse_masks <- import_cluster_masks(files = mouse_cluster_mask_files,
                                      df = df_mouse_meta_cluster)
  
  df_mouse_cluster_fractions <- map_dfr(.x = mouse_masks,
                                        .f = calc_cluster_region_fractions,
                                        .id = "cluster_id",
                                        labels = mouse_labels_reduced,
                                        defs = mouse_defs_reduced,
                                        mask_type = mask_type)
  
  df_spider <- prepare_spider_data(mouse = df_mouse_cluster_fractions,
                                   human = df_human_cluster_fractions,
                                   homologues = df_neuro_pairs)
  
  p_spider <- ggplot(data = df_spider, 
                     mapping = aes(x = name,
                                   y = fvoxels_expr, 
                                   col = species,
                                   group = group_id)) + 
    geom_line() + 
    coord_polar() + 
    labs(title = "Neuroanatomical weighting of clusters",
         subtitle = str_c("Number of mouse meta-clusters: ", mouse_kcut, "\n",
                          "Number of human meta-clusters: ", human_kcut, "\n",
                          "Mouse meta-cluster: ", mouse_meta_k, "\n",
                          "Human meta-cluster: ", human_meta_k, "\n"),
         x = NULL,
         y = "Fraction of expressing voxels in cluster",
         col = "Species") + 
    theme_bw() 
  
  if (save) {
    outfile <- str_c("SpiderPlot",
                     gene_space,
                     "jacobians", jacobians,
                     "threshold", threshold_method,
                     threshold, mask_type,
                     "kcut", kcut,
                     "mouse_meta_k", mouse_meta_k,
                     "human_meta_k", human_meta_k,
                     spokes,
                     sep = "_")
    outfile <- str_c(outfile, ".pdf")
    outfile <- file.path(outdir, outfile)
    pdf(file = outfile,
        width = unit(fig_width, "inch"),
        height = unit(fig_height, "inch"))
    print(p_spider)   
    dev.off()
  }
  
  return(p_spider)
  
}
