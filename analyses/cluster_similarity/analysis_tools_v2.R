#' Calculate cluster region fractions
#'
#' @param cluster (mincSingleDim)
#' @param labels (mincSingleDim)
#' @param defs (data.frame)
#'
#' @return (data.frame)
calc_cluster_region_fractions <- function(cluster, labels, defs) {
  voxels_cluster <- cluster == 1
  voxels_nonzero <- labels != 0
  labels_cluster <- labels[voxels_cluster & voxels_nonzero]
  if (length(labels_cluster) == 0){
    out <- tibble(label = defs[["label"]],
                  nvoxels = 0)
  } else {
    out <- table(labels_cluster) %>% 
      as_tibble() %>% 
      rename(label = labels_cluster,
             nvoxels = n) %>% 
      mutate(label = as.integer(label))
  }
  out <- out%>% 
    right_join(defs, by = "label") %>%
    mutate(nvoxels = ifelse(is.na(nvoxels), 0, nvoxels),
           nvoxels_cluster = sum(voxels_cluster),
           nvoxels_expr = length(labels_cluster),
           fvoxels_cluster = nvoxels/nvoxels_cluster,
           fvoxels_expr = nvoxels/nvoxels_expr,
           fvoxels_cluster = ifelse(is.nan(fvoxels_cluster), 0, fvoxels_cluster),
           fvoxels_expr = ifelse(is.nan(fvoxels_expr), 0, fvoxels_expr)) %>% 
    select(name,
           label,
           nvoxels,
           nvoxels_cluster,
           nvoxels_expr,
           fvoxels_cluster,
           fvoxels_expr)
  return(out)
}

#' Import cluster map
#'
#' @param species (character scalar) Species to use. Options: {mouse, 
#' human}.
#' @param jacobians (character scalar) Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param nk (numeric scalar) Number of clusters in the solution.
#' @param k (numeric scalar) Cluster number. 
#' @param mask (character scalar) Type of mask to apply to the cluster
#' maps. Ignored if NULL. Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) Thresholding method used 
#' to create the mask. Ignored if mask = NULL. Options: {intensity, 
#' topn}.
#' @param threshold (numeric scalar) Thresholding value applied to 
#' create the mask. Ignored if mask = NULL.
#'
#' @return (mincSingleDim) Cluster map
import_cluster_map <- function(species, nk, k, jacobians = NULL, mask = NULL, 
                               threshold_method = NULL, threshold = NULL, datadir = "../../data/") {
  
  #Argument options: jacobians
  opts_jacobians <- c("absolute", "relative")
  msg_jacobians <- str_c("jacobians must be one of {", 
                         str_flatten(opts_jacobians, collapse = ", "),
                         "}")
  if (is.null(jacobians)){
    stop(msg_jacobians)
  } else if (!(jacobians %in% opts_jacobians)) {
    stop(msg_jacobians)    
  }
  
  #Cluster directory
  if (species == "mouse") {
    cluster_dir <- str_c(file.path(datadir, "mouse/clustering/cluster_maps/"), 
                         jacobians, "/resolution_200/mean/")
  } else if (species == "human") {
    cluster_dir <- str_c(file.path(datadir, "human/clustering/cluster_maps/"), 
                         jacobians, "/resolution_1.0/mean/")
  } else {
    stop("species must be one of {mouse, human}")
  }
  
  #Import cluster map
  cluster_map_file <- list.files(cluster_dir, full.names = TRUE) %>% 
    str_subset(str_c("Clusternum_", nk)) %>%
    str_subset(str_c("Group_", k))
  cluster_map <- mincGetVolume(cluster_map_file)
  
  #Option to apply mask
  if (!is.null(mask)) {
    cluster_mask <- import_cluster_mask(datadir = datadir,
                                        species = species,
                                        jacobians = jacobians,
                                        nk = nk, k = k,
                                        type = mask,
                                        threshold_method = threshold_method,
                                        threshold = threshold)
    cluster_map[cluster_mask == 0] <- 0
  }
  
  return(cluster_map)    
  
}


#' Import cluster mask
#'
#' @param species (character scalar) Species to use. Options: {mouse, 
#' human}.
#' @param jacobians (character scalar) Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param nk (numeric scalar) Number of clusters in the solution.
#' @param k (numeric scalar) Cluster number. 
#' @param type (character scalar) Type of mask to return. 
#' Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) Thresholding method used 
#' to create the mask. Options: {intensity, topn}.
#' @param threshold (numeric scalar) Thresholding value applied to 
#' create the mask.
#'
#' @return (mincSingleDim) Cluster mask
import_cluster_mask <- function(species, nk, k, jacobians = NULL, 
                                type = "symmetric", threshold_method = NULL, 
                                threshold = NULL, datadir = "../../data/") {
  
  #Argument options: jacobians
  opts_jacobians <- c("absolute", "relative")
  msg_jacobians <- str_c("jacobians must be one of {", 
                         str_flatten(opts_jacobians, collapse = ", "),
                         "}")
  if (is.null(jacobians)){
    stop(msg_jacobians)
  } else if (!(jacobians %in% opts_jacobians)) {
    stop(msg_jacobians)    
  } 
  
  #Argument options: threshold_method  
  opts_threshold_method <- c("intensity", "topn")
  msg_threshold_method <- str_c("threshold_method must be one of {", 
                                str_flatten(opts_threshold_method, collapse = ", "),
                                "}")
  if (is.null(threshold_method)){
    stop(msg_threshold_method)
  } else if (!(threshold_method %in% opts_threshold_method)) {
    stop(msg_threshold_method)    
  }
  
  #Cluster mask directory
  if (species == "mouse") {
    maskdir <- str_c(file.path(datadir, "mouse/clustering/cluster_masks/"), 
                     jacobians, "/resolution_200/mean/threshold_", 
                     threshold_method, "/symmetric/")
    
  } else if (species == "human") {
    maskdir <- str_c(file.path(datadir, "human/clustering/cluster_masks/"), jacobians,
                     "/resolution_1.0/mean/threshold_", threshold_method,
                     "/symmetric/")
  } else {
    stop("species must be one of {mouse, human}")
  }
  
  if (is.null(threshold)) {
    opt_threshold <- str_remove(list.files(maskdir), "threshold_")
    stop(str_c("No threshold specified. Available thresholds are:\n", 
               str_flatten(opt_threshold, collapse = "\n")))
  } else {
    maskdir <- str_c(maskdir, "threshold_", threshold)
    if(!(dir.exists(maskdir))) {
      stop(str_c("Directory not found:", maskdir, "\n",
                 "Check threshold value:", threshold,
                 sep = " "))
    }
  }
  
  #Import cluster mask
  infile <- list.files(maskdir, full.names = TRUE) %>% 
    str_subset(str_c("Clusternum_", nk)) %>%
    str_subset(str_c("Group_", k))
  mask <- floor(mincGetVolume(infile))
  
  #Mask options
  if (type == "positive") {
    mask[mask < 0] <- 0
  } else if (type == "negative") {
    mask[mask > 0] <- 0
    mask <- abs(mask)
  } else if (type == "symmetric"){
    mask <- abs(mask)
  } else {
    opts_mask <- c("symmetric", "positive", "negative")
    msg_mask <- str_c("mask must be one of {", 
                      str_flatten(opts_mask, collapse = ", "),
                      "}")
    stop(msg_mask)
  }
  
  return(mask)
}


#' Import cluster similarity matrix
#'
#' @param gene_space (character scalar) Gene expression space used to 
#' compute the similarity matrix. Options: {average_latent_space, 
#' latent_space, input_space}.
#' @param jacobians (character scalar) Jacobians used to compute 
#' the expression signatures and similarity matrix. Options: {absolute,
#' relative}.
#' @param threshold_method (character scalar) Thresholding method used
#' to determine which voxels where included in the expression signatures.
#' Options: {intensity, topn}.
#' @param threshold (numeric scalar) Threshold used.
#' @param mask (character scalar) Type of mask used to compute the 
#' expression signatures and similarity matrx. Options: {symmetric, 
#' positive, negative}.
#' @param latent_space_id (integer scalar) Latent space to use. Ignored 
#' if gene_space is not "latent_space".
#'
#' @return (matrix) Mouse-human cluster similarity matrix.
import_similarity_matrix <- function(gene_space = "average_latent_space", 
                                     jacobians = NULL, 
                                     threshold_method = NULL, 
                                     threshold = NULL, mask = NULL, 
                                     latent_space_id = 1) {
  
  #Argument options: jacobians
  opts_jacobians <- c("absolute", "relative")
  msg_jacobians <- str_c("jacobians must be one of {", 
                         str_flatten(opts_jacobians, collapse = ", "),
                         "}")
  if (is.null(jacobians)){
    stop(msg_jacobians)
  } else if (!(jacobians %in% opts_jacobians)) {
    stop(msg_jacobians)    
  } else {
    if (jacobians == "absolute") {
      jacobians <- "abs"
    } else {
      jacobians <- "rel"
    }
  }
  
  #Argument options: threshold_method  
  opts_threshold_method <- c("intensity", "topn")
  msg_threshold_method <- str_c("threshold_method must be one of {", 
                                str_flatten(opts_threshold_method, collapse = ", "),
                                "}")
  if (is.null(threshold_method)){
    stop(msg_threshold_method)
  } else if (!(threshold_method %in% opts_threshold_method)) {
    stop(msg_threshold_method)    
  }
  
  #Argument options: mask
  opts_mask <- c("symmetric", "positive", "negative")
  msg_mask <- str_c("mask must be one of {", 
                    str_flatten(opts_mask, collapse = ", "),
                    "}")
  if (is.null(mask)){
    stop(msg_mask)
  } else if (!(mask %in% opts_mask)) {
    stop(msg_mask)    
  }
  
  msg_threshold <- str_c("threshold is NULL. Specify a threshold.")
  if (is.null(threshold)){
    stop(msg_threshold)
  }
  
  if (gene_space == "input_space") {
    stop("Feature currently unavailable")    
  } else if (gene_space == "latent_space") {
    stop("Feature currently unavailable")    
  } else if (gene_space == "average_latent_space") {
    
    indir <- "../../data/similarity_matrix/latent_space_100/"
    infile <- str_c("similarity_hm_", jacobians, "_mean_threshold_", 
                    threshold_method, "_", threshold, "_", mask,
                    "_latentspace100.csv")
    infile <- file.path(indir, infile)
    
    if (!(file.exists(infile))){
      msg_infile <- str_c("Input file not found: ", infile, "\n",
                          "Check threshold value: ", threshold)
      stop(msg_infile)
    }
    
    similarity <- read_csv(file = infile, show_col_types = FALSE) %>% 
      column_to_rownames("cluster_id") %>% 
      as.matrix()
    
  } else {
    stop(str_c("gene_space must be one of ",
               "{average_latent_space, latent_space, input_space}"))
  }
  return(similarity)
}



#' Plot cluster map slice series
#'
#' @param species (character scalar) Species to use. Options: {mouse, 
#' human}.
#' @param nk (numeric scalar) Number of clusters in the solution.
#' @param k (numeric scalar) Cluster number. 
#' @param nrow (numeric scalar) Number of rows in slice series.
#' @param ncol (numeric scalar) Number of columns in slice series.
#' @param jacobians (character scalar) Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param mask (character scalar) Type of mask to apply to the cluster
#' maps. Ignored if NULL. Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) Thresholding method used 
#' to create the mask. Ignored if mask = NULL. Options: {intensity, 
#' topn}.
#' @param threshold (numeric scalar) Thresholding value applied to 
#' create the mask. Ignored if mask = NULL.
#' @param overlay_low (numeric scalar) Lower threshold for slice series 
#' overlay. Generated programmatically if NULL.
#' @param overlay_high (numeric scalar) Upper threshold for slice series 
#' overlay. Generated programmatically if NULL.
#'
#' @return 
plot_cluster_map <- function(species, nk, k, nrow = 5, ncol = 5, jacobians = NULL,
                             mask = NULL, threshold_method = NULL, threshold = NULL,
                             overlay_low = NULL, overlay_high = NULL) {
  
  if (species == "mouse") {
    anat_file <- "../../data/mouse/atlas/DSURQE_CCFv3_average_200um.mnc"
    slice_begin <- 10
    slice_end <- 60
    anat_low <- 700
    anat_high <- 1400 
  } else if (species == "human") {
    anat_file <- "../../data/human/registration/reference_files/model_1.0mm.mnc"
    slice_begin <- 50
    slice_end <- 200
    anat_low <- 3
    anat_high <- 7
  } else {
    stop("species must be one of {mouse, human}")
  }
  
  anat_vol <- mincArray(mincGetVolume(anat_file))
  
  #Cluster maps and masks
  cluster_map <- import_cluster_map(species = species,
                                    nk = nk, k = k,
                                    jacobians = jacobians,
                                    mask = mask,
                                    threshold_method = threshold_method,
                                    threshold = threshold)
  
  if (is.null(overlay_low)) {
    overlay_low <- min(abs(cluster_map[cluster_map != 0]))
    overlay_low <- round(overlay_low, 1)
  }
  
  if (is.null(overlay_high)) {
    overlay_high <- max(abs(cluster_map[cluster_map != 0]))
    overlay_high <- round(0.6*overlay_high, 1)
  }
  
  #Base slice series anatomy
  ss_cluster_map <- sliceSeries(nrow = nrow, 
                                ncol = ncol, 
                                begin = slice_begin, 
                                end = slice_end) %>% 
    anatomy(anat_vol, 
            low = anat_low, 
            high = anat_high) %>% 
    overlay(mincArray(cluster_map),
            low = overlay_low,
            high = overlay_high,
            symmetric = TRUE) %>% 
    legend("Effect size") 
  
  return(draw(ss_cluster_map))
  
}


reduce_atlas <- function(species, tree, labels, defs, structs, remove = TRUE) {
  
  tree_reduced <- Clone(tree)
  pruneAnatTree(tree_reduced, 
                nodes = structs,
                method = "BelowNode")
  
  if (remove) {
    flag <- TRUE
    leafs <- tree_reduced$Get("name", filterFun = isLeaf)
    while (flag) {
      leafs_to_prune <- leafs[!(leafs %in% structs)]
      pruneAnatTree(tree_reduced,
                    nodes = leafs_to_prune,
                    method = "AtNode")
      leafs <- tree_reduced$Get("name", filterFun = isLeaf)
      if (all(leafs %in% structs)) {
        flag <- FALSE
      }
    } 
  }
  
  defs <- defs[defs$label %in% labels,]
  
  if (species == "mouse") {
    labels_reduced <- hanatToAtlas(tree_reduced, mincArray(labels))
    defs_reduced <- hanatToAtlasDefs(tree_reduced)
    defs_reduced <- defs_reduced %>% 
      rename(name = Structure,
             label = Label)
  } else if (species == "human") {
    labels_reduced <- reduce_human_labels(tree = tree_reduced,
                                          labels = labels,
                                          defs = defs)
    defs_reduced <- reduce_human_defs(tree = tree_reduced,
                                      defs = defs,
                                      simplify = TRUE)
  } else {
    stop("species must be one of {mouse, human}")
  }
  
  return(list(labels = labels_reduced,
              defs = defs_reduced))
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


relabel_matrix_axis <- function(x, tree, structs, axis = "columns"){
  
  if (axis == "columns") {
    names_old <- colnames(x)
  } else if (axis == "rows") {
    names_old <- rownames(x)
  }  
  
  names_new <- character(length(names_old))
  paths <- tree$Get("path", filterFun = isLeaf)
  for(i in 1:length(paths)){
    ind <- which(names_old == names(paths)[[i]])
    names_new[[ind]] <- structs[structs %in% paths[[i]]]
  }
  
  if (axis == "columns") {
    colnames(x) <- names_new
  } else if (axis == "rows") {
    rownames(x) <- names_new
  } 
  
  return(x)
  
}



#' Title
#'
#' @param x 
#' @param tree 
#' @param labels 
#' @param mask 
#'
#' @return
label_mouse_data <- function(x, tree, labels, mask) {
  
  labels_pruned <- hanatToAtlas(tree, mincArray(labels))
  defs_pruned <- hanatToAtlasDefs(tree)
  
  labels_pruned <- labels_pruned[mask == 1]
  
  x <- x %>% 
    mutate(label = labels_pruned) %>% 
    left_join(defs_pruned, by = c('label' = 'Label')) %>% 
    select(-label) %>% 
    rename(name = Structure)
  
  return(x)
  
}

#' Title
#'
#' @param x 
#' @param tree 
#' @param samples 
#'
#' @return
label_human_data <- function(x, tree, samples) {
  
  tree <- Clone(tree)
  
  tree$Do(function(node){
    node$struct <- rep(node$name, length(node$samples))
  }, traversal = 'post-order')
  
  tree_structs <- unlist(tree$Get('struct', filterFun = isLeaf))
  names(tree_structs) <- NULL
  
  tree_samples <- unlist(tree$Get('samples', filterFun = isLeaf))
  names(tree_samples) <- NULL
  
  sample_defs <- tibble(name = tree_structs,
                        sample_id = tree_samples)
  
  x <- x %>% 
    mutate(sample_id = samples) %>% 
    left_join(sample_defs, by = 'sample_id') %>% 
    select(-sample_id) 
  
  return(x)
  
}

aggregate_data <- function(x, groupby, method = 'mean', output = 'tibble'){
  
  if (method == 'mean') {
    aggfunc <- mean
  } else if (method == 'median') {
    aggfunc <- median
  } else if (method == 'sd') {
    aggfunc <- sd
  } else {
    stop()
  }
  
  if (!(groupby %in% colnames(x))) {
    stop()
  }
  
  if (any(is.na(x[[groupby]]))){
    warning(paste("Grouping variable", groupby, "contains NA.",
                  "These observations will be ommitted when aggregating."))
  }
  
  x <- x %>% 
    filter(!is.na(.data[[groupby]])) %>% 
    group_by_at(.vars = vars(groupby)) %>% 
    summarise_all(.funs = aggfunc) %>% 
    ungroup() 
  
  if (output == 'tibble') {
    x %>% 
      as_tibble() %>% 
      return()
  } else if (output == 'data.frame') {
    x %>% 
      as.data.frame(stringsAsFactors = FALSE) %>% 
      return()
  } else if (output == 'matrix') {
    x %>% 
      column_to_rownames(groupby) %>% 
      as.matrix() %>% 
      return()
  } else {
    stop()
  }
}

prepare_cluster_fractions <- function(species, structs, nk, k, jacobians, 
                                      mask, threshold_method, threshold, 
                                      tree, labels, defs, remove = FALSE, datadir = "../../data/") {
  
  atlas_reduced <- reduce_atlas(species = species,
                                tree = tree,
                                labels = labels,
                                defs = defs,
                                structs = structs,
                                remove = remove)
  
  cluster_mask <- import_cluster_mask(datadir = datadir,
                                      species = species,
                                      nk = nk, k = k,
                                      jacobians = jacobians,
                                      type = mask,
                                      threshold_method = threshold_method,
                                      threshold = threshold)
  
  df_fractions <- calc_cluster_region_fractions(cluster = cluster_mask,
                                                labels = atlas_reduced[["labels"]],
                                                defs = atlas_reduced[["defs"]]) %>% 
    mutate(cluster_id = str_c(nk, k, sep = "-")) 
  
  return(df_fractions)
  
}

prepare_histogram_data <- function(x) {
  
  regions <- c("cortex", "cerebellum")
  df_regions <- expand_grid(mouse = regions, 
                            human = regions)
  
  mouse_labels <- list("Isocortex",
                       "Cerebellar cortex")
  names(mouse_labels) <- regions
  
  human_labels <- list(c("frontal lobe", 
                         "insula", 
                         "occipital lobe", 
                         "parietal lobe", 
                         "temporal lobe"),
                       "cerebellar cortex")
  names(human_labels) <- regions
  
  list_regions <- vector(mode = "list", length = nrow(df_regions))
  for(i in 1:nrow(df_regions)) {
    
    mouse_region <- df_regions[[i, "mouse"]]
    human_region <- df_regions[[i, "human"]]
    
    ind_mouse <- colnames(x) %in% mouse_labels[[mouse_region]]
    ind_human <- rownames(x) %in% human_labels[[human_region]]
    
    x_region <- x[ind_human, ind_mouse]
    
    list_regions[[i]] <- tibble(similarity = as.numeric(x_region),
                                mouse = mouse_region,
                                human = human_region)
  }
  df <- reduce(.x = list_regions,
               .f = bind_rows)
  
  df <- df %>% 
    unite(col = "pair",
          mouse, human,
          sep = "-",
          remove = FALSE) %>% 
    mutate(pair = ifelse(pair %in% c("cortex-cerebellum", 
                                     "cerebellum-cortex"),
                         "cortex-cerebellum", 
                         pair))
  return(df)
}


prepare_spider_data <- function(mouse, human, homologues) {
  
  human <- human %>% 
    rename(human = name)
  
  mouse <- mouse %>% 
    rename(mouse = name)
  
  #Filter for regions in homologues
  human <- semi_join(human, homologues, by = "human")
  mouse <- semi_join(mouse, homologues, by = "mouse")
  
  #Rename human regions with proper name
  human <- human %>% 
    inner_join(homologues, by = "human") %>% 
    select(-human, -mouse)
  
  #Rename mouse regions with proper name
  mouse <- mouse %>% 
    inner_join(homologues, by = "mouse") %>% 
    select(-human, -mouse)
  
  #Convert names to factor
  human <- human %>% 
    mutate(name = factor(name, levels = homologues[["name"]]),
           species = "human")
  mouse <- mouse %>% 
    mutate(name = factor(name, levels = homologues[["name"]]),
           species = "mouse")
  
  #Combine mouse and human data
  out <- bind_rows(human, mouse)
  
  #Create a unique group ID
  out <- out %>% 
    unite(col = group_id, 
          species, cluster_id, 
          sep = "-", remove = FALSE)
  
  return(out)
}

# plot_heatmap <- function(x, clustering = FALSE, cutree_rows = 1, cutree_cols = 1, 
#                          gene_space, jacobians, mask_type, threshold_method,
#                          threshold, metric, header = TRUE, padding = 0.2, 
#                          save = FALSE, outdir = "plots/", 
#                          fig_width = 10, fig_height = 10) {

plot_heatmap <- function(x, clustering = FALSE, cutree_rows = 1, 
                         cutree_cols = 1, annotations = "nk-k", 
                         padding = 0.2, color = NULL) {
  
  
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
    
    palette <- mako(n = 10, direction = -1, begin = 0.3)
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
    
    if (annotations == "nk") {
      df_annotations <- select(df_annotations, nk)
      annotation_colours <- list(nk = nk_colours)
    } else if (annotations == "k") {
      df_annotations <- select(df_annotations, k)
      annotation_colours <- list(k = k_colours)
    } else {
      annotation_colours <- list(nk = nk_colours,
                                 k = k_colours)
    }
    
    if(is.null(color)) {
      color <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                                 "RdYlBu")))(100)
    }
    p_heatmap <- as.ggplot(pheatmap(mat = x,
                                    color = color,
                                    cluster_cols = F, cluster_rows = F,
                                    annotation_row = df_annotations,
                                    annotation_col = df_annotations,
                                    annotation_colors = annotation_colours,
                                    silent = TRUE,
                                    na_col = 'black'))
    
  }
  
  # header <- FALSE
  # if (header) {
  #   p_heatmap <- p_heatmap #+
  #     # labs(title = "Cluster similarity",
  #     #      subtitle = str_c("Gene space: ", gene_space, "\n",
  #     #                       "Jacobians: ", jacobians, "\n",
  #     #                       "Mask type: ", mask_type, "\n",
  #     #                       "Thresholding method: ", threshold_method, "\n",
  #     #                       "Threshold: ", threshold, "\n",
  #     #                       "Metric: ", metric, "\n",
  #     #                       "Rows: Human", "\n",
  #     #                       "Columns: Mouse", "\n"))
  #
  # }
  #
  p_heatmap <- p_heatmap +
    theme(plot.margin = margin(t = padding,
                               r = padding,
                               b = padding,
                               l = padding,
                               unit = "in"))
  #
  # if (save) {
  #   outfile <- str_c("ClusterSimilarity",
  #                    gene_space,
  #                    "jacobians", jacobians,
  #                    "threshold", threshold_method,
  #                    threshold,
  #                    mask_type,
  #                    sep = "_")
  #   if (clustering) {
  #     outfile <- str_c(outfile,
  #                      "clustered",
  #                      "rowcut", cutree_rows,
  #                      "colcut", cutree_cols,
  #                      sep = "_")
  #   }
  #   outfile <- str_c(outfile, ".pdf")
  #   outfile <- file.path(outdir, outfile)
  #   pdf(file = outfile,
  #       width = unit(fig_width, "inch"),
  #       height = unit(fig_height, "inch"))
  #   print(p_heatmap)
  #   dev.off()
  # }
  
  return(p_heatmap)
  
}



plot_spider_chart <- function(nk, k, spokes, jacobians, threshold_method, threshold, step = 0.2, datadir = "../../data/", save = FALSE, outdir = "./", fig_width = 10, fig_height = 10) {
  
  if (spokes == "neuro-coarse") {
    infile <- file.path(datadir, "MouseHumanMatches_H10M09.csv")
  } else if (spokes == "neuro-mid") {
    infile <- file.path(datadir, "MouseHumanMatches_H23M23.csv")
  } else {
    stop("spokes must be one of {'neuro-coarse', 'neuro-mid'}")
  }
  df_spokes <- read_csv(infile, show_col_types = FALSE)
  colnames(df_spokes) <- str_to_lower(colnames(df_spokes))
  df_spokes <- df_spokes %>% 
    filter(!(name %in% c("White matter", "Ventricles")))
  
  species <- c("mouse", "human")
  cluster_fractions <- vector(mode = "list", length = length(species))
  names(cluster_fractions) <- species
  for (i in 1:length(species)) {
    
    if (species[[i]] == "mouse") {
      
      #Mouse DSURQE labels
      label_file <- file.path(datadir, "mouse/atlas/DSURQE_CCFv3_labels_200um.mnc")
      labels <- mincGetVolume(label_file)
      
      #Mouse definitions for DSURQE labels
      defs_file <- file.path(datadir, "mouse/atlas/DSURQE_40micron_R_mapping_long.csv")
      defs <- read_csv(defs_file, show_col_types = FALSE) %>% 
        select(name = Structure, label = Label)
      
      #Mouse hierarchy
      tree_file <- file.path(datadir, "mouse/expression/MouseExpressionTree_DSURQE.RData")
      load(tree_file)
      tree <- Clone(treeMouseExpr)
      rm(treeMouseExpr)
      
    } else {
      
      #Human labels for microarray sample
      label_file <- file.path(datadir, "human/expression/AHBA_microarray_labels_studyspace_1.0mm.mnc")
      labels <- mincGetVolume(label_file)
      
      #Human definitions for microarray sample labels
      defs_file <- file.path(datadir, "human/expression/AHBA_microarray_coordinates_mni_defs.csv")
      defs <- read_csv(defs_file, show_col_types = FALSE)
      
      #Human hierarchy
      tree_file <- file.path(datadir, "human/expression/HumanExpressionTree.RData")
      load(tree_file)
      tree <- Clone(treeHumanExpr)
      rm(treeHumanExpr)
      
    }
    
    df_species_frac_pos <- prepare_cluster_fractions(datadir = datadir, 
                                                     species = species[[i]],
                                                     structs = df_spokes[[species[[i]]]],
                                                     nk = nk[[species[[i]]]],
                                                     k = k[[species[[i]]]],
                                                     jacobians = jacobians,
                                                     mask = "positive",
                                                     threshold_method = threshold_method,
                                                     threshold = threshold,
                                                     tree = tree,
                                                     labels = labels,
                                                     defs = defs) %>% 
      mutate(sign = "positive")
    
    df_species_frac_neg <- prepare_cluster_fractions(datadir = datadir,
                                                     species = species[[i]],
                                                     structs = df_spokes[[species[[i]]]],
                                                     nk = nk[[species[[i]]]],
                                                     k = k[[species[[i]]]],
                                                     jacobians = jacobians,
                                                     mask = "negative",
                                                     threshold_method = threshold_method,
                                                     threshold = threshold,
                                                     tree = tree,
                                                     labels = labels,
                                                     defs = defs) %>% 
      mutate(sign = "negative",
             fvoxels_expr = -1*fvoxels_expr)
    
    cluster_fractions[[species[[i]]]] <- bind_rows(df_species_frac_pos,
                                                   df_species_frac_neg) %>% 
      mutate(species = species[[i]])
    
  }
  
  df_cluster_fractions <- prepare_spider_data(mouse = cluster_fractions[["mouse"]],
                                              human = cluster_fractions[["human"]],
                                              homologues = df_spokes)
  
  df_cluster_fractions <- df_cluster_fractions %>% 
    unite(col = "group_id",
          species, cluster_id, sign,
          sep = "-",
          remove = FALSE)
  
  breaks_end <- round(max(abs(df_cluster_fractions$fvoxels_expr)), digits = 1)+0.1
  breaks_end <- ifelse(breaks_end > 1, 1, breaks_end)
  breaks_start <- -1*breaks_end
  radial_breaks <- seq(breaks_start, breaks_end, by = step)
  radial_labels <- tibble(name = levels(df_cluster_fractions$name)[length(levels(df_cluster_fractions$name))],
                          breaks = radial_breaks)
  
  p_spider <- ggplot(data = df_cluster_fractions, 
                     mapping = aes(x = name,
                                   y = fvoxels_expr, 
                                   col = group_id,
                                   group = group_id)) + 
    # geom_line(size = 1) + 
    geom_line(size = 2.0) + 
    geom_hline(yintercept = 0,
               size = 1.2,
               linetype = "dashed") + 
    geom_text(data = radial_labels,
              inherit.aes = FALSE,
              mapping = aes(x = name, 
                            y = breaks,
                            label = breaks),
              size = 6,
              # size = 4,
              nudge_x = 0.5,
              nudge_y = 0.25*step) + 
    coord_polar() +
    scale_y_continuous(breaks = radial_breaks) +
    scale_color_manual(values = c("deepskyblue", "deepskyblue4", "darkorange", "orange4"), 
                       labels = c("Human -- negative", 
                                  "Human -- positive", 
                                  "Mouse -- negative", 
                                  "Mouse -- positive")) + 
    labs(x = NULL,
         y = "Fraction of expressing voxels in cluster",
         col = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  if (save) {
    outfile <- str_c("spider", spokes, sep = "_")
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



plot_cluster_map_array <- function(species, nk, k, nrow = 8, jacobians = NULL,
                                   mask = NULL, threshold_method = NULL, threshold = NULL,
                                   overlay_low = NULL, overlay_high = NULL, datadir = "../../data/") {
  
  if (species == "mouse") {
    anat_file <- file.path(datadir, "mouse/atlas/DSURQE_CCFv3_average_200um.mnc")
    slice_begin <- 10
    slice_end <- 60
    anat_low <- 700
    anat_high <- 1400 
  } else if (species == "human") {
    anat_file <- file.path(datadir, "human/registration/reference_files/model_1.0mm.mnc")
    slice_begin <- 50
    slice_end <- 180
    anat_low <- 3
    anat_high <- 7
  } else {
    stop("species must be one of {mouse, human}")
  }
  
  anat_vol <- mincArray(mincGetVolume(anat_file))
  
  cluster_maps <- vector(mode = "list", length = length(nk))
  names(cluster_maps) <- str_c(nk, k, sep = "-")
  for (i in 1:length(cluster_maps)) {
    cluster_maps[[i]] <- import_cluster_map(datadir = datadir,
                                            species = species,
                                            nk = nk[i], 
                                            k = k[i],
                                            jacobians = jacobians,
                                            mask = mask,
                                            threshold_method = threshold_method,
                                            threshold = threshold)
  }
  
  # if (is.null(overlay_low)) {
  #   overlay_low <- min(abs(cluster_map[cluster_map != 0]))
  #   overlay_low <- round(overlay_low, 1)
  # }
  # 
  # if (is.null(overlay_high)) {
  #   overlay_high <- max(abs(cluster_map[cluster_map != 0]))
  #   overlay_high <- round(0.6*overlay_high, 1)
  # }
  
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
  #Base slice series anatomy
  ss_cluster_maps <- sliceSeries(nrow = nrow, 
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
  
  return(draw(ss_cluster_maps))
  
}

