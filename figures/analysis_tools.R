
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
  out <- out %>% 
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


#Radar coordinate system for ggplot2
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  
  #dirty
  rename_data <- function(coord, data) {
    if (coord$theta == "y") {
      plyr::rename(data, c("y" = "theta", "x" = "r"), warn_missing = FALSE)
    } else {
      plyr::rename(data, c("y" = "r", "x" = "theta"), warn_missing = FALSE)
    }
  }
  
  theta_rescale <- function(coord, x, scale_details) {
    rotate <- function(x) (x + coord$start) %% (2 * pi) * coord$direction
    rotate(scales::rescale(x, c(0, 2 * pi), scale_details$theta.range))
  }
  
  r_rescale <- function(coord, x, scale_details) {
    scales::rescale(x, c(0, 0.4), scale_details$r.range)
  }
  
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE,
          render_bg = function(self, scale_details, theme) {
            scale_details <- rename_data(self, scale_details)
            
            theta <- if (length(scale_details$theta.major) > 0)
              theta_rescale(self, scale_details$theta.major, scale_details)
            thetamin <- if (length(scale_details$theta.minor) > 0)
              theta_rescale(self, scale_details$theta.minor, scale_details)
            # thetafine <- seq(0, 2 * pi, length.out = 100)
            thetafine <- theta
            
            rfine <- c(r_rescale(self, scale_details$r.major, scale_details))
            
            # This gets the proper theme element for theta and r grid lines:
            #   panel.grid.major.x or .y
            majortheta <- paste("panel.grid.major.", self$theta, sep = "")
            minortheta <- paste("panel.grid.minor.", self$theta, sep = "")
            majorr     <- paste("panel.grid.major.", self$r,     sep = "")
            
            ggplot2:::ggname("grill", grid::grobTree(
              ggplot2:::element_render(theme, "panel.background"),
              if (length(theta) > 0) ggplot2:::element_render(
                theme, majortheta, name = "angle",
                x = c(rbind(0, 0.45 * sin(theta))) + 0.5,
                y = c(rbind(0, 0.45 * cos(theta))) + 0.5,
                id.lengths = rep(2, length(theta)),
                default.units = "native"
              ),
              if (length(thetamin) > 0) ggplot2:::element_render(
                theme, minortheta, name = "angle",
                x = c(rbind(0, 0.45 * sin(thetamin))) + 0.5,
                y = c(rbind(0, 0.45 * cos(thetamin))) + 0.5,
                id.lengths = rep(2, length(thetamin)),
                default.units = "native"
              ),
              
              ggplot2:::element_render(
                theme, majorr, name = "radius",
                x = rep(rfine, each = length(thetafine)) * sin(thetafine) + 0.5,
                y = rep(rfine, each = length(thetafine)) * cos(thetafine) + 0.5,
                id.lengths = rep(length(thetafine), length(rfine)),
                default.units = "native"
              )
            ))
          })
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
                               threshold_method = NULL, threshold = NULL, datadir = "data/") {
  
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
                                threshold = NULL, datadir = "data/") {
  
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
    
    indir <- "../data/similarity_matrix/latent_space_100/"
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


#' Intersect data with neuroanatomical homologues
#'
#' @param x (list) A list with named elements "mouse" and "human" 
#' containing data to intersect. 
#' @param homologues (data.frame) A data frame containing the set
#' of mouse and human neuroanatomical homologues.
#'
#' @return (list) A list with named elements "mouse" and "human" 
#' containing the intersected data.
intersect_neuro_homologues <- function(x, homologues) {
  for (i in 1:length(x)) {
    species <- names(x)[[i]]
    cols_init <- colnames(x[[i]])
    colnames(x[[i]])[cols_init == "name"] <- species
    x[[i]] <- inner_join(x[[i]], homologues, by = species)
    x[[i]][["name"]] <- factor(x[[i]][["name"]], levels = homologues[["name"]])
    x[[i]][["species"]] <- species
    x[[i]] <- x[[i]][,match(cols_init, colnames(x[[i]]))]
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


#' Prepare cluster regional fractions data
#'
#' @param datadir (character scalar) Path to data directory.
#' @param species (character scalar) Species to use {"mouse", "human"}.
#' @param structs (character vector) Vector of structures to compute 
#' voxel fractions.
#' @param nk (integer scalar) Number of clusters in the solution.
#' @param k (integer scalar) Cluster number. 
#' @param jacobians (character scalar) Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param mask (character scalar) Type of mask to return. 
#' Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) Thresholding method used
#' to create the mask. Options: {intensity, topn}.
#' @param threshold (character scalar). Threshold value applied to
#' create the mask.
#' @param tree (data.tree) A data tree containing the neuroanatomical
#' hierarchy. Leaf nodes must correspond to the atlas regions.
#' @param labels (mincSingleDim) Atlas labels.
#' @param defs (data.frame) Atlas definitions.
#' @param remove (logical scalar) Option to remove extra nodes not 
#' specified in argument structs.
#'
#' @return
prepare_cluster_fractions <- function(datadir = "data/", species, structs, nk, k, 
                                      jacobians, mask, threshold_method, 
                                      threshold, tree, labels, defs, 
                                      remove = TRUE) {
  
  if (species == "mouse") {
    atlas_reduced <- reduce_mouse_atlas(tree = tree,
                                        labels = labels,
                                        defs = defs,
                                        nodes = structs,
                                        remove = remove)
  } else if (species == "human") {
    atlas_reduced <- reduce_human_atlas(tree = tree,
                                        labels = labels,
                                        defs = defs,
                                        nodes = structs,
                                        remove = remove)
  }
  
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


#' Prepare cluster data for radar chart visualization
#'
#' @param datadir (character scalar) Path to data directory. 
#' @param spokes (data.frame) Data frame containing spokes for radar 
#' chart.
#' @param nk (integer scalar) Number of clusters in the solution.
#' @param k (integer scalar) Cluster number. 
#' @param jacobians (character scalar) Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param mask (character scalar) Type of mask to return. 
#' Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) Thresholding method used
#' to create the mask. Options: {intensity, topn}.
#' @param threshold (character scalar). Threshold value applied to
#' create the mask.
#' @param trees (list)  A list with named elements "mouse" and "human" 
#' containing data trees.
#' @param labels (list)  A list with named elements "mouse" and "human" 
#' containing atlas labels.
#' @param defs (list)  A list with named elements "mouse" and "human" 
#' containing atlas definitions.
#'
#' @return (data.frame) Radar chart data.
prepare_radar_chart <- function(datadir = "data/", spokes, nk, k, jacobians, 
                                mask, threshold_method, threshold, trees, 
                                labels, defs) {
  
  #Compute cluster fractions for mouse and human clusters
  species <- c("mouse", "human")
  fractions <- vector(mode = "list", length = length(species))
  names(fractions) <- species
  for (i in 1:length(species)) {
    
    positive_fractions <- prepare_cluster_fractions(datadir = datadir, 
                                                    species = species[[i]],
                                                    structs = spokes[[species[[i]]]],
                                                    nk = nk[[species[[i]]]],
                                                    k = k[[species[[i]]]],
                                                    jacobians = jacobians,
                                                    mask = "positive",
                                                    threshold_method = threshold_method,
                                                    threshold = threshold,
                                                    tree = trees[[species[[i]]]],
                                                    labels = labels[[species[[i]]]],
                                                    defs = defs[[species[[i]]]]) %>% 
      select(name, positive = fvoxels_expr)
    
    negative_fractions <- prepare_cluster_fractions(datadir = datadir, 
                                                    species = species[[i]],
                                                    structs = spokes[[species[[i]]]],
                                                    nk = nk[[species[[i]]]],
                                                    k = k[[species[[i]]]],
                                                    jacobians = jacobians,
                                                    mask = "negative",
                                                    threshold_method = threshold_method,
                                                    threshold = threshold,
                                                    tree = trees[[species[[i]]]],
                                                    labels = labels[[species[[i]]]],
                                                    defs = defs[[species[[i]]]]) %>% 
      select(name, negative = fvoxels_expr) 
    
    fractions[[species[[i]]]] <- inner_join(positive_fractions,
                                            negative_fractions,
                                            by = "name") %>% 
      mutate(negative = -1*negative,
             species = species[[i]])
    
    
  }
  
  #Intersect regions with neuroanatomical homologues
  fractions <- intersect_neuro_homologues(x = fractions, 
                                          homologues = spokes)
  
  #Reduce list to a data frame
  radar <- reduce(.x = fractions, .f = bind_rows)
  
  #Relevel regions so that the radar chart closes on itself
  lvls <- levels(radar$name)
  lvls <- c(" ", lvls)
  radar_dummy <- radar %>% 
    mutate(name = as.character(name)) %>% 
    filter(name == lvls[length(lvls)]) %>% 
    mutate(name = " ")
  radar <- bind_rows(radar, radar_dummy)
  radar <- mutate(radar, name = factor(name, levels = lvls))
  
  return(radar)
  
}


#' Create a reduced version of a human atlas. 
#'
#' @param tree (data.tree) A data tree containing the neuroanatomical
#' hierarchy. Leaf nodes must correspond to the atlas regions.
#' @param labels (mincSingleDim) Atlas labels.
#' @param defs (data.frame) Atlas definitions.
#' @param nodes (character vector) Tree nodes at which to aggregate the
#' atlas.
#' @param remove (logical scalar) Option to remove extra nodes not 
#' specified in argument nodes.
#'
#' @return (list) A list with fields "labels" and "defs" containing the 
#' reduced atlas labels and definitions.
reduce_human_atlas <- function(tree, labels, defs, nodes, remove = TRUE) {
  
  #Prune the tree to the specified nodes
  tree_reduced <- prune_tree(tree = tree,
                             nodes = nodes,
                             method = "below",
                             remove = remove)
  
  #Filter definitions labels for those in labels
  defs <- defs[defs$label %in% labels,]
  
  #Reduce atlas labels
  labels_reduced <- reduce_human_labels(tree = tree_reduced,
                                        labels = labels,
                                        defs = defs)
  
  #Reduce atlas definitions
  defs_reduced <- reduce_human_defs(tree = tree_reduced,
                                    defs = defs,
                                    simplify = TRUE)
  
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


#' Create a reduced version of a mouse atlas. 
#'
#' @param tree (data.tree) A data tree containing the neuroanatomical
#' hierarchy. Leaf nodes must correspond to the atlas regions.
#' @param labels (mincSingleDim) Atlas labels.
#' @param defs (data.frame) Atlas definitions.
#' @param nodes (character vector) Tree nodes at which to aggregate the
#' atlas.
#' @param remove (logical scalar) Option to remove extra nodes not 
#' specified in argument nodes.
#'
#' @return (list) A list with fields "labels" and "defs" containing the 
#' reduced atlas labels and definitions.
reduce_mouse_atlas <- function(tree, labels, defs, nodes, remove = TRUE) {
  
  #Prune the tree to the specified nodes
  tree_reduced <- prune_tree(tree = tree,
                             nodes = nodes,
                             method = "below",
                             remove = remove)
  
  #Reduce atlas labels
  labels_reduced <- hanatToAtlas(tree_reduced, mincArray(labels))
  
  #Reduce atlas definitions
  defs <- defs[defs$label %in% labels,]
  defs_reduced <- hanatToAtlasDefs_new(tree_reduced)
  
  return(list(labels = labels_reduced,
              defs = defs_reduced))
  
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





# prepare_histogram_data <- function(x) {
#   
#   regions <- c("cortex", "cerebellum")
#   df_regions <- expand_grid(mouse = regions, 
#                             human = regions)
#   
#   mouse_labels <- list("Isocortex",
#                        "Cerebellar cortex")
#   names(mouse_labels) <- regions
#   
#   human_labels <- list(c("frontal lobe", 
#                          "insula", 
#                          "occipital lobe", 
#                          "parietal lobe", 
#                          "temporal lobe"),
#                        "cerebellar cortex")
#   names(human_labels) <- regions
#   
#   list_regions <- vector(mode = "list", length = nrow(df_regions))
#   for(i in 1:nrow(df_regions)) {
#     
#     mouse_region <- df_regions[[i, "mouse"]]
#     human_region <- df_regions[[i, "human"]]
#     
#     ind_mouse <- colnames(x) %in% mouse_labels[[mouse_region]]
#     ind_human <- rownames(x) %in% human_labels[[human_region]]
#     
#     x_region <- x[ind_human, ind_mouse]
#     
#     list_regions[[i]] <- tibble(similarity = as.numeric(x_region),
#                                 mouse = mouse_region,
#                                 human = human_region)
#   }
#   df <- reduce(.x = list_regions,
#                .f = bind_rows)
#   
#   df <- df %>% 
#     unite(col = "pair",
#           mouse, human,
#           sep = "-",
#           remove = FALSE) %>% 
#     mutate(pair = ifelse(pair %in% c("cortex-cerebellum", 
#                                      "cerebellum-cortex"),
#                          "cortex-cerebellum", 
#                          pair))
#   return(df)
# }





