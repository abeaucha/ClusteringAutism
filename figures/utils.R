suppressPackageStartupMessages(library(tidyverse))


#' Import the representative map for a cluster.
#'
#' @param datadir (character scalar) Path to the clustering directory.
#' @param nk (numeric scalar) The number of clusters in the solution.
#' @param k (numeric scalar) The cluster ID number.
#' @param jacobians (character scalar) The Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param mask (character scalar) The type of mask to apply to the 
#' cluster maps. Ignored if NULL. 
#' Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) The thresholding method 
#' used to create the mask. Ignored if mask = NULL. 
#' Options: {intensity, topn}.
#' @param threshold (numeric scalar) The threshold value applied to 
#' create the mask. Ignored if mask = NULL.
#'
#' @return (mincSingleDim) The cluster map.
import_cluster_map <- function(datadir, species, nk, k, jacobians = NULL, mask = NULL, 
                               threshold_method = NULL, threshold = NULL) {
  
  # TO-DO
  # -Fix directories so that I don't have to have this species 
  #  argument in order for this to work. 
  
  #Check if data directory exists
  mapdir <- file.path(datadir, "cluster_maps", "")
  msg_mapdir <- paste("Data directory not found:", mapdir)
  if (!dir.exists(mapdir)) {stop(msg_mapdir)}
  
  #Argument options: jacobians
  opts_jacobians <- c("absolute", "relative")
  msg_jacobians <- paste("Argument jacobians must be one of",
                         "{\"absolute\", \"relative\"}.")
  if (is.null(jacobians)){
    stop(msg_jacobians)
  } else if (!(jacobians %in% opts_jacobians)) {
    stop(msg_jacobians)    
  }
  
  #Extend path to maps directory
  if (species == "mouse") {
    mapdir <- file.path(mapdir, jacobians, "resolution_200", "mean", "")
  } else if (species == "human") {
    mapdir <- file.path(mapdir, jacobians, "resolution_1.0", "mean", "")
  } else {
    stop("species must be one of {mouse, human}")
  }
  if (!dir.exists(mapdir)) {stop(msg_mapdir)}
  
  #Import cluster map
  infile <- list.files(mapdir, full.names = TRUE) %>% 
    str_subset(str_c("Clusternum_", nk)) %>%
    str_subset(str_c("Group_", k))
  cluster <- mincGetVolume(infile)
  
  #Option to apply mask
  if (!is.null(mask)) {
    opts_mask <- c("symmetric", "positive", "negative")
    msg_mask <- paste("Argument mask must be one of",
                      "{\"symmetric\", \"positive\", \"negative\"}.")
    if (!(mask %in% opts_mask)) {stop(msg_mask)}
    mask <- import_cluster_mask(datadir = datadir,
                                species = species,
                                nk = nk, k = k,
                                jacobians = jacobians,
                                type = mask,
                                threshold_method = threshold_method,
                                threshold = threshold)
    cluster[mask == 0] <- 0
  }
  
  return(cluster)
  
}


#' Import a cluster mask
#'
#' @param datadir (character scalar) Path to the clustering directory.
#' @param nk (numeric scalar) The number of clusters in the solution.
#' @param k (numeric scalar) The cluster ID number.
#' @param jacobians (character scalar) The Jacobians used to compute 
#' the cluster maps. Options: {absolute, relative}.
#' @param type (character scalar) The type of mask to return. 
#' Options: {symmetric, positive, negative}.
#' @param threshold_method (character scalar) The thresholding method 
#' used to create the mask. Options: {intensity, topn}.
#' @param threshold (numeric scalar) The threshold value applied to 
#' create the mask.
#'
#' @return (mincSingleDim) The cluster mask.
import_cluster_mask <- function(datadir, species, nk, k, jacobians = NULL, 
                                type = "symmetric", threshold_method = NULL, 
                                threshold = NULL) {
  
  # TO-DO
  # -Fix directories so that I don't have to have this species 
  #  argument in order for this to work. 
  
  #Check if data directory exists
  maskdir <- file.path(datadir, "cluster_masks", "")
  msg_maskdir <- paste("Data directory not found:", maskdir)
  if (!dir.exists(maskdir)) {stop(msg_maskdir)}
  
  #Argument options: jacobians
  opts_jacobians <- c("absolute", "relative")
  msg_jacobians <- paste("Argument jacobians must be one of",
                         "{\"absolute\", \"relative\"}.")
  if (is.null(jacobians)){
    stop(msg_jacobians)
  } else if (!(jacobians %in% opts_jacobians)) {
    stop(msg_jacobians)    
  }
  
  #Argument options: threshold_method  
  opts_threshold_method <- c("intensity", "topn")
  msg_threshold_method <- paste("Argument threshold_method must be one of",
                                "{\"intensity\", \"topn\"}.")
  if (is.null(threshold_method)){
    stop(msg_threshold_method)
  } else if (!(threshold_method %in% opts_threshold_method)) {
    stop(msg_threshold_method)    
  }
  
  #Extend path to mask directory
  if (species == "mouse") {
    maskdir <- str_c(maskdir, jacobians, "/resolution_200/mean/threshold_", 
                     threshold_method, "/symmetric/")
  } else if (species == "human") {
    maskdir <- str_c(maskdir, jacobians, "/resolution_1.0/mean/threshold_", 
                     threshold_method, "/symmetric/")
  } else {
    stop("species must be one of {mouse, human}")
  }
  if (!dir.exists(maskdir)) {stop(msg_maskdir)}
  
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
    opts_type <- c("symmetric", "positive", "negative")
    msg_type <- paste("Argument type must be one of",
                      "{\"symmetric\", \"positive\", \"negative\"}.")
    stop(msg_type)
  }
  
  return(mask)
  
}


#' Import cluster similarity matrix
#'
#' @param datadir (character scalar) Path to the directory containing 
#' the similarity matrix CSV file.
#' @param gene_space (character scalar) Gene expression space used to 
#' compute the similarity matrix. Options: {average_latent_space, 
#' latent_space, input_space}.
#' @param jacobians (character scalar) Jacobians used to compute 
#' the expression signatures and similarity matrix. 
#' Options: {absolute, relative}.
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
import_similarity_matrix <- function(datadir, 
                                     gene_space = "average_latent_space", 
                                     jacobians = NULL,
                                     threshold_method = NULL,
                                     threshold = NULL,
                                     mask = NULL,
                                     latent_space_id = 1) {
  
  #Argument options: jacobians
  opts_jacobians <- c("absolute", "relative")
  msg_jacobians <- "Argument jacobians must be one of {\"absolute\", \"relative\"}."
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
  msg_threshold_method <- "Argument threshold_method must be one of {\"intensity\", \"topn\"}."
  if (is.null(threshold_method)){
    stop(msg_threshold_method)
  } else if (!(threshold_method %in% opts_threshold_method)) {
    stop(msg_threshold_method)    
  }
  
  msg_threshold <- "Argument threshold is NULL. Specify a threshold."
  if (is.null(threshold)){
    stop(msg_threshold)
  }
  
  #Argument options: mask
  opts_mask <- c("symmetric", "positive", "negative")
  msg_mask <- paste("Argument mask must be one of {\"symmetric\", \"positive\", \"negative\"}.")
  if (is.null(mask)){
    stop(msg_mask)
  } else if (!(mask %in% opts_mask)) {
    stop(msg_mask)    
  }
  
  if (gene_space == "input_space") {
    #Implement this feature
    stop("Feature currently unavailable")    
  } else if (gene_space == "latent_space") {
    #Implement this feature
    stop("Feature currently unavailable")    
  } else if (gene_space == "average_latent_space") {
    
    infile <- str_c("similarity_hm_", jacobians, "_mean_threshold_", 
                    threshold_method, "_", threshold, "_", mask,
                    "_latentspace100.csv")
    infile <- file.path(datadir, infile)
    
    if (!(file.exists(infile))){
      msg_infile <- paste("Input file not found:", infile)
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