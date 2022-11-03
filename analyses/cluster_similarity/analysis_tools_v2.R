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
  
  %>% %>% %>% %>% %>% 