library(tools)

# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")

# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "processing.R"))


compare_centroid_anatomy <- function(centroid_dirs, species, nk, k, nodes, trees, 
                                     labels, defs, masks, threshold, threshold_value,
                                     threshold_symmetric, threshold_comparison = NULL) {
  
  freq <- vector(mode = "list", length = 2)
  for (i in 1:2) {
    
    if (species[i] == "human") {
      reduce_atlas <- reduce_human_atlas
    } else {
      reduce_atlas <- reduce_mouse_atlas
    }
    
    atlas <- reduce_atlas(tree = trees[[i]],
                          labels = labels[[i]],
                          defs = defs[[i]],
                          nodes = nodes[[species[i]]],
                          remove = TRUE)
    
    freq_signed <- vector(mode = "list", length = 2)
    names(freq_signed) <- c("positive", "negative")
    for (s in names(freq_signed)) {
      freq_signed[[s]] <- compute_cluster_fractions(cluster_dir = centroid_dirs[[i]],
                                                    nk = nk[[i]],
                                                    k = k[[i]],
                                                    labels = atlas[["labels"]],
                                                    defs = atlas[["defs"]],
                                                    mask = masks[[i]],
                                                    sign = s,
                                                    threshold = threshold,
                                                    threshold_value = threshold_value,
                                                    threshold_symmetric = threshold_symmetric,
                                                    threshold_comparison = threshold_comparison)
    }
    
    freq[[i]] <- bind_rows(freq_signed, .id = "sign")
    
  }
  
  freq <- intersect_neuro_homologues(x = freq, 
                                     species = species,
                                     homologues = nodes)
  
  return(freq)
  
}


#' Compute the fraction of voxels
#'
#' @param cluster_dir (character scalar) Path to the directory 
#' containing cluster images.
#' @param nk (integer scalar) Number of clusters in the cluster 
#' solution. 
#' @param k (integer scalar) Cluster ID in the solution set.
#' @param labels (mincSingleDim, numeric) Atlas labels.
#' @param defs (data.frame) Atlas definitions.
#' @param mask (character scalar or NULL) Path to a mask image.
#' @param sign (character scalar or NULL) Option to tabulate using 
#' positive or negative intensity voxels only. 
#' @param threshold (character scalar or NULL) Method used to threshold 
#' the cluster centroid image. No thresholding is applied if NULL.
#' @param threshold_value (numeric scalar or NULL) Value at which to 
#' threshold the image.
#' @param threshold_symmetric (logical scalar or NULL) Option to apply 
#' the threshold symmetrically. 
#' @param threshold_comparison (character scalar or NULL) String indicating 
#' how to apply the threshold
#'
#' @return (data.frame) Cluster fractions
compute_cluster_fractions <- function(cluster_dir, nk, k, labels, defs, mask,
                                      sign = "both", threshold = NULL, 
                                      threshold_value = NULL, threshold_symmetric = NULL,
                                      threshold_comparison = NULL) {
  
  #Import cluster map
  cluster_map <- import_cluster_map(imgdir = cluster_dir,
                                    mask = mask,
                                    nk = nk, k = k,
                                    threshold = threshold,
                                    threshold_value = threshold_value, 
                                    threshold_symmetric = threshold_symmetric, 
                                    threshold_comparison = threshold_comparison,
                                    flatten = TRUE)
  
  #Compute cluster anatomical fractions
  cluster_fractions <- tabulate_labels_in_img(img = cluster_map,
                                              mask = mask, 
                                              labels = labels,
                                              defs = defs,
                                              sign = sign)
  
  #Include cluster ID
  cluster_fractions[["cluster_id"]] <- str_c(nk, k, sep = "-")
  
  return(cluster_fractions)
  
}


compute_enrichment_HGtest <- function(mouse_enrichment, human_enrichment, alpha_mouse, alpha_human = NULL) {
  
  if (is.null(alpha_human)) {
    alpha_human <- alpha_mouse
  } 
  
  N <- nrow(mouse_enrichment)
  
  df_hgtest <- expand_grid(alpha_mouse = alpha_mouse,
                           alpha_human = alpha_human) %>% 
    mutate(b = 0, B = 0, n = 0, N = N, E = 0, pval = 0)
  for (i in 1:nrow(df_hgtest)) {
    mouse_pass <- mouse_enrichment[["adj.P.Val"]] < df_hgtest[[i, "alpha_mouse"]]
    human_pass <- human_enrichment[["adj.P.Val"]] < df_hgtest[[i, "alpha_human"]]
    B <- sum(mouse_pass)
    n <- sum(human_pass)
    b <- sum(mouse_pass & human_pass)
    pval <- phyper(b-1, B, N-B, n, lower.tail = FALSE)
    df_hgtest[[i, "b"]] <- b
    df_hgtest[[i, "B"]] <- B
    df_hgtest[[i, "n"]] <- n
    df_hgtest[[i, "E"]] <- (b/B)/(n/N)
    df_hgtest[[i, "pval"]] <- pval
  }
  return(df_hgtest)
}

compute_enrichment_PR <- function(mouse_enrichment, human_enrichment, alpha = 0.05) {
  out <- vector(mode = "list", length = length(alpha))
  names(out) <- alpha
  for (i in 1:length(alpha)) {
    response <- mouse_enrichment[["adj.P.Val"]] < alpha[i]
    if (sum(response) > 0) {
      predictor <- 1 - human_enrichment[["adj.P.Val"]]
      out[[i]] <- PRROC::pr.curve(scores.class0 = predictor[response],
                                  scores.class1 = predictor[!response],
                                  curve = TRUE)
    } else {
      out[[i]] <- list(type = "PR", auc.integral = NA, 
                       auc.davis.goadrich = NA, curve = NA)
    }
    out[[i]] <- c(out[[i]], true.proportion = sum(response)/length(response))
  }
  return(out)
}


compute_enrichment_ROC <- function(mouse_enrichment, human_enrichment, alpha = 0.05) {
  out <- vector(mode = "list", length = length(alpha))
  names(out) <- alpha
  for (i in 1:length(alpha)) {
    response <- mouse_enrichment[["adj.P.Val"]] < alpha[i]
    if (sum(response) > 0) {
      predictor <- 1 - human_enrichment[["adj.P.Val"]]
      out[[i]] <- pROC::roc(response = response, predictor = predictor, quiet = TRUE)
    } else {
      out[[i]] <- NA
    }
  }
  return(out)
}


simulate_enrichment_null <- function(enrichment, n_true, method = "ROC", n_sim = 500){
  
  if (!(method %in% c("ROC", "PR", "HG"))) {stop()}
  
  # Iterate over number of pathways enriched
  out <- vector(mode = "list", length = length(n_true))
  for (i in 1:length(n_true)) {
    
    # Iterate over simulations
    out_i <- vector(mode = "list", length = n_sim)
    for (s in 1:n_sim) {
      
      # Randomly sample which pathways are enriched for the given number 
      target <- mutate(enrichment, adj.P.Val = 1)
      set.seed(s)
      idx_flip <- sample(1:nrow(target), size = n_true[i], replace = FALSE)
      target$adj.P.Val[idx_flip] <- 0
      
      if (method == "ROC") {
        out_i[[s]] <- compute_enrichment_ROC(mouse_enrichment = target, 
                                             human_enrichment = enrichment,
                                             alpha = 0.05) %>% 
          map_dbl(pROC::auc) %>% 
          enframe(name = "alpha", value = "AUC") %>% 
          mutate(n_enriched = n_true[i]) %>% 
          select(-alpha)
      } else if (method == "PR") {
        out_i[[s]] <- compute_enrichment_PR(mouse_enrichment = target, 
                                            human_enrichment = enrichment,
                                            alpha = 0.05) %>% 
          map_dfr(.f = function(x){tibble(AUC = x[["auc.integral"]], 
                                          Prop = x[["true.proportion"]])},
                  .id = "alpha") %>% 
          mutate(n_enriched = n_true[i]) %>% 
          select(-alpha)
      } else if (method == "HG") {
        out_i[[s]] <- compute_enrichment_HGtest(mouse_enrichment = target, 
                                                human_enrichment = enrichment,
                                                alpha_mouse = 0.05,
                                                alpha_human = c(0.05, 1e-3, 1e-5, 1e-10))
      } else {
        stop()
      }
    }
    
    out[[i]] <- bind_rows(out_i, .id = "sample")
    
  }
  
  return(bind_rows(out))
  
}


compute_similarity_significance <- function(similarity, permutations, off_diag = 1) {
  
  nk_max_1 <- max(similarity[["img1_nk"]])
  nk_max_2 <- max(similarity[["img2_nk"]])
  
  significance <- tibble()
  for (nk_1 in 2:nk_max_1) {
    for (nk_2 in (nk_1 - off_diag):(nk_1 + off_diag)) {
      
      if ((nk_2 > 1) & (nk_2 <= nk_max_2)) {
        
        sim_nk <- similarity %>% 
          select(img1_cluster_id, img1_nk, img1_k,
                 img2_cluster_id, img2_nk, img2_k,
                 similarity) %>% 
          filter(img1_nk == nk_1,
                 img2_nk == nk_2) %>% 
          mutate(pval = 0)
        
        sim_perm_nk <- permutations %>% 
          filter(img1_nk == nk_1,
                 img2_nk == nk_2) %>% 
          pull(similarity) %>% 
          sort()
        
        for (i in 1:nrow(sim_nk)) {
          ntail <- sum(sim_perm_nk >= sim_nk[[i, "similarity"]])
          sim_nk[[i, "pval"]] <- ntail/length(sim_perm_nk)
        }
        
        significance <- bind_rows(significance, sim_nk)
        
      }
    } 
  }
  
  return(significance)
  
}


#' Export grid grob to PDF
#'
#' @param x (grob) Grob to export
#' @param file (character scalar) Output file
#' @param width (numeric scalar) Width of PDF
#' @param height (numeric scalar) Height of PDF 
#' @param units (character scalar) Dimension units
#'
#' @return NULL
export_pdf <- function(x, file, width, height, units = "in") {
  
  if (!("grob" %in% class(x))) {
    stop("Argument \"x\" must be a grob")
  }
  
  if (units == "in") {
    width <- width
    height <- height
  } else if (units == "inches") {
    width <- width
    height <- height
  } else if (units == "bigpts") {
    pt_per_in <- 72
    width <- width/pt_per_in
    height <- height/pt_per_in
  } else {
    stop("Unknown units: ", units)
  }
  
  if (file_ext(file) != "pdf") {
    stop("Argument \"file\" must have extension .pdf")
  }
  
  tryCatch(
    {pdf(file = file,
         width = unit(width, "in"),
         height = unit(height, "in"))
      grid.draw(x)
      dev.off()}, 
    error = function(c) {
      message("Error: ", conditionMessage(c))
    }
  )
}


#' Compute within-cluster sum of squared distances
#'
#' @param x (matrix) A matrix to cluster according to columns
#' @param method (character scalar) Distance metric used for 
#' hierarchical clustering
#'
#' @return (numeric vector) Within-cluster sum of squared distances
hclust_wcss <- function(x, method = "euclidean") {
  hc <- hclust(dist(x = t(x), method = method))
  wcss <- numeric(ncol(x))
  for (nk in 1:length(wcss)) {
    labs <- cutree(hc, k = nk)
    wcss_nk <- numeric(nk)
    for (k in 1:nk) {
      xk <- as.matrix(x[,labs == k])
      wcss_nk[k] <- sum(colSums((xk - rowMeans(xk))^2))/ncol(xk)
    }
    wcss[nk] <- mean(wcss_nk)
  }
  return(wcss)
}


import_cluster_map <- function(imgdir, nk, k, mask = NULL, flatten = TRUE, threshold = NULL, threshold_value = NULL, threshold_symmetric = NULL, threshold_comparison = NULL) {
  
  pattern = paste("nk", nk, "k", k, sep = "_")
  img <- list.files(imgdir, full.names = TRUE, pattern = pattern)
  img <- import_image(img = img, mask = mask, flatten = FALSE)
  
  if (!is.null(threshold)) {
    img <- threshold_image(img = img, 
                           method = threshold, 
                           threshold = threshold_value,
                           symmetric = threshold_symmetric,
                           comparison = threshold_comparison)
  }
  
  if (flatten) {
    img_attr <- attributes(img)
    img_attr[["dim"]] <- NULL
    img <- as.numeric(img)
    attributes(img) <- img_attr
  }
  
  return(img)
}

import_enrichment_human <- function(params_id, pipeline_dir = "data/human/derivatives/v3/",
                                    nk, k, gene_score = 950, stringdb_version = "12.0",
                                    bader_version = 2023, file_prefix = "cluster_pathway_enrichment"){
  
  enrichment_dir <- file.path(pipeline_dir, params_id, "enrichment", 
                              paste("StringDB", stringdb_version, 
                                    "Bader", bader_version, sep = "_"), 
                              gene_score)
  
  enrichment_file <- paste(file_prefix, nk, k, gene_score, sep = "_")
  enrichment_file <- paste0(enrichment_file, ".csv")
  enrichment_file <- file.path(enrichment_dir, enrichment_file)
  
  enrichment <- read_csv(enrichment_file, show_col_types = FALSE) %>% 
    arrange(rank)
  
  return(enrichment)
}

import_enrichment_mouse <- function(params_id, pipeline_dir = "data/human/derivatives/v3/",
                                    nk, k, gene_score = 950, stringdb_version = "12.0",
                                    bader_version = 2023, file_prefix = "NewBader_enrichment_clusterneighbourhood_vs_brain_all"){
  
  enrichment_dir <- file.path(pipeline_dir, params_id, "enrichment", 
                              paste("StringDB", stringdb_version, 
                                    "Bader", bader_version, sep = "_"), 
                              "NeighbourhoodEnrichment", gene_score)
  
  enrichment_file <- paste(file_prefix, nk, k, gene_score, sep = "_")
  enrichment_file <- paste0(enrichment_file, ".csv")
  enrichment_file <- file.path(enrichment_dir, enrichment_file)
  
  enrichment <- read_csv(enrichment_file, show_col_types = FALSE) %>% 
    arrange(rank)
  
  return(enrichment)
  
}


import_similarity <- function(param_id, pipeline_dir = "data/cross_species/v3/", combine_jacobians = TRUE) {
  
  input_dir <- file.path(pipeline_dir, param_id, "similarity")
  input_file <- file.path(input_dir, "similarity.csv")
  
  # Import similarity data
  similarity <- read_csv(input_file, show_col_types = FALSE) %>% 
    mutate(img1_nk = img1 %>% 
             basename() %>% 
             str_extract("_nk_[0-9]+") %>% 
             str_extract("[0-9]+") %>% 
             as.numeric(),
           img1_k = img1 %>% 
             basename() %>% 
             str_extract("_k_[0-9]+") %>% 
             str_extract("[0-9]+") %>% 
             as.numeric(),
           img1_jacobians = img1 %>% 
             str_extract("absolute|relative"),
           img2_nk = img2 %>% 
             basename() %>% 
             str_extract("_nk_[0-9]+") %>% 
             str_extract("[0-9]+") %>% 
             as.numeric(),
           img2_k = img2 %>% 
             basename() %>% 
             str_extract("_k_[0-9]+") %>% 
             str_extract("[0-9]+") %>% 
             as.numeric(),
           img2_jacobians = img2 %>% 
             str_extract("absolute|relative")) %>% 
    unite(col = "img1_cluster_id", img1_nk, img1_k, 
          sep = "-", remove = FALSE) %>% 
    unite(col = "img2_cluster_id", img2_nk, img2_k, 
          sep = "-", remove = FALSE)
  
  # Filter similarity data for desired cluster numbers
  # and combine Jacobians
  if (combine_jacobians) {
    similarity <- similarity %>% 
      group_by(img1_cluster_id, img1_nk, img1_k, 
               img2_cluster_id, img2_nk, img2_k) %>% 
      summarise(similarity = mean(similarity),
                .groups = "drop")
  }
  
  return(similarity)
  
}


import_similarity_permutations <- function(param_id, pipeline_dir, combine_jacobians = TRUE) {
  
  input_dir <- file.path(pipeline_dir, param_id, "permutations", "similarity")
  input_files <- list.files(input_dir)
  
  np <- length(input_files)
  list_permutations <- vector(mode = "list", length = np)
  for (i in 1:length(list_permutations)) {
    
    # Input file  
    input_file <- input_files[i]
    
    # Permutation number
    p <- input_file %>% 
      str_extract("[0-9]+") %>% 
      as.numeric()
    
    input_file <- file.path(input_dir, input_file)
    
    # Import permutation i
    list_permutations[[i]] <- read_csv(input_file, show_col_types = FALSE) %>% 
      mutate(img1_nk = img1 %>% 
               basename() %>% 
               str_extract("_nk_[0-9]+") %>% 
               str_extract("[0-9]+") %>% 
               as.numeric(),
             img1_k = img1 %>% 
               basename() %>% 
               str_extract("_k_[0-9]+") %>% 
               str_extract("[0-9]+") %>% 
               as.numeric(),
             img1_jacobians = img1 %>% 
               str_extract("absolute|relative"),
             img2_nk = img2 %>% 
               basename() %>% 
               str_extract("_nk_[0-9]+") %>% 
               str_extract("[0-9]+") %>% 
               as.numeric(),
             img2_k = img2 %>% 
               basename() %>% 
               str_extract("_k_[0-9]+") %>% 
               str_extract("[0-9]+") %>% 
               as.numeric(),
             img2_jacobians = img2 %>% 
               str_extract("absolute|relative")) %>% 
      unite(col = "img1_cluster_id", img1_nk, img1_k, 
            sep = "-", remove = FALSE) %>% 
      unite(col = "img2_cluster_id", img2_nk, img2_k, 
            sep = "-", remove = FALSE) %>% 
      mutate(permutation = p)
    
  }
  
  df_permutations <- bind_rows(list_permutations)
  
  if (combine_jacobians) {
    df_permutations <- df_permutations %>% 
      group_by(permutation,
               img1_cluster_id, img1_nk, img1_k, 
               img2_cluster_id, img2_nk, img2_k) %>% 
      summarise(similarity = mean(similarity),
                .groups = "drop")
  }
  
  return(df_permutations)
}


#' Intersect data with neuroanatomical homologues
#'
#' @param x (list) A list containing the data to intersect. 
#' @param species (character vector) 
#' @param homologues (data.frame) A data frame containing the set
#' of mouse and human neuroanatomical homologues.
#'
#' @return (list) A list containing the intersected data.
intersect_neuro_homologues <- function(x, species, homologues) {
  for (i in 1:length(x)) {
    cols_init <- colnames(x[[i]])
    colnames(x[[i]])[cols_init == "name"] <- species[i]
    x[[i]] <- inner_join(x[[i]], homologues, by = species[i])
    x[[i]][["name"]] <- factor(x[[i]][["name"]], levels = homologues[["name"]])
    x[[i]][["species"]] <- species[i]
    x[[i]] <- x[[i]][,match(cols_init, colnames(x[[i]]))]
  }
  return(x)
}

# intersect_neuro_homologues_old <- function(x, homologues) {
#   for (i in 1:length(x)) {
#     species <- names(x)[[i]]
#     cols_init <- colnames(x[[i]])
#     colnames(x[[i]])[cols_init == "name"] <- species
#     x[[i]] <- inner_join(x[[i]], homologues, by = species)
#     x[[i]][["name"]] <- factor(x[[i]][["name"]], levels = homologues[["name"]])
#     x[[i]][["species"]] <- species
#     x[[i]] <- x[[i]][,match(cols_init, colnames(x[[i]]))]
#   }
#   return(x)
# }


prepare_radar_chart <- function(cluster_dirs, nk, k, spokes, trees, labels, 
                                defs, masks, threshold, threshold_value, 
                                threshold_symmetric, threshold_comparison = NULL) {
  
  species <- c("human", "mouse")
  fractions <- vector(mode = "list", length = length(species))
  names(fractions) <- species
  for (s in species) {
    
    if (s == "human") {
      reduce_atlas <- reduce_human_atlas
    } else {
      reduce_atlas <- reduce_mouse_atlas
    }
    
    atlas <- reduce_atlas(tree = trees[[s]],
                          labels = labels[[s]],
                          defs = defs[[s]],
                          nodes = spokes[[s]],
                          remove = TRUE)
    
    fractions_positive <- compute_cluster_fractions(cluster_dir = cluster_dirs[[s]],
                                                    nk = nk[[s]],
                                                    k = k[[s]],
                                                    labels = atlas[["labels"]],
                                                    defs = atlas[["defs"]],
                                                    mask = masks[[s]],
                                                    sign = "positive",
                                                    threshold = threshold,
                                                    threshold_value = threshold_value,
                                                    threshold_symmetric = threshold_symmetric,
                                                    threshold_comparison = threshold_comparison)
    
    fractions_negative <- compute_cluster_fractions(cluster_dir = cluster_dirs[[s]],
                                                    nk = nk[[s]],
                                                    k = k[[s]],
                                                    labels = atlas[["labels"]],
                                                    defs = atlas[["defs"]],
                                                    mask = masks[[s]],
                                                    sign = "negative",
                                                    threshold = threshold,
                                                    threshold_value = threshold_value,
                                                    threshold_symmetric = threshold_symmetric,
                                                    threshold_comparison = threshold_comparison)
    
    fractions_positive <- select(fractions_positive, name, positive = f_per_label)
    fractions_negative <- select(fractions_negative, name, negative = f_per_label)
    
    fractions[[s]] <- inner_join(fractions_positive, 
                                 fractions_negative, 
                                 by = "name") %>% 
      mutate(negative = -1*negative, 
             species = s)
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

#' Create reduced human atlas definitions
#'
#' @param tree (data.tree) Tree whose leaves will be new labels
#' @param defs (data.frame) Definitions containing labels for all human microarray samples
#' @param simplify (logical scalar) Option to simplify returned data frame
#'
#' @return (data.frame) Reduced atlas label definitions
reduce_human_defs <- function(tree, defs, simplify = TRUE) {
  
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
  labels_reduced <- defs_reduced[ind_match, "label"][[1]]
  labels_reduced[is.na(labels_reduced)] <- 0
  attributes(labels_reduced) <- attributes(labels)
  
  return(labels_reduced)
}


#' Create a reduced version of a human atlas. 
#'
#' @param tree (data.tree) A data tree containing the neuroanatomical
#' hierarchy. Leaf nodes must correspond to the atlas regions.
#' @param labels (mincSingleDim, numeric) Atlas labels.
#' @param defs (data.frame) Atlas definitions.
#' @param nodes (character vector) Tree nodes at which to aggregate the
#' atlas.
#' @param remove (logical scalar) Option to remove extra nodes not 
#' specified in `nodes`.
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


#' Tabulate atlas labels in an image
#'
#' This function will tabulate the atlas labels using non-zero voxels
#' in the input image. 
#'
#' @param img (mincSingleDim, numeric) Input image.
#' @param labels (mincSingleDim, numeric) Atlas labels.
#' @param defs (data.frame) Atlas definitions.
#' @param mask (character scalar or NULL) Path to a mask image.
#' @param sign (character scalar or NULL) Option to tabulate using 
#' positive or negative intensity voxels only. 
#'
#' @return (data.frame) Number and fraction of atlas voxels found in 
#' the image.
tabulate_labels_in_img <- function(img, labels, defs, mask = NULL, sign = NULL) {
  
  #Create a mask from the image
  img_mask <- mask_from_image(img = img, signed = TRUE)
  
  #Filter for signed voxels if specified
  if (is.null(sign)) {
    img_mask <- abs(img_mask)
  } else {
    if (sign == "positive") {
      img_mask[img_mask < 1] <- 0
    } else if (sign == "negative") {
      img_mask[img_mask > -1] <- 0
      img_mask <- abs(img_mask)
    } else {
      stop("`sign` must be NULL or one of {'positive', 'negative'}")
    }
  }
  
  #Voxels in the image
  voxels_in_img <- img_mask > 0.5
  
  #Mask if specified
  if (!is.null(mask)) {
    mask <- round(mincGetVolume(mask))
    voxels_in_mask <- mask > 0.5
    voxels_in_img <- voxels_in_img & voxels_in_mask
  }
  
  #Nonzero labels
  labels_nonzero <- labels != 0
  
  #Nonzero labels in the image
  labels_in_img <- labels[voxels_in_img & labels_nonzero]
  
  #Tabulate number of labels in the image
  if (length(labels_in_img) != 0) {
    out <- table(labels_in_img) %>% 
      as_tibble() %>% 
      rename(label = labels_in_img, 
             n_per_label = n) %>% 
      mutate(label = as.integer(label))
  } else {
    out <- tibble(label = defs[["label"]],
                  n_per_label = 0)
  }
  
  #Compute label fractions
  out <- out %>% 
    right_join(defs, by = "label") %>% 
    mutate(n_per_label = ifelse(is.na(n_per_label), 0, n_per_label),
           n_labels_in_img = length(labels_in_img),
           f_per_label = n_per_label/n_labels_in_img) %>% 
    select(name, label, n_per_label, f_per_label, n_labels_in_img)
  
  return(out)
}


#Radar coordinate system for ggplot2
coord_radar <- function (theta = "x", start = 0, direction = 1, clip = "on") {
  
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
          clip = clip,
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