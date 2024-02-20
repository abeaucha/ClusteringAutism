# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(SNFtool))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")

# REMOVE THESE LINES WHEN FINISHED
# SRCPATH <- "src"


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "minc_parallel.R"))


#' Fit and predict normative model
#'
#' @param y (numeric vector) Voxel values across study participants.
#' @param demographics (data.frame) Demographics information for study 
#' participants.
#' @param group (character scalar) Group of participants for which to 
#' compute effect sizes.
#' @param batch (character scalar) Batch variable to residualize.
#' @param df (numeric scalar) Degrees of freedom in natural spline 
#' model.
#'
#' @return (data.frame) Model predictions for test participants.
fit_predict_model <- function(y, demographics, group = "patients",
                              batch = NULL, df = 3) {
  
  if (length(y) != nrow(demographics)) {stop()}
  
  # Residualize using batch variable if specified
  if (!is.null(batch)) {
    batch <- demographics %>%
      select(all_of(batch)) %>%
      unite(col = batch) %>%
      pull(batch)
    y <- residuals(lm(y ~ batch))
    names(y) <- NULL
  }
  
  # Filters for train and test sets
  ind_fit <- demographics[["DX"]] == "Control"
  if (group == "patients") {
    ind_pred <- !ind_fit
  } else if (group == "controls") {
    ind_pred <- ind_fit
  } else if (group == "all") {
    ind_pred <- !logical(nrow(demographics))
  }
  
  # Training data frame
  df_fit <- demographics[ind_fit, c("Age", "Sex")]
  df_fit[["y"]] <- y[ind_fit]
  
  # Test data frame
  df_pred <- demographics[ind_pred, c("Age", "Sex")]
  df_pred[["y"]] <- y[ind_pred]
  
  # Fit model and predict on test set
  model_fit <- lm(y ~ Sex + ns(Age, df = df), data = df_fit)
  model_pred <- predict(model_fit,
                        newdata = df_pred,
                        interval = "prediction",
                        level = pnorm(q = 1) - pnorm(q = -1))
  
  # Extract model parameters of interest
  df_pred <- df_pred %>%
    mutate(y_pred = model_pred[,"fit"],
           y_lwr = model_pred[,"lwr"],
           y_upr = model_pred[,"upr"],
           y_sd = y_pred - y_lwr)
  
  return(df_pred)
  
}


#' Compute z-score
#'
#' @param x (data.frame) Data frame containing normative growth 
#' model outputs. 
#'
#' @return (data.frame) Input data frame with new column containing 
#' z-scores.
zscore <- function(x){
  cols_check <- c("y", "y_pred", "y_sd")
  if (any(!(cols_check %in% colnames(x)))){stop()}
  x <- mutate(x, z = (y - y_pred)/y_sd)
  return(x)
}


#' Compute normative z-score for a voxel
#'
#' @param y (numeric vector) Voxel values across study participants.
#' @param demographics (data.frame) Demographics information for study 
#' participants.
#' @param group (character scalar) Group of participants for which to 
#' compute effect sizes.
#' @param batch (character scalar) Batch variable to residualize.
#' @param df (numeric scalar) Degrees of freedom in natural spline 
#' model.
#'
#' @return (numeric vector) Voxel normative z-scores
compute_normative_zscore <- function(y, demographics, group = "patients",
                                     batch = NULL, df = 3) {
  y_pred <- fit_predict_model(y = y, demographics = demographics,
                              group = group, batch = batch, df = df)
  z <- pull(zscore(y_pred), "z")
  return(z)
}



#' Calculate human effect sizes using normative growth modelling.
#'
#' @param imgdir (character scalar) Path to the directory containing the 
#' images (.mnc) to use to compute the effect sizes.
#' @param demographics (character scalar) Path to the file (.csv) 
#' containing the demographics data.
#' @param mask (character scalar) Path to the mask file (.mnc).
#' @param outdir (character scalar) Path to the directory in which to 
#' save the effect size images.
#' @param key (character scalar) Primary key between demographics data 
#' and constructed voxel matrix.
#' @param group (character scalar) Group of participants for which to 
#' compute effect sizes.
#' @param df (numeric scalar) Degrees of freedom to use in normative 
#' model natural splines.
#' @param batch (character scalar) Variables to use in normalization 
#' prior to modelling.
#' @param nbatches (numeric scalar) Number of voxel batches to use for 
#' computation of voxel-wise normative models.
#' @param nproc (numeric scalar) Number of processors to use.
#'
#' @return (character vector) Paths to the effect size images.
normative_growth_norm <- function(imgdir, demographics, mask, outdir,
                                  key = "file", group = "patients",
                                  df = 3, batch = NULL, nbatches = 1,
                                  execution = "local", nproc = 1, 
                                  resources = list()) {
  
  # Import demographics data
  if (verbose) {message("Importing demographics information...")}
  demographics <- as_tibble(data.table::fread(demographics, header = TRUE))
  
  # Check existence of key column in demographics
  if (!(key %in% colnames(demographics))) {
    stop(paste("demographics data is missing key column:", key))
  }
  
  # Remove entries with missing diagnosis, age, or sex
  demographics <- demographics %>%
    filter(!is.na(DX),
           !is.na(Age),
           !is.na(Sex),
           !is.na(Site),
           !is.na(Scanner))
  
  # Check existence of batch columns
  if (!is.null(batch)) {
    batch <- str_split(batch, pattern = "-")[[1]]
    batch_check <- batch %in% colnames(demographics)
    if (!all(batch_check)) {
      stop("Batch columns not found in demographics:\n",
           str_flatten(batch, collapse = "\n"))
    }
  }
  
  # Image files
  imgfiles <- list.files(imgdir, full.names = TRUE)
  
  # Match image files to demographics
  if (verbose) {message("Matching image files to demographics...")}
  imgs_in_demographics <- basename(imgfiles) %in% demographics[[key]]
  imgfiles <- imgfiles[imgs_in_demographics]
  row_match <- match(basename(imgfiles), demographics[[key]])
  demographics <- demographics[row_match,]
  
  # Run normative growth modelling
  if (verbose) {message("Evaluating normative growth models...")}
  if (execution == "local") {
    print("Executing locally...")
    voxels <- mcMincApply(filenames = imgfiles,
                          fun = compute_normative_zscore,
                          demographics = demographics,
                          group = group,
                          batch = batch,
                          df = df,
                          mask = mask,
                          cores = nproc,
                          return_raw = TRUE)
  } else if (execution == "slurm") {
    print("Executing on HPC...")
    voxels <- qMincApply(filenames = imgfiles,
                         fun = compute_normative_zscore,
                         demographics = demographics,
                         group = group,
                         batch = batch,
                         df = df,
                         mask = mask,
                         batches = nproc,
                         source = c(file.path(SRCPATH, "pipelines/processing.R")),
                         cleanup = FALSE,
                         return_raw = TRUE,
                         resources = resources)
  } else {
    stop()
  }
  
  # Convert voxel list into matrix
  voxels <- simplify_masked(voxels[["vals"]])
  
  # Export images
  if (verbose) {message("Exporting normalized images...")}
  if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
  if (group == "patients") {
    outfiles <- demographics[demographics[["DX"]] != "Control", key][[1]]
  } else if (group == "controls") {
    outfiles <- demographics[demographics[["DX"]] == "Control", key][[1]]
  } else if (group == "all") {
    outfiles <- demographics[[key]]
  }
  outfiles <- file.path(outdir, outfiles)
  matrix_to_images(x = voxels, outfiles = outfiles, mask = mask,
                   margin = 2, nproc = nproc)
  
  return(outfiles)
}


#' Calculate human effect sizes using propensity-matching.
#'
#' @param imgdir (character scalar) Path to the directory containing the 
#' images (.mnc) to use to compute the effect sizes.
#' @param demographics (character scalar) Path to the file (.csv) 
#' containing the demographics data.
#' @param mask (character scalar) Path to the mask file (.mnc).
#' @param outdir (character scalar) Path to the directory in which to 
#' save the effect size images.
#' @param ncontrols (numeric scalar) Number of propensity-matched 
#' controls to use when computing the effect sizes.
#' @param nproc (numeric scalar) Number of processors to use.
#'
#' @return (character vector) Paths to the effect size images.
propensity_matching_norm <- function(imgdir, demographics, mask, outdir,
                                     ncontrols = 10, nproc = 1) {
  return(outfiles)
}


#' Run similarity network fusion (SNF)
#'
#' @param x1 (matrix) Input matrix
#' @param x2 (matrix) Input matrix
#' @param metric (character scalar) Distance metric used to compute the 
#' affinity matrices.
#' @param K (numeric scalar) Number of nearest-neighbours used to 
#' compute the SNF affinity matrices. 
#' @param sigma (numeric scalar) Variance for the local model in the 
#' SNF affinity matrices.
#' @param t (numeric scalar) Number of iterations for the diffusion
#' process in SNF.
#' @param outfile (character scalar) Path to file in which to save
#' affinity matrix.
#'
#' @return (matrix) SNF affinity matrix.
similarity_network <- function(x1, x2, metric = "correlation", K = 10,
                               sigma = 0.5, t = 20, outfile = NULL){
  
  if (metric == "correlation") {
    d1 <- (1-cor(t(x1)))
    d2 <- (1-cor(t(x2)))
  } else if (metric == "euclidean") {
    d1 <- (dist2(as.matrix(x1), as.matrix(x1)))^(1/2)
    d2 <- (dist2(as.matrix(x2), as.matrix(x2)))^(1/2)
  } else {
    stop(paste("Argument metric must be one of {correlation, euclidean}:", metric))
  }
  
  W1 <- affinityMatrix(d1, K = K, sigma = sigma)
  W2 <- affinityMatrix(d2, K = K, sigma = sigma)
  
  W <- SNF(list(W1, W2), K = K, t = t)
  
  if (!is.null(outfile)){
    data.table::fwrite(x = as_tibble(W), file = outfile)
  }
  
  return(W)
  
}


#' Create clusters from SNF affinity matrix
#'
#' @param W (matrix) SNF affinity matrix.
#' @param nk (numeric scalar) Maximum number of clusters to use in 
#' clustering.
#' @param outfile (character scalar) Optional path to .csv file in 
#' which to save cluster assignments.
#'
#' @return (data.frame) Cluster assignments.
create_clusters <- function(W, nk = 10, outfile = NULL) {
  
  for(k in 2:nk) {
    group <- spectralClustering(affinity = W, K = k)
    group_name <- paste0("nk", k)
    assign(group_name, group)
    if (k == 2) {
      if (is.null(rownames(W))) {
        ids <- as.character(1:nrow(W))
      } else {
        ids <- rownames(W)
      }
      all_clusters <- data.frame(ids, group, stringsAsFactors = F)
      colnames(all_clusters) <- c("ID", group_name)
    } else {
      group <- data.frame(group)
      colnames(group) <- group_name
      all_clusters <- cbind(all_clusters, group)
    }
  }
  
  if (!is.null(outfile)) {
    write.csv(x = all_clusters, file = outfile, row.names = FALSE)
  }
  
  return(all_clusters)
  
}