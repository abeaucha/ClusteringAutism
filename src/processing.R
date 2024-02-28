# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(SNFtool))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doSNOW))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "minc_parallel.R"))

import_image <- function(img, mask = NULL, flatten = TRUE) {
  
  # Import image
  img <- mincGetVolume(img)
  
  # Convert to 3D if specified
  if (!flatten) {
    img <- mincArray(img)
  }
  
  # Apply mask if specified
  if (!is.null(mask)) {
    mask <- mincGetVolume(mask)
    if (length(img) != length(mask)) {
      stop("Input image and mask contain a different number of voxels.")
    }
    if (flatten) {
      img <- img[mask > 0.5]
    } else {
      mask <- mincArray(mask)
      img[mask < 0.5] <- 0
    }
  }
  return(img)
}



import_images <- function(imgfiles, mask = NULL, output_format = "list", 
                          flatten = TRUE, margin = 1, inparallel = FALSE, nproc = NULL) {
  
  #Check output format
  format_opts <- c("list", "matrix", "tibble")
  if (!(output_format %in% format_opts)) {
    format_err <- str_c("Argument output_format must be one of (", 
                        str_flatten(format_opts, collapse = ", "), 
                        "): ", output_format)
    stop(format_err)
  }
  
  #Warning when flatten FALSE
  if (!flatten) {
    if (output_format != "list") {
      message(paste("flatten = FALSE is only valid when output_format = 'list'.",
                    "Proceeding with flattened images."))
      flatten <- TRUE
    }
  }
  
  #Import images
  pb <- txtProgressBar(max = length(imgfiles), style = 3)
  progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
  if (inparallel) {
    if (is.null(nproc)) {stop("Specific nproc when running in parallel.")}
    cl <- makeSOCKcluster(nproc)
    registerDoSNOW(cl)
    opts <- list(progress=progress)
    imgs <- foreach(i = 1:length(imgfiles),
                    .packages = "RMINC",
                    .export = c("import_image"),
                    .options.snow = opts) %dopar% {
                      import_image(img = imgfiles[i],
                                   mask = mask,
                                   flatten = flatten)
                    }
    close(pb)
    stopCluster(cl)
  } else {
    imgs <- foreach(i = 1:length(imgfiles),
                    .packages = "RMINC") %do% {
                      progress(n = i)
                      import_image(img = imgfiles[i],
                                   mask = mask,
                                   flatten = flatten)
                    }
  }
  
  #Check image sizes
  imgsize <- unique(map_dbl(imgs, length))
  imgsize_test <- length(imgsize)
  if (imgsize_test != 1) {
    stop("Images provided contain different numbers of voxels.")
  }
  
  #Convert to output format
  if (output_format != "list") {
    if (margin == 1) {
      nrow = length(imgs)
      ncol = imgsize
    } else if (margin == 2) {
      nrow = imgsize
      ncol = length(imgs)
    } else {
      stop()
    }
    out <- matrix(data = 0, nrow = nrow, ncol = ncol)
    for (i in 1:length(imgs)) {
      if (margin == 1) {
        out[i,] <- imgs[[i]]
      } else if (margin == 2) {
        out[,i] <- imgs[[i]]
      }
    }
    if (output_format == "tibble") {
      colnames(out) <- as.character(1:ncol(out))
      out <- as_tibble(out)
    }
  } else {
    out <- imgs
  }
  
  return(out)
  
}


build_voxel_matrix <- function(imgfiles, mask = NULL, file_col = FALSE, 
                               sort = FALSE, save = FALSE, 
                               outfile = "voxel_matrix.csv", inparallel = FALSE,
                               nproc = NULL) {
  
  #Import images as tibble
  df_imgs <- import_images(imgfiles = imgfiles, 
                           mask = mask, 
                           output_format = "tibble",
                           flatten = TRUE,
                           inparallel = inparallel,
                           nproc = nproc)
  
  #Save input files in a column
  df_imgs[["file"]] = imgfiles
  
  #Sort if desired
  if (sort) {df_imgs <- arrange(df_imgs, file)}
  
  #Remove input file column if desired
  if (!file_col) {df_imgs <- select(df_imgs, -file)}
  
  #Save data frame to file if desired
  if (save) {data.table::fwrite(x = df_imgs, file = outfile)}
  
  return(df_imgs)
  
}


vector_to_image <- function(x, outfile, mask) {
  
  #Check output file
  if (is.null(outfile)) {
    stop("Specify output file.")
  }
  
  #Check mask
  if (is.null(mask)) {
    stop("Specify mask file.")
  }
  
  #Import mask
  mask <- mincGetVolume(mask)
  
  #Check that x and mask match
  if (length(x) != sum(mask > 0.5)) {
    stop(paste("Number of elements in x does not match the number of non-zero",
               "voxels in the mask."))
  }
  
  # Export vector as image
  img <- numeric(length(mask))
  img[mask == 1] <- x
  attributes(x) <- attributes(mask)
  sink(nullfile(), type = "output")
  mincWriteVolume(buffer = img,
                  output.filename = outfile,
                  clobber = TRUE,
                  like.filename = attr(mask, "likeVolume"))
  sink(NULL)
  
}


matrix_to_images <- function(x, outfiles, mask, margin = 1, nproc = NULL) {
  
  # Check output files
  if (is.null(outfiles)) {
    stop("Specify output files.")
  }
  
  # Check mask
  if (is.null(mask)) {
    stop("Specify mask file.")
  }
  
  # Check that x is a matrix
  if (!is.matrix(x)) {
    stop("x must be a matrix.")  
  }
  
  # Split matrix into array along margin
  x <- asplit(x, MARGIN = margin)

  # Check that number of output files matches the number of images
  if (length(x) != length(outfiles)) {
    stop("Number of entries in x along margin ", margin, 
         " must be equal to the number of entries in outfiles")  
  }
  
  # Export to images
  if (!is.null(nproc)) {
    # out <- mcmapply(vector_to_image, x, outfiles, 
    #                 MoreArgs = list(mask = mask),
    #                 SIMPLIFY = TRUE, mc.cores = nproc)
    cl <- makeCluster(nproc)
    snk <- clusterEvalQ(cl, library(RMINC))
    out <- clusterMap(cl = cl, 
                      fun = vector_to_image, x, outfiles, 
                      MoreArgs = list(mask = mask))
    stopCluster(cl)
  } else {
    out <- mapply(vector_to_image, x, outfiles,
                  MoreArgs = list(mask = mask),
                  SIMPLIFY = TRUE)
  }
  
}



threshold_intensity <- function(img, threshold = 0.5, symmetric = TRUE, comparison = "gt") {
  
  out <- img
  if (symmetric) {
    if (comparison == "gt") {
      ind <- abs(img) > abs(threshold)
    } else if (comparison == "lt") {
      ind <- abs(img) < abs(threshold)
    } else {
      stop(paste("Argument comparison must be one of ['gt', 'lt']:", comparison))
    }
  } else {
    if (comparison == "gt") {
      ind <- img > threshold
    } else if (comparison == "lt") {
      ind <- img < threshold
    } else {
      stop(paste("Argument comparison must be one of ['gt', 'lt']:", comparison))
    }
  }
  
  if (sum(ind) == 0) {
    stop("No voxels survive the threshold. Select another threshold.")
  }
  
  out[!ind] = 0
  
  return(out)
  
}

threshold_top_n <- function(img, n = 0.2, symmetric = TRUE, tolerance = 1e-5) {
  
  # Raise error if symmetric is True and n < 0
  if (symmetric & (n < 0)) {
    stop(paste("Setting n < 0 while symmetric = True",
               "will return an empty mask."))
  }
  
  # Flatten image
  values <- as.numeric(img)
  
  # If symmetric, use absolute values
  if (symmetric) {
    values <- abs(values)
  }
  
  # Sort values and corresponding indices
  sorted_index <- order(values)
  sorted_values <- values[sorted_index]
  
  # Tolerance filter
  tolerance_filter <- abs(sorted_values) > tolerance
  
  # Compute top n values
  if (n > 0) {
    positive_filter <- sorted_values > 0
    sorted_values <- sorted_values[positive_filter & tolerance_filter]
    sorted_index <- sorted_index[positive_filter & tolerance_filter]
    if (n < 1) {n <- as.integer(floor(n*length(sorted_values)))}
    top_n_index <- sorted_index[length(sorted_index):(length(sorted_index)-n+1)]
  } else if (n < 0) {
    negative_filter <- sorted_values < 0
    sorted_values <- sorted_values[negative_filter & tolerance_filter]
    sorted_index <- sorted_index[negative_filter & tolerance_filter]
    n <- abs(n)
    if (n < 1) {n <- as.integer(floor(n*length(sorted_values)))}
    top_n_index <- sorted_index[1:n]
  } else {
    stop("Argument n cannot be 0.")
  }
  
  # Threshold the image
  out <- as.numeric(img)
  index <- 1:length(out)
  out[!(index %in% top_n_index)] <- 0
  attributes(out) <- attributes(img)
  out <- mincArray(out)
  
  return(out)
  
}


threshold_image <- function(img, method = "top_n", threshold = 0.2, symmetric = TRUE, comparison = "gt") {
  
  if (method == "intensity") {
    img <- threshold_intensity(img = img, 
                               threshold = threshold,
                               symmetric = symmetric,
                               comparison = comparison)
  } else if (method == "top_n") {
    img <- threshold_top_n(img = img,
                           n = threshold,
                           symmetric = symmetric)
  } else {
    stop(paste("Argument method must be one of ",
               "['intensity', 'top_n']"))
  }
  
  return(img)
  
}


mask_from_image <- function(img, signed = FALSE) {
  
  mask <- img
  if (signed) {
    mask[img > 0] = 1
    mask[img < 0] = -1
  } else {
    mask[abs(img) > 0] = 1
  }
  
  return(mask)
  
}




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
#' @param nproc (numeric scalar) Number of processors to use.
#'
#' @return (character vector) Paths to the effect size images.
normative_growth_norm <- function(imgdir, demographics, mask, outdir,
                                  key = "file", group = "patients",
                                  df = 3, batch = NULL,
                                  execution = "local", nproc = 1,
                                  njobs = NULL, resources = list()) {

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
    voxels <- qMincApply(filenames = imgfiles,
                         fun = compute_normative_zscore,
                         demographics = demographics,
                         group = group,
                         batch = batch,
                         df = df,
                         mask = mask,
                         batches = njobs,
                         source = file.path(SRCPATH, "processing.R"),
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



#' Compute cluster centroids
#'
#' @param i 
#' @param clusters 
#' @param mask 
#' @param outdir 
#' @param method 
#' @param execution 
#' @param nproc 
#' @param njobs 
#' @param resources 
#'
#' @return
#' @export
#'
#' @examples
compute_cluster_centroids <- function(i, clusters, mask, outdir, method = "mean",
                                      execution = "local", nproc = 1, 
                                      njobs = NULL, resources = list()){
  
  labels <- clusters[,i]
  files <- rownames(clusters)
  
  # Iterate over clusters
  krange <- sort(unique(labels))  
  centroids <- character(length(krange))
  for (k in krange) {
    
    message(paste("Cluster", k, "of", max(krange)))
    
    # Centroid function
    if (method == "mean") {
      centroid_fun <- mean
    } else if (method == "median") {
      centroid_fun <- median
    } else {
      stop("method must be one of {mean, median}.")
    }
    
    # Images for cluster k
    files_k <- files[labels == k]
    
    # Create centroid image
    centroid <- mcMincApply(filenames = files_k,
                            fun = centroid_fun,
                            mask = mask,
                            cores = nproc)
    
    # Export image
    outfile <- paste0("centroid_nk_", max(krange), "_k_", k, ".mnc")
    outfile <- file.path(outdir, outfile)
    mincWriteVolume(centroid,
                    output.filename = outfile,
                    clobber = TRUE)
    
    centroids[[k]] <- outfile
    
  }
  
  return(centroids)
  
}
