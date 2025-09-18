# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(SNFtool))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doSNOW))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

# Import modified minc_parallel.R functions from RMINC
# (I forget why right now)
source(file.path(SRCPATH, "minc_parallel.R"))

#' Import a MINC image
#'
#' @param img (character scalar) Path to image file (.mnc)
#' @param mask (character scalar) Path to mask file (.mnc)
#' @param flatten (logical scalar) Return image as flattened vector or array
#' @param version 
#'
#' @returns (mincSingleDim or mincArray) Image
import_image <- function(img, mask = NULL, flatten = TRUE, version = "v1") {
  
  # Import image
  img <- mincGetVolume(img)
  
  # Convert to 3D if specified
  if (!flatten) {img <- mincArray(img)}
  
  # Apply mask if specified
  if (!is.null(mask)) {
    mask <- mincGetVolume(mask)
    if (length(img) != length(mask)) {
      stop("Input image and mask contain a different number of voxels.")
    }
    if (flatten) {
      if (version == "v2") {
        img <- img[mask > 0.5]
      } else {
        img <- img[mask == 1]
      }
    } else {
      mask <- mincArray(mask)
      if (version == "v2") {
        img[mask < 0.5] <- 0
      } else {
        img[mask != 1] <- 0
      }
    }
  }
  return(img)
}

#' Import a set of MINC images
#'
#' @param imgfiles (character vector) Path to image files (.mnc)
#' @param mask (character scalar) Path to mask file (.mnc)
#' @param output_format (character scalar) Output format for imported images
#' @param flatten (logical scalar) Flatten images into vectors
#' @param margin (numeric scalar) Margin along which to store images in matrix
#' @param version 
#' @param nproc (numeric scalar) Number of processors to use. 
#'
#' @returns (list, matrix, tbl) Images
import_images <- function(imgfiles, mask = NULL, output_format = "list",
                          flatten = TRUE, margin = 1, version = "v1",
                          nproc = 1) {
  
  # Check output format
  format_opts <- c("list", "matrix", "tibble")
  if (!(output_format %in% format_opts)) {
    format_err <- str_c("Argument output_format must be one of (",
                        str_flatten(format_opts, collapse = ", "),
                        "): ", output_format)
    stop(format_err)
  }
  
  # Warning when flatten FALSE
  if (!flatten) {
    if (output_format != "list") {
      message(paste("flatten = FALSE is only valid when output_format = 'list'.",
                    "Proceeding with flattened images."))
      flatten <- TRUE
    }
  }
  
  # Import images
  pb <- txtProgressBar(max = length(imgfiles), style = 3)
  progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
  if (nproc > 1) {
    cl <- makeSOCKcluster(nproc)
    registerDoSNOW(cl)
    opts <- list(progress=progress)
    imgs <- foreach(i = 1:length(imgfiles),
                    .packages = "RMINC",
                    .export = c("import_image"),
                    .options.snow = opts) %dopar% {
                      import_image(img = imgfiles[i],
                                   mask = mask,
                                   flatten = flatten,
                                   version = version)
                    }
    close(pb)
    stopCluster(cl)
  } else {
    imgs <- foreach(i = 1:length(imgfiles),
                    .packages = "RMINC") %do% {
                      progress(n = i)
                      import_image(img = imgfiles[i],
                                   mask = mask,
                                   flatten = flatten,
                                   version = version)
                    }
  }
  
  # Check image sizes
  imgsize <- unique(map_dbl(imgs, length))
  imgsize_test <- length(imgsize)
  if (imgsize_test != 1) {
    stop("Images provided contain different numbers of voxels.")
  }
  
  # Convert to output format
  if (output_format != "list") {
    if (margin == 1) {
      nrow <- length(imgs)
      ncol <- imgsize
    } else if (margin == 2) {
      nrow <- imgsize
      ncol <- length(imgs)
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


#' Build a voxel matrix (data frame) from a set of images
#'
#' @param imgfiles (character vector) Path to image files (.mnc)
#' @param mask (character scalar) Path to mask file (.mnc)
#' @param file_col (logical scalar) Keep image files in a column named "file".
#' @param sort (logical scalar) Option to sort rows by image file names.
#' @param save (logical scalar) Option to export file to CSV
#' @param outfile (character scalar) Name of file (.csv) in which to export
#' @param version (character scalar)
#' @param nproc (numeric scalar) Number of processors to use.
#'
#' @returns (tbl, data.frame) Data frame of image voxels
build_voxel_matrix <- function(imgfiles, mask = NULL, file_col = FALSE,
                               sort = FALSE, save = FALSE,
                               outfile = "voxel_matrix.csv", version = "v1",
                               nproc = 1) {
  
  # Import images as tibble
  df_imgs <- import_images(imgfiles = imgfiles,
                           mask = mask,
                           output_format = "tibble",
                           flatten = TRUE,
                           version = version,
                           nproc = nproc)
  
  # Save input files in a column
  df_imgs[["file"]] <- imgfiles
  
  # Sort if desired
  if (sort) {df_imgs <- arrange(df_imgs, file)}
  
  # Remove input file column if desired
  if (!file_col) {df_imgs <- select(df_imgs, -file)}
  
  # Save data frame to file if desired
  if (save) {data.table::fwrite(x = df_imgs, file = outfile)}
  
  return(df_imgs)
  
}


#' Export a vector as MINC image
#'
#' @param x (numeric vector) Vector to export
#' @param outfile (character scalar) Path to output image file (.mnc)
#' @param mask (character scalar) Path to mask file (.mnc). This is used to 
#' define the image.
#' @param version 
#'
#' @returns NULL
vector_to_image <- function(x, outfile, mask, version = "v1") {
  
  # Check output file
  if (is.null(outfile)) {stop("Specify output file.")}
  
  # Check mask
  if (is.null(mask)) {stop("Specify mask file.")}
  
  # Import mask
  mask <- mincGetVolume(mask)
  
  if (version == "v2") {
    ind_mask <- mask > 0.5
  } else {
    ind_mask <- mask == 1
  }
  
  # Check that x and mask match
  if (length(x) != sum(ind_mask)) {
    stop(paste("Number of elements in x does not match the number of non-zero",
               "voxels in the mask."))
  }
  
  # Export vector as image
  img <- numeric(length(mask))
  img[ind_mask] <- x
  attributes(x) <- attributes(mask)
  sink(nullfile(), type = "output")
  mincWriteVolume(buffer = img,
                  output.filename = outfile,
                  clobber = TRUE,
                  like.filename = attr(mask, "likeVolume"))
  sink(NULL)
  
}


#' Export a matrix to a set of images
#'
#' @param x (matrix) Matrix to export
#' @param outfiles (character vector) Vector of paths to output files (.mnc)
#' @param mask (character scalar) Path to mask file (.mnc). This is used to 
#' define the image.
#' @param margin (numeric scalar) Margin specifying the images. 
#' @param version 
#' @param nproc (numeric scalar) Number of processors to use.
#'
#' @returns
matrix_to_images <- function(x, outfiles, mask, margin = 1, version = "v1",
                             nproc = 1) {
  
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
  if (nproc > 1) {
    cl <- makeCluster(nproc)
    snk <- clusterEvalQ(cl, library(RMINC))
    out <- clusterMap(cl = cl,
                      fun = vector_to_image, x, outfiles,
                      MoreArgs = list(mask = mask,
                                      version = version))
    stopCluster(cl)
  } else {
    out <- mapply(vector_to_image, x, outfiles,
                  MoreArgs = list(mask = mask,
                                  version = version),
                  SIMPLIFY = TRUE)
  }
  
}


#' Threshold an image based on intensity values
#'
#' @param img (mincSingleDim or mincArray) Image to threshold.
#' @param threshold (numeric scalar) Intensity threshold to apply.
#' @param symmetric (logical scalar) Apply threshold symmetrically. 
#' @param comparison (character scalar) Method of comparison used to apply the 
#' threshold.
#'
#' @returns (mincSingleDim or mincArray) Thresholded image.
threshold_intensity <- function(img, threshold = 0.5, symmetric = TRUE,
                                comparison = "gt") {
  
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
  
  out[!ind] <- 0
  
  return(out)
  
}

#' Threshold an image based on rank
#'
#' @param img (mincSingleDim or mincArray) Image to threshold.
#' @param n (numeric scalar) Top number or fraction of voxels to return. 
#' If negative, the function will return the top negative values. 
#' @param symmetric (logical scalar) Apply threshold symmetrically. 
#' @param tolerance (numeric scalar)
#'
#' @returns (mincSingleDim or mincArray) Thresholded image.
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


#' Threshold an image
#'
#' @param img (mincSingleDim or mincArray) Image to threshold.
#' @param method (character scalar) Method used for thresholding.
#' @param threshold (numeric scalar) Threshold value.
#' @param symmetric (logical scalar) Apply threshold symmetrically. 
#' @param comparison (character scalar) Method of comparison used to apply the 
#' threshold.
#'
#' @returns (mincSingleDim or mincArray) Thresholded image.
threshold_image <- function(img, method = "top_n", threshold = 0.2,
                            symmetric = TRUE, comparison = "gt") {
  
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


#' Create a mask from nonzero voxels
#'
#' @param img (mincSingleDim or mincArray) Image to threshold.
#' @param signed (logical scalar) Return a mask with -1 for negative voxels 
#' and 1 for positive voxels
#'
#' @returns (mincSingleDim or mincArray) Mask
mask_from_image <- function(img, signed = FALSE) {
  
  mask <- img
  if (signed) {
    mask[img > 0] <- 1
    mask[img < 0] <- -1
  } else {
    mask[abs(img) > 0] <- 1
  }
  return(mask)
}


#' Compute normative z-score for a voxel
#'
#' @param y (numeric vector) Voxel values across study participants.
#' @param idx (list) Train and test indices
#' @param demographics (data.frame) Demographics information for study
#' participants.
#' @param df (numeric scalar) Degrees of freedom in natural spline
#' model.
#'
#' @return (numeric vector) Voxel normative z-scores
compute_normative_zscore <- function(y, idx, demographics, df = 3) {
  
  if ("batch" %in% colnames(demographics)) {
    y <- residuals(lm(y ~ batch, data = demographics))
    names(y) <- NULL
  }
  
  # Training data frame
  df_fit <- demographics[idx[[1]], c("Age", "Sex")]
  df_fit[["y"]] <- y[idx[[1]]]
  
  # Test data frame
  df_pred <- demographics[idx[[2]], c("Age", "Sex")]
  df_pred[["y"]] <- y[idx[[2]]]
  
  # Fit model and predict on test set
  model_fit <- lm(y ~ Sex + ns(Age, df = df), data = df_fit)
  model_pred <- predict(model_fit,
                        newdata = df_pred,
                        interval = "prediction",
                        level = pnorm(q = 1) - pnorm(q = -1))
  
  # Compute z-score
  y_pred <- model_pred[,"fit"]
  y_sd <- y_pred - model_pred[,"lwr"]
  z <- (df_pred[["y"]] - y_pred)/y_sd
  
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
#' @param execution (character scalar) Flag indicating whether to run
#' locally or using Slurm.
#' @param nproc (numeric scalar) Number of processors to use. Executed
#' in parallel if > 1
#' @param registry_name (character scalar) Name of the registry directory
#' for batched jobs.
#' @param registry_cleanup (logical scalar) Option to clean up registry
#' after completion of batched jobs.
#' @param njobs (numeric scalar) Number of jobs to deploy on Slurm.
#' @param resources (list) List of resources for Slurm jobs.
#'
#' @return (character vector) Paths to the effect size images.
normative_growth_norm <- function(imgdir, demographics, mask, outdir,
                                  cv_seed = cv_seed,
                                  key = "file", group = "patients",
                                  df = 3, batch = NULL,
                                  execution = "local", nproc = 1,
                                  registry_name = NULL,
                                  registry_cleanup = TRUE,
                                  njobs = NULL, resources = list(), 
                                  verbose = TRUE) {

  # Raise cross-validation error if group not controls
  if (!is.null(cv_seed) & group != "controls") {
    stop("Cross-validation only possible for controls.")
  }
  
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
           !is.na(Sex))
  
  # TODO: Remove rows where batch columns are NA
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
  
  
  ind_ctrl <- demographics[["DX"]] == "Control"
  ind_dx <- !ind_ctrl
  if (group == "controls") {
    demographics <- demographics[ind_ctrl,]
    imgfiles <- imgfiles[ind_ctrl]
    
    # If cross-validation, sample train-test split
    if (!is.null(cv_seed)) {
      ind_fit <- logical(nrow(demographics))
      set.seed(cv_seed)
      ind_sample <- sample(1:nrow(demographics), size = nrow(demographics), replace = FALSE)
      ind_fit[ind_sample[1:(floor(0.8*nrow(demographics)))]] <- TRUE
      ind_pred <- !ind_fit
    } else {
      ind_fit <- !logical(nrow(demographics))
      ind_pred <- ind_fit
    }
  } else if (group == "patients") {
    ind_fit <- ind_ctrl
    ind_pred <- ind_dx
  } else if (group == "all") {
    ind_fit <- ind_ctrl
    ind_pred <- ind_ctrl | ind_dx
  }
  
  # 
  idx <- list(ind_fit, ind_pred)
  
  if (!is.null(batch)) {
    demographics <- demographics %>% 
      unite(col = "batch", all_of(batch)) %>% 
      select(all_of(key), Age, Sex, batch)
  } else {
    demographics <- demographics %>% 
      select(all_of(key), Age, Sex)
  }
  
  # Run normative growth modelling
  if (verbose) {message("Evaluating normative growth models...")}
  if (execution == "local") {
    voxels <- mcMincApply(filenames = imgfiles,
                          fun = compute_normative_zscore,
                          idx = idx,
                          demographics = demographics,
                          df = df,
                          mask = mask,
                          cores = nproc,
                          return_raw = TRUE)
  } else if (execution == "slurm") {
    voxels <- qMincApply(filenames = imgfiles,
                         fun = compute_normative_zscore,
                         idx = idx,
                         demographics = demographics,
                         df = df,
                         mask = mask,
                         batches = njobs,
                         source = file.path(SRCPATH, "processing.R"),
                         registry_name = ifelse(is.null(registry_name),
                                                "registry_normative_growth",
                                                registry_name),
                         cleanup = registry_cleanup,
                         return_raw = TRUE,
                         resources = resources)
  } else {
    stop()
  }
  
  # Convert voxel list into matrix
  # This matrix has number of voxels consistent with mask > 0.5
  voxels <- simplify_masked(voxels[["vals"]])
  
  # Export images
  if (verbose) {message("Exporting normalized images...")}
  if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
  outfiles <- demographics[idx[[2]], key][[1]]
  outfiles <- file.path(outdir, outfiles)
  matrix_to_images(x = voxels, outfiles = outfiles, mask = mask,
                   margin = 2, version = "v2", nproc = nproc)
  
  
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
  #TODO: Fill out this function
  return(outfiles)
}


#' Run similarity network fusion (SNF)
#'
#' @param x (list) List of input matrices
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
similarity_network <- function(x, metric = "correlation", K = 10,
                               sigma = 0.5, t = 20, outfile = NULL){

  if (metric == "correlation") {
    d <- map(x, function(x) {1-cor(t(x))})
  } else if (metric == "euclidean") {
    d <- map(x, function(x){(dist2(as.matrix(x), as.matrix(x)))^(1/2)})
  } else {
    stop(paste("Argument metric must be one of {correlation, euclidean}:", metric))
  }

  W_list <- map(d, affinityMatrix, K = K, sigma = sigma)

  W <- SNF(W_list, K = K, t = t)

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
#' @param i (numeric scalar) Column index.
#' @param clusters (data.frame, tbl) Data frame containing cluster assignments.
#' @param mask (character scalar) Path to a mask image (.mnc).
#' @param outdir (character scalar) Path to the output directory.
#' @param method (character scalar) Method used to compute centroids.
#' @param execution (character scalar) Method of execution.
#' @param nproc (numeric scalar) Number of processors to use.
#' @param njobs (numeric scalar) Number of jobs to deploy.
#' @param resources (list) Resources to use on a Slurm cluster.
#'
#' @return (character vector) Paths to the centroid images.
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


scaler <- function(data, scale = TRUE, axis = "columns"){
  
  if (axis == "columns"){
    
    colMeansMat <- matrix(colMeans(data), nrow = nrow(data), ncol = ncol(data), byrow = T)
    out <- data - colMeansMat
    
    if (scale == TRUE){
      
      sigma <- sqrt(colSums(out^2)/(nrow(data)-1))
      sigmaMat <- matrix(sigma, nrow = nrow(data), ncol = ncol(data), byrow = T)
      out <- out/sigmaMat
      
    }
    
  } else if (axis == "rows") {
    
    rowMeansMat <- matrix(rowMeans(data), nrow = nrow(data), ncol = ncol(data), byrow = F)
    out <- data - rowMeansMat
    
    if (scale == TRUE){
      
      sigma <- sqrt(rowSums(out^2)/(ncol(data)-1))
      sigmaMat <- matrix(sigma, nrow = nrow(data), ncol = ncol(data), byrow = F)
      out <- out/sigmaMat
      
    }
  }
  
  return(out)
  
}


MakeModelList<-function(scanbase_scans_file="Data/Resources/scanbase_40um - Scans_22July19.csv",
                        scanbase_studies_file="Data/Resources/scanbase_40um - Studies_Feb2023.csv",
                        scanbase_genotypes_file="Data/Resources/scanbase_40um - Genotypes_Feb2023.csv",
                        scanbase_focus="Autism",
                        scanbase_sample_size_threshold=6,
                        base_directory="Data/Outputs/ModelInfo"){

  dir.create(base_directory, recursive = TRUE, showWarnings = FALSE)

  scans <- read.csv(scanbase_scans_file)
  studies <- read.csv(scanbase_studies_file)
  genotypes_file <- read.csv(scanbase_genotypes_file)

  # Use only autism studies
  studies %<>% filter(Focus==scanbase_focus)

  # Keep only scans in autism studies
  scans %<>% filter(Study_Name %in% studies$Study_Name)
  genotypes_file %<>% filter(Study_Name %in% studies$Study_Name)

  scans %<>% filter(Genotype_Code %in% genotypes_file$Genotype_Code)
  scans %<>% filter(Study_Name %in% genotypes_file$Study_Name)


  # Make genotype code a character
  scans %<>% mutate(Study_Name=as.character(Study_Name),
                    Genotype_Code=as.character(Genotype_Code))


  # Remove sex tags in scan genotype code
  # Why is this even here? There is a variable associated with sex
  scans %<>% mutate(Genotype_Code=(Genotype_Code %>%
                                     gsub("_M", "", .) %>%
                                     gsub("_F", "", .)))

  # Fix Pten_Kim study
  scans$Study_Name[which(scans$Study_Name=="PTEN_Kim" & endsWith(scans$Genotype_Code, "_131"))] <- "PTEN_Kim_131"
  scans$Study_Name[which(scans$Study_Name=="PTEN_Kim" & endsWith(scans$Genotype_Code, "_167"))] <- "PTEN_Kim_167"

  # Relabel genotypes for consistency
  wt_codes <- c("B6", "Wt", "WT_131", "WT_167", "Chd7++", "Chd7ff", "CON","Kctd13Wt_LatWt","WND","WNA","Snf2L_WT","Snf2Hcko-nestincre_WT")
  scans$Genotype_Code[which(scans$Genotype_Code %in% wt_codes)] <- "WT"

  het_codes <- c("HT", "HETero", "HT_131", "HT_167", "Het", "Hetero")
  scans$Genotype_Code[which(scans$Genotype_Code %in% het_codes)] <- "HET"

  hom_codes <- c("Hom", "HOMo","Homo")
  scans$Genotype_Code[which(scans$Genotype_Code %in% hom_codes)] <- "HOM"

  hem_codes <- c("Hemi","Hem","HMZ","Hmz")
  scans$Genotype_Code[which(scans$Genotype_Code %in% hem_codes)] <- "HEM"

  mut_codes <- c("Mut")
  scans$Genotype_Code[which(scans$Genotype_Code %in% mut_codes)] <- "MUT"

  # Sample size data per genotype
  table(scans$Study_Name, scans$Genotype_Code)

  # Remove studies with fewer than 8 wildtypes
  scans %<>% filter(Study_Name %in% names(which(table(scans$Study_Name, scans$Genotype_Code)[,"WT"] >= scanbase_sample_size_threshold)))

  #Removes scans with bad label files based on NA values in later steps.
  #scans<-scans[-c(1674,2402,2417,3097,3357,3387,3675),]

  write.csv(scans, file = glue("{base_directory}/scanbase_scans_autism_filtered_Feb2023.csv"))

  unique_studies <- unique(scans$Study_Name)
  model_list <- list()
  for (i in seq_along(unique_studies)) {
    study <- unique_studies[i]
    study_scans <- scans %>% filter(Study_Name==study)
    unique_genotypes <- setdiff(unique(study_scans$Genotype_Code), "WT")
    genotypes <- c()
    for (j in seq_along(unique_genotypes)) {
      genotype <- unique_genotypes[j]
      N <- nrow(study_scans %>% filter(Genotype_Code==genotype))
      if (N >= scanbase_sample_size_threshold) {
        genotypes <- c(genotypes, genotype)
      }
    }
    if (length(genotypes) >= 1) {
      model_list[[study]] <- list(wildtype="WT",
                                  genotypes=genotypes)
    }
  }

  save(list = "model_list", file = glue("{base_directory}/model_list_Feb2023.RData"))
  return(model_list)
}



MakeEffectSizeMaps<-function(jdtype,model_list,resolution="200",
                             dir_determinants="Data/Raw/jacobian_determinants",
                             output_phenotype_dir="Data/Outputs/EffectSizeMaps",
                             base_directory="Data/Outputs/ModelInfo",
                             boot="Y",
                             num_boot="100"){

  scans<-read.csv(glue("{base_directory}/scanbase_scans_autism_filtered_Feb2023.csv"))

  if(boot=="Y"){
    output_phenotype_dir<-glue("{output_phenotype_dir}/{resolution}/Boot")
  } else {
    output_phenotype_dir<-glue("{output_phenotype_dir}/{resolution}")
  }

  dir.create(output_phenotype_dir,recursive = TRUE,showWarnings = FALSE)
  mvol <- mincGetVolume(glue("{dir_determinants}/scanbase_second_level-nlin-3_mask_{resolution}um.mnc"))
  nvoxels <- length(which(mvol > 0.5))
  scans<-read.csv(file = glue("{base_directory}/scanbase_scans_autism_filtered_Feb2023.csv"))

  if(startsWith(tolower(jdtype),"a")){
    jddir <- glue("abs_{resolution}")
    jdtype <- "Absolute"
  }else{
    jddir <- glue("rel_{resolution}")
    jdtype <- "Relative"
  }

  effect_size_data_matrix <- matrix(nrow=length(as.character(unlist(sapply(model_list, "[[", 2)))), ncol=nvoxels)
  effect_size_data_matrix_rownames <- c()

  if (boot=="Y"){
    for (nboot in 1:num_boot) {
      for (model in names(model_list)) {
        genotypes <- model_list[[model]]$genotypes

        study_scans_wt <- scans %>% filter(Study_Name==model, Genotype_Code %in% c("WT"))
        wt_data_matrix <- matrix(nrow=nrow(study_scans_wt), ncol = nvoxels)

        for (i in 1:nrow(study_scans_wt)) {
          filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_wt$Mouse_ID)[i]}_{jddir}.mnc")
          wt_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
        }

        rowboot<-sample(nrow(wt_data_matrix),size=nrow(wt_data_matrix),replace=TRUE)
        wt_data_matrix<-wt_data_matrix[rowboot,]

        wt_mean <- colMeans(wt_data_matrix)
        # wt_sd <- colSd(wt_data_matrix)
        # Equivalent to: apply(wt_data_matrix, MARGIN=2, FUN=function(x) {mean(x)})
        wt_sd<-apply(wt_data_matrix, MARGIN=2, sd)

        for (genotype in genotypes) {

          print(paste("Working on model:", model, ",", "genotype:", genotype, ",", "num_boot:", nboot))

          # Make a dataframe with only wildtypes and genotype of interest
          study_scans_mut <- scans %>% filter(Study_Name==model, Genotype_Code %in% c(genotype))
          mut_data_matrix <- matrix(nrow=nrow(study_scans_mut), ncol = nvoxels)

          rowboot<-sample(nrow(mut_data_matrix),size=nrow(mut_data_matrix),replace=TRUE)
          mut_data_matrix<-mut_data_matrix[rowboot,]

          for (i in 1:nrow(study_scans_mut)) {
            filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_mut$Mouse_ID)[i]}_{jddir}.mnc")
            mut_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
          }

          mut_mean <- colMeans(mut_data_matrix)
          effect_size <- (mut_mean - wt_mean) / wt_sd

          filenameES <- glue("{output_phenotype_dir}/{as.character(model)}_{as.character(genotype)}_ES_{jdtype}_{resolution}_{nboot}.mnc")

          ovol <- mvol
          ovol[] <- 0
          ovol[mvol > 0.5] <- effect_size
          #ovol[is.infinite()]<-0
          ovol[is.na(ovol)]<-0
          mincWriteVolume(ovol, filenameES)

        }
      }
    }
  } else {
    k <- 1
    for (model in names(model_list)) {
      genotypes <- model_list[[model]]$genotypes

      study_scans_wt <- scans %>% filter(Study_Name==model, Genotype_Code %in% c("WT"))
      wt_data_matrix <- matrix(nrow=nrow(study_scans_wt), ncol = nvoxels)

      for (i in 1:nrow(study_scans_wt)) {
        filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_wt$Mouse_ID)[i]}_{jddir}.mnc")
        wt_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
      }

      wt_mean <- colMeans(wt_data_matrix)
      # wt_sd <- colSd(wt_data_matrix)
      # Equivalent to: apply(wt_data_matrix, MARGIN=2, FUN=function(x) {mean(x)})
      wt_sd<-apply(wt_data_matrix, MARGIN=2, sd)

      for (genotype in genotypes) {

        message(paste("Working on model:", model, ",", "genotype:", genotype))

        # Make a dataframe with only wildtypes and genotype of interest
        study_scans_mut <- scans %>% filter(Study_Name==model, Genotype_Code %in% c(genotype))
        mut_data_matrix <- matrix(nrow=nrow(study_scans_mut), ncol = nvoxels)

        for (i in 1:nrow(study_scans_mut)) {
          filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_mut$Mouse_ID)[i]}_{jddir}.mnc")
          mut_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
        }

        mut_mean <- colMeans(mut_data_matrix)
        effect_size <- (mut_mean - wt_mean) / wt_sd

        filenameES <- glue("{output_phenotype_dir}/{as.character(model)}_{as.character(genotype)}_ES_{jdtype}_{resolution}.mnc")

        ovol <- mvol
        ovol[] <- 0
        ovol[mvol > 0.5] <- effect_size
        #ovol[is.infinite()]<-0
        ovol[is.na(ovol)]<-0
        mincWriteVolume(ovol, filenameES, clobber = TRUE)

        effect_size_data_matrix[k,] <- effect_size
        effect_size_data_matrix_rownames <- c(effect_size_data_matrix_rownames, paste(model, genotype, sep="_"))
        k <- k + 1


      }
    }

  }

  if(boot!="Y"){
    rownames(effect_size_data_matrix)<-effect_size_data_matrix_rownames
    save(file=glue("{output_phenotype_dir}/ES_Matricies_{jdtype}.RData"),effect_size_data_matrix)
    return(effect_size_data_matrix)
  }
}