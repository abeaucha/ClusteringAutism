# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doSNOW))


# Functions ------------------------------------------------------------------

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



