suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(tcltk))


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
      img <- img[mask == 1]
    } else {
      mask <- mincArray(mask)
      img[mask == 0] <- 0
    }
  }
  return(img)
}



import_images <- function(imgfiles, mask = NULL, output_format = "list", 
                          flatten = TRUE, inparallel = FALSE, nproc = NULL) {
  
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
  imgsize_test <- length(unique(map_dbl(imgs, length)))
  if (imgsize_test != 1) {
    stop("Images provided contain different numbers of voxels.")
  }
  
  #Convert to output format
  if (output_format != "list") {
    imgs <- reduce(imgs, rbind)
    rownames(imgs) <- NULL
    colnames(imgs) <- NULL
    if (output_format == "tibble") {
      colnames(imgs) <- as.character(1:ncol(imgs))
      imgs <- as_tibble(imgs)
    }
  }
  
  return(imgs)
  
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
  
  if (is.null(outfile)) {
    stop("Specify output file.")
  }
  
  if (is.null(mask)) {
    stop("Specify mask file.")
  }
  
  mask <- mincGetVolume(mask)
  
  if (length(x) != sum(mask)) {
    stop(paste("Number of elements in x does not match the number of non-zero",
               "voxels in the mask."))
  }
  
  img <- numeric(length(mask))
  img[mask == 1] <- x
  attributes(x) <- attributes(mask)
  mincWriteVolume(buffer = img,
                  output.filename = outfile,
                  clobber = TRUE,
                  like.filename = attr(mask, "likeVolume"))
  
}


matrix_to_images <- function(x, outfiles, mask) {
  
  #Check output files
  if (is.null(outfiles)) {
    stop("Specify output files.")
  }
  
  #Check mask
  if (is.null(mask)) {
    stop("Specify mask file.")
  }
  
  #Check that x is a matrix
  if (!is.matrix(x)) {
    stop("x must be a matrix.")  
  }
  
  #Check that number of output files matches the number of images
  if (nrow(x) != length(outfiles)) {
    stop(paste("Number of rows in x must be equal to the number of entries in",
               "outfiles"))  
  }
  
  #Export matrix rows
  sink(file = nullfile(), type = "output")
  for (i in 1:nrow(x)) {
    vector_to_image(x = x[i,], outfile = outfiles[i], mask = mask)
  }
  sink(file = NULL)
  
  return(outfiles)
  
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



