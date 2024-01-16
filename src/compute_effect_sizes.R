#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# compute_effect_sizes.R
# Author: Antoine Beauchamp
# Created: December 23rd, 2023
#
# Calculate voxel-wise effect sizes for human patients

# Description
# -----------


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--imgdir",
              type = "character",
              help = paste("Path to directory containing the images (.mnc) to",
                           "use to compute the effect sizes.")),
  make_option("--demographics",
              type = "character",
              help = "Path to file (.csv) containing the demographics data."),
  make_option("--mask",
              type = "character",
              help = "Path to the mask file (.mnc)."),
  make_option("--outdir",
              type = "character",
              help = paste("Path to directory in which to save the effect",
                           "size images.")),
  make_option("--method",
              type = "character",
              default = "normative-growth",
              help = paste("Method used to compute effect sizes.",
                           "[default %default]")),
  make_option("--key",
              type = "character",
              default = "file",
              help = paste("Primary key between demographics data and",
                           "constructed voxel matrix. [default %default]")),
  make_option("--group",
              type = "character",
              default = "patients",
              help = paste("Group of participants for which to compute",
                           "effect sizes. [default %default]")),
  make_option("--nbatches",
              type = "numeric",
              default = 1),
  make_option("--df",
              type = "numeric",
              default = 3, 
              help = paste("Degrees of freedom to use in normative model",
                           "natural splines. [default %default]")),
  make_option("--batch",
              type = "character",
              help = paste("Variables to use in normalization prior to",
                           "modelling.")),
  make_option("--ncontrols",
              type = "numeric",
              help = paste("Number of propensity-matched controls to use when",
                           "computing the effect sizes. [default %default]")),
  make_option("--matrix-file",
              type = "character",
              help = paste("File in which to export effect size matrix.",
                           "Ignored if NULL.")),
  make_option("--matrix-res",
              type = "numeric",
              help = paste("Resolution of the effect size matrix.",
                           "If specified, effect size images will be",
                           "resampled to this resolution before being",
                           "converted to a matrix.",
                           "Ignored if --matrix-file is NULL.")),
  make_option("--nproc",
              type = "numeric",
              default = 1,
              help = paste("Number of processors to use in parallel.",
                           "Executed serially if 1.",
                           "[default %default]")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = "Verbosity [default %default]")
)


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "pipelines/processing.R"))

# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

# REMOVE THESE LINES WHEN FINISHED
# args[["imgdir"]] <- "data/test/human/derivatives/v2/310/jacobians/absolute/"
# args[["demographics"]] <- "data/human/registration/v2/subject_info/demographics.csv"
# args[["mask"]] <- "data/human/registration/v2/reference_files/mask_3.0mm.mnc"
# args[["outdir"]] <- "data/test/human/derivatives/v2/310/effect_sizes/resolution_3.0/absolute/"
# args[["nbatches"]] <- 1
# args[["nproc"]] <- 8

imgdir <- args[["imgdir"]]
demographics <- args[["demographics"]]
mask <- args[["mask"]]
outdir <- args[["outdir"]]
method <- args[["method"]]
key <- args[["key"]]
group <- args[["group"]]
nbatches <- args[["nbatches"]]
df <- args[["df"]]
batch <- args[["batch"]]
matrix_file <- args[["matrix-file"]]
matrix_res <- args[["matrix-res"]]
nproc <- args[["nproc"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# Check required arguments
args_req <- c("imgdir", "demographics", "mask", "outdir")
for (arg in args_req) {
  if (is.null(args[[arg]])) {
    arg <- paste0("--", arg)
    stop("Argument ", arg, " must be specified.")
  }
}

# Generate effect size images
if (method == "normative-growth") {
  if (nbatches > 1) {
    
    # Create voxel batches
    mask_array <- import_image(img = mask, mask = mask)
    batches <- parallel::splitIndices(length(mask_array), ncl = nbatches)
    
    # Execute script in batches
    batch_dirs <- character(nbatches)
    for (i in 1:nbatches) {
      batch_dir <- file.path(outdir, paste("batch", i, sep = "_"))
      batch_mask <- numeric(length(mask_array))
      batch_mask[batches[[i]]] <- 1
      
      if (!dir.exists(batch_dir)) {
        dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)
      }
      batch_maskfile <- file.path(batch_dir, "batch_mask.mnc")
      vector_to_image(x = batch_mask, 
                      outfile = batch_maskfile,
                      mask = mask)
      
      files <- normative_growth_norm(imgdir = imgdir, 
                                     demographics = demographics,
                                     mask = batch_maskfile,
                                     outdir = batch_dir,
                                     key = key,
                                     group = group,
                                     df = df,
                                     batch = batch,
                                     nproc = nproc)
      
    }
    
    # Collate batch images
    # print("Collating batched images...")
    # outfiles = os.listdir(batch_dirs[0])
    # outfiles = [file for file in outfiles if file != 'batch_mask.mnc']
    # for outfile in outfiles:
    #   
    #   img = np.zeros_like(mask_array)
    # for b, batch in enumerate(batches):
    #   batch_img = os.path.join(batch_dirs[b], outfile)
    # batch_mask = os.path.join(batch_dirs[b], 'batch_mask.mnc')
    # img[batch] = import_image(img = batch_img, mask = batch_mask)
    # 
    # outfile = os.path.join(outdir, outfile)
    # vector_to_image(x = img, outfile = outfile, maskfile = mask)
    # 
    # print("Cleaning up...")
    # for batch_dir in batch_dirs:
    #   rmtree(batch_dir)
    
    
    quit()
  } else {
    print("In batch = 1 condition")
    files <- normative_growth_norm(imgdir = imgdir,
                                   demographics = demographics,
                                   mask = mask,
                                   outdir = outdir,
                                   key = key,
                                   group = group,
                                   df = df,
                                   batch = batch,
                                   nproc = nproc)
  }
} else if (method == "propensity-matching") {
  files <- propensity_matching_norm(imgdir = imgdir, 
                                    demographics = demographics,
                                    mask = mask,
                                    outdir = outdir,
                                    ncontrols = ncontrols,
                                    nproc = nproc)
} else {
  stop("Argument --method must be one of ",
       "{normative-growth, propensity-matching")
}

# Generate effect size matrix
if (!is.null(matrix_file)) {
  if (!is.null(matrix_res)) {
    # Resample images if alternate resolution is specified
    files <- resample_images(files = files,
                             outdir = ...,
                             isostep = matrix_res,
                             nproc = nproc)
  } 
  
  # Build effect size matrix
  build_voxel_matrix(imgfiles = files,
                     mask = ...,
                     file_col = TRUE,
                     sort = TRUE,
                     outfile = matrix_file,
                     nproc = nproc)
}

