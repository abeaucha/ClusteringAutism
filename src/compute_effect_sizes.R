#!Rscript
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
  make_option("--df",
              type = "numeric",
              default = 3, 
              help = paste("Degrees of freedom to use in normative model",
                           "natural splines. [default %default]")),
  make_option("--batch",
              type = "character",
              help = paste("Variables to use in normalization prior to",
                           "modelling.")),
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

# Processing functions
source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "pipelines/processing.R"))


# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
imgdir <- args[["imgdir"]]
demographics <- args[["demographics"]]
mask <- args[["mask"]]
outdir <- args[["outdir"]]
method <- args[["method"]]
key <- args[["key"]]
group <- args[["group"]]
df <- args[["df"]]
batch <- args[["batch"]]
matrix_file <- args[["matrix-file"]]
matrix_res <- args[["matrix-res"]]
nproc <- args[["nproc"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# Check necessary input args
input_check <- c("imgdir", "demographics", "mask", "outdir")
for (input in input_check) {
  if (is.null(args[[input]])) {
    input_arg <- paste0("--", input)
    stop("Argument ", input_arg, " must be specified.")
  }
}

# Check nproc
if (is.null(nproc)) {
  stop("Specify the number of processors to use in parallel.")
}

# Generate effect size images
if (method == "normative-growth") {
  files <- normative_growth_norm(imgdir = imgdir, 
                                 demographics = demographics,
                                 mask = mask,
                                 outdir = outdir,
                                 key = key,
                                 group = group,
                                 df = df,
                                 batch = batch,
                                 nbatches = nbatches,
                                 nproc = nproc)
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

