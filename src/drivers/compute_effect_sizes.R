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
  make_option("--matrix-resolution",
              type = "numeric",
              help = paste("Resolution of the effect size matrix.",
                           "If specified, effect size images will be",
                           "resampled to this resolution before being",
                           "converted to a matrix.",
                           "Ignored if --matrix-file is NULL.")),
  make_option("--execution",
              type = "character",
              default = "local",
              help = paste("Flag indicating whether to run locally or",
                           "using Slurm. [default %default]")),
  make_option("--nproc",
              type = "numeric",
              default = 1,
              help = paste("Number of processors to use.",
                           "Executed in parallel if > 1.",
                           "[default %default]")),
  make_option("--registry-name",
              type = "character",
              help = "Name of the registry directory for batched jobs."),
  make_option("--registry-cleanup",
              type = "character",
              default = "true",
              help = paste("Option to clean up registry after completion",
                           "of batched jobs. [default %default]" )),
  make_option("--slurm-njobs",
              type = "numeric",
              help = "Number of jobs to deploy on Slurm."),
  make_option("--slurm-mem",
              type = "character",
              help = paste("Memory per CPU  when --execution slurm.",
                           "Example: '16G'")),
  make_option("--slurm-time",
              type = "numeric",
              help = paste("Walltime in minutes for Slurm jobs",
                           "when --execution slurm.")),
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
execution <- args[["execution"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)


print(demographics)
quit()

# Check required arguments
args_req <- c("imgdir", "demographics", "mask", "outdir")
for (arg in args_req) {
  if (is.null(args[[arg]])) {
    arg <- paste0("--", arg)
    stop("Argument ", arg, " must be specified.")
  }
}

# Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Check execution option
if (execution == "local") {
  resources <- list()
  registry_name <- NULL
  registry_cleanup <- NULL
} else if (execution == "slurm") {
  njobs <- args[["slurm-njobs"]]
  registry_name <- args[["registry-name"]]
  registry_cleanup <- ifelse(args[["registry-cleanup"]] == "true",
                             TRUE, FALSE)
  resources <- list(memory = args[["slurm-mem"]],
                    walltime = args[["slurm-time"]]*60)
} else {
  stop()
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
                                   execution = execution,
                                   nproc = nproc,
                                   registry_name = registry_name,
                                   registry_cleanup = registry_cleanup,
                                   njobs = njobs,
                                   resources = resources)
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
