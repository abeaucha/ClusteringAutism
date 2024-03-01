
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
              help = paste("[default %default]")),
  make_option("--nproc",
              type = "numeric",
              default = 1,
              help = paste("Number of processors to use in parallel.",
                           "Executed serially if 1.",
                           "[default %default]")),
  make_option("--slurm-njobs",
              type = "numeric",
              help = "Number of jobs to deploy on Slurm."),
  make_option("--slurm-mem",
              type = "character",
              help = paste("Memory per CPU core")),
  make_option("--slurm-time",
              type = "numeric",
              help = paste("Walltime in minutes")),
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
#TODO remove line
# source(file.path(SRCPATH, "pipelines/processing.R"))


# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

#TODO remove lines when script works
# REMOVE THESE LINES WHEN FINISHED
# args[["imgdir"]] <- "data/test/human/derivatives/v2/700/jacobians/absolute/"
# args[["demographics"]] <- "data/human/registration/v2/subject_info/demographics.csv"
# args[["mask"]] <- "data/human/registration/v2/reference_files/mask_0.8mm.mnc"
# args[["outdir"]] <- "data/test/human/derivatives/v2/700/effect_sizes/resolution_0.8/absolute/"
# args[["slurm-njobs"]] <- 300
# args[["slurm-time"]] <- 60
# args[["slurm-mem"]] <- "32G"

args[["imgdir"]] <- "data/test/human/derivatives/v2/547/jacobians/absolute/"
args[["demographics"]] <- "data/human/registration/v2/subject_info/demographics.csv"
# args[["mask"]] <- "data/human/registration/v2/reference_files/mask_3.0mm.mnc"
args[["mask"]] <- "data/human/registration/v2/reference_files/mask_0.8mm_3.0mm.mnc"
args[["outdir"]] <- "data/test/human/derivatives/v2/547/effect_sizes/resolution_3.0/absolute/"
args[["slurm-njobs"]] <- 50
args[["slurm-time"]] <- 30
args[["slurm-mem"]] <- "8G"

args[["batch"]] <- "Site-Scanner"
args[["nbatches"]] <- 1
args[["matrix-file"]] <- "effect_sizes.csv"
args[["matrix-res"]] <- 3.0
args[["nproc"]] <- 8
args[["execution"]] <- "slurm"


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
execution <- args[["execution"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

if (execution == "local") {
  resources <- list() 
} else if (execution == "slurm") {
  njobs <- args[["slurm-njobs"]]
  resources <- list(memory = args[["slurm-mem"]],
                    walltime = args[["slurm-time"]]*60)
} else {
  stop()
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


filenames = imgfiles
fun = compute_normative_zscore
demographics = demographics
group = group
batch = batch
df = df
mask = mask
batches = njobs
source = file.path(SRCPATH, "processing.R")
cleanup = FALSE
return_raw = TRUE
resources = resources
cores = nproc
slab_sizes <- NULL
temp_dir <- getwd()
mask_vals <- NULL
tinyMask = FALSE



# Convert voxel list into matrix
voxels_simplify_masked <- simplify_masked(voxels[["vals"]])


voxels_simplify2minc <- qMincApply(filenames = imgfiles,
                                   fun = compute_normative_zscore,
                                   demographics = demographics,
                                   group = group,
                                   batch = batch,
                                   df = df,
                                   mask = mask,
                                   batches = njobs,
                                   source = file.path(SRCPATH, "processing.R"),
                                   cleanup = FALSE,
                                   resources = resources)


mask_vol <- mincGetVolume(mask)
ind_mask <- mask_vol > 0.5

voxels_simplify2minc_masked <- voxels_simplify2minc[ind_mask,]

j <- sample(1:ncol(voxels_simplify_masked), size = 1)
print(j)
img_simplify_masked <- voxels_simplify_masked[,j]
img_simplify2minc_masked <- voxels_simplify2minc_masked[,j]
sum(img_simplify_masked - img_simplify2minc_masked)

# 
# imgdir_true <- "data/human/derivatives/v2/700/effect_sizes/resolution_0.8/absolute/"
# imgfiles_true <- list.files(imgdir_true, full.names = TRUE, pattern = "*.mnc")
# 
# img_true <- import_image(imgfiles_true[1], mask = mask, flatten = TRUE)

