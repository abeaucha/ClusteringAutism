suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(MRIcrotome))

SRCPATH <- Sys.getenv("SRCPATH")

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))

fit_predict_mean <- function(y, demographics, newdata, batch = NULL, df = 3) {
  
  # Residualize using batch variable if specified
  if (!is.null(batch)) {
    batch <- demographics %>%
      select(all_of(batch)) %>%
      unite(col = batch) %>%
      pull(batch)
    y <- residuals(lm(y ~ batch))
    names(y) <- NULL
  }
  
  ind_fit <- demographics[["DX"]] == "Control"
  
  # Training data frame
  df_fit <- demographics[ind_fit, c("Age", "Sex")]
  df_fit[["y"]] <- y[ind_fit]
  
  model_fit <- lm(y ~ Sex + ns(Age, df = df), data = df_fit)
  
  y_pred <- predict(model_fit, newdata = newdata)
  
  return(y_pred)
  
}


# Output directory
output_dir <- "outputs/human_normative_models/"

# Pipeline directory
version <- "v3"
pipeline_dir <- "../../data/human/derivatives/"
pipeline_dir <- file.path(pipeline_dir, version)

# Registration directory
registration_dir <- "../../data/human/registration/v3/"


# Identify parameter set ID
metadata <- file.path(pipeline_dir, "metadata.csv")
params <- fetch_params_metadata(metadata, 
                                dataset = "POND-SickKids",
                                resolution = 3.0,
                                es_group = "patients")


# Parameter set ID
param_id <- "547"

# Image resolution
resolution <- params %>% 
  filter(id == param_id) %>% 
  pull(resolution)
resolution <- sprintf("%.1f", resolution)

# Effect size batch variables
batch <- params %>% 
  filter(id == param_id) %>% 
  pull(es_batch)

df <- params %>% 
  filter(id == param_id) %>% 
  pull(es_df)

# Jacobians
jacobians <- c("absolute", "relative")

# Directory to Jacobian images
jacobians_dir <- file.path(registration_dir, "jacobians_resampled/resolution_3.0/")

# Image mask
mask <- paste0("mask_", resolution, "mm.mnc")
mask <- file.path(registration_dir, "reference_files", mask)

# Path to demographics
demographics <- file.path(pipeline_dir, param_id, "demographics.csv")


demographics <- read_csv(demographics, show_col_types = FALSE)

demographics <- demographics %>%
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex))

if (!is.null(batch)) {
  batch <- str_split(batch, pattern = "-")[[1]]
  batch_check <- batch %in% colnames(demographics)
  if (!all(batch_check)) {
    stop("Batch columns not found in demographics:\n",
         str_flatten(batch, collapse = "\n"))
  }
}


j <- 1

# Image files
imgdir <- file.path(jacobians_dir, jacobians[j])
imgfiles <- list.files(imgdir, full.names = TRUE)

imgs_in_demographics <- basename(imgfiles) %in% demographics[["file"]]
imgfiles <- imgfiles[imgs_in_demographics]
row_match <- match(basename(imgfiles), demographics[["file"]])
demographics <- demographics[row_match,]


ages <- seq(5, 25, by = 5)
sex <- c("Male", "Female")
df_pred <- expand_grid(Age = ages,
                       Sex = sex)

voxels <- qMincApply(filenames = imgfiles,
                     fun = fit_predict_mean,
                     newdata = df_pred,
                     demographics = demographics,
                     batch = batch,
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

