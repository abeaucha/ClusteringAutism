# Packages
suppressPackageStartupMessages(library(tidyverse))

# Functions
filter_latent_space <- function(infile, outfile, ind) {
  data.table::fread(infile, header = TRUE) %>%
    as_tibble() %>%
    filter(ind) %>%
    data.table::fwrite(file = outfile)
}

# Directories
expr_dir <- "data/human/expression/"
output_dir <- "vae_latent_space"
input_dir <- paste0(output_dir, "_all_samples")

input_dir <- file.path(expr_dir, input_dir)
output_dir <- file.path(expr_dir, output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Microarray sample annotations filter
annotations <- "AHBA_microarray_sample_annotations.csv"
annotations <- file.path(expr_dir, annotations)
annotations <- read_csv(annotations, show_col_types = FALSE)

# Apply function over files
input_files <- list.files(input_dir, full.names = TRUE)
output_files <- file.path(output_dir, basename(input_files))
# result <- parallel::mcmapply(filter_latent_space, input_files, output_files,
#                              MoreArgs = list(ind = annotations[["keep"]]),
#                              mc.cores = 8)
result <- mapply(
  filter_latent_space,
  input_files,
  output_files,
  MoreArgs = list(ind = annotations[["keep"]])
)
