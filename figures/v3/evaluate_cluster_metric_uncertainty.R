# Packages -------------------------------------------------------------------

library(doParallel)

# Environment variables ------------------------------------------------------

PROJECTPATH <- Sys.getenv("PROJECTPATH")
SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "analysis.R"))


import_mouse_effect_sizes <- function(path, clusters, mask) {
  list_es <- vector(mode = "list", length = 2)
  names(list_es) <- c("absolute", "relative")
  for (i in 1:length(list_es)) {
    jacobians <- str_to_title(names(list_es)[i])

    es_files <- clusters %>%
      pull(ID) %>%
      str_replace(".mnc", paste0("_ES_", jacobians, "_200.mnc"))

    es_files <- file.path(path, es_files)

    es_mat <- import_images(
      imgfiles = es_files,
      mask = mask,
      output_format = "matrix",
      version = "v2"
    )

    es_mat[is.nan(es_mat)] <- 0

    rownames(es_mat) <- clusters$ID

    list_es[[i]] <- es_mat
  }

  return(list_es)
}


import_human_effect_sizes <- function(path) {
  list_es <- vector(mode = "list", length = 2)
  names(list_es) <- c("absolute", "relative")
  for (i in 1:length(list_es)) {
    jacobians <- names(list_es)[i]

    es_file <- file.path(path, jacobians, "effect_sizes.csv")

    df_es <- as_tibble(data.table::fread(es_file, header = TRUE))

    es_mat <- df_es %>%
      column_to_rownames("file") %>%
      as.matrix()

    list_es[[i]] <- es_mat
  }

  return(list_es)
}


sample_effect_sizes <- function(x, size = 0.8, replace = FALSE, seed = NULL) {
  N <- nrow(x[[1]])
  n <- floor(size * N)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  idx <- sample(1:N, size = n, replace = replace)
  x_sampled <- map(.x = x, .f = function(x) {
    x[idx, ]
  })

  return(x_sampled)
}


sample_snf_metrics <- function(x, size = 0.8, replace = FALSE, seed = NULL) {
  x_sampled <- sample_effect_sizes(
    x = x,
    size = size,
    replace = replace,
    seed = seed
  )
  W <- similarity_network(
    x = x_sampled,
    metric = "correlation",
    K = 10,
    sigma = 0.5,
    t = 20,
    outfile = NULL
  )

  metrics <- estimate_cluster_metrics(W = W, NUMC = 2:10)

  return(metrics)
}


sample_ARI <- function(x, clusters, size = 0.8, replace = FALSE, seed = NULL) {
  x_sampled <- sample_effect_sizes(
    x = x,
    size = size,
    replace = replace,
    seed = seed
  )
  W_sampled <- similarity_network(
    x = x_sampled,
    metric = "correlation",
    K = 10,
    sigma = 0.5,
    t = 20,
    outfile = NULL
  )

  clusters_sampled <- create_clusters(W = W_sampled, nk = 10)
  clusters_sampled <- clusters_sampled %>%
    arrange(ID)

  clusters <- clusters %>%
    semi_join(clusters_sampled, by = "ID") %>%
    arrange(ID)

  clusters_sampled <- clusters_sampled %>%
    column_to_rownames("ID")

  clusters <- clusters %>%
    column_to_rownames("ID")

  df_ARI <- tibble(nk = 2:10, ARI = 0)
  for (i in 1:nrow(df_ARI)) {
    nki <- df_ARI[[i, "nk"]]
    nk_col <- paste0("nk", nki)

    ARI <- mclust::adjustedRandIndex(
      x = clusters[[nk_col]],
      y = clusters_sampled[[nk_col]]
    )
    df_ARI[[i, "ARI"]] <- ARI
  }

  return(df_ARI)
}


# Main -----------------------------------------------------------------------

# dataset <- "MICe"
dataset <- "POND"
# dataset <- "HBN"

B <- 500
replace <- FALSE
nproc <- 8
size <- 0.8


if (dataset == "MICe") {
  params_id <- "107"

  pipeline_dir <- file.path(
    PROJECTPATH,
    "data/mouse/derivatives/v3/",
    params_id
  )

  es_dir <- file.path(pipeline_dir, "effect_sizes", "200")
  cluster_dir <- file.path(pipeline_dir, "clusters", "resolution_0.2")

  mask_file <- file.path(
    PROJECTPATH,
    "data/mouse/registration/reference_files",
    "scanbase_second_level-nlin-3_mask_200um.mnc"
  )

  clusters_file <- file.path(cluster_dir, "clusters.csv")
  clusters <- read_csv(clusters_file, show_col_types = FALSE)

  import_effect_sizes <- import_mouse_effect_sizes
} else if (dataset == "POND") {
  params_id <- "700"

  pipeline_dir <- file.path(
    PROJECTPATH,
    "data/human/derivatives/v3/",
    params_id
  )

  es_dir <- file.path(pipeline_dir, "effect_sizes", "resolution_3.0")
  cluster_dir <- file.path(pipeline_dir, "clusters", "resolution_3.0")

  clusters_file <- file.path(cluster_dir, "clusters.csv")
  clusters <- read_csv(clusters_file, show_col_types = FALSE)

  import_effect_sizes <- import_human_effect_sizes
} else if (dataset == "HBN") {
  params_id <- "013"

  pipeline_dir <- file.path(
    PROJECTPATH,
    "data/human/derivatives/v3/",
    params_id
  )

  es_dir <- file.path(pipeline_dir, "effect_sizes", "resolution_3.0")
  cluster_dir <- file.path(pipeline_dir, "clusters", "resolution_3.0")

  clusters_file <- file.path(cluster_dir, "clusters.csv")
  clusters <- read_csv(clusters_file, show_col_types = FALSE)

  import_effect_sizes <- import_human_effect_sizes
} else {
  stop()
}


#
if (dataset == "MICe") {
  list_es <- import_effect_sizes(
    path = es_dir,
    clusters = clusters,
    mask = mask_file
  )
} else {
  list_es <- import_effect_sizes(path = es_dir)
}


# pb <- txtProgressBar(max = B, style = 3)
# progress <- function(n) {
#   setTxtProgressBar(pb = pb, value = n)
# }

cl <- makePSOCKcluster(nproc)
registerDoParallel(cl)
# cl <- makeSOCKcluster(nproc)
# registerDoSNOW(cl)
# opts <- list(progress = progress)

list_snf_metrics <- foreach(
  i = 1:B,
  .packages = c("tidyverse", "SNFtool")
) %dopar%
  {
    sample_snf_metrics(x = list_es, seed = i, size = size, replace = replace)
  }

list_ARI <- foreach(
  i = 1:B,
  .packages = c("tidyverse", "SNFtool")
) %dopar%
  {
    sample_ARI(
      x = list_es,
      clusters = clusters,
      seed = i,
      size = size,
      replace = replace
    )
  }

# close(pb)
stopCluster(cl)

df_snf_metrics <- bind_rows(list_snf_metrics, .id = "seed")
df_ARI <- bind_rows(list_ARI, .id = "seed")
df_metrics <- inner_join(df_snf_metrics, df_ARI, by = c("seed", "nk"))

outfile <- paste0("cluster_metrics_uncertainty_", dataset, ".csv")
outfile <- file.path(PROJECTPATH, "figures/v3/resources", outfile)
write_csv(df_metrics, outfile)
