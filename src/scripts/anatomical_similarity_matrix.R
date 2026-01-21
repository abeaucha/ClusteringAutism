# ----------------------------------------------------------------------------
# anatomical_similarity_matrix.R
# Author: Antoine Beauchamp
# Created: November 1st, 2022
#
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(tcltk))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option(
    '--mouse-dir',
    type = 'character',
    help = paste("Path to direcotry")
  ),
  make_option("--mouse-labels", type = "character", help = paste("")),
  make_option("--mouse-mask", type = "character", help = paste("")),
  make_option("--mouse-tree", type = "character", help = paste("")),
  make_option(
    '--human-dir',
    type = 'character',
    help = paste("Path to direcotry")
  ),
  make_option("--human-tree", type = "character", help = paste("")),
  make_option("--human-metadata", type = "character", help = paste("")),
  make_option("--tree-labels", type = "character", help = paste("")),
  make_option(
    '--metric',
    type = 'character',
    default = 'correlation',
    help = paste(
      "Metric used to compute the similarity matrix.",
      "[default %default]"
    )
  ),
  make_option(
    '--outfile',
    type = 'character',
    help = paste("Path to .csv file containing similarity matrix.")
  ),
  make_option(
    '--save-intermediate',
    type = 'character',
    default = 'false',
    help = paste(
      "Option to save intermediate latent space",
      "similarity matrices. File names will be ",
      "adapted from --outfile. [default %default]"
    )
  ),
  make_option(
    '--verbose',
    type = 'character',
    default = 'true',
    help = paste("Verbosity option. [default %default]")
  ),
  make_option(
    "--parallel",
    type = "character",
    default = 'false',
    help = "Option to run in parallel. [default %default]"
  ),
  make_option(
    "--nproc",
    type = "numeric",
    help = paste(
      "Number of processors to use in parallel.",
      "Ignored if --parallel is false."
    )
  )
)


# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>%
  str_subset("--file=") %>%
  str_remove("--file=") %>%
  dirname()

path_func <- file.path(
  working_dir,
  script_dir,
  "functions",
  "buildSimilarityMatrix.R"
)

path_func_tree <- file.path(
  working_dir,
  script_dir,
  "functions",
  "tree_tools.R"
)
source(path_func)
source(path_func_tree)


#' Title
#'
#' @param x
#' @param tree
#' @param labels
#' @param mask
#'
#' @return
label_mouse_data <- function(x, tree, labels, mask) {
  labels_pruned <- hanatToAtlas(tree, mincArray(labels))
  defs_pruned <- hanatToAtlasDefs(tree)

  labels_pruned <- labels_pruned[mask == 1]

  x <- x %>%
    mutate(label = labels_pruned) %>%
    left_join(defs_pruned, by = c('label' = 'Label')) %>%
    select(-label) %>%
    rename(name = Structure)

  return(x)
}

#' Title
#'
#' @param x
#' @param tree
#' @param samples
#'
#' @return
label_human_data <- function(x, tree, samples) {
  tree <- Clone(tree)

  tree$Do(
    function(node) {
      node$struct <- rep(node$name, length(node$samples))
    },
    traversal = 'post-order'
  )

  tree_structs <- unlist(tree$Get('struct', filterFun = isLeaf))
  names(tree_structs) <- NULL

  tree_samples <- unlist(tree$Get('samples', filterFun = isLeaf))
  names(tree_samples) <- NULL

  sample_defs <- tibble(name = tree_structs, sample_id = tree_samples)

  x <- x %>%
    mutate(sample_id = samples) %>%
    left_join(sample_defs, by = 'sample_id') %>%
    select(-sample_id)

  return(x)
}


aggregate_data <- function(x, groupby, method = 'mean', output = 'tibble') {
  if (method == 'mean') {
    aggfunc <- mean
  } else if (method == 'median') {
    aggfunc <- median
  } else if (method == 'sd') {
    aggfunc <- sd
  } else {
    stop()
  }

  if (!(groupby %in% colnames(x))) {
    stop()
  }

  if (any(is.na(x[[groupby]]))) {
    warning(paste(
      "Grouping variable",
      groupby,
      "contains NA.",
      "These observations will be ommitted when aggregating."
    ))
  }

  x <- x %>%
    filter(!is.na(.data[[groupby]])) %>%
    group_by_at(.vars = vars(all_of(groupby))) %>%
    summarise_all(.funs = aggfunc) %>%
    ungroup()

  if (output == 'tibble') {
    x %>%
      as_tibble() %>%
      return()
  } else if (output == 'data.frame') {
    x %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      return()
  } else if (output == 'matrix') {
    x %>%
      column_to_rownames(groupby) %>%
      as.matrix() %>%
      return()
  } else {
    stop()
  }
}


compute_anatomical_similarity <- function(
  latent_space,
  mouse_files,
  human_files,
  mouse_tree,
  mouse_labels,
  mouse_mask,
  human_tree,
  human_samples,
  metric,
  save_intermediate = FALSE,
  outfile = "similarity.csv"
) {
  mouse_file <- str_subset(
    mouse_files,
    str_c("mousetransform_", latent_space, ".csv")
  )
  human_file <- str_subset(
    human_files,
    str_c("humantransform_", latent_space, ".csv")
  )

  mouse <- data.table::fread(mouse_file, header = TRUE) %>%
    as_tibble()
  human <- data.table::fread(human_file, header = TRUE) %>%
    as_tibble()

  mouse <- mouse %>%
    label_mouse_data(
      tree = tree_mouse_pruned,
      labels = mouse_labels,
      mask = mouse_mask
    ) %>%
    filter(!is.na(name)) %>%
    aggregate_data(groupby = "name", output = "matrix") %>%
    t()

  human <- human %>%
    label_human_data(tree = tree_human_pruned, samples = human_sample_ids) %>%
    filter(!is.na(name)) %>%
    aggregate_data(groupby = "name", output = "matrix") %>%
    t()

  similarity <- buildSimilarityMatrix(x1 = human, x2 = mouse, method = metric)

  similarity <- similarity[, match(mouse_regions, colnames(similarity))]
  similarity <- similarity[match(human_regions, rownames(similarity)), ]

  if (save_intermediate) {
    data.table::fwrite(
      x = as_tibble(similarity, rownames = "Human"),
      file = outfile
    )
  }

  return(similarity)
}

# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
mouse_dir <- args[["mouse-dir"]]
mouse_labels <- args[["mouse-labels"]]
mouse_mask <- args[["mouse-mask"]]
mouse_tree <- args[["mouse-tree"]]
human_dir <- args[["human-dir"]]
human_tree <- args[["human-tree"]]
human_metadata <- args[["human-metadata"]]
tree_labels <- args[["tree-labels"]]
metric <- args[["metric"]]
outfile <- args[["outfile"]]
save_intermediate <- ifelse(args[["save-intermediate"]] == "true", TRUE, FALSE)
inparallel <- ifelse(args[["parallel"]] == "true", TRUE, FALSE)
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# mouse_dir <- "data/mouse/expression/latent_space/"
# mouse_labels <- "data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc"
# mouse_mask <- "data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc"
# mouse_tree <- "data/mouse/expression/MouseExpressionTree_DSURQE.RData"
# human_dir <- "data/human/expression/latent_space/"
# human_tree <- "data/human/expression/HumanExpressionTree.RData"
# human_metadata <- "data/human/expression/SampleInformation_pipeline_abagen.csv"
# tree_labels <- "data/TreeLabelsReordered.RData"
# metric <- "correlation"
# outfile <- "tmp.csv"
# save_intermediate <- FALSE

#Metric error catch
metric_choices <- c("correlation", "cosine", "euclidean")
if (!(metric %in% metric_choices)) {
  stop(paste(
    "Argument --metric must be one of",
    "[correlation, cosine, euclidean]"
  ))
}

#Create outdir if needed
outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Trees
load(mouse_tree)
tree_mouse <- Clone(treeMouseExpr)
rm(treeMouseExpr)

load(human_tree)
tree_human <- Clone(treeHumanExpr)
rm(treeHumanExpr)

#Tree labels
load(tree_labels)

#Prune mouse tree to desired level
mouse_regions <- c(
  listLabelsMouseReordered$Region67,
  "fiber tracts",
  "ventricular systems"
)
tree_mouse_pruned <- Clone(tree_mouse)
pruneAnatTree(tree_mouse_pruned, nodes = mouse_regions, method = "BelowNode")

#Pruned human tree to desired level
human_regions <- c(
  listLabelsHumanReordered$Region88,
  "white matter",
  "sulci & spaces"
)
tree_human_pruned <- Clone(tree_human)
pruneAnatTree(tree_human_pruned, nodes = human_regions, method = "BelowNode")

#Mouse labels and mask
mouse_labels <- mincGetVolume(mouse_labels)
mouse_mask <- mincGetVolume(mouse_mask)

#Human microarray sample IDs
human_sample_ids <- data.table::fread(human_metadata, header = TRUE) %>%
  as_tibble() %>%
  pull(SampleID)

mouse_files <- list.files(mouse_dir, full.names = TRUE) %>%
  str_subset("mousetransform")
human_files <- list.files(human_dir, full.names = TRUE) %>%
  str_subset("humantransform")

mouse_latent_space_ids <- mouse_files %>%
  str_extract("mousetransform_[0-9]+.csv") %>%
  str_extract("[0-9]+") %>%
  as.numeric() %>%
  sort()

human_latent_space_ids <- human_files %>%
  str_extract("humantransform_[0-9]+.csv") %>%
  str_extract("[0-9]+") %>%
  as.numeric() %>%
  sort()

if (length(mouse_latent_space_ids) != length(human_latent_space_ids)) {
  stop("Different number of latent spaces")
} else {
  if (sum(mouse_latent_space_ids != human_latent_space_ids) != 0) {
    stop("Inconsistent latent spaces")
  }
}

latent_space_ids <- mouse_latent_space_ids

pb <- txtProgressBar(max = length(latent_space_ids), style = 3)
progress <- function(n) {
  setTxtProgressBar(pb = pb, value = n)
}
if (inparallel) {
  nproc <- args[["nproc"]]
  cl <- makeSOCKcluster(nproc)
  registerDoSNOW(cl)
  opts <- list(progress = progress)
  similarity <- foreach(
    i = 1:length(latent_space_ids),
    .packages = c("tidyverse", "RMINC", "data.tree"),
    .options.snow = opts
  ) %dopar%
    {
      compute_anatomical_similarity(
        latent_space = latent_space_ids[i],
        mouse_files = mouse_files,
        human_files = human_files,
        mouse_tree = tree_mouse_pruned,
        mouse_labels = mouse_labels,
        mouse_mask = mouse_mask,
        human_tree = tree_human_pruned,
        human_samples = human_samples,
        metric = metric,
        save_intermediate = save_intermediate,
        outfile = str_replace(
          outfile,
          ".csv",
          str_c("_", latent_space_ids[i], ".csv")
        )
      )
    }
  close(pb)
  stopCluster(cl)
} else {
  similarity <- foreach(
    i = 1:length(latent_space_ids),
    .packages = c("tidyverse", "RMINC", "data.tree")
  ) %do%
    {
      progress(n = i)
      compute_anatomical_similarity(
        latent_space = latent_space_ids[i],
        mouse_files = mouse_files,
        human_files = human_files,
        mouse_tree = tree_mouse_pruned,
        mouse_labels = mouse_labels,
        mouse_mask = mouse_mask,
        human_tree = tree_human_pruned,
        human_samples = human_samples,
        metric = metric,
        save_intermediate = save_intermediate,
        outfile = str_replace(
          outfile,
          ".csv",
          str_c("_", latent_space_ids[i], ".csv")
        )
      )
    }
  close(pb)
}

similarity_array <- array(
  unlist(similarity),
  dim = c(nrow(similarity[[1]]), ncol(similarity[[1]]), length(similarity))
)
similarity_avg <- rowMeans(similarity_array, dims = 2)
rownames(similarity_avg) <- rownames(similarity[[1]])
colnames(similarity_avg) <- colnames(similarity[[1]])

data.table::fwrite(
  x = as_tibble(similarity_avg, rownames = "Human"),
  file = outfile
)
