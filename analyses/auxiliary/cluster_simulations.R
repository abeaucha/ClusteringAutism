suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(SNFtool))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doSNOW))

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")



#' Create clusters from SNF affinity matrix
#'
#' @param W (matrix) SNF affinity matrix.
#' @param nk (numeric scalar) Maximum number of clusters to use in 
#' clustering.
#' @param outfile (character scalar) Optional path to .csv file in 
#' which to save cluster assignments.
#'
#' @return (data.frame) Cluster assignments.
create_clusters <- function(W, nk = 10, outfile = NULL) {
  
  for(k in 2:nk) {
    group <- spectralClustering(affinity = W, K = k)
    group_name <- paste0("nk", k)
    assign(group_name, group)
    if (k == 2) {
      if (is.null(rownames(W))) {
        ids <- as.character(1:nrow(W))
      } else {
        ids <- rownames(W)
      }
      all_clusters <- data.frame(ids, group, stringsAsFactors = F)
      colnames(all_clusters) <- c("ID", group_name)
    } else {
      group <- data.frame(group)
      colnames(group) <- group_name
      all_clusters <- cbind(all_clusters, group)
    }
  }
  
  if (!is.null(outfile)) {
    write.csv(x = all_clusters, file = outfile, row.names = FALSE)
  }
  
  return(all_clusters)
  
}

estimate_cluster_metrics <- function (W, NUMC = 2:5){
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC <- NUMC[NUMC > 1]
  }
  W <- (W + t(W))/2
  diag(W) <- 0
  if (length(NUMC) <= 0) {
    warning(paste("Invalid NUMC provided, must be an integer vector",
                  "with atleast one other number than 1.", "Using default NUMC=c(2,3,4,5)",
                  sep = ""))
    NUMC <- 2:5
  }
  degs <- rowSums(W)
  degs[degs == 0] <- .Machine$double.eps
  D <- diag(degs)
  L <- D - W
  Di <- diag(1/sqrt(degs))
  L <- Di %*% L %*% Di
  eigs <- eigen(L)
  eigs_order <- sort(eigs$values, index.return = T)$ix
  eigs$values <- eigs$values[eigs_order]
  eigs$vectors <- eigs$vectors[, eigs_order]
  eigengap <- abs(diff(eigs$values))
  quality <- list()
  for (c_index in 1:length(NUMC)) {
    ck <- NUMC[c_index]
    UU <- eigs$vectors[, 1:ck]
    EigenvectorsDiscrete <- SNFtool:::.discretisation(UU)[[1]]
    EigenVectors <- EigenvectorsDiscrete^2
    temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors),
                                                function(i) EigenVectors[, i])), ]
    temp1 <- t(apply(temp1, 1, sort, TRUE))
    quality[[c_index]] <- (1 - eigs$values[ck + 1])/(1 -
                                                       eigs$values[ck]) * sum(sum(diag(1/(temp1[, 1] + .Machine$double.eps)) %*%
                                                                                    temp1[, 1:max(2, ck - 1)]))
  }
  
  out <- tibble(nk = NUMC,
                eigengap = eigengap[NUMC],
                rotation = unlist(quality))
  
  return(out)
}



simulate_clusters <- function(i, grid, nobs, outdir = "./", run_umap = FALSE) {
  
  # print(paste("i:", i, "of", length(results)))
  # 
  nvars <- grid[[i, "nvars"]]
  delta <- grid[[i, "delta"]]
  seed <- grid[[i, "seed"]]
  
  # Generate normal distributions
  set.seed(seed)
  
  x1 <- matrix(rnorm(nobs[1]*nvars), nrow = nobs[1], ncol = nvars)
  x2_1 <- matrix(rnorm(nobs[2]*(nvars-1)), nrow = nobs[2], ncol = nvars-1)
  x2_2 <- rnorm(nobs[2], mean = delta)
  x2 <- cbind(x2_1, x2_2)
  x <- rbind(x1, x2)
  rownames(x) <- paste0("p", 1:sum(nobs))
  
  # x1 <- rmvnorm(n = nobs_nk2[1], mean = rep(0, nvars), sigma = diag(nvars))
  # x2 <- rmvnorm(n = nobs_nk2[2], mean = c(delta, rep(0, nvars-1)), sigma = diag(nvars))
  # x <- rbind(x1, x2)
  
  # UMAP projection  
  if (run_umap) {
    umap_out <- umap(x, n_components = 2, 
                     random_state = seed)
  }
  
  # Distance matrix
  d <- 1 - cor(t(x))
  
  # Free up memory
  rm(list = c("x1", "x2", "x")); gc()
  
  # Affinity matrix
  W <- affinityMatrix(d, K = 10, sigma = 0.5)
  
  # Run SNF because it seems to apply some normalization?
  W <- SNF(list(W, W), K = 10, t = 20)
  
  # Estimate clustering metrics
  metrics <- estimate_cluster_metrics(W = W, NUMC = 2:10)
  
  
  if (run_umap) {
    clusters <- create_clusters(W = W, nk = 2)
    
    umap_x <- umap_out$layout
    colnames(umap_x) <- paste0("x", 1:2)
    p_umap <- umap_x %>% 
      as_tibble() %>% 
      bind_cols(clusters) %>% 
      ggplot(aes(x = x1, y = x2, col = factor(nk2))) + 
      geom_point() + 
      coord_equal() +
      labs(color = "Cluster",
           title = paste0("Number of variables: ", nvars, "; Delta: ", delta, "; Seed: ", seed)) + 
      theme_bw()
    
    outfile_umap <- paste("umap_nvars", nvars, "delta", delta, "seed", seed, sep = "_")
    outfile_umap <- paste0(outfile_umap, ".pdf")
    outfile_umap <- file.path(outdir, outfile_umap)
    pdf(file = outfile_umap,
        width = unit(8, "inch"),
        height = unit(8, "inch"))
    print(p_umap)
    dev.off()
    
  }
  
  results <- metrics %>% 
    mutate(nvars = nvars, delta = delta, seed = seed)
  
  outfile_metrics <- paste("metrics_nvars", nvars, "delta", delta, "seed", seed, sep = "_")
  outfile_metrics <- paste0(outfile_metrics, ".csv")
  outfile_metrics <- file.path(outdir, outfile_metrics)
  write_csv(x = results, file = outfile_metrics)
  
  return(1)
  
}


output_dir <- file.path(PROJECTPATH, "analyses", "auxiliary", "outputs", "cluster_simulations")
if (!file.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

human_pipeline_dir <- file.path(PROJECTPATH, "data/human/derivatives/v3/700/")
human_clusters_dir <- file.path(human_pipeline_dir, "clusters", "resolution_3.0")
human_es_dir <- file.path(human_pipeline_dir, "effect_sizes", "resolution_3.0")

human_clusters <- "clusters.csv"
human_clusters <- file.path(human_clusters_dir, human_clusters)
df_human_clusters <- read_csv(human_clusters, show_col_types = FALSE)

human_affinity <- "affinity.csv"
human_affinity <- file.path(human_clusters_dir, human_affinity)
df_human_affinity <- read_csv(human_affinity, show_col_types = FALSE)

human_es <- "effect_sizes.csv"
human_es_rel <- file.path(human_es_dir, "relative", human_es)
df_human_es_rel <- as_tibble(data.table::fread(human_es_rel, header = TRUE))
nvoxels <- ncol(df_human_es_rel)-1
rm("df_human_es_rel"); gc()

df_human_clusters_nk <- df_human_clusters %>% 
  select(ID, k = nk2)

# Number of observations per cluster
nobs <- df_human_clusters_nk %>% 
  group_by(k) %>% 
  count() %>% 
  pull(n)

# Number of variables
# nvars <- nvoxels

niter <- 50
df_grid <- expand_grid(nvars = nvoxels, 
                       delta = seq(0, 20, by = 1),
                       seed = 1:niter)
# results <- vector(mode = "list", length = nrow(df_grid))


pb <- txtProgressBar(max = nrow(df_grid), style = 3)
progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
cl <- makeSOCKcluster(6)
registerDoSNOW(cl)
opts <- list(progress=progress)
results <- foreach(i = 1:nrow(df_grid),
                   .packages = c("tidyverse", "SNFtool", "umap"),
                   .options.snow = opts) %dopar% {
                     simulate_clusters(i = i,
                                       grid = df_grid, 
                                       nobs = nobs, 
                                       outdir = output_dir,
                                       run_umap = FALSE)
                   }
close(pb)
stopCluster(cl)

# df_metrics <- bind_rows(results)
# 
# outfile <- "clustering_metrics.csv"
# outfile <- file.path(output_dir, outfile)
# write_csv(x = df_metrics, file = outfile)

