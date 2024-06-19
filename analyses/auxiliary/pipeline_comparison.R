suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(parallel))

source("src/processing.R")

image_cor <- function(imgs, masks) {
  imgs <- map2(.x = imgs, .y = masks, .f = import_image)
  correlation <- reduce(.x = imgs, .f = cor)
  return(correlation)
}

# Pipeline v2 700 ------------------------------------------------------------

registration_dir_old <- "data/human/registration/v2/"
registration_dir_new <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/human/registration/v2/"

pipeline_dir_old <- "data/human/derivatives/v2/700/"
pipeline_dir_new <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/human/derivatives/v2/700/"

registration_dir <- list(old = registration_dir_old,
                         new = registration_dir_new)

pipeline_dir <- list(old = pipeline_dir_old,
                     new = pipeline_dir_new)


version <- c("old", "new")
es_dir = list(old="", new="")
cluster_dir = list(old="", new="")
centroid_dir = list(old="", new="")
for (i in version) {
  es_dir[[i]] = file.path(pipeline_dir[[i]], "effect_sizes")
  cluster_dir[[i]] = file.path(pipeline_dir[[i]], "clusters")
  centroid_dir[[i]] = file.path(pipeline_dir[[i]], "centroids")
}
centroid_dir[["old"]] <- file.path(pipeline_dir[["old"]], "cluster_maps")

## Effect sizes at 0.8mm -----------------------------------------------------
resolution <- "resolution_0.8"
masks <- map(.x = registration_dir, .f = function(x){file.path(x, "reference_files", "mask_0.8mm.mnc")})

es_dir_abs <- map(.x = es_dir, .f = function(x){file.path(x, resolution, "absolute")}) 
es_dir_rel <- map(.x = es_dir, .f = function(x){file.path(x, resolution, "relative")}) 

es_files_abs <- map(.x = es_dir_abs, .f = function(x){sort(list.files(x, full.names = TRUE, pattern = "*.mnc"))})
es_files_rel <- map(.x = es_dir_rel, .f = function(x){sort(list.files(x, full.names = TRUE, pattern = "*.mnc"))})

### Absolute effect sizes ----------------------------------------------------
length(es_files_abs[["old"]]) == length(es_files_abs[["new"]])

nfiles <- length(es_files_abs[["old"]])

es_files_abs_reorder <- vector(mode = "list", length = nfiles)
for (i in 1:nfiles) {
  es_files_abs_reorder[[i]][["old"]] <- es_files_abs[["old"]][[i]]
  es_files_abs_reorder[[i]][["new"]] <- es_files_abs[["new"]][[i]]
}

ti <- Sys.time()
es_abs_cor <- mclapply(X = es_files_abs_reorder, 
                       FUN = image_cor, 
                       masks = masks, 
                       mc.cores = 8)
es_abs_cor <- reduce(es_abs_cor, c)
tf <- Sys.time()
print(tf - ti)

message(paste("Mean correlation:", mean(es_abs_cor)))
message(paste("Minimum correlation:", min(es_abs_cor)))

hist(es_abs_cor)

# Outcomes: 
# Mean correlation of 1.0
# Min correlation of 1.0
# Images with correlation < 0.95: 0
# Images with correlation < 0.9: 0


### Relative effect sizes ----------------------------------------------------
length(es_files_rel[["old"]]) == length(es_files_rel[["new"]])

es_files_rel_reorder <- vector(mode = "list", length = nfiles)
for (i in 1:nfiles) {
  es_files_rel_reorder[[i]][["old"]] <- es_files_rel[["old"]][[i]]
  es_files_rel_reorder[[i]][["new"]] <- es_files_rel[["new"]][[i]]
}

ti <- Sys.time()
es_rel_cor <- mclapply(X = es_files_rel_reorder, 
                       FUN = image_cor, 
                       masks = masks, 
                       mc.cores = 12)
es_rel_cor <- reduce(es_rel_cor, c)
tf <- Sys.time()
print(tf - ti)

message(paste("Mean correlation:", mean(es_rel_cor)))
message(paste("Minimum correlation:", min(es_rel_cor)))

hist(es_rel_cor)

# Outcomes: 
# Mean correlation of 0.993
# Min correlation of 0.981
# Images with correlation < 0.95: 0
# Images with correlation < 0.9: 0


## Effect sizes resampled to 3.0mm ------------------------------------------------------
resolution <- "resolution_3.0"
masks <- list(old = file.path(es_dir[["old"]], resolution, "mask_0.8mm_3.0mm.mnc"),
              new = file.path(es_dir[["new"]], resolution, "mask_0.8mm_autocrop_3.0mm.mnc"))

es_dir_abs <- map(.x = es_dir, .f = function(x){file.path(x, resolution, "absolute")}) 
es_dir_rel <- map(.x = es_dir, .f = function(x){file.path(x, resolution, "relative")}) 

es_files_abs <- map(.x = es_dir_abs, .f = function(x){sort(list.files(x, full.names = TRUE, pattern = "*.mnc"))})
es_files_rel <- map(.x = es_dir_rel, .f = function(x){sort(list.files(x, full.names = TRUE, pattern = "*.mnc"))})


### Absolute effect sizes ----------------------------------------------------
length(es_files_abs[["old"]]) == length(es_files_abs[["new"]])

nfiles <- length(es_files_abs[["old"]])

es_files_abs_reorder <- vector(mode = "list", length = nfiles)
for (i in 1:nfiles) {
  es_files_abs_reorder[[i]][["old"]] <- es_files_abs[["old"]][[i]]
  es_files_abs_reorder[[i]][["new"]] <- es_files_abs[["new"]][[i]]
}

ti <- Sys.time()
es_abs_cor <- mclapply(X = es_files_abs_reorder, 
                       FUN = image_cor, 
                       masks = masks, 
                       mc.cores = 12)
es_abs_cor <- reduce(es_abs_cor, c)
tf <- Sys.time()
print(tf - ti)

message(paste("Mean correlation:", mean(es_abs_cor)))
message(paste("Minimum correlation:", min(es_abs_cor)))
print(which(es_abs_cor < 0.95))
print(which(es_abs_cor < 0.90))

hist(es_abs_cor)

# Outcomes: 
# Mean correlation of 1.0
# Min correlation of 1.0
# Images with correlation < 0.95: None
# Images with correlation < 0.9: None


### Relative effect sizes ----------------------------------------------------
length(es_files_rel[["old"]]) == length(es_files_rel[["new"]])

es_files_rel_reorder <- vector(mode = "list", length = nfiles)
for (i in 1:nfiles) {
  es_files_rel_reorder[[i]][["old"]] <- es_files_rel[["old"]][[i]]
  es_files_rel_reorder[[i]][["new"]] <- es_files_rel[["new"]][[i]]
}

ti <- Sys.time()
es_rel_cor <- mclapply(X = es_files_rel_reorder, 
                       FUN = image_cor, 
                       masks = masks, 
                       mc.cores = 12)
es_rel_cor <- reduce(es_rel_cor, c)
tf <- Sys.time()
print(tf - ti)

message(paste("Mean correlation:", mean(es_rel_cor)))
message(paste("Minimum correlation:", min(es_rel_cor)))
print(which(es_rel_cor < 0.95))
print(which(es_rel_cor < 0.90))

hist(es_rel_cor)

# Outcomes: 
# Mean correlation of 1.0
# Min correlation of 1.0
# Images with correlation < 0.95: 0
# Images with correlation < 0.9: 0

## Effect size matrices ------------------------------------------------------

### Absolute effect sizes ----------------------------------------------------

es_files_abs <- map(.x = es_dir_abs, .f = function(x){file.path(x, "effect_sizes.csv")})
es_matrices_abs <- map(.x = es_files_abs, .f = function(x){as_tibble(data.table::fread(x, header = TRUE))})
map(es_matrices_abs, dim)

es_matrices_abs <- map(.x = es_matrices_abs, .f = function(x){x %>% select(-file)})

es_matrices_abs_cor <- map2_dbl(.x = es_matrices_abs[["new"]],
                                .y = es_matrices_abs[["old"]], 
                                .f = cor)

mean(es_matrices_abs_cor)
median(es_matrices_abs_cor)

es_matrices_abs_dist <- map2_dbl(.x = es_matrices_abs[["new"]],
                                 .y = es_matrices_abs[["old"]], 
                                 .f = function(x,y){sum((x-y)^2)})

mean(es_matrices_abs_dist)
median(es_matrices_abs_dist)
hist(es_matrices_abs_dist)


### Relative effect sizes ----------------------------------------------------

es_files_rel <- map(.x = es_dir_rel, .f = function(x){file.path(x, "effect_sizes.csv")})
es_matrices_rel <- map(.x = es_files_rel, .f = function(x){as_tibble(data.table::fread(x, header = TRUE))})
map(es_matrices_rel, dim)

es_matrices_rel <- map(.x = es_matrices_rel, .f = function(x){x %>% select(-file)})

es_matrices_rel_cor <- map2_dbl(.x = es_matrices_rel[["new"]],
                                .y = es_matrices_rel[["old"]], 
                                .f = cor)

mean(es_matrices_rel_cor)
median(es_matrices_rel_cor)

es_matrices_rel_dist <- map2_dbl(.x = es_matrices_rel[["new"]],
                                 .y = es_matrices_rel[["old"]], 
                                 .f = function(x,y){sum((x-y)^2)})

mean(es_matrices_rel_dist)
median(es_matrices_rel_dist)
hist(es_matrices_rel_dist)


## Affinity matrices 

affinity <- cluster_dir %>% 
  map(.f = function(x){file.path(x, resolution, "affinity.csv")}) %>% 
  map(.f = read_csv, show_col_types = FALSE)

map(affinity, dim)

affinity_cor <- map2_dbl(.x = affinity[["new"]],
                         .y = affinity[["old"]], 
                         .f = cor)

mean(affinity_cor)
median(affinity_cor)

affinity_dist <- map2_dbl(.x = affinity[["new"]],
                                 .y = affinity[["old"]], 
                                 .f = function(x,y){sum((x-y)^2)})

mean(affinity_dist)
median(affinity_dist)
hist(affinity_dist)


## Clusters ------------------------------------------------------------------
resolution <- "resolution_3.0"
clusters <- cluster_dir %>% 
  map(.f = function(x){file.path(x, resolution, "clusters.csv")}) %>% 
  map(.f = read_csv, show_col_types = FALSE) %>% 
  map(.f = function(x){column_to_rownames(x, "ID")})


files <- rownames(clusters[["old"]])

nk_range <- 2:10
cluster_comparison <- vector(mode = "list", length = length(nk_range))
names(cluster_comparison) <- nk_range
for (nk in nk_range) {
  
  clusters_nk <- map(clusters, .f = function(x){x[,paste0("nk", nk)]})
  
  nk_grid <- expand_grid(k_old = 1:nk,
                         k_new = 1:nk,
                         overlap = 0)
  
  for (i in 1:nrow(nk_grid)) {
    
    k_old <- nk_grid[[i, "k_old"]]
    k_new <- nk_grid[[i, "k_new"]]
    
    ind_k_old <- clusters_nk[["old"]] == k_old
    ind_k_new <- clusters_nk[["new"]] == k_new
    
    files_k_old <- files[ind_k_old]
    files_k_new <- files[ind_k_new]
  
    nk_grid[[i, "overlap"]] <- length(intersect(files_k_old, files_k_new))/length(union(files_k_old, files_k_new))
    
  }
  
  cluster_comparison[[as.character(nk)]] <- nk_grid
  
}

df_cluster_comparison <- bind_rows(cluster_comparison, .id = "nk")

df_cluster_comparison %>% 
  group_by(nk, k_new) %>% 
  filter(overlap == max(overlap)) %>% 
  pull(overlap)

## Centroids -----------------------------------------------------------------
resolution <- "resolution_0.8"
masks <- map(.x = registration_dir, .f = function(x){file.path(x, "reference_files", "mask_0.8mm.mnc")})

centroid_dir_abs <- map(.x = centroid_dir, .f = function(x){file.path(x, resolution, "absolute")}) 
centroid_dir_rel <- map(.x = centroid_dir, .f = function(x){file.path(x, resolution, "relative")}) 

centroid_files_abs <- map(.x = centroid_dir_abs, .f = function(x){sort(list.files(x, full.names = TRUE, pattern = "*.mnc"))})
centroid_files_rel <- map(.x = centroid_dir_rel, .f = function(x){sort(list.files(x, full.names = TRUE, pattern = "*.mnc"))})

### Absolute centroids -------------------------------------------------------
length(centroid_files_abs[["old"]]) == length(centroid_files_abs[["new"]])

nfiles <- length(centroid_files_abs[["old"]])

centroid_files_abs_reorder <- vector(mode = "list", length = nfiles)
for (i in 1:nfiles) {
  centroid_files_abs_reorder[[i]][["old"]] <- centroid_files_abs[["old"]][[i]]
  centroid_files_abs_reorder[[i]][["new"]] <- centroid_files_abs[["new"]][[i]]
}

ti <- Sys.time()
centroid_abs_cor <- mclapply(X = centroid_files_abs_reorder, 
                       FUN = image_cor, 
                       masks = masks, 
                       mc.cores = 12)
centroid_abs_cor <- reduce(centroid_abs_cor, c)
tf <- Sys.time()
print(tf - ti)

message(paste("Mean correlation:", mean(centroid_abs_cor)))
message(paste("Minimum correlation:", min(centroid_abs_cor)))
print(which(centroid_abs_cor < 0.95))
print(which(centroid_abs_cor < 0.90))

hist(centroid_abs_cor)


### Relative centroids -------------------------------------------------------

centroid_files_rel_reorder <- vector(mode = "list", length = nfiles)
for (i in 1:nfiles) {
  centroid_files_rel_reorder[[i]][["old"]] <- centroid_files_rel[["old"]][[i]]
  centroid_files_rel_reorder[[i]][["new"]] <- centroid_files_rel[["new"]][[i]]
}

ti <- Sys.time()
centroid_rel_cor <- mclapply(X = centroid_files_rel_reorder, 
                             FUN = image_cor, 
                             masks = masks, 
                             mc.cores = 12)
centroid_rel_cor <- reduce(centroid_rel_cor, c)
tf <- Sys.time()
print(tf - ti)

message(paste("Mean correlation:", mean(centroid_rel_cor)))
message(paste("Minimum correlation:", min(centroid_rel_cor)))
print(which(centroid_rel_cor < 0.95))
print(which(centroid_rel_cor < 0.90))

hist(centroid_rel_cor)
