
library(tidyverse)
library(RMINC)

sample_metadata <- read_csv('data/human/SampleInformation_pipeline_v1.csv')
template <- 'data/human/registration/reference_files/model_downsampled_3.0mm.mnc'

donors <- unique(sample_metadata[['Donor']])
for (i in 1:length(donors)) {
  url <- paste0("https://raw.githubusercontent.com/gdevenyi/AllenHumanGeneMNI/master/transformed-points/recombine/", 
                donors[i], 
                "_SampleAnnot.csv")
  df_tmp <- suppressMessages(read_csv(url))
  df_tmp <- df_tmp %>% 
    mutate(Donor = donors[i]) %>% 
    unite(SampleID, structure_id, slab_num, well_id, sep = '-', remove = FALSE) %>% 
    column_to_rownames('SampleID') %>% 
    select(contains('mni_nlin')) %>% 
    as.matrix()
  if (i == 1) {
    sample_coordinates_world <- df_tmp
  } else {
    sample_coordinates_world <- rbind(sample_coordinates_world,
                                      df_tmp)
  }
}

ind_match_metadata <- match(rownames(sample_coordinates_world),
                            sample_metadata[['SampleID']])
sample_coordinates_world <- sample_coordinates_world[ind_match_metadata,]

sample_coordinates_voxel <- mincConvertWorldMatrix(world_matrix = t(sample_coordinates_world),
                                                   file = template,
                                                   nearest_voxel = TRUE)
sample_coordinates_voxel <- t(sample_coordinates_voxel)
colnames(sample_coordinates_voxel) <- c('x', 'y', 'z')
sample_coordinates <- sample_coordinates_voxel %>% 
  as_tibble() %>% 
  unite(coords, x, y, z, sep = '-') %>% 
  pull(coords)

expr_dir <- 'data/'
# expr_dir <- 'data/MLP_outcomes/'
expr_glob <- 'HumanExpressionMatrix_samples_pipeline_v1_homologs_scaled.csv'
# expr_glob <- '*humantransform*'
expr_files <- Sys.glob(file.path(expr_dir, expr_glob))


cluster_dir <- 'data/human/clustering/cluster_masks/resolution_3.0/'
cluster_files <- list.files(cluster_dir, full.names = TRUE)

infiles <- expand_grid(clusterfile = cluster_files, 
                       exprfile = expr_files)

signatures <- tibble()
for (i in 1:nrow(infiles)){

  cluster <- mincArray(mincGetVolume(infiles[[i, 1]]))
  
  cluster_voxels <- which(cluster == 1, arr.ind = TRUE)
  
  cluster_coordinates <- cluster_voxels %>% 
    as_tibble() %>% 
    unite(coords, 1:3, sep = '-') %>% 
    pull(coords)
  
  expression <- data.table::fread(infiles[[i, 2]], header = TRUE) %>% 
    as_tibble()
  
  ind <- sample_coordinates %in% cluster_coordinates
  
  signature <- expression[ind,] %>% 
    summarise_all(.funs = mean)
  
  signatures <- bind_rows(signatures,
                          signature)  
}

tmp <- bind_cols(infiles, signatures)


