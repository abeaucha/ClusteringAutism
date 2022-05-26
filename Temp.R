library(tidyverse)
library(RMINC)
library(MRIcrotome)
library(RColorBrewer)
library(pheatmap)

source('functions/buildSimilarityMatrix.R')
source('../Paper_TranscriptomicSimilarity/functions/tree_tools.R')
source('../Paper_TranscriptomicSimilarity/functions/metrics_tools.R')

# 
# cluster_ids <- read_csv('data/human/clustering/clusters_groups10_3.0mm.csv')
# 
# library(ggalluvial)
# 
# cluster_ids_long <- cluster_ids %>% 
#   pivot_longer(cols = -ID, names_to = 'nk_name', values_to = 'k') %>% 
#   mutate(nk = str_extract(nk_name, '[0-9]+'),
#          nk = as.numeric(nk),
#          nk = factor(nk, levels = 2:10),
#          k = factor(k, levels = 1:10))
# 
# ggplot(data = cluster_ids_long,
#        aes(x = nk,
#            stratum = k,
#            alluvium = ID,
#            fill = k,
#            label = k)) + 
#   geom_flow(stat = 'alluvium', colour = 'darkgrey', aes.flow = 'forward') + 
#   geom_stratum(alpha = 0.5) + 
#   theme_bw()
#   
# ---

imgdir <- 'data/human/clustering/cluster_maps/resolution_3.0/'

template <- mincArray(mincGetVolume('data/human/registration/reference_files/model_downsampled_3.0mm.mnc'))
mask <- mincArray(mincGetVolume('data/human/registration/reference_files/mask_downsampled_3.0mm.mnc'))

clusternum <- 4
groupnum <- 4
jacobians <- 'absolute'
method <- 'mean'

file_abs <- str_c('Group_', groupnum, '_Clusternum_', clusternum, '_ES_absolute_3_', method, '.mnc')
file_rel <- str_c('Group_', groupnum, '_Clusternum_', clusternum, '_ES_relative_3_', method, '.mnc')
file_abs <- file.path(imgdir, file_abs)
file_rel <- file.path(imgdir, file_rel)

cluster_map_abs <- mincArray(mincGetVolume(file_abs))
cluster_map_rel <- mincArray(mincGetVolume(file_rel))

sliceSeries(nrow = 5, ncol = 1, begin = 15, end = 60) %>% 
  anatomy(template, low = 3, high = 7) %>% 
  addtitle(str_c('nk: ', clusternum, ', k: ', groupnum)) %>% 
  sliceSeries() %>% anatomy() %>% 
  overlay(cluster_map_abs, symmetric = T, low = 0.1, high = 0.9) %>% 
  addtitle('Absolute') %>% 
  sliceSeries() %>% anatomy() %>% 
  overlay(cluster_map_rel, symmetric = T, low = 0.1, high = 0.9) %>% 
  addtitle('Relative') %>% 
  legend('mean z-score') %>% 
  draw()


# ---------

file <- file_human
species <- 'human'
jacobians <- 'absolute'
method <- 'mean'
threshold <- 0.1

process_signatures <- function(file, species, jacobians, method, threshold) {
  
  signatures <- suppressMessages(read_csv(file))
  
  if (species == 'mouse') {
    if (jacobians == 'absolute') {
      jacobians <- 'abs'
    } else {
      jacobians <- 'rel'
    }
  }
  
  signatures <- signatures %>% 
    filter(str_detect(clusterfile, jacobians),
           str_detect(clusterfile, method),
           str_detect(clusterfile, str_c('threshold', threshold))) %>% 
    mutate(nk = str_extract(clusterfile, 'Clusternum_[0-9]+'),
           nk = str_extract(nk, '[0-9]+'),
           nk = as.numeric(nk)) %>% 
    mutate(k = str_extract(clusterfile, 'Group_[0-9]+'),
           k = str_extract(k, '[0-9]+'),
           k = as.numeric(k)) %>% 
    unite(ID, nk, k, sep = '-', remove = FALSE) %>% 
    arrange(nk, k)
  
  cluster_expr <- signatures %>% 
    select(-nk, -k, -clusterfile, -exprfile) %>% 
    column_to_rownames('ID') %>% 
    as.matrix() %>% 
    t()
  
  cluster_ids <- signatures %>% 
    select(ID, nk, k) %>% 
    column_to_rownames('ID') %>% 
    mutate(nk = factor(nk),
           k = factor(k))
  
  params <- str_c('Jacobians: ', jacobians, '; Effect sizes: ', method, '; Threshold: ', threshold)
  
  return(list(expr = cluster_expr, ids = cluster_ids, params = params))
}


file_human <- 'data/human/human_cluster_signatures.csv'
file_mouse <- 'data/mouse/mouse_cluster_signatures.csv'

human_abs_mean_0.1 <- process_signatures(file_human, 
                                         species = 'human',
                                         jacobians = 'absolute',
                                         method = 'mean',
                                         threshold = 0.1)

mouse_abs_mean_0.1 <- process_signatures(file_mouse, 
                                         species = 'mouse',
                                         jacobians = 'absolute',
                                         method = 'mean',
                                         threshold = 0.1)

sim_cor_abs_mean_0.1 <- buildSimilarityMatrix(x1 = mouse_abs_mean_0.1$expr,
                                              x2 = human_abs_mean_0.1$expr, 
                                              method = "correlation")

colpalette <- tibble(vals = factor(1:10),
                     colour = brewer.pal(n = 10, name = 'Set3'))

human_annotations <- human_abs_mean_0.1$ids %>% 
  left_join(colpalette, by = c('nk' = 'vals')) %>% 
  rename(nk_col = colour) %>% 
  left_join(colpalette, by = c('k' = 'vals')) %>% 
  rename(k_col = colour) 

nk_colours <- human_annotations$nk_col
names(nk_colours) <- human_annotations$nk

k_colours <- human_annotations$k_col
names(k_colours) <- human_annotations$k

annotation_colours <- list(nk = nk_colours,
                           k = k_colours)

imgfile <- 'plots_temp/absolute_mean_0.1.png'
pheatmap(sim_cor_abs_mean_0.1, cluster_rows = F, cluster_cols = F,
         annotation_row = mouse_abs_mean_0.1$ids,
         annotation_col = human_abs_mean_0.1$ids,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         main = human_abs_mean_0.1$params,
         width = 10,
         height = 10,
         filename = imgfile)



#--

human_rel_mean_0.1 <- process_signatures(file_human, 
                                         species = 'human',
                                         jacobians = 'relative',
                                         method = 'mean',
                                         threshold = 0.1)

mouse_rel_mean_0.1 <- process_signatures(file_mouse, 
                                         species = 'mouse',
                                         jacobians = 'relative',
                                         method = 'mean',
                                         threshold = 0.1)

sim_cor_rel_mean_0.1 <- buildSimilarityMatrix(x1 = mouse_rel_mean_0.1$expr,
                                              x2 = human_rel_mean_0.1$expr, 
                                              method = "correlation")

colpalette <- tibble(vals = factor(1:10),
                     colour = brewer.pal(n = 10, name = 'Set3'))

human_annotations <- human_rel_mean_0.1$ids %>% 
  left_join(colpalette, by = c('nk' = 'vals')) %>% 
  rename(nk_col = colour) %>% 
  left_join(colpalette, by = c('k' = 'vals')) %>% 
  rename(k_col = colour) 

nk_colours <- human_annotations$nk_col
names(nk_colours) <- human_annotations$nk

k_colours <- human_annotations$k_col
names(k_colours) <- human_annotations$k

annotation_colours <- list(nk = nk_colours,
                           k = k_colours)

imgfile <- 'plots_temp/relative_mean_0.1.png'
pheatmap(sim_cor_rel_mean_0.1, cluster_rows = F, cluster_cols = F,
         annotation_row = mouse_rel_mean_0.1$ids,
         annotation_col = human_rel_mean_0.1$ids,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         main = human_rel_mean_0.1$params,
         width = 10,
         height = 10,
         filename = imgfile)

#-- 

human_abs_mean_0.5 <- process_signatures(file_human, 
                                         species = 'human',
                                         jacobians = 'absolute',
                                         method = 'mean',
                                         threshold = 0.5)

mouse_abs_mean_0.5 <- process_signatures(file_mouse, 
                                         species = 'mouse',
                                         jacobians = 'absolute',
                                         method = 'mean',
                                         threshold = 0.5)

sim_cor_abs_mean_0.5 <- buildSimilarityMatrix(x1 = mouse_abs_mean_0.5$expr,
                                              x2 = human_abs_mean_0.5$expr, 
                                              method = "correlation")

human_annotations <- human_abs_mean_0.5$ids %>% 
  left_join(colpalette, by = c('nk' = 'vals')) %>% 
  rename(nk_col = colour) %>% 
  left_join(colpalette, by = c('k' = 'vals')) %>% 
  rename(k_col = colour) 

nk_colours <- human_annotations$nk_col
names(nk_colours) <- human_annotations$nk

k_colours <- human_annotations$k_col
names(k_colours) <- human_annotations$k

annotation_colours <- list(nk = nk_colours,
                           k = k_colours)

imgfile <- 'plots_temp/absolute_mean_0.5.png'
pheatmap(sim_cor_abs_mean_0.5, cluster_rows = F, cluster_cols = F,
         annotation_row = mouse_abs_mean_0.5$ids,
         annotation_col = human_abs_mean_0.5$ids,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         na_col = 'black',
         main = human_abs_mean_0.5$params,
         width = 10,
         height = 10,
         filename = imgfile)

imgfile <- 'plots_temp/absolute_mean_0.5_hc.png'
pheatmap(sim_cor_abs_mean_0.5[,c(-1,-3,-9, -15)], cluster_rows = T, cluster_cols = T,
         annotation_row = mouse_abs_mean_0.5$ids,
         annotation_col = human_abs_mean_0.5$ids,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         na_col = 'black',
         main = human_abs_mean_0.5$params,
         width = 10,
         height = 10,
         filename = imgfile)


# --

human_rel_mean_0.5 <- process_signatures(file_human, 
                                         species = 'human',
                                         jacobians = 'relative',
                                         method = 'mean',
                                         threshold = 0.5)

mouse_rel_mean_0.5 <- process_signatures(file_mouse, 
                                         species = 'mouse',
                                         jacobians = 'relative',
                                         method = 'mean',
                                         threshold = 0.5)

sim_cor_rel_mean_0.5 <- buildSimilarityMatrix(x1 = mouse_rel_mean_0.5$expr,
                                              x2 = human_rel_mean_0.5$expr, 
                                              method = "correlation")

human_annotations <- human_rel_mean_0.5$ids %>% 
  left_join(colpalette, by = c('nk' = 'vals')) %>% 
  rename(nk_col = colour) %>% 
  left_join(colpalette, by = c('k' = 'vals')) %>% 
  rename(k_col = colour) 

nk_colours <- human_annotations$nk_col
names(nk_colours) <- human_annotations$nk

k_colours <- human_annotations$k_col
names(k_colours) <- human_annotations$k

annotation_colours <- list(nk = nk_colours,
                           k = k_colours)

imgfile <- 'plots_temp/relative_mean_0.5.png'
pheatmap(sim_cor_rel_mean_0.5, cluster_rows = F, cluster_cols = F,
         annotation_row = mouse_rel_mean_0.5$ids,
         annotation_col = human_rel_mean_0.5$ids,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         na_col = 'black',
         main = human_rel_mean_0.5$params,
         width = 10,
         height = 10,
         filename = imgfile)

imgfile <- 'plots_temp/relative_mean_0.5_hc.png'
pheatmap(sim_cor_rel_mean_0.5[,c(-1,-3)], cluster_rows = T, cluster_cols = T,
         annotation_row = mouse_rel_mean_0.5$ids,
         annotation_col = human_rel_mean_0.5$ids,
         annotation_colors = annotation_colours,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         na_col = 'black',
         main = human_rel_mean_0.5$params,
         width = 10,
         height = 10,
         filename = imgfile)

# sim_mat_euc <- buildSimilarityMatrix(x1 = mouse_cluster_expr,
#                                      x2 = human_cluster_expr, 
#                                      method = "euclidean")
# sim_mat_euc = 1 - (sim_mat_euc/max(sim_mat_euc))
# 
# pheatmap(sim_mat_euc, cluster_rows = F, cluster_cols = F,
#          annotation_row = mouse_cluster_ids,
#          annotation_col = human_cluster_ids,
#          main = 'Rows: Human clusters; Cols: Mouse clusters; Metric: 1 - normalized distance')


# ------

mouse_files <- list.files('data/mouse/cluster_maps/', full.names = T)
human_files <- list.files('data/human/clustering/cluster_maps/resolution_3.0/', full.names = T)

# human_img_1 <- human_files %>% 
#   str_subset('relative') %>% 
#   str_subset('mean') %>% 
#   str_subset('Group_6_Clusternum_7')

human_img_1 <- human_files %>% 
  str_subset('relative') %>% 
  str_subset('mean') %>% 
  str_subset('Group_6_Clusternum_7')

# human_img_2 <- human_files %>% 
#   str_subset('relative') %>% 
#   str_subset('mean') %>% 
#   str_subset('Group_7_Clusternum_8')

human_img_2 <- human_files %>% 
  str_subset('relative') %>% 
  str_subset('mean') %>% 
  str_subset('Group_2_Clusternum_9')

human_img_3 <- human_files %>%
  str_subset('relative') %>%
  str_subset('mean') %>%
  str_subset('Group_7_Clusternum_10')

human_template <- mincArray(mincGetVolume('data/human/registration/reference_files/model_downsampled_3.0mm.mnc'))
human_mask <- mincArray(mincGetVolume('data/human/registration/reference_files/mask_downsampled_3.0mm.mnc'))
human_vol_1 <- mincArray(mincGetVolume(human_img_1))
human_vol_2 <- mincArray(mincGetVolume(human_img_2))
human_vol_3 <- mincArray(mincGetVolume(human_img_3))

mouse_img_1 <- mouse_files %>%
  str_subset('rel') %>%
  str_subset('mean') %>%
  str_subset('Group_1_Clusternum_4')

# mouse_img_1 <- mouse_files %>% 
#   str_subset('rel') %>% 
#   str_subset('mean') %>% 
#   str_subset('Group_3_Clusternum_4')

mouse_img_2 <- mouse_files %>%
  str_subset('rel') %>%
  str_subset('mean') %>%
  str_subset('Group_1_Clusternum_3')

# mouse_img_2 <- mouse_files %>% 
#   str_subset('rel') %>% 
#   str_subset('mean') %>% 
#   str_subset('Group_4_Clusternum_6')

mouse_template <- mincArray(mincGetVolume('data/mouse/atlas/DSURQE_CCFv3_average_200um.mnc'))
mouse_vol_1 <- mincArray(mincGetVolume(mouse_img_1))
mouse_vol_2 <- mincArray(mincGetVolume(mouse_img_2))

sliceSeries(nrow = 8, ncol = 1, begin = 8, end = 54) %>% 
  anatomy(mouse_template, low = 700, high = 1400) %>% 
  overlay(mouse_vol_1, low = 0.5, high = 0.8, symmetric = T) %>% 
  addtitle('nk = 4, k = 3') %>% 
  sliceSeries() %>% anatomy() %>% 
  overlay(mouse_vol_2, low = 0.5, high = 0.8, symmetric = T) %>% 
  addtitle('nk = 6, k = 4') %>% 
  sliceSeries(nrow = 8, ncol = 1, begin = 15, end = 60) %>% 
  anatomy(human_template, low = 2, high = 7) %>% 
  overlay(human_vol_1, low = 0.5, high = 0.8, symmetric = T) %>% 
  addtitle('nk = 7, k = 6') %>% 
  sliceSeries() %>% anatomy() %>% 
  overlay(human_vol_2, low = 0.5, high = 0.8, symmetric = T) %>% 
  addtitle('nk = 9, k = 2') %>% 
  sliceSeries() %>% anatomy() %>% 
  overlay(human_vol_3, low = 0.5, high = 0.8, symmetric = T) %>% 
  addtitle('nk = 10, k = 7') %>% 
  legend("Effect size") %>% 
  draw()
