library(tidyverse)
library(RMINC)
library(MRIcrotome)

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




