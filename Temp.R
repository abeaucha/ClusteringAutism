library(tidyverse)
library(RMINC)
library(MRIcrotome)

avg <- mincGetVolume(filename = "average_template_200um.mnc")
effect <- mincGetVolume(filename = "Group_1_Clusternum_2_ES_abs_200_mean.mnc")

avg_vol <- mincArray(avg)
effect_vol <- mincArray(effect)

dim(avg)

sliceSeries(nrow = 5, ncol = 5, begin = 10, end = 60) %>% 
  anatomy(avg_vol, low = 60, high = 150) %>% 
  overlay(effect_vol, low = 0.1, high = 0.5, symmetric = T) %>% 
  draw()

length(avg)
length(effect)
