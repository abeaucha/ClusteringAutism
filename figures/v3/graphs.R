setwd("/Users/abeauchamp/Documents/Projects/ClusteringAutism/repository/main")


library(tidyverse)
# library(ggraph)
# library(tidygraph)
library(igraph)


file <- "data/mouse/derivatives/v3/107/clusters/resolution_0.2/clusters.csv"
clusters <- read_csv(file)


cluster_ids <- c()
for (nk in 2:10) {
  for (k in 1:nk) {
    cluster_ids <- c(cluster_ids, paste(nk, k, sep = "-"))
  }
}

df_cluster_hierarchy <- expand_grid(cluster_id = cluster_ids, 
                                    cluster_id_child = cluster_ids) %>% 
  separate(col = "cluster_id", into = c("nk", "k"), remove = FALSE) %>% 
  separate(col = "cluster_id_child", into = c("nk_child", "k_child"), remove = FALSE) %>% 
  mutate(nk = as.numeric(nk), k = as.numeric(k),
         nk_child = as.numeric(nk_child), k_child = as.numeric(k_child)) %>% 
  filter(nk_child == nk + 1)

df_cluster_graph <- df_cluster_hierarchy %>% 
  select(from = cluster_id, to = cluster_id_child)

g <- graph_from_data_frame(df_cluster_graph, directed = TRUE)

plot(g)
  
  
