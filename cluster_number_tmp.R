library(tidyverse)


cluster_file <- "data/human/derivatives/v2/700/clusters/resolution_3.0/clusters.csv"

df_clusters <- read_csv(cluster_file)

df_clusters <- df_clusters %>% 
  column_to_rownames("ID")


for (j in 1:ncol(df_clusters)) {

  col <- colnames(df_clusters)[j]
  
  x <- df_clusters %>% 
    group_by_at(.vars = all_of(col)) %>% 
    count() %>% 
    ungroup() 
  
  colnames(x) <- c("k", col)
  
  if (j == 1) {
    df_cluster_n <- x
  } else {
    df_cluster_n <- right_join(df_cluster_n, x, by = "k")
  }
}

df_cluster_n_long <- df_cluster_n %>% 
  pivot_longer(cols = -k, names_to = "nk", values_to = "n") %>% 
  mutate(nk = str_remove(nk, "nk"),
         nk = as.numeric(nk),
         k = factor(k),
         nk = factor(nk))

ggplot(df_cluster_n_long, aes(x = nk, y = fct_rev(k), fill = n)) +
  geom_tile(col = "black") +
  geom_text(aes(label = n)) + 
  scale_x_discrete(expand = expansion()) + 
  scale_y_discrete(expand = expansion()) + 
  scale_fill_distiller(na.value = "white", palette = "Reds", direction = 1) +
  labs(x = "nk",
       y = "k") + 
  theme_bw()

