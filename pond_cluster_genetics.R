
library(tidyverse)


df_clusters <- read_csv("data/human/derivatives/v3/700/clusters/resolution_3.0/clusters.csv")

demographics <- read_csv("data/human/derivatives/v3/700/demographics.csv")

pond_cluster_ids <- demographics %>% 
  semi_join(df_clusters, by = c("file" = "ID")) %>% 
  filter(Dataset == "POND") %>% 
  pull(Subject_ID) %>% 
  str_remove("sub-")

pond_cluster_ids

df_genetics <- readxl::read_excel("data/human/registration/v3/subject_info/POND_variants.xlsx")
pond_genetics_ids <- df_genetics %>% 
  filter(str_detect(POND_ID, "POND", negate = TRUE)) %>% 
  pull(POND_ID)

sum(pond_cluster_ids %in% pond_genetics_ids)
