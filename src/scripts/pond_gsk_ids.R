
library(tidyverse)

source("src/analysis.R")


file <- "data/mouse/derivatives/v3/107/clusters/resolution_0.2/clusters.csv"
df_mouse_clusters <- read_csv(file)

file <- "data/human/derivatives/v3/700/clusters/resolution_3.0/clusters.csv"
df_human_clusters <- read_csv(file)

df_similarity <- compute_similarity_significance(
  similarity = import_similarity(param_id = "375", 
                                 pipeline_dir = "data/cross_species/v3/"), 
  permutations = import_similarity_permutations(param_id = "375", 
                                                pipeline_dir = "data/cross_species/v3/")
)

df_similarity <- df_similarity %>% 
  mutate(match = pval < 0.05)

df_mouse_clusters_long <- df_mouse_clusters %>% 
  pivot_longer(cols = -ID, names_to = "nk_name", values_to = "k") %>% 
  mutate(nk = str_remove(nk_name, "nk")) %>% 
  unite(col = "cluster_id", nk, k, sep = "-", remove = FALSE) %>% 
  select(-nk_name)

df_human_clusters_long <- df_human_clusters %>% 
  pivot_longer(cols = -ID, names_to = "nk_name", values_to = "k") %>% 
  mutate(nk = str_remove(nk_name, "nk")) %>% 
  unite(col = "cluster_id", nk, k, sep = "-", remove = FALSE) %>% 
  select(-nk_name)


df_mouse_gsk <- df_mouse_clusters_long %>% 
  filter(str_detect(ID, "GSK3")) 

df_similarity %>% 
  semi_join(df_mouse_gsk, by = c("img2_cluster_id" = "cluster_id")) %>% 
  filter(match)

# Import updated mouse cluster labels
df_mouse_cluster_labels <- readxl::read_excel("figures/v3/resources/MICe_cluster_labels_colours.xlsx")

# Clean up mouse cluster labels
df_mouse_cluster_labels <- df_mouse_cluster_labels %>% 
  separate(col = "cluster_id", into = c("nk_old", "k_old"), remove = FALSE) %>% 
  separate(col = "cluster_id_new", into = c("nk_new", "k_new"), remove = FALSE) %>%
  mutate(nk_old = as.numeric(nk_old),
         k_old = as.numeric(k_old),
         nk_new = as.numeric(nk_new),
         k_new = as.numeric(k_new))

df_human_cluster_labels <- readxl::read_excel("figures/v3/resources/POND_cluster_labels_colours.xlsx")

# Clean up human cluster labels
df_human_cluster_labels <- df_human_cluster_labels %>% 
  separate(col = "cluster_id", into = c("nk_old", "k_old"), remove = FALSE) %>% 
  separate(col = "cluster_id_new", into = c("nk_new", "k_new"), remove = FALSE) %>%
  mutate(nk_old = as.numeric(nk_old),
         k_old = as.numeric(k_old),
         nk_new = as.numeric(nk_new),
         k_new = as.numeric(k_new))

df_gsk_match <- df_similarity %>% 
  semi_join(df_mouse_gsk, by = c("img2_cluster_id" = "cluster_id")) %>% 
  filter(match) %>% 
  select(img1_cluster_id, img2_cluster_id, pval) %>% 
  left_join(df_mouse_cluster_labels %>% 
              select(cluster_id, mouse_cluster_id = cluster_id_new, mouse_colour = colour), by = c("img2_cluster_id" = "cluster_id")) %>% 
  left_join(df_human_cluster_labels %>% 
              select(cluster_id, human_cluster_id = cluster_id_new, human_colour = colour), by = c("img1_cluster_id" = "cluster_id"))

df_gsk_match

df_gsk_human_match <- df_human_clusters_long %>% 
  semi_join(df_gsk_match, by = c("cluster_id" = "img1_cluster_id"))

df_gsk_human_match %>% 
  group_by(nk) %>% count()

df_gsk_human_match %>% 
  filter(nk == 6) 

file <- "data/human/registration/v3/subject_info/demographics.csv"
df_demographics <- read_csv(file)


  

file <- "data/human/registration/v3/subject_info/POND/POND_database_export_20251210.csv"
df_database <- read_csv(file)

df_site_codes <- df_database %>%
  select(subject_id, sub_id, site, redcap_event_name, redcap_repeat_instance) %>% 
  filter(str_detect(subject_id, "deprecated", negate = TRUE)) %>% 
  select(subject_id, site) %>% 
  filter(!is.na(site)) %>% 
  mutate(site_code = subject_id %>% 
           str_remove("PND01_") %>% 
           str_remove("_[0-9]+")) %>% 
  select(site, site_code) %>% 
  distinct() %>% 
    filter(!(site == 394 & site_code == "MCU")) %>% 
    mutate(site = as.character(site),
           site = str_pad(site, width = 3, side = "left", pad = "0"))

  
df_gsk_ids <- df_demographics %>% 
  semi_join(df_gsk_human_match %>% 
              filter(nk == 6), by = c("file" = "ID")) %>% 
  filter(Dataset == "POND", DX == "ASD") %>% 
  select(POND_ID = Subject_ID) %>% 
  mutate(POND_ID = str_remove(POND_ID, "sub-"),
         site = str_trunc(POND_ID, width = 3, side = "right", ellipsis = "")) %>% 
  left_join(df_site_codes, by = "site") %>% 
  select(POND_ID, site = site_code) 

outfile <- "POND_GSK_IDs.csv"
write_csv(df_gsk_ids, file = outfile)