
library(tidyverse)


jacob_dir <- "/projects/jacob/ClusteringAutism_125Models_Mar2020/"
es_dirs_in <- file.path(jacob_dir, "Data", "Outputs", "EffectSizeMaps_Paper/", c("200", "50"))

registration_dir <- "data/mouse/registration/"
pipeline_dir <- "data/mouse/derivatives/v3/107/"
es_dirs_out <- file.path(pipeline_dir, "effect_sizes")

res_um <- c(200, 50)
res_mm <- res_um/1000

jacobians <- c("absolute", "relative")

for (i in 1:2) {
  for (j in 1:2) {
    
    res_str <- paste0("resolution_", res_mm[i])
    es_dir_out <- file.path(es_dirs_out, res_str, jacobians[j])
    dir.create(path = es_dir_out, recursive = TRUE, showWarnings = FALSE)
    
    pattern <- paste0(str_to_title(jacobians[j]), "*.mnc")
    es_files_in <- es_dirs_in[i] %>%
      list.files(pattern = "*.mnc") %>% 
      str_subset(str_to_title(jacobians[j])) 
    
    es_files_out <- es_files_in %>%
      str_remove(paste("_ES", str_to_title(jacobians[j]), res_um[i], sep = "_"))
    
    es_files_in <- file.path(es_dirs_in[i], es_files_in)
    es_files_out <- file.path(es_dir_out, es_files_out)
    
    for (f in 1:length(es_files_in)) {
      out <- file.copy(from = es_files_in[f], to = es_files_out[f])
    }
  }
}

names_file_in <- file.path(jacob_dir, "Data", "Raw", "Names_Paper.csv")
names_file_out <- "model_names.csv"
names_file_out <- file.path(pipeline_dir, names_file_out)
df_names <- read_csv(names_file_in, show_col_types = FALSE) %>% 
  select(file = Name, ID = NewName) %>% 
  mutate(file = paste0(file, ".mnc"))
write_csv(x = df_names, file = names_file_out)

clusters <- read_csv("data/mouse/derivatives/v3/107/clusters/resolution_0.2/clusters_init.csv")
colnames(clusters) <- c("ID", paste0("nk", 2:ncol(clusters)))
clusters <- clusters %>% 
  inner_join(df_names, by = "ID") %>% 
  select(-ID) %>% 
  rename(ID = file) %>% 
  select(ID, contains("nk"))
outfile <- "data/mouse/derivatives/v3/107/clusters/resolution_0.2/clusters.csv"
write_csv(x = clusters, file = outfile)
