library(tidyverse)


# Version 3 comparison ---- 

# Absolute
sim_old_v3_abs <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_old_coords/984/similarity/similarity.csv"
sim_old_v3_abs <- read_csv(sim_old_v3_abs, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "absolute"),
         str_detect(img2, "absolute")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_new_v3_abs <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_new_coords/984/similarity/similarity.csv"
sim_new_v3_abs <- read_csv(sim_new_v3_abs, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "absolute"),
         str_detect(img2, "absolute")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_old_v3_abs
sim_new_v3_abs


# Relative 
sim_old_v3_rel <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_old_coords/984/similarity/similarity.csv"
sim_old_v3_rel <- read_csv(sim_old_v3_rel, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "relative"),
         str_detect(img2, "relative")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_new_v3_rel <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_new_coords/984/similarity/similarity.csv"
sim_new_v3_rel <- read_csv(sim_new_v3_rel, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "relative"),
         str_detect(img2, "relative")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_old_v3_rel
sim_new_v3_rel




# Version 3 HBN -----

# Absolute
sim_old_v3_abs <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_old_coords/394/similarity/similarity.csv"
sim_old_v3_abs <- read_csv(sim_old_v3_abs, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "absolute"),
         str_detect(img2, "absolute")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_new_v3_abs <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_new_coords/394/similarity/similarity.csv"
sim_new_v3_abs <- read_csv(sim_new_v3_abs, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "absolute"),
         str_detect(img2, "absolute")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_old_v3_abs
sim_new_v3_abs


# Relative 
sim_old_v3_rel <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_old_coords/394/similarity/similarity.csv"
sim_old_v3_rel <- read_csv(sim_old_v3_rel, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "relative"),
         str_detect(img2, "relative")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_new_v3_rel <- "/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/test/similarity_troubleshoot/v3_new_coords/394/similarity/similarity.csv"
sim_new_v3_rel <- read_csv(sim_new_v3_rel, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "relative"),
         str_detect(img2, "relative")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric())

sim_old_v3_rel
sim_new_v3_rel



# Microarray samples -----------------

library(RMINC)

source("src/analysis.R")

human_pipeline_dir <- "data/human/derivatives/v3/700/"
mouse_pipeline_dir <- "data/mouse/derivatives/v3/107/"


# Mouse anatomy
mouse_anat_file <- "data/mouse/atlas/DSURQE_CCFv3_average_200um.mnc"
mouse_anat <- mincGetVolume(mouse_anat_file)
mouse_anat_vol <- mincArray(mouse_anat)

# Human anatomy
human_anat_file <- "data/human/registration/v3/reference_files/model_0.8mm.mnc"
human_anat <- mincGetVolume(human_anat_file)
human_anat_vol <- mincArray(human_anat)

# Cropped human images along sagittal and transverse planes
human_slices_dim_1 <- 25:200
human_slices_dim_3 <- 25:220
human_anat_vol_cropped <- human_anat_vol[human_slices_dim_1,,human_slices_dim_3]


# Human mask
human_mask <- "data/human/registration/v2/reference_files/mask_0.8mm.mnc"

# Mouse mask
mouse_mask <- "data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc"

jacobians <- c("absolute", "relative")

# Human cluster map directories
human_centroid_dirs <- file.path(human_pipeline_dir, "centroids")
human_centroid_dirs <- file.path(human_centroid_dirs, "resolution_0.8")
human_centroid_dirs <- file.path(human_centroid_dirs, jacobians)
names(human_centroid_dirs) <- jacobians

# Mouse cluster map directories
mouse_centroid_dirs <- file.path(mouse_pipeline_dir, "centroids")
mouse_centroid_dirs <- file.path(mouse_centroid_dirs, "resolution_0.2")
mouse_centroid_dirs <- file.path(mouse_centroid_dirs, jacobians)
names(mouse_centroid_dirs) <- jacobians  

# Import mouse atlas labels
mouse_labels_file <- "../../data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc"
mouse_labels <- mincGetVolume(mouse_labels_file)

# Import mouse atlas definitions
mouse_defs_file <- "../../data/mouse/atlas/DSURQE_40micron_R_mapping_long.csv"
mouse_defs <- read_csv(mouse_defs_file, show_col_types = FALSE) %>% 
  select(name = Structure, label = Label)

# Import mouse neuroanatomical tree
mouse_tree_file <- "../../data/mouse/expression/MouseExpressionTree_DSURQE.RData"
load(mouse_tree_file)
mouse_tree <- Clone(treeMouseExpr)
rm(treeMouseExpr)



centroid_dirs_abs <- c("human" = human_centroid_dirs[[1]], 
                       "mouse" = mouse_centroid_dirs[[1]])

centroid_dirs_rel <- c("human" = human_centroid_dirs[[2]], 
                       "mouse" = mouse_centroid_dirs[[2]])

# Combine cluster map directories
centroid_dirs <- list(centroid_dirs_abs, centroid_dirs_rel)
names(centroid_dirs) <- jacobians

# Human and mouse trees
trees <- list("human" = human_tree,
              "mouse" = mouse_tree)

# Human and mouse labels
labels <- list("human" = human_labels,
               "mouse" = mouse_labels)

# Human and mouse definitions
defs <- list("human" = human_defs,
             "mouse" = mouse_defs)

# Human and mouse masks
masks <- list("human" = human_mask,
              "mouse" = mouse_mask)
