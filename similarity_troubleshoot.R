library(tidyverse)

# Absolute
sim_v2_abs <- "data/cross_species/v2/405/similarity/similarity_absolute.csv"
sim_v2_abs <- read_csv(sim_v2_abs, show_col_types = FALSE) %>% 
  mutate(human_nk = human_img %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         human_k = human_img %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         mouse_nk = mouse_img %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         mouse_k = mouse_img %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric()) %>% 
  filter(human_nk == 2, mouse_nk == 2) 


sim_v3_abs <- "data/cross_species/v3/394/similarity/similarity.csv"
sim_v3_abs <- read_csv(sim_v3_abs, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "absolute"),
         str_detect(img2, "absolute")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric()) %>% 
  filter(img_1_nk == 2, img_2_nk == 2)

sim_v2_abs$human_img
sim_v3_abs$img1

# Relative
sim_v2_rel <- "data/cross_species/v2/405/similarity/similarity_relative.csv"
sim_v2_rel <- read_csv(sim_v2_rel, show_col_types = FALSE) %>% 
  mutate(human_nk = human_img %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         human_k = human_img %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         mouse_nk = mouse_img %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         mouse_k = mouse_img %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric()) %>% 
  filter(human_nk == 2, mouse_nk == 2)


sim_v3_rel <- "data/cross_species/v3/394/similarity/similarity.csv"
sim_v3_rel <- read_csv(sim_v3_rel, show_col_types = FALSE) %>% 
  filter(str_detect(img1, "relative"),
         str_detect(img2, "relative")) %>% 
  mutate(img_1_nk = img1 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_1_k = img1 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_nk = img2 %>% basename() %>% str_extract("_nk_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric(),
         img_2_k = img2 %>% basename() %>% str_extract("_k_[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric()) %>% 
  filter(img_1_nk == 2, img_2_nk == 2)

sim_v2_abs
sim_v3_abs