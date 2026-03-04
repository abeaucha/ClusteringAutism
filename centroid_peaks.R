library(tidyverse)
library(RMINC)

PROJECTPATH <- Sys.getenv("PROJECTPATH")
SRCPATH <- Sys.getenv("SRCPATH")

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "analysis.R"))
source(file.path(SRCPATH, "enrichment.R"))

pipeline_dir <- file.path(PROJECTPATH, "data/cross_species/v3")

human_datasets <- c("POND", "HBN")

params_ids <- c(POND = "375", HBN = "861")



# Centroid peaks -------

mouse_centroid_dir <- "data/mouse/derivatives/v3/107/centroids/scanbase/resolution_0.2/"

nk_max <- 10
jacobians <- "relative"

df_peaks_all <- tibble()
for (nk in 2:nk_max) {
  for (k in 1:nk) {
    
    centroid_file <- paste0("centroid_nk_", nk, "_k_", k, ".mnc")
    centroid_file <- file.path(mouse_centroid_dir, jacobians, centroid_file)
    
    df_peaks <- mincFindPeaks(inputStats = centroid_file)
    
    df_peaks <- df_peaks %>% 
      as_tibble() %>% 
      mutate(nk = nk, k = k)
    
    df_peaks_all <- bind_rows(df_peaks_all, df_peaks)
    
  }
}  


# Atlas labels ------

labels <- "data/mouse/registration/reference_files/scanbase_second_level-nlin-3_labels_200um.mnc"
defs <- "data/mouse/registration/reference_files/DSURQE_40micron_R_mapping.csv"
defs <- read_csv(defs, show_col_types = FALSE)

defs_unilat <- defs %>% 
  filter(`right label` != `left label`)

defs_right <- defs_unilat %>% 
  select(structure = Structure, label = `right label`) %>% 
  mutate(structure = paste0("right ", structure))

defs_left <- defs_unilat %>% 
  select(structure = Structure, label = `left label`) %>% 
  mutate(structure = paste0("left ", structure))

defs_bilat <- defs %>% 
  filter(`right label` == `left label`) %>% 
  select(structure = Structure, label = `right label`) 

defs <- bind_rows(defs_right, defs_left, defs_bilat)

df_peaks_all <- df_peaks_all %>% 
  mutate(label = 0, structure = "")
for (i in 1:nrow(df_peaks_all)) {
  
  v1 <- df_peaks_all[[i, "x"]]
  v2 <- df_peaks_all[[i, "y"]]
  v3 <- df_peaks_all[[i, "z"]]
  label_i <- mincGetWorldVoxel(filenames = labels, v1 = v1, v2 = v2, v3 = v3)
  label_i <- round(label_i)
  
  struct <- defs %>% 
    filter(label == label_i) %>% 
    pull(structure)

  if (length(struct) == 0){
    struct <- NA
  }
    
  df_peaks_all[[i, "label"]] <- label_i[1]
  df_peaks_all[[i, "structure"]] <- struct
  
}


# Effect sizes -------

es_dir <- "data/mouse/derivatives/v3/107/effect_sizes/200/"
es_files <- list.files(es_dir, pattern = "Relative_200.mnc")
model_names <- es_files %>%
  str_remove("_ES_Relative_200.mnc")

es_mat <- matrix(data = 0, nrow = nrow(df_peaks_all), ncol = length(es_files))
colnames(es_mat) <- model_names

for (i in 1:nrow(df_peaks_all)) {

  if (i %% 200 == 0) {message(paste0("i = ", i))}
  
  v1 <- df_peaks_all[[i, "x"]]
  v2 <- df_peaks_all[[i, "y"]]
  v3 <- df_peaks_all[[i, "z"]]

  peaks <- mincGetWorldVoxel(filenames = file.path(es_dir, es_files),
                             v1 = v1, v2 = v2, v3 = v3)
  es_mat[i,] <- as.numeric(peaks)
}

df_peaks_all <- bind_cols(df_peaks_all, as_tibble(es_mat))


# Similarity ----


list_similarity <- vector(mode = "list", length = length(params_ids))
names(list_similarity) <- names(params_ids)
for (i in 1:length(list_similarity)) {
  
  df_similarity <- compute_similarity_significance(
    similarity = import_similarity(param_id = params_ids[i], 
                                   pipeline_dir = pipeline_dir), 
    permutations = import_similarity_permutations(param_id = params_ids[i], 
                                                  pipeline_dir = pipeline_dir)
  )
  
  cluster_ids <- df_similarity %>% 
    select(img1_cluster_id, img1_nk, img1_k) %>% 
    arrange(img1_nk, img1_k) %>% 
    distinct() %>% 
    pull(img1_cluster_id)
  
  df_cluster_grid <- expand_grid(img1_cluster_id = cluster_ids,
                                 img2_cluster_id = cluster_ids) 
  
  
  df_matches <- df_similarity %>% 
    mutate(match = pval < 0.05) %>% 
    right_join(df_cluster_grid,
               by = c("img1_cluster_id", "img2_cluster_id")) %>% 
    mutate(match = ifelse(is.na(match), FALSE, match))
  
  list_similarity[[i]] <- df_matches
  
}  

df_similarity_POND <- list_similarity[["POND"]]
colnames(df_similarity_POND) <- str_replace(colnames(df_similarity_POND), "img1", "POND")
colnames(df_similarity_POND) <- str_replace(colnames(df_similarity_POND), "img2", "MICe")
df_match_POND <- df_similarity_POND %>% 
  select(nk = MICe_nk, k = MICe_k, match) %>% 
  group_by(nk, k) %>% 
  summarise(POND_match = any(match), .groups = "drop")

outfile <- "MICe_POND_similarity.csv"
write_csv(x = df_similarity_POND, file = outfile)

df_similarity_HBN <- list_similarity[["HBN"]]
colnames(df_similarity_HBN) <- str_replace(colnames(df_similarity_HBN), "img1", "HBN")
colnames(df_similarity_HBN) <- str_replace(colnames(df_similarity_HBN), "img2", "MICe")
df_match_HBN <- df_similarity_HBN %>% 
  select(nk = MICe_nk, k = MICe_k, match) %>% 
  group_by(nk, k) %>% 
  summarise(HBN_match = any(match), .groups = "drop")

outfile <- "MICe_HBN_similarity.csv"
write_csv(x = df_similarity_HBN, file = outfile)


df_peaks_all <- df_peaks_all %>% 
  left_join(df_match_POND, by = c("nk", "k")) %>% 
  left_join(df_match_HBN, by = c("nk", "k")) 

df_peaks_all <- df_peaks_all %>%
  select(nk, k, POND_match, HBN_match, d1, d2, d3, x, y, z, label, structure, peak = value, `15q_pDp`:Wdfy3_HET)

outfile <- "MICe_peaks.csv"
write_csv(df_peaks_all, file = outfile)


# Enrichment ------

# Path to base enrichment dir
enrichment_dir <- file.path(PROJECTPATH, "data", "enrichment")

# Path to Reactome hierarchy file
reactome_hierarchy <- file.path(enrichment_dir, "reactome_hierarchy.csv")

# Import Reactome hierarchy
df_reactome_hierarchy <- read_csv(reactome_hierarchy, show_col_types = FALSE)

# Initialize tree levels for root pathways
df_reactome_hierarchy <- df_reactome_hierarchy %>% 
  mutate(level = ifelse(is.na(ParentID), 0, NA))

# Iterate lvl down the tree
nmissing <- sum(is.na(df_reactome_hierarchy$level))
lvl <- 0
while (nmissing > 0) {
  
  df_reactome_lvl <- df_reactome_hierarchy %>% 
    filter(level == lvl)
  
  # Iterate over nodes at given level
  for (i in 1:nrow(df_reactome_lvl)){
    
    node <- df_reactome_lvl[[i, "Name"]]
    
    df_reactome_hierarchy <- df_reactome_hierarchy %>% 
      mutate(level = ifelse(is.na(Parent), 0, 
                            ifelse(Parent == node, lvl+1, level)))
    
  }
  
  nmissing <- sum(is.na(df_reactome_hierarchy$level))
  lvl <- lvl + 1
  
}

# Clean up names
df_reactome_hierarchy <- df_reactome_hierarchy %>% 
  mutate(Name = str_replace_all(Name, "  ", " "),
         Name = str_replace_all(Name, "/", " "))

# Bader Reactome gene sets
module_file <- file.path(enrichment_dir, "Human_Reactome_June_01_2025_symbol.gmt")

# Import Bader sets
df_modules <- get_module_sizes(modules = module_file)

# Clean up IDs
df_modules <- df_modules %>% 
  rename(IDBader = ID) %>% 
  mutate(ID = IDBader %>% 
           str_split_i("%", i = -1) %>% 
           str_remove("R-HSA-") %>% 
           str_remove("\\.[0-9]+"),
         ID = paste0("R-HSA-", ID)) %>% 
  mutate(Title = str_replace_all(Title, "  ", " "))

# Include hierarchy information
df_reactome <- inner_join(df_reactome_hierarchy,
                          df_modules,
                          by = "ID")

# Pathways to keep
pathway_ids_keep <- df_reactome %>% 
  filter(B > 10) %>% 
  pull(IDBader)


mouse_enrichment_dir <- "data/mouse/derivatives/v3/107/enrichment/StringDB_12.0_Bader_2025/950/NeighbourhoodEnrichment/"


# Prefix for pathway data files
pathways_file_prefix <- "cluster_pathway_enrichment"

nk_max <- 10
stringdb_threshold <- 950

# Iterate over cluster solutions and import mouse pathway enrichment files
list_pathways <- vector(mode = "list", length = nk_max-1)
names(list_pathways) <- 2:nk_max
for (nk in 2:nk_max) {
  
  # Iterate over cluster number
  list_pathways[[nk-1]] <- vector(mode = "list", length = nk)
  for (k in 1:nk) {
    pathways_file <- paste(pathways_file_prefix, nk, k, stringdb_threshold, sep = "_")
    pathways_file <- paste0(pathways_file, ".csv")
    pathways_file <- file.path(mouse_enrichment_dir, pathways_file)
    list_pathways[[nk-1]][[k]] <- read_csv(pathways_file, show_col_types = FALSE)  %>% 
      semi_join(df_reactome, by = c("ID" = "IDBader")) %>% 
      left_join(df_reactome %>% 
                  select(ID = IDBader, level, Parent, Root), by = "ID") %>% 
      filter(ID %in% pathway_ids_keep) %>% 
      filter(!(Root == "Signal Transduction" & Title == "Integrin signaling"))
  }
  
  # Combine clusters per solution
  list_pathways[[nk-1]] <- bind_rows(list_pathways[[nk-1]]) %>% 
    mutate(nk = nk) 
  
}

# Reduce all pathway data frames into one
df_pathways_all <- bind_rows(list_pathways)

df_pathways_all <- df_pathways_all %>% 
  rename(pathway = Title) %>% 
  mutate(NLQ = -log10(adj.P.Val)) %>% 
  mutate(pathway = ifelse(pathway == "Signaling Pathways", "Signal Transduction", pathway)) %>% 
  unite(col = "cluster_id", nk, k, sep = "-", remove = FALSE)

df_pathways_all <- df_pathways_all %>%
  filter(Root != "Disease") %>% 
  filter(!(Root == "Reproduction" & pathway == "Meiosis")) %>% 
  filter(!(Root == "DNA Replication" & pathway == "Synthesis of DNA"))  

outfile <- "MICe_pathways.csv"
write_csv(df_pathways_all, file = outfile)
