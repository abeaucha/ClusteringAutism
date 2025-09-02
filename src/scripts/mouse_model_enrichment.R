library(tidyverse)
suppressPackageStartupMessages(library(ggalluvial))

# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "enrichment.R"))
source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "analysis.R"))
source(file.path(SRCPATH, "Clustering_Functions_AB.R"))

# Pipeline directory
pipeline_dir <- "data/mouse/derivatives/model_pathways/"
if (!dir.exists(pipeline_dir)) {
  dir.create(pipeline_dir, recursive = TRUE)
}

# Registration directory
enrichment_dir <- "data/enrichment/"
registration_dir <- "data/mouse/registration/"



# List of scanbase files
scanbase_files <- list(scans = "scanbase_40um - Scans_31Jan22.csv",
                       studies = "scanbase_40um - Studies_Feb2023.csv",
                       genotypes = "scanbase_40um - Genotypes_Feb2023.csv")

# Prepend path to scanbase files
scanbase_files <- map(.x = scanbase_files, 
                      .f = function(x) {
                        file.path(registration_dir, "resources", x)
                      })

# Image resolution 
resolution_um <- 200
resolution_mm <- resolution_um/1000

# List of paths
paths <- list(
  jacobians = file.path(registration_dir, "jacobians"),
  effect_sizes = file.path(pipeline_dir, "effect_sizes"),
  clusters = file.path(pipeline_dir, "clusters", paste0("resolution_", resolution_mm)),
  centroids = file.path(pipeline_dir, "centroids", paste0("resolution_", resolution_mm))
)

# Create output directories
for (path in paths[-1]) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}


# ANATOMY -----------

# Generate list of models
model_list <- MakeModelList(
  scanbase_scans_file = scanbase_files[["scans"]],
  scanbase_studies_file = scanbase_files[["studies"]],
  scanbase_genotypes_file = scanbase_files[["genotypes"]],
  scanbase_focus = "Autism",
  scanbase_sample_size_threshold = 6,
  base_directory = pipeline_dir
)

effect_size_data_matrix_abs <- MakeEffectSizeMaps(
  jdtype = "Absolute", model_list = model_list, 
  resolution = as.character(resolution_um),
  dir_determinants = paths[["jacobians"]],
  output_phenotype_dir = paths[["effect_sizes"]],
  base_directory = pipeline_dir, boot = "N"
)

effect_size_data_matrix_rel <- MakeEffectSizeMaps(
  jdtype = "Relative", model_list = model_list,
  resolution = as.character(resolution_um),
  dir_determinants = paths[["jacobians"]],
  output_phenotype_dir = paths[["effect_sizes"]],
  base_directory = pipeline_dir, boot = "N"
)

# Replace NaN with 0
effect_size_data_matrix_abs[is.nan(effect_size_data_matrix_abs)] <- 0
effect_size_data_matrix_rel[is.nan(effect_size_data_matrix_rel)] <- 0



# PATHWAYS ---------

gene_score <- 950
stringdb_version <- "12.0"
bader_modules <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  

modules_list <- importMsigDBGMT(bader_modules)
df_modules_size <- map_dbl(modules_list[["MODULES2GENES"]], length) %>% 
  enframe(name = "ID", value = "B") %>% 
  inner_join(modules_list[["MODULES"]], by = "ID")

B_threshold <- 10

# Pathways to keep
pathway_ids_keep <- df_modules_size %>% 
  filter(B >= B_threshold,
         Title != "Formation of ATP by chemiosmotic coupling") %>% 
  pull(ID)

n_modules <- length(pathway_ids_keep)

background_set <- file.path(enrichment_dir, "sagittal_gene_table_normalized_filtered.csv")
background_set <- read_csv(background_set, show_col_types = FALSE) %>% 
  pull(msg.genes.acronym) %>% 
  as.character()

# Copy model names file into pipeline directory
models_file_in <- file.path(registration_dir, "Names_Paper.csv")
models_file_out <- "model_names.csv"
models_file_out <- file.path(pipeline_dir, models_file_out)
models <- read_csv(models_file_in, show_col_types = FALSE) %>% 
  select(file = Name, ID = NewName) %>% 
  mutate(file = paste0(file, ".mnc"))
write_csv(x = models, file = models_file_out)

models <- get_model_genes(models_file_out)

# models <- models %>% 
#   filter(gene != "Ar")

n_models <- nrow(models)

enrichment <- matrix(data = 0, nrow = n_models, ncol = n_modules)

genes_unique <- unique(models[["gene"]])
for (gene in genes_unique) {
  
  print(gene)
  
  neighbourhood <- get_gene_neighbourhood(genes = gene, 
                                          score = gene_score, 
                                          stringdb_version = stringdb_version)
  
  if (nrow(neighbourhood) == 0){
    target_set <- gene
  } else {
    target_set <- unique(c(neighbourhood[["gene_A"]], neighbourhood[["gene_B"]]))
  }
  
  enrichment_gene <- get_neighbourhood_enrichment(target = target_set,
                                                  background = background_set, 
                                                  modules = bader_modules) %>% 
    filter(ID %in% pathway_ids_keep) %>% 
    arrange(ID)
  
  idx_gene <- which(models[["gene"]] == gene)
  for (i in idx_gene) {
    enrichment[i,] <- enrichment_gene$NLQ
  }
}

rownames(enrichment) <- str_remove(models$file, ".mnc")
  
  


# CLUSTERING ---------------

similarity_network_new <- function(x, metric = "correlation", K = 10,
                                   sigma = 0.5, t = 20, outfile = NULL){
  
  if (metric == "correlation") {
    d <- map(x, function(x) {1-cor(t(x))})
  } else if (metric == "euclidean") {
    d <- map(x, function(x){(dist2(as.matrix(x), as.matrix(x)))^(1/2)})
  } else {
    stop(paste("Argument metric must be one of {correlation, euclidean}:", metric))
  }
  
  W_list <- map(d, affinityMatrix, K = K, sigma = sigma)
  
  W <- SNF(W_list, K = K, t = t)
  
  if (!is.null(outfile)){
    data.table::fwrite(x = as_tibble(W), file = outfile)
  }
  
  return(W)
  
}


idx_single_genes <- rownames(effect_size_data_matrix_abs) %in% str_remove(models[["file"]], ".mnc")
es_matrix_abs <- effect_size_data_matrix_abs[idx_single_genes,]
es_matrix_rel <- effect_size_data_matrix_rel[idx_single_genes,]

outfile <- file.path(paths$clusters, "affinity.csv")
affinity <- similarity_network_new(x = list(es_matrix_abs, es_matrix_rel, enrichment),
                                   metric = "euclidean", outfile = outfile)

outfile <- file.path(paths$clusters, "clusters.csv")
clusters <- create_clusters(W = affinity, nk = 10, outfile = outfile)


mask <- file.path(paths[["jacobians"]], "scanbase_second_level-nlin-3_mask_200um.mnc")
jacobians <- c("absolute", "relative")
for (jtype in jacobians) {
  
  clusters_jtype <- clusters %>%
    as_tibble() %>%
    mutate(ID = file.path(paths[["effect_sizes"]], resolution_um, ID),
           ID = paste0(ID, "_ES_", str_to_title(jtype), "_", resolution_um, ".mnc")) %>%
    column_to_rownames("ID")
  
  outdir <- file.path(paths[["centroids"]], jtype)
  if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
  for (j in 1:ncol(clusters_jtype)) {
    compute_cluster_centroids(i = j, clusters = clusters_jtype, mask = mask,
                              method = "mean", outdir = outdir)
  }
}



heatmap_palette_cols <- RColorBrewer::brewer.pal(n = 9, name = "OrRd")
heatmap_palette <- colorRampPalette(colors = heatmap_palette_cols)(255)


nk <- 7

clusters_nk <- clusters %>% 
  select(ID, k = paste0("nk", nk)) %>% 
  arrange(k, ID)

idx_match <- match(clusters_nk$ID, rownames(enrichment))
mat <- enrichment[idx_match,]
mat[mat > 20] <- 20
mat[mat == 0] <- NA

annotation_row <- clusters_nk %>% 
  column_to_rownames("ID") %>% 
  mutate(k = factor(k))

outfile <- file.path(pipeline_dir, "enrichment_heatmap_nk7.pdf")
pheatmap::pheatmap(mat = mat, 
                   cluster_cols = FALSE, cluster_rows = FALSE, 
                   annotation_row = annotation_row,
                   color = heatmap_palette, border_color = "grey85", na_col = "grey85",
                   fontsize_row = 6, filename = outfile, width = unit(16, "in"), height = unit(10, "in"))
