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
pipeline_dir <- "data/mouse/derivatives/model_pathways_test/"
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

enrichment_E <- enrichment_NLQ <- enrichment_binary <- matrix(data = 0, 
                                                              nrow = n_models, 
                                                              ncol = n_modules)

gene_scores <- c(500, 700, 900, 950)
genes_unique <- unique(models[["gene"]])
list_enrichment <- vector(mode = "list", length = length(gene_scores))
names(list_enrichment) <- gene_scores
for (l in 1:length(list_enrichment)) {
  
  print(paste("Gene score", gene_scores[l]))
  
  for (gene in genes_unique) {
    
    print(gene)
    
    neighbourhood <- get_gene_neighbourhood(genes = gene, 
                                            score = gene_scores[l], 
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
      enrichment_E[i,] <- enrichment_gene$E
      enrichment_NLQ[i,] <- enrichment_gene$NLQ
      enrichment_binary[i,] <- as.numeric(enrichment_gene$adj.P.Val < 0.05)
    }
  }
  
  rownames(enrichment_E) <- str_remove(models$file, ".mnc")
  rownames(enrichment_NLQ) <- str_remove(models$file, ".mnc")
  rownames(enrichment_binary) <- str_remove(models$file, ".mnc")
  
  list_enrichment[[l]][["E"]] <- enrichment_E
  list_enrichment[[l]][["NLQ"]] <- enrichment_NLQ
  list_enrichment[[l]][["binary"]] <- enrichment_binary
  
}
  


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

list_clusters <- vector(mode = "list", length = length(list_enrichment))
list_metrics <- vector(mode = "list", length = length(list_enrichment))
names(list_clusters) <- names(list_enrichment)
names(list_metrics) <- names(list_enrichment)
for (l in 1:length(list_clusters)) {
  
  affinity <- similarity_network_new(x = list(es_matrix_abs, es_matrix_rel, list_enrichment[[l]][["E"]]),
                                                      metric = "euclidean")
  # affinity <- similarity_network_new(x = list(es_matrix_rel, list_enrichment[[l]][["E"]]),
  #                                    metric = "euclidean")
  list_metrics[[l]][["E"]] <- estimate_cluster_metrics(W = affinity, NUMC = 2:10)
  list_clusters[[l]][["E"]] <- create_clusters(W = affinity, nk = 10)
  
  affinity <- similarity_network_new(x = list(es_matrix_abs,es_matrix_rel, list_enrichment[[l]][["NLQ"]]),
                                                      metric = "euclidean")
  # affinity <- similarity_network_new(x = list(es_matrix_rel, list_enrichment[[l]][["NLQ"]]),
  #                                    metric = "euclidean")
  list_metrics[[l]][["NLQ"]] <- estimate_cluster_metrics(W = affinity, NUMC = 2:10)
  list_clusters[[l]][["NLQ"]] <- create_clusters(W = affinity, nk = 10)
  
  affinity <- similarity_network_new(x = list(es_matrix_abs, es_matrix_rel, list_enrichment[[l]][["binary"]]),
                                                      metric = "euclidean")
  # affinity <- similarity_network_new(x = list(es_matrix_rel, list_enrichment[[l]][["binary"]]),
  #                                    metric = "euclidean")
  list_metrics[[l]][["binary"]] <- estimate_cluster_metrics(W = affinity, NUMC = 2:10)
  list_clusters[[l]][["binary"]] <- create_clusters(W = affinity, nk = 10)
}


df_metrics <- list_metrics %>% 
  map(bind_rows, .id = "measure") %>% 
  bind_rows(.id = "score") %>% 
  mutate(measure = factor(measure, levels = c("E", "NLQ", "binary")))

plt_metrics <- ggplot(df_metrics, aes(x = nk, y = eigengap, group = measure, col = measure)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~score, ncol = 4) + 
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  labs(y = "Eigengap") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

outfile <- "cluster_metrics.pdf"
outfile <- file.path(pipeline_dir, outfile)
pdf(file = outfile, 
    width = unit(10, "in"),
    height = unit(5, "in"))
print(plt_metrics)
dev.off()

outfile <- "cluster_metrics_NLQ.pdf"
outfile <- file.path(pipeline_dir, outfile)
pdf(file = outfile, 
    width = unit(10, "in"),
    height = unit(5, "in"))
df_metrics %>% 
  filter(measure == "NLQ") %>% 
  ggplot(aes(x = nk, y = eigengap, group = factor(score), col = factor(score))) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  labs(y = "Eigengap",
       col = "StringDB score") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())
dev.off()

# Convert cluster information to long format
measures <- names(list_clusters[[1]])
for (score in gene_scores) {
  for (measure in measures) {
    
    score <- as.character(score)
    
    clusters_long <- list_clusters[[score]][[measure]] %>% 
      pivot_longer(cols = -ID, names_to = "nk_name", values_to = "k") %>% 
      mutate(nk = str_remove(nk_name, "nk"),
             nk = as.numeric(nk),
             nk = factor(nk),
             k = factor(k))
    
    # Mouse cluster Sankey plot
    sankey <- ggplot(clusters_long,
                     aes(x = nk,
                         stratum = k,
                         alluvium = ID,
                         fill = k, 
                         label = k)) + 
      geom_flow(stat = "alluvium", aes.flow = "forward") + 
      geom_stratum(alpha = 0.5) + 
      labs(x = 'Number of clusters',
           y = 'Number of models',
           title = paste("Gene score", score, "; Measure", measure)) + 
      theme_bw() + 
      theme(panel.grid.major.x = element_blank())
    
    outfile <- paste("sankey", score, measure, sep = "_")
    outfile <- paste0(outfile, ".pdf")
    outfile <- file.path(pipeline_dir, outfile)
    pdf(file = outfile, 
        width = unit(10, "in"),
        height = unit(5, "in"))
    print(sankey)
    dev.off()
    
  }
}


affinity_euclidean <- similarity_network_new(x = list(es_matrix_abs, es_matrix_rel),
                                             metric = "euclidean")

affinity_cor <- similarity_network_new(x = list(es_matrix_abs, es_matrix_rel),
                                      metric = "correlation")
df_metrics <- bind_rows(estimate_cluster_metrics(W = affinity_euclidean, NUMC = 2:10) %>% 
            mutate(metric = "euclidean"),
          estimate_cluster_metrics(W = affinity_cor, NUMC = 2:10) %>% 
            mutate(metric = "correlation"))

 
outfile <- "cluster_metrics_anat_only.pdf"
outfile <- file.path(pipeline_dir, outfile)
plt_metrics <- ggplot(df_metrics, aes(x = nk, y = eigengap, group = metric, col = metric)) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  labs(y = "Eigengap") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

pdf(file = outfile, 
    width = unit(10, "in"),
    height = unit(5, "in"))
print(plt_metrics)
dev.off()



affinity_euclidean <- similarity_network_new(x = list(effect_size_data_matrix_abs, effect_size_data_matrix_rel), 
                                             metric = "euclidean")

affinity_cor <- similarity_network_new(x = list(effect_size_data_matrix_abs, effect_size_data_matrix_rel),
                                       metric = "correlation")
df_metrics <- bind_rows(estimate_cluster_metrics(W = affinity_euclidean, NUMC = 2:10) %>% 
                          mutate(metric = "euclidean"),
                        estimate_cluster_metrics(W = affinity_cor, NUMC = 2:10) %>% 
                          mutate(metric = "correlation"))

outfile <- "cluster_metrics_anat_all_models.pdf"
outfile <- file.path(pipeline_dir, outfile)
plt_metrics <- ggplot(df_metrics, aes(x = nk, y = eigengap, group = metric, col = metric)) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  labs(y = "Eigengap") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

pdf(file = outfile, 
    width = unit(10, "in"),
    height = unit(5, "in"))
print(plt_metrics)
dev.off()


# Mixing correlation and euclidean ---- 


list_metrics <- vector(mode = "list", length = length(list_enrichment))
names(list_metrics) <- names(list_enrichment)
for (l in 1:length(list_metrics)) {
  for(measure in measures) {
    x <- list_enrichment[[l]][[measure]]
    d1 <- 1-cor(t(es_matrix_abs))
    d2 <- 1-cor(t(es_matrix_rel))
    d3 <- (dist2(as.matrix(x), as.matrix(x)))^(1/2)
    # W_list <- map(list(d1, d2, d3), affinityMatrix, K = 10, sigma = 0.5)
    W_list <- map(list(d2, d3), affinityMatrix, K = 10, sigma = 0.5)
    W <- SNF(W_list, K = 10, t = 20)
    list_metrics[[l]][[measure]] <- estimate_cluster_metrics(W = W, NUMC = 2:10)
    
  }
}

df_metrics <- list_metrics %>% 
  map(bind_rows, .id = "measure") %>% 
  bind_rows(.id = "score") %>% 
  mutate(measure = factor(measure, levels = c("E", "NLQ", "binary")))

outfile <- "cluster_metrics_mixed_rel_only.pdf"
outfile <- file.path(pipeline_dir, outfile)
plt_metrics <- ggplot(df_metrics, aes(x = nk, y = eigengap, group = measure, col = measure)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~score, ncol = 4) + 
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  labs(y = "Eigengap") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

pdf(file = outfile,
    width = unit(10, "in"),
    height = unit(5, "in"))
print(plt_metrics)
dev.off()


# ----

list_clusters$`950`$NLQ

mat <- list_enrichment$`950`$NLQ
sum(!(rowSums(mat) > 0))

length(rowSums(mat))


