library(tidyverse)

# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "enrichment.R"))

# Pipeline directory
pipeline_dir <- "data/mouse/derivatives/test/"
if (!dir.exists(pipeline_dir)) {
  dir.create(pipeline_dir, recursive = TRUE)
}

# Registration directory
enrichment_dir <- "data/enrichment/"
registration_dir <- "data/mouse/registration/"

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
  filter(B >= B_threshold) %>% 
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
models <- read_csv(names_file_in, show_col_types = FALSE) %>% 
  select(file = Name, ID = NewName) %>% 
  mutate(file = paste0(file, ".mnc"))
write_csv(x = models, file = models_file_out)

models <- get_model_genes(models_file_out)

n_models <- nrow(models)

enrichment_NLQ <- matrix(data = 0, nrow = n_models, ncol = n_modules)
enrichment_E <- matrix(data = 0, nrow = n_models, ncol = n_modules)

genes_unique <- unique(models[["gene"]])
for (gene in genes_unique) {
  
  print(gene)
  
  neighbourhood <- get_gene_neighbourhood(genes = gene, 
                                          score = gene_score, 
                                          stringdb_version = stringdb_version)
  
  target_set <- unique(c(neighbourhood[["gene_A"]], neighbourhood[["gene_B"]]))
  
  enrichment_gene <- get_neighbourhood_enrichment(target = target_set,
                                                  background = background_set, 
                                                  modules = bader_modules) %>% 
    filter(ID %in% pathway_ids_keep) %>% 
    arrange(ID)
  
  idx_gene <- which(models[["gene"]] == gene)
  print(idx_gene)
  for (i in idx_gene) {
    enrichment_NLQ[i,] <- enrichment_gene$NLQ
    enrichment_E[i,] <- enrichment_gene$E
  }
}


hist(enrichment)
