
library(tidyverse)

list_permutations <- readRDS(file = "data/mouse/derivatives/v3/107/enrichment/StringDB_12.0_Bader_2025/950/cluster_pathway_enrichment_permutations_1_100.rds")

B <- 100
for (i in 1:100) {

file <- paste0("cluster_pathway_enrichment_seed_", i, ".csv")

# file_full <- "permutations_test/cluster_pathway_enrichment_seed_1.csv"
file_full <- file.path("data/mouse/derivatives/v3/107/enrichment/StringDB_12.0_Bader_2025/950/permutations/", file)

df_1 <- read_csv(file_full) 

df_1 <- df_1 |> 
  unite(col = "cluster_id", c("nk", "k"), sep = "-", remove = FALSE)

cluster_lvls <- df_1 |> 
  select(cluster_id, nk, k) |> 
  distinct() |> 
  arrange(nk, k) |> 
  pull(cluster_id)

mat_enrichment_csv <- df_1 |> 
  select(ID, cluster_id, E) |> 
  mutate(cluster_id = factor(cluster_id, levels = cluster_lvls)) |> 
  pivot_wider(id_cols = "ID", names_from = "cluster_id", values_from = "E") |> 
  arrange(ID) |> 
  column_to_rownames("ID") |> 
  as.matrix()


mat_enrichment_rds <- list_permutations[[i]]$result$enrichment

# dim(mat_enrichment_csv)
# dim(mat_enrichment_rds)

# sum(rownames(mat_enrichment_csv) != rownames(mat_enrichment_rds))

# sum(mat_enrichment_csv != mat_enrichment_rds, na.rm = TRUE)

# idx_which <- which(mat_enrichment_csv != mat_enrichment_rds)

sumcount <- sum(abs(mat_enrichment_csv - mat_enrichment_rds) > 1e-8, na.rm = TRUE)

if (sumcount > 0) {print(i)}

}
# ------------------------------

