library(tidyverse)

source(file.path(SRCPATH, "enrichment.R"))

# Registration directory
pipeline_dir <- "data/mouse/derivatives/test/"
enrichment_dir <- "data/enrichment/"
registration_dir <- "data/mouse/registration/"

gene_score <- 950
stringdb_version <- "12.0"
bader_modules <- file.path(
  enrichment_dir,
  "Human_Reactome_October_01_2023_symbol.gmt"
)

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


background_set <- file.path(
  enrichment_dir,
  "sagittal_gene_table_normalized_filtered.csv"
)
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

genes_unique <- unique(models[["gene"]])

gene_scores <- c(seq(500, 900, by = 100), 950)

# df_neighbourhood <- expand_grid(gene = genes_unique,
#                                 score = gene_scores,
#                                 size = 0)
#
# for (i in 1:nrow(df_neighbourhood)) {
#
#   if (i %% 50 == 0){
#     print(i)
#   }
#
#   neighbourhood <- get_gene_neighbourhood(genes = df_neighbourhood[[i, "gene"]],
#                                           score = df_neighbourhood[[i, "score"]],
#                                           stringdb_version = stringdb_version)
#   df_neighbourhood[[i, "size"]] <- nrow(neighbourhood)
#
# }

df_neighbourhood
plt_neighbourhood <- ggplot(
  df_neighbourhood,
  aes(x = score, y = size, group = gene)
) +
  geom_line(alpha = 0.5) +
  geom_point(col = "grey30") +
  labs(x = "StringDB score", y = "Neighbourhood size") +
  scale_x_continuous(breaks = gene_scores) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

outfile <- "~/Downloads/neighbourhood_size.pdf"
pdf(file = outfile, width = unit(8, "in"), height = unit(6, "in"))
print(plt_neighbourhood)
dev.off()

score <- gene_scores[1]
size_seq <- c(0, 5, seq(10, 100, by = 10))
df_size_props <- expand_grid(score = gene_scores, size = size_seq, prop = 0)

for (i in 1:nrow(df_size_props)) {
  df_size_props[[i, "prop"]] <- df_neighbourhood %>%
    filter(
      score == df_size_props[[i, "score"]],
      size <= df_size_props[[i, "size"]]
    ) %>%
    nrow() /
    length(genes_unique)
}

plt_neighbourhood_prop <- ggplot(
  df_size_props,
  aes(x = size, y = prop, group = factor(score), col = factor(score))
) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = size_seq) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x = "Neighbourhood size",
    y = "Proportion of genes",
    col = "StringDB score"
  ) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

outfile <- "~/Downloads/neighbourhood_size_props.pdf"
pdf(file = outfile, width = unit(8, "in"), height = unit(6, "in"))
print(plt_neighbourhood_prop)
dev.off()

df_neighbourhood %>%
  filter(score == 500) %>%
  arrange(size)
