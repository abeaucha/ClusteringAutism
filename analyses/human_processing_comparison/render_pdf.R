#Packages
suppressPackageStartupMessages(library(tidyverse))

#Data sets
datasets <- c("POND", "SickKids")

#Parameter set 1
params_1 = list("es_method" = "normative-growth",
                "es_df" = 3,
                "es_combat" = TRUE,
                "es_combat_batch" = "Site-Scanner",
                "es_ncontrols" = NA,
                "cluster_nk_max" = 10,
                "cluster_metric" = "correlation",
                "cluster_K" = 30,
                "cluster_sigma" = 2.0,
                "cluster_t" = 20,
                "cluster_map_method" = "mean")
# params_1 = list("es_method" = "propensity-matching",
#                 "es_df" = NA,
#                 "es_combat" = FALSE,
#                 "es_combat_batch" = NA,
#                 "es_ncontrols" = 10,
#                 "cluster_nk_max" = 10,
#                 "cluster_metric" = "correlation",
#                 "cluster_K" = 10,
#                 "cluster_sigma" = 2.0,
#                 "cluster_t" = 20,
#                 "cluster_map_method" = "mean")

#Parameter set 2
# params_2 = list("es_method" = "normative-growth",
#                 "es_df" = 3,
#                 "es_combat" = TRUE,
#                 "es_combat_batch" = "Site-Scanner",
#                 "es_ncontrols" = NA,
#                 "cluster_nk_max" = 10,
#                 "cluster_metric" = "correlation",
#                 "cluster_K" = 20,
#                 "cluster_sigma" = 2.0,
#                 "cluster_t" = 20,
#                 "cluster_map_method" = "mean")
params_2 = list("es_method" = "propensity-matching",
                "es_df" = NA,
                "es_combat" = FALSE,
                "es_combat_batch" = NA,
                "es_ncontrols" = 20,
                "cluster_nk_max" = 10,
                "cluster_metric" = "correlation",
                "cluster_K" = 30,
                "cluster_sigma" = 2.0,
                "cluster_t" = 20,
                "cluster_map_method" = "mean")

#Path to pipeline directory
pipeline_dir <- "../../data/human/derivatives/"
pipeline_dir <- file.path(pipeline_dir, str_flatten(datasets, collapse = "_"), "")

#Path to complete parameter metadata
metadata <- file.path(pipeline_dir, "cluster_maps", "metadata.csv")
df_metadata <- read_csv(metadata, show_col_types = FALSE)
# View(df_metadata)

#ID for parameter set 1
id_params_1 <- inner_join(as_tibble(params_1),
                          df_metadata) %>% 
  pull(id)

if (length(id_params_1) == 0) {
  stop("Parameter set 1 not found.")
}

#ID for parameter set 2
id_params_2 <- inner_join(as_tibble(params_2),
                          df_metadata) %>% 
  pull(id)

if (length(id_params_2) == 0) {
  stop("Parameter set 2 not found.")
}

#Output file
outfile <- str_c("human_processing_comparison",
                 str_flatten(datasets, collapse = "-"),
                 id_params_1, 
                 id_params_2,
                 sep = "_")
outfile <- str_c(outfile, ".pdf")
outfile <- file.path("outputs", outfile)

#Render PDF
rmarkdown::render(input = "human_processing_comparison.Rmd",
                  params = list("datasets" = datasets,
                                "params_1" = params_1,
                                "params_2" = params_2),
                  output_file = outfile)
