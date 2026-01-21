library(tidyverse)

# Get list of mouse-human homologues from JAX
hom <- read_tsv(
  "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
) %>%
  select(
    DBClassKey = `DB Class Key`,
    Species = `Common Organism Name`,
    Gene = Symbol,
    EntrezID = `EntrezGene ID`
  ) %>%
  mutate(Species = ifelse(Species == "human", "human", "mouse"))

hom_human <- hom %>%
  filter(Species == "human") %>%
  select(-Species, -EntrezID) %>%
  rename(Human = Gene)

hom_mouse <- hom %>%
  filter(Species == "mouse") %>%
  select(-Species, -EntrezID) %>%
  rename(Mouse = Gene)

hom_joined <- inner_join(hom_human, hom_mouse, by = "DBClassKey")

expr_mouse <- "data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed.csv"
genes_mouse_amba <- as_tibble(data.table::fread(expr_mouse, header = TRUE)) %>%
  pull(Gene)

expr_human <- "data/human/expression/HumanExpressionMatrix_samples_pipeline_abagen.csv"
genes_human_ahba <- as_tibble(data.table::fread(expr_human, header = TRUE)) %>%
  pull(Gene)

# 3724 of 3958 with human homologues and in AMBA coronal
hom_joined %>%
  filter(Mouse %in% genes_mouse_amba) %>%
  nrow()

# 15,078 of 15,627 with mouse homologues and in AHBA
hom_joined %>%
  filter(Human %in% genes_human_ahba) %>%
  nrow()

hom_joined <- hom_joined %>%
  filter(Mouse %in% genes_mouse_amba, Human %in% genes_human_ahba)

# 3,252 homologous genes present in both AMBA and AHBA
nrow(hom_joined)
