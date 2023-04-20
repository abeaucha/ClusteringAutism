#Packages
suppressPackageStartupMessages(library(tidyverse))

#Paths
expr_dir = "data/human/expression"
metadata <- file.path(expr_dir, "SampleInformation_pipeline_abagen.csv")
annotations <- file.path(expr_dir, "all_samples_coarse_structure_assignments.csv")

#Import sample metadata
metadata <- read_csv(metadata, show_col_types = FALSE)

#Import annotations
annotations <- read_csv(annotations, show_col_types = FALSE) %>% 
  filter(Coordinates == "New") %>% 
  unite(col = "SampleID", structure_id, slab_num, well_id, 
        sep = "-", remove = FALSE) %>% 
  select(SampleID, ahba_coarse_structure, atlas_coarse_structure, Correct)

#Match order of annotated samples to metadata
ind_match <- match(metadata[["SampleID"]], annotations[["SampleID"]])
annotations <- annotations[ind_match,] %>% 
  select(sample_id = SampleID,
         ahba_coarse_structure,
         atlas_coarse_structure,
         keep = Correct)

#Export
annotations_out <- "AHBA_microarray_sample_annotations.csv"
annotations_out <- file.path(expr_dir, annotations_out)
write_csv(x = annotations, file = annotations_out)
