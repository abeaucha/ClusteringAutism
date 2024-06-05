# ----------------------------------------------------------------------------
# template.R
# Author: Antoine Beauchamp
# Created:
#
# Brief description
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(magrittr))
# suppressPackageStartupMessages(library(glue))
# suppressPackageStartupMessages(library(tmod))


# Command line arguments -----------------------------------------------------

# option_list <- list(
#   make_option("--arg",
#               type = "character",
#               help = paste("Help message")),
#   make_option("--verbose",
#               type = "character",
#               default = "true",
#               help = paste("Verbosity option. [default %default]"))
# ) 


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")


# Functions ------------------------------------------------------------------

source(file.path(SRCPATH, "utils.R"))


#' Function description
#'
#' @param x (class) Parameter description.
#'
#' @return (class) Return description.
template <- function(x){
  return()
}


# Main -----------------------------------------------------------------------

#Parse command line args
# args <- parse_args(OptionParser(option_list = option_list))
# arg <- args[["arg"]]
# verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# if (verbose) {message("")}

# Directory containing mouse cluster files
Cluster_Dir <- "/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters_Paper/"
cluster_dir <- file.path(PROJECTPATH, "data/mouse/derivatives/v3/107/clusters/resolution_0.2/")


# Command line args? 
pipeline_dir <- file.path(PROJECTPATH, "data/mouse/derivatives/v3/")
params_id <- 107
models <- "model_names.csv"

# Get pipeline parameters
metadata <- file.path(pipeline_dir, "metadata.csv")
params <- fetch_params_metadata(metadata = metadata, id = params_id)

# Update pipeline directory with parameter set ID
pipeline_dir <- file.path(pipeline_dir, params_id)

# Get cluster resolution
cluster_resolution <- params[["cluster_resolution"]]
cluster_resolution <- sprintf("%.1f", cluster_resolution)

# Import model names dictionary
models <- file.path(pipeline_dir, models)
models <- read_csv(models, show_col_types = FALSE)

# Update cluster directory
cluster_dir <- file.path(pipeline_dir, "clusters", 
                         paste0("resolution_", cluster_resolution))

clusters <- file.path(cluster_dir, "clusters.csv")
clusters <- read_csv(clusters, show_col_types = FALSE) %>% 
  rename(file = ID) %>% 
  left_join(models, by = "file") %>% 
  select(-file) %>% 
  rename(Model = ID)

# Enrichment output directory
enrichment_dir <- file.path(PROJECTPATH, "data/mouse/enrichment/")
if (!dir.exists(enrichment_dir)) {
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
}

clusters %>% 
  mutate(Gene = Model %>%
           str_split("\\(") %>% 
           map_chr(.f = function(x){x[[1]]}))



cl <- read.csv(glue("{Cluster_Dir}/Clusters.csv"))

cl <- cl %>% dplyr::mutate(Model=X)%>% dplyr::select(-X)

test1 <- clusters %>% 
  pull(Model) %>% 
  str_split("\\(") %>% 
  map_chr(.f = function(x){x[[1]]})

test2 <- sapply(strsplit(clusters$Model, "\\("), "[[", 1)

sapply(strsplit(as.character(cl$Model),"\\("), "[[", 1)

# Get gene names
cl$Gene <- sapply(strsplit(as.character(cl$Model),"\\("), "[[", 1)
cl$Gene[cl$Gene=="Andr"] <- "Ar"
cl$Gene[cl$Gene=="Caspr2"] <- "Cntnap2"
cl$Gene[cl$Gene=="Dat"] <- "Slc6a3"
cl$Gene[cl$Gene=="Mor"] <- "Oprm1"
cl$Gene[cl$Gene=="Nl1"] <- "Nlgn1"
cl$Gene[cl$Gene=="Nl3"] <- "Nlgn3"
cl$Gene[cl$Gene=="Nrxn1a"] <- "Nrxn1"
cl$Gene[cl$Gene=="Sert"] <- "Slc6a4"
cl$Gene[cl$Gene=="Pcdh"] <- "Pcdhga3"
cl$Gene[cl$Gene=="Chd7;En1Cre"] <- "Chd7"
cl$Gene[cl$Gene=="Ube3a.2"] <- "Ube3a"
cl$Gene[cl$Gene=="FusDelta14"] <- "Fus"
cl$Gene[cl$Gene=="Nr1a"] <- "Nmdar1"
cl$Gene[cl$Model=="itsn1(+/+);itsn2(-/-)"]<- "itsn2"
cl$Gene[cl$Model=="Snf2H(+/+);Snf2L(-/-);emxcre"]<-"Snf2l"
cl$Gene[cl$Gene=="Snf2L"] <- "Smarca1"
cl$Gene[cl$Gene=="Snf2H"] <- "Smarca5"
cl$Gene[cl$Model=="Gsk3(a)"]<- "Gsk3A"
cl$Gene[cl$Model=="Gsk3(B)"]<- "Gsk3B"

# Remove CNVs, chromosomal, and behavioural models
valid_indices <- (! (cl$Gene %in% c("15q11-13", "16p11.2", "22q11.2", "XO", "Btbr", "Balbc", "MAR","15q25","TCDD","VPA","BtbrTT")))
cl_filt <- cl[valid_indices, ]
return(cl_filt)
}

