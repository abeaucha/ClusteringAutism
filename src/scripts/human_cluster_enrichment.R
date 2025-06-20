# ----------------------------------------------------------------------------
# human_cluster_enrichment.R
# Authors: Antoine Beauchamp
# Created: November 15th, 2023


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(tmod))


# Functions ------------------------------------------------------------------

# NOTE: These functions were copied from 
# /projects/jacob/ClusteringAutism_125Models_Mar2020/Code/Clustering_Functions.R

get_homologs <- function(genes, species, ordered=T) {
  if (!exists("hom", envir = globalenv())) {
    hom <<- read_tsv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
  }
  species <- switch (species,
                     mouse={"mouse, laboratory"},
                     human={"human"}
  )
  
  if (ordered) {
    genes_human <- genes_mouse <- c()
    prog <- txtProgressBar(max=length(genes), style = 3)
    for (i in 1:length(genes)) {
      g <- genes[i]
      homologene_ids <- hom$`DB Class Key`[which((hom$Symbol==g) & hom$`Common Organism Name`==species)]
      genes_mouse <- c(genes_mouse, hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="mouse, laboratory"))])
      genes_human <- c(genes_human, hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="human"))])
      setTxtProgressBar(prog, i)
    }  
    close(prog)
  } else {
    # Faster, but unordered
    homologene_ids <- hom$`DB Class Key`[which((hom$Symbol %in% genes) & hom$`Common Organism Name`==species)]
    genes_mouse <- hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="mouse, laboratory"))]
    genes_human <- hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="human"))]
  }
  
  return(list(input_genes=genes, mouse_genes=genes_mouse, human_genes=genes_human))  
}


importMsigDBGMT <- function(file) {
  # stop("This does not work at the present.")
  msig <- list()
  con <- file(file, open = "r")
  lines <- readLines(con)
  close(con)
  ids <- gsub("\t.*", "", lines)
  desc <- gsub("^[^\t]*\t([^\t]*)\t.*", "\\1", lines)
  genes <- gsub("^[^\t]*\t[^\t]*\t(.*)", "\\1", lines)
  msig$MODULES <- data.frame(ID = ids, Title = desc, stringsAsFactors = FALSE)
  if (any(duplicated(msig$MODULES$ID))) {
    warning("Duplicated IDs found; automatic IDs will be generated")
    msig$MODULES$oldID <- msig$MODULES$ID
    msig$MODULES$ID <- make.unique(as.character(msig$MODULES$ID))
  }
  rownames(msig$MODULES) <- msig$MODULES[, "ID"]
  msig$MODULES2GENES <- strsplit(genes, "\t")
  names(msig$MODULES2GENES) <- ids
  msig$GENES <- data.frame(ID = unique(unlist(msig$MODULES2GENES)))
  # msig <- new("tmod", msig)
  msig
}


HGtest <- function (fg, bg, mset, qval = 0.05, cols = "Title", nodups = TRUE) {
  
  if (nodups) {
    fg <- fg[!duplicated(fg)]
    bg <- bg[!duplicated(bg)]
  }
  
  if (is.null(bg)) {
    warning("No genes in bg match any of the genes in the GENES")
  }
  
  if (is.null(fg)) {
    warning("No genes in fg match any of the genes in the GENES")
    return(NULL)
  }
  
  bg <- setdiff(bg, fg)
  if (length(bg) == 0){stop("All features from bg in fg.")}
  tot <- unique(c(fg, bg))
  n <- length(tot)
  k <- length(fg)
  
  n_modules <- length(mset[["MODULES2GENES"]])
  enrichment <- matrix(data = 0, nrow = n_modules, ncol = 6)
  for (i in 1:n_modules) {
    
    mg <- mset[["MODULES2GENES"]][[i]]
    q <- sum(fg %in% mg)
    m <- sum(tot %in% mg)
    if (m == 0) {
      E <- NA
    } else {
      E <- (q/k)/(m/n)
    }
    
    pv <- phyper(q - 1, m, n - m, k, lower.tail = FALSE)
    
    enrichment[i,] <- c(q, m, k, n, E, pv)
  }
  
  colnames(enrichment) <- c("b", "B", "n", "N", "E", "P.Value")
  enrichment <- as_tibble(enrichment)
  enrichment[["ID"]] <- mset[["MODULES"]][["ID"]]
  enrichment[["Title"]] <- mset[["MODULES"]][["Title"]]
  
  return(enrichment)
  
}


# Main -----------------------------------------------------------------------

# Parameters

# Human pipeline version
pipeline_version <- "v3"

# Human pipeline parameter set ID
params_id <- 700

# Cluster solution
# nk <- 6
nk_max <- 10

# Path to pipeline directory
pipeline_dir <- "data/human/derivatives/"

# Path to enrichment directory
enrichment_dir <- "data/human/enrichment/"

# Background set data file
background_set <- "sagittal_gene_table_normalized_filtered.csv"

# Gene score for StringDB
# gene_scores <- c(850, 900, 950)
gene_scores <- 950
# gene_scores <- seq(400, 950, by = 50)

# Version of StringDB
# stringdb_versions <- c("11.5", "12.0")
stringdb_versions <- "12.0"

# Version of Bader pathways
# bader_versions <- c(2020, 2023)
bader_versions <- 2023


# Path to pipeline version directory
pipeline_dir <- file.path(pipeline_dir, pipeline_version, params_id)
if (!dir.exists(pipeline_dir)) {stop("Pipeline version not found.")}

# Directory for enrichment outputs
output_dir <- file.path(pipeline_dir, "enrichment")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# File containing cluster variants
# variants <- "analyses/primary/outputs/human_cluster_genetics/"
variants <- file.path(pipeline_dir, "cluster_variants.csv")

# Import variants
df_variants <- read_csv(variants, show_col_types = FALSE)

df_variants <- df_variants %>% 
  filter(!is.na(gene))

# Identify genes for those patients with single gene mutations
# Version 2
# df_variants[["gene"]] <- c(NA, NA, 
#                            "Kcnq2", "Gria2",
#                            "Rai1", "Fmr1",
#                            "Tnpo3", "Stxbp1",
#                            NA, "Rbm10", 
#                            NA, NA, 
#                            "Fgf14", "Tsc1",
#                            "Meis2", "Rfx3",
#                            "Pten", "Rtl9",
#                            "Nlgn1")

# Version 3
# df_variants[["gene"]] <- c("Pten", "Rtl9", 
#                            "Nlgn1", NA,
#                            "Rbm10", NA,
#                            NA, "Fgf14",
#                            "Tsc1", "Meis2",
#                            NA, "Kcnq2",
#                            "Gria2", "Rfx3",
#                            "Rai1", "Fmr1",
#                            "Tnpo3", "Stxbp1",
#                            NA) 

# Import background gene set
background_set <- file.path(enrichment_dir, background_set)
background_set <- read_csv(background_set, show_col_types = FALSE) %>% 
  pull(msg.genes.acronym) %>% 
  as.character()

# Iterate over database versions and gene scores
params <- expand_grid(stringdb = stringdb_versions,
                      bader = bader_versions,
                      score = gene_scores)
for (i in 1:nrow(params)) {
  
  # Extract params
  bader_version <- params[[i, "bader"]]
  stringdb_version <- params[[i, "stringdb"]]
  gene_score <- params[[i, "score"]]
  
  message(
    paste(paste("StringDB version:", stringdb_version),
          paste("Bader version:", bader_version),
          paste("Gene score:", gene_score),
          sep = "\n")
  )
  
  # Create output directory for parameter set
  output_dir_i <- paste("StringDB", stringdb_version, "Bader", 
                        bader_version, sep = "_")
  output_dir_i <- file.path(output_dir, output_dir_i, gene_score)
  if (!dir.exists(output_dir_i)) {
    dir.create(output_dir_i, showWarnings = FALSE, recursive = TRUE)
  }
  
  # StringDB version URL
  if (stringdb_version == "11.5") {
    stringdb_url <- "https://version-11-5.string-db.org/"
  } else if (stringdb_version == "12.0") {
    stringdb_url <- "https://string-db.org/"
  } else {
    stop()
  }
  
  # Get Bader data version
  if (bader_version == 2020) {
    bader_modules <- file.path(enrichment_dir, "Human_Reactome_March_01_2020_symbol.gmt")
  } else if (bader_version == 2023) {
    bader_modules <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  
  } else {
    stop()
  }
  # bader_modules <- my_tmodImportMSigDB(bader_modules, format = "gmt")
  bader_modules <- importMsigDBGMT(bader_modules)
  
  # Iterate over clusters
  # list_cluster_interactions <- vector(mode = "list", length = nk)
  for (nk in 2:nk_max) {
    for (k in 1:nk) {
      
      message(
        paste("\tCluster", k, "of", nk)
      )
      
      message(
        "\t\tBuilding gene neighbourhoods..."
      )
      
      # Filter variants for cluster
      cols <- c(paste0("nk", nk), "gene")
      ind_k <- df_variants[[paste0("nk", nk)]] == k
      if (sum(ind_k) == 0) {
        message("\t\tNo genes in cluster.")
        next
      }
      df_variants_k <- df_variants[ind_k, cols]
      df_variants_k[["queryIndex"]] <- 0:(nrow(df_variants_k)-1)
      
      # StringDB API query to get gene IDs
      api_query_string_ids <- paste0(stringdb_url, "api/tsv/get_string_ids?identifiers=", 
                                     paste(df_variants_k[["gene"]], collapse = "%0D"), 
                                     "&species=10090&limit=1")
      
      df_string_ids <- read_tsv(api_query_string_ids, show_col_types = FALSE) %>% 
        select(queryIndex, stringID = stringId, preferredName) 
      
      # Join StringDB IDs to variants
      df_variants_k <- df_variants_k %>% 
        left_join(df_string_ids, by = "queryIndex")
      
      # StringDB API query for gene neighbourhoods
      api_query_interactions <- paste0(stringdb_url, "api/tsv/interaction_partners?identifiers=", 
                                       paste(df_variants_k[["stringID"]], collapse = "%0D"),
                                       "&species=10090&required_score=", gene_score)
      
      df_interactions_k <- read_tsv(api_query_interactions, show_col_types = FALSE) %>% 
        rename(stringID_A = stringId_A, stringID_B = stringId_B) %>% 
        left_join(df_variants_k[, c(paste0("nk", nk), "stringID")],
                  by = c("stringID_A" = "stringID"))  
      
      message(
        "\t\tRunning hypergeometric tests..."
      )
      
      # Target gene set for cluster k
      target_set <- unique(c(df_interactions_k[["preferredName_A"]],
                             df_interactions_k[["preferredName_B"]]))
      
      # Convert mouse gene names to human
      fg <- get_homologs(target_set, "mouse", ordered = F)[["human_genes"]]
      bg <- get_homologs(background_set, "mouse", ordered = F)[["human_genes"]]
      
      # Run hypergeometric test against pathway modules
      hgtest_out <- HGtest(fg = fg, bg = bg, mset = bader_modules)
      
      # Clean up output
      hgtest_out <- hgtest_out %>% 
        as_tibble() %>% 
        arrange(P.Value) %>% 
        mutate(rank = 1:nrow(.),
               adj.P.Val = p.adjust(P.Value, method = "fdr"),
               NLQ = -log(adj.P.Val),
               k = k) %>% 
        select(rank, Title, NLQ, P.Value, adj.P.Val, everything()) %>% 
        arrange(adj.P.Val)
      
      # Export the enrichment data to file
      outfile <- paste("cluster_pathway_enrichment", nk, k, gene_score, sep = "_")
      outfile <- paste0(outfile, ".csv")
      outfile <- file.path(output_dir_i, outfile)
      write_csv(x = hgtest_out, file = outfile)
      
    }
  }
}