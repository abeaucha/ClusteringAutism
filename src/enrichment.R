# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

get_model_genes <- function(models) {
  
  # Genes to rename
  genes_rename <- tibble(gene = c("Andr", "Caspr2", "Dat", "Mor", "Nl1", "Nl3",
                                  "Nrxn1a", "Sert", "Pcdh", "Chd7;En1Cre",
                                  "Ube3a.2", "FusDelta14", "Nr1a", "Snf2L", "Snf2H"),
                         gene_new = c("Ar", "Cntnap2", "Slc6a3", "Oprm1", "Nlgn1",
                                      "Nlgn3", "Nrxn1", "Slc6a4", "Pcdhga3", "Chd7",
                                      "Ube3a", "Fus", "Nmdar1", "Smarca1", "Smarca5"))
  
  # Genes to exclude
  genes_exclude <- c("15q11-13", "16p11.2", "22q11.2", "XO", "Btbr", "Balbc", 
                     "MAR", "15q25", "TCDD", "VPA", "BtbrTT")
  
  # Import models 
  models <- read_csv(models, show_col_types = FALSE)
  
  # Create genes column
  models <- models %>% 
    mutate(gene = ID %>%
             str_split("\\(") %>% 
             map_chr(.f = function(x){x[[1]]}))
  
  models <- models %>% 
    left_join(genes_rename, by = "gene") %>% 
    mutate(gene_new = ifelse(is.na(gene_new), gene, gene_new),
           gene_new = ifelse(ID == "itsn1(+/+);itsn2(-/-)", "itsn2", gene_new),
           gene_new = ifelse(ID == "Snf2H(+/+);Snf2L(-/-);emxcre", "Snf2l", gene_new),
           gene_new = ifelse(ID == "Gsk3(a)", "Gsk3A", gene_new),
           gene_new = ifelse(ID == "Gsk3(B)", "Gsk3B", gene_new)) %>% 
    filter(!(gene_new %in% genes_exclude)) %>% 
    select(-gene) %>% 
    rename(gene = gene_new)  
  
  return(models)
  
}

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



get_gene_neighbourhood <- function(genes, score = 950, stringdb_version = "12.0", limit = 1000){
  
  # StringDB version URL
  if (stringdb_version == "11.5") {
    stringdb_url <- "https://version-11-5.string-db.org/"
  } else if (stringdb_version == "12.0") {
    stringdb_url <- "https://string-db.org/"
  } else {
    stop()
  }
  
  # StringDB API query to get gene IDs
  api_query_string_ids <- paste0(stringdb_url, "api/tsv/get_string_ids?identifiers=", 
                                 paste(genes, collapse = "%0D"), 
                                 "&species=10090&limit=1")
  
  # Fetch gene IDs
  string_ids <- read_tsv(api_query_string_ids, show_col_types = FALSE) %>% 
    select(queryIndex, stringID = stringId, preferredName) %>% 
    mutate(gene = genes)
  
  # StringDB API query for gene neighbourhoods
  api_query_interactions <- paste0(stringdb_url, "api/tsv/interaction_partners?identifiers=", 
                                   paste(string_ids[["stringID"]], collapse = "%0D"),
                                   "&species=10090&required_score=", score, "&limit=", limit)
  
  # Fetch gene neighbourhoods
  neighbourhood <- read_tsv(api_query_interactions, show_col_types = FALSE) %>% 
    rename(stringID_A = stringId_A, stringID_B = stringId_B, 
           gene_A = preferredName_A, gene_B = preferredName_B)  
  
  return(neighbourhood)
  
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


get_neighbourhood_enrichment <- function(target, background, modules = "Human_Reactome_October_01_2023_symbol.gmt") {
  
  modules <- importMsigDBGMT(modules)
  
  # Convert mouse gene names to human
  fg <- get_homologs(target, "mouse", ordered = F)[["human_genes"]]
  bg <- get_homologs(background, "mouse", ordered = F)[["human_genes"]]
  
  # Run hypergeometric test against pathway modules
  out <- HGtest(fg = fg, bg = bg, mset = modules)
  
  out <- out %>% 
    arrange(P.Value) %>% 
    mutate(rank = 1:nrow(.),
           adj.P.Val = p.adjust(P.Value, method = "fdr"),
           NLQ = -log10(adj.P.Val)) %>% 
    select(rank, ID, Title, P.Value, adj.P.Val, NLQ, b, B, n, N, E)
  
  return(out)
  
}


get_module_sizes <- function(modules = "Human_Reactome_October_01_2023_symbol.gmt") {
  modules <- importMsigDBGMT(modules)
  module_sizes <- map_dbl(modules[["MODULES2GENES"]], length) %>% 
    enframe(name = "ID", value = "B") %>% 
    inner_join(modules[["MODULES"]], by = "ID")
  return(module_sizes)
}