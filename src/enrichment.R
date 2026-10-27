# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Environment variables ------------------------------------------------------

SRCPATH <- Sys.getenv("SRCPATH")


# Functions ------------------------------------------------------------------

#' Filter mouse models for enrichment
#'
#' @param models (character scalar) Path to the file (.csv) containing models
#'
#' @returns (data.frame, tbl) Filtered set of models.
get_model_genes <- function(models) {
  # Genes to rename
  genes_rename <- tibble(
    gene = c(
      "Andr",
      "Caspr2",
      "Dat",
      "Mor",
      "Nl1",
      "Nl3",
      "Nrxn1a",
      "Sert",
      "Pcdh",
      "Chd7;En1Cre",
      "Ube3a.2",
      "FusDelta14",
      "Nr1a",
      "Snf2L",
      "Snf2H"
    ),
    gene_new = c(
      "Ar",
      "Cntnap2",
      "Slc6a3",
      "Oprm1",
      "Nlgn1",
      "Nlgn3",
      "Nrxn1",
      "Slc6a4",
      "Pcdhga3",
      "Chd7",
      "Ube3a",
      "Fus",
      "Nmdar1",
      "Smarca1",
      "Smarca5"
    )
  )

  # Genes to exclude
  genes_exclude <- c(
    "15q11-13",
    "16p11.2",
    "22q11.2",
    "XO",
    "Btbr",
    "Balbc",
    "MAR",
    "15q25",
    "TCDD",
    "VPA",
    "BtbrTT"
  )

  # Import models
  models <- read_csv(models, show_col_types = FALSE)

  # Create genes column
  models <- models %>%
    mutate(
      gene = ID %>%
        str_split("\\(") %>%
        map_chr(.f = function(x) {
          x[[1]]
        })
    )

  models <- models %>%
    left_join(genes_rename, by = "gene") %>%
    mutate(
      gene_new = ifelse(is.na(gene_new), gene, gene_new),
      gene_new = ifelse(ID == "itsn1(+/+);itsn2(-/-)", "itsn2", gene_new),
      gene_new = ifelse(
        ID == "Snf2H(+/+);Snf2L(-/-);emxcre",
        "Snf2l",
        gene_new
      ),
      gene_new = ifelse(ID == "Gsk3(a)", "Gsk3A", gene_new),
      gene_new = ifelse(ID == "Gsk3(B)", "Gsk3B", gene_new)
    ) %>%
    filter(!(gene_new %in% genes_exclude)) %>%
    select(-gene) %>%
    rename(gene = gene_new)

  return(models)
}



#' Fetch gene homologs
#'
#' @param genes (character vector) Vector of gene names
#' @param species (character scalar) One of "mouse" or "human"
#'
#' @returns (tibble)
fetch_gene_homologs <- function(genes, species) {
  # Update species name
  species <- switch(
    species,
    mouse = {
      "mouse, laboratory"
    },
    human = {
      "human"
    }
  )

  # Get unique genes
  genes <- unique(genes)

  # Fetch homologs database from JAX
  hom <- read_tsv(
    "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
    show_col_types = FALSE
  )

  # Subset columns and get distinct entries
  hom <- hom %>%
    select(`DB Class Key`, Symbol, `Common Organism Name`) %>%
    distinct()

  # Build data frame from input genes
  hom_input <- tibble(Symbol = genes, `Common Organism Name` = species)

  # Get homolog information for input genes
  hom_output <- hom_input %>%
    left_join(
      select(hom, `DB Class Key`, Symbol, `Common Organism Name`),
      by = c("Symbol", "Common Organism Name")
    ) %>%
    select(`DB Class Key`, SymbolInput = Symbol)

  # Get homolog information for mouse
  hom_mouse <- hom %>%
    filter(`Common Organism Name` == "mouse, laboratory") %>%
    select(`DB Class Key`, SymbolMouse = Symbol)

  # Get homolog information for human
  hom_human <- hom %>%
    filter(`Common Organism Name` == "human") %>%
    select(`DB Class Key`, SymbolHuman = Symbol)

  # Join mouse and human homologs to input
  hom_output <- hom_output %>%
    left_join(hom_mouse, by = "DB Class Key") %>%
    left_join(hom_human, by = "DB Class Key")

  # Rename columns
  hom_output <- hom_output %>%
    rename(
      HomologeneID = `DB Class Key`,
      GeneInput = SymbolInput,
      GeneMouse = SymbolMouse,
      GeneHuman = SymbolHuman
    )

  return(hom_output)
}


#' Get the STRING interaction network for a set of genes
#'
#' @param genes (character vector) Genes whose interaction network to build.
#' @param score (numeric scalar) StringDB score for gene network confidence.
#' @param stringdb_version (character scalar) StringDB version.
#' @param limit (numeric scalar) Max number of records to return per gene from
#' the StringDB API.
#'
#' @returns (data.frame, tbl) Neighbourhood information for the input genes.
build_gene_interaction_network <- function(
  genes,
  score = 950,
  stringdb_version = "12.0",
  limit = 1000
) {
  # StringDB version URL
  if (stringdb_version == "11.5") {
    stringdb_url <- "https://version-11-5.string-db.org/"
  } else if (stringdb_version == "12.0") {
    stringdb_url <- "https://string-db.org/"
  } else {
    stop()
  }

  # StringDB API query to get gene IDs
  api_query_string_ids <- paste0(
    stringdb_url,
    "api/tsv/get_string_ids?identifiers=",
    paste(genes, collapse = "%0D"),
    "&species=10090&limit=1"
  )

  # Fetch gene IDs
  string_ids <- read_tsv(api_query_string_ids, show_col_types = FALSE) %>%
    select(queryIndex, stringID = stringId, preferredName) %>%
    mutate(gene = genes)

  # StringDB API query for gene neighbourhoods
  api_query_interactions <- paste0(
    stringdb_url,
    "api/tsv/interaction_partners?identifiers=",
    paste(string_ids[["stringID"]], collapse = "%0D"),
    "&species=10090&required_score=",
    score,
    "&limit=",
    limit
  )

  # Fetch gene neighbourhoods
  network <- read_tsv(api_query_interactions, show_col_types = FALSE) %>%
    rename(
      stringID_A = stringId_A,
      stringID_B = stringId_B,
      gene_A = preferredName_A,
      gene_B = preferredName_B
    ) %>%
    left_join(
      string_ids,
      by = c("stringID_A" = "stringID", "gene_A" = "preferredName")
    ) %>%
    rename(gene_in = gene)

  return(network)
}


#' Import GMT file
#'
#' Modified version of the function from tmod, which currently doesn't work
#' in the package.
#'
#' @param file (character scalar) File to import (.gmt)
#'
#' @returns (list) Modules and genes.
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


HGtest <- function(fg, bg, mset, qval = 0.05, cols = "Title", nodups = TRUE) {
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
  if (length(bg) == 0) {
    stop("All features from bg in fg.")
  }
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
      E <- (q / k) / (m / n)
    }

    pv <- phyper(q - 1, m, n - m, k, lower.tail = FALSE)

    enrichment[i, ] <- c(q, m, k, n, E, pv)
  }

  colnames(enrichment) <- c("b", "B", "n", "N", "E", "P.Value")
  enrichment <- as_tibble(enrichment)
  enrichment[["ID"]] <- mset[["MODULES"]][["ID"]]
  enrichment[["Title"]] <- mset[["MODULES"]][["Title"]]

  return(enrichment)
}



#' Evaluate the enrichment for a target set
#'
#' @param target (character vector) Target set of genes
#' @param background (character vector) Background set of genes
#' @param modules (character scalar) File (.gmt) containing the pathway modules
#'
#' @returns (data.frame, tbl) Enrichment data frame
gene_set_enrichment <- function(
  target,
  background,
  modules
) {
  # Run hypergeometric test against pathway modules
  enrichment <- HGtest(fg = target, bg = background, mset = modules)

  enrichment <- enrichment %>%
    arrange(P.Value) %>%
    mutate(
      rank = 1:nrow(.),
      adj.P.Val = p.adjust(P.Value, method = "fdr"),
      NLQ = -log10(adj.P.Val)
    ) %>%
    select(rank, ID, Title, P.Value, adj.P.Val, NLQ, b, B, n, N, E)

  return(enrichment)
}

get_module_sizes <- function(
  modules = "Human_Reactome_June_01_2025_symbol.gmt"
) {
  modules <- importMsigDBGMT(modules)
  module_sizes <- map_dbl(modules[["MODULES2GENES"]], length) %>%
    enframe(name = "ID", value = "B") %>%
    inner_join(modules[["MODULES"]], by = "ID")
  return(module_sizes)
}



compute_cluster_enrichment <- function(
    clusters,
    gene_network,
    background_set,
    modules,
    string_score = 950,
    nk_max = 10,
    permute = FALSE,
    seed = NULL
) {
  
  # If permutation, reshuffle cluster assignments
  if (permute) {
    
    if (is.null(seed)) {
      seed <- sample(1E5, size = 1)
    }
    
    set.seed(seed)
    idx_sample <- sample(
      seq_len(nrow(clusters)),
      size = nrow(clusters),
      replace = FALSE
    )
    
    clusters$model <- clusters$model[idx_sample]
    clusters$gene <- clusters$gene[idx_sample]
  }
  
  # Iterate over cluster solutions
  nk_range <- 2:nk_max
  list_enrichment <- vector(mode = "list", length = length(nk_range))
  for (l in 1:length(nk_range)) {
    nk <- nk_range[l]
    
    # Identify unique genes in each cluster
    cluster_genes_nk <- clusters |>
      select(cluster = paste0("nk", nk), gene) |>
      unite(dup_select, "gene", "cluster", remove = FALSE) |>
      filter(!duplicated(dup_select)) |>
      select(cluster, gene)
    
    # Iterate over clusters in solution
    list_enrichment[[l]] <- vector(mode = "list", length = nk)
    for (k in 1:nk) {
      message(paste0("On cluster ", nk, "-", k))
      
      # Pull genes for cluster
      cluster_genes_k <- cluster_genes_nk %>%
        filter(cluster == k) %>%
        pull(gene)
      
      # Extract gene neighbourhood for cluster
      neighbourhood_k <- gene_network %>%
        filter(gene_in %in% cluster_genes_k)
      
      # Construct target set for cluster
      target_set <- unique(c(
        neighbourhood_k[["gene_A"]],
        neighbourhood_k[["gene_B"]]
      ))
      
      # Evaluate pathway enrichment for cluster
      enrichment <- gene_set_enrichment(
        target = target_set,
        background = background_set,
        modules = modules
      )
      
      enrichment[["nk"]] <- nk
      enrichment[["k"]] <- k
      
      list_enrichment[[l]][[k]] <- enrichment
    }
    
    list_enrichment[[l]] <- bind_rows(list_enrichment[[l]])
  }
  
  df_enrichment <- bind_rows(list_enrichment) |> 
    unite(col = "cluster_id", c("nk", "k"), sep = "-", remove = FALSE)
  
  return(df_enrichment)
  
}


enrichment_to_matrices <- function(df) {
  
  # Get cluster ID levels
  cluster_lvls <- df |> 
    select(cluster_id, nk, k) |> 
    distinct() |> 
    arrange(nk, k) |> 
    pull(cluster_id)
  
  # Create enrichment matrix
  mat_enrichment <- df |> 
    select(ID, cluster_id, E) |> 
    mutate(cluster_id = factor(cluster_id, levels = cluster_lvls)) |> 
    pivot_wider(id_cols = "ID", names_from = "cluster_id", values_from = "E") |> 
    arrange(ID) |> 
    column_to_rownames("ID") |> 
    as.matrix()
  
  # Create p-value matrix
  mat_pvals <- df |> 
    select(ID, cluster_id, P.Value) |> 
    mutate(cluster_id = factor(cluster_id, levels = cluster_lvls)) |> 
    pivot_wider(id_cols = "ID", names_from = "cluster_id", values_from = "P.Value") |> 
    arrange(ID) |> 
    column_to_rownames("ID") |> 
    as.matrix()
  
  return(list("enrichment" = mat_enrichment,
              "p-values" = mat_pvals))
  
}

#' Run cluster-enrichment permutations in parallel
#'
#' @param clusters
#' @param gene_network
#' @param background_set
#' @param modules
#' @param string_score
#' @param nk_max
#' @param output_dir
#' @param B Positive integer number of permutations; seeds are seq_len(B).
#' @param n_workers Positive integer maximum number of worker R processes.
#'
#' @return A data frame with seed, success, and error columns, in seed order.
#'   Enrichment results are written to the seed-specific output directories.
run_cluster_enrichment_permutations <- function(
  clusters,
  gene_network,
  background_set,
  modules,
  string_score = 950,
  nk_max = 10,
  simplify = TRUE,
  output_dir = "data/enrichment/StringDB_12.0_Bader_2025",
  B = 1000L,
  seed_start = 1L,
  n_workers = 4L
) {
  # Reject invalid counts before starting any worker processes.
  stopifnot(
    length(B) == 1L,
    is.finite(B),
    B >= 1,
    B == floor(B),
    length(n_workers) == 1L,
    is.finite(n_workers),
    n_workers >= 1,
    n_workers == floor(n_workers)
  )

  # Start separate R processes, with no more workers than permutations.
  workers <- parallel::makePSOCKcluster(min(n_workers, B), outfile = "")
  # Release workers and reset foreach to sequential execution on exit,
  # including when setup or execution raises an error.
  on.exit(
    {
      parallel::stopCluster(workers)
      foreach::registerDoSEQ()
    },
    add = TRUE
  )
  # doSNOW supplies the per-result progress callback used below.
  doSNOW::registerDoSNOW(workers)

  # Only the parent process updates the bar; workers return their results.
  progress <- utils::txtProgressBar(min = 0, max = B, style = 3)
  on.exit(close(progress), add = TRUE)
  snow_options <- list(
    # Dispatch tasks dynamically so progress is reported for each seed.
    preschedule = FALSE,
    progress = function(n) utils::setTxtProgressBar(progress, n)
  )

  # Keep relative input/output paths consistent with the parent R session.
  parallel::clusterCall(workers, setwd, getwd())

  # Export shared data once per worker, including indirectly called helpers.
  shared_objects <- c(
    "compute_cluster_enrichment",
    "enrichment_to_matrices",
    "gene_set_enrichment",
    "HGtest",
    "clusters",
    "gene_network",
    "background_set",
    "modules",
    "string_score",
    "nk_max",
    "output_dir"
  )
  # Resolve names here and in enclosing environments. Each worker receives
  # its own copy, reused across its assigned permutations.
  parallel::clusterExport(workers, shared_objects, envir = environment())

  # Bind the operator locally without attaching foreach to the search path.
  `%dopar%` <- foreach::`%dopar%`
  foreach::foreach(
    b = seq.int(from = seed_start, length.out = B),
    # Collect one small status row per seed instead of enrichment tables.
    # .combine = rbind,
    .packages = "tidyverse",
    # These objects were already exported during worker initialization.
    .noexport = shared_objects,
    .options.snow = snow_options
  ) %dopar%
    {
      # Record an analysis failure without abandoning the remaining seeds.
      # Failed seeds also advance the bar, which counts completed attempts.
      tryCatch(
        {
          enrichment <- compute_cluster_enrichment(
            clusters = clusters,
            gene_network = gene_network,
            background_set = background_set,
            modules = modules,
            string_score = string_score,
            nk_max = nk_max,
            permute = TRUE,
            # The analysis sets its RNG seed and output directory from b.
            seed = b
          )
          
          if (simplify) {
            result <- enrichment_to_matrices(enrichment)
          } else {
            outfile <- "cluster_pathway_enrichment"
            outfile <- paste(outfile, "seed", b, sep = "_")
            outfile <- paste0(outfile, ".csv")
            outfile <- file.path(output_dir, outfile)
            write_csv(enrichment, file = outfile)
            result <- NULL
          }
          
          # result <- enrichment_to_matrices(enrichment)
          # 
          # outfile <- "cluster_pathway_enrichment"
          # outfile <- paste(outfile, "seed", b, sep = "_")
          # outfile <- paste0(outfile, ".csv")
          # outfile <- file.path(output_dir, outfile)
          # write_csv(enrichment, file = outfile)
          # 
          # # TEST - REMOVE
          # cluster_lvls <- enrichment |> 
          #   select(cluster_id, nk, k) |> 
          #   distinct() |> 
          #   arrange(nk, k) |> 
          #   pull(cluster_id)
          # 
          # mat_enrichment_test <- enrichment |> 
          #   select(ID, cluster_id, E) |> 
          #   mutate(cluster_id = factor(cluster_id, levels = cluster_lvls)) |> 
          #   pivot_wider(id_cols = "ID", names_from = "cluster_id", values_from = "E") |> 
          #   arrange(ID) |> 
          #   column_to_rownames("ID") |> 
          #   as.matrix()
          # 
          # testsum <- sum(mat_enrichment_test != result[["enrichment"]], na.rm = TRUE)
          # print(paste0("run_cluster_enrichment_permutations Testsum: ", testsum))
          
          list(
            seed = b,
            success = TRUE,
            result = result,
            test = enrichment,
            error = NA_character_
          )
          
        },
        error = function(e) {
            list(
              seed = b,
              success = FALSE,
              result = NULL,
              test = NULL,
              error = conditionMessage(e)
            )
        }
      )
    }
}
