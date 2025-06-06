file <- bader_modules

my_importMsigDBGMT <- function(file) {
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

my_tmodImportMSigDB <- function (file = NULL, format = "xml", organism = "Homo sapiens", 
                                 fields = c("STANDARD_NAME", "CATEGORY_CODE", "SUB_CATEGORY_CODE", 
                                            "EXACT_SOURCE", "EXTERNAL_DETAILS_URL")) 
{
  if (length(file) != 1) 
    stop("Incorrect file parameter")
  if (!file.exists(file)) 
    stop(sprintf("File %s does not exist", file))
  format <- match.arg(format, c("xml", "gmt"))
  msig <- switch(format, xml = .importMsigDBXML(file, fields, 
                                                organism), gmt = my_importMsigDBGMT(file))
  s <- msig$gs$Title
  msig$gs$Title <- paste0(toupper(substring(s, 1, 1)), tolower(substring(s, 
                                                                         2)))
  msig$gs$Title <- gsub("^Gse([0-9])", "GSE\\1", msig$gs$Title)
  msig$gs$Title <- gsub("_", " ", msig$gs$Title)
  msig$gs$B <- sapply(msig$gs2gv, length)
  msig
}


tmodHGtest(fg = fg, bg = bg, mset = bader_modules, 
           qval = 1.1, filter = FALSE, order.by = "pval")


# Target gene set for cluster k
target_set <- unique(c(df_interactions_k[["preferredName_A"]],
                       df_interactions_k[["preferredName_B"]]))

# Convert mouse gene names to human
fg <- get_homologs(target_set, "mouse", ordered = F)[["human_genes"]]
bg <- get_homologs(background_set, "mouse", ordered = F)[["human_genes"]]
mset <- bader_modules
modules <- NULL
qval = 1.1
filter = FALSE
cols = "Title"
order.by <- "pval"
nodups <- TRUE

my_tmodHGtest <- function (fg, bg, modules = NULL, qval = 0.05, order.by = "pval", 
          filter = FALSE, mset = "all", cols = "Title", nodups = TRUE) 
{
  # mset <- tmod:::.getmodules_gs(modules, mset)
  # fg <- .prep_list(fg, mset, filter = filter, nodups = nodups)
  # bg <- .prep_list(bg, mset, filter = filter, nodups = nodups)
  if (is.null(bg)) {
    warning("No genes in bg match any of the genes in the GENES")
  }
  if (is.null(fg)) {
    warning("No genes in fg match any of the genes in the GENES")
    return(NULL)
  }
  bg <- setdiff(bg, fg)
  if (length(bg) == 0) 
    stop("All features from bg in fg.")
  tot <- unique(c(fg, bg))
  n <- length(tot)
  k <- length(fg)
  mod.test <- function(n_id) {
    mg <- mset$gs2gv[[n_id]]
    q <- sum(fg %in% mg)
    m <- sum(tot %in% mg)
    if (m == 0) {
      E <- NA
    }
    else {
      E <- (q/k)/(m/n)
    }
    if (q == 0 || m == 0) 
      return(c(n_id = n_id, b = q, B = m, n = k, N = n, 
               E = E, P.Value = 1))
    pv <- phyper(q - 1, m, n - m, k, lower.tail = FALSE)
    c(n_id = n_id, b = q, B = m, n = k, N = n, E = E, P.Value = pv)
  }
  ret <- .tmodTest(mod.test, NULL, qval = qval, order.by = order.by, 
                   mset = mset, cols = cols)
  attr(ret, "effect_size") <- "E"
  attr(ret, "pvalue") <- "adj.P.Val"
  ret
}


function (modules = NULL, mset = "all", known.only = FALSE, skipcheck = FALSE) 
{
  if (is(mset, "tmod")) {
    warning("You are loading an obsolete version of tmod R object.\nThe class `tmod` has been retired.\nThe data will still work, but it will incur a penalty\non the computational time. Please use the `tmod2tmodGS`\nfunction to convert your object to the new tmodGS class.")
    mset <- tmod2tmodGS(mset)
  }
  if (is(mset, "list")) {
    mset <- .mset_sanity_check_gs(mset, modules)
    if (!is.null(modules)) {
      mset <- mset[modules]
    }
  }
  else {
    tmod <- .gettmod()
    mset <- match.arg(mset, c("all", unique(tmod$gs$SourceID)))
    if (mset != "all") {
      if (is.null(modules)) 
        modules <- tmod$gs$ID
      sel <- tmod$gs$SourceID[match(modules, tmod$gs$ID)] == 
        mset
      modules <- modules[sel]
    }
    if (!is.null(modules)) {
      modules <- modules[modules %in% tmod$gs$ID]
      mset <- tmod[modules]
    }
    else {
      mset <- tmod
    }
  }
  if (known.only && "Title" %in% colnames(mset$gs)) {
    mset <- mset[!is.na(mset$gs$Title) & !mset$gs$Title %in% 
                   c("TBA", "Undetermined", ""), ]
  }
  mset
}


# Import Bader modules
file <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  
bader_modules <- my_tmodImportMSigDB(file, format = "gmt")

out <- my_importMsigDBGMT(file)

# Convert mouse gene names to human
fg <- get_homologs(target_set, "mouse", ordered = F)[["human_genes"]]
bg <- get_homologs(background_set, "mouse", ordered = F)[["human_genes"]]

mset <- bader_modules
modules <- NULL
qval = 1.1
filter = FALSE
cols = "Title"
order.by <- "pval"
nodups <- TRUE


bg <- setdiff(bg, fg)
if (length(bg) == 0){stop("All features from bg in fg.")}
tot <- unique(c(fg, bg))
n <- length(tot)
k <- length(fg)

n_modules <- length(bader_modules$MODULES2GENES)
enrichment <- matrix(data = 0, nrow = n_modules, ncol = 6)
for (i in 1:n_modules) {
  
  # n_id <- names(bader_modules$MODULES2GENES)[i]
  # Not sure if this is right
  mg <- bader_modules$MODULES2GENES[[i]]
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
enrichment[["n_id"]] <- names(bader_modules$MODULES2GENES)
