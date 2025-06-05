file <- bader_modules

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
msig <- new("tmod", msig)
msig