# ----------------------------------------------------------------------------
# label_expression_matrix.R 
# Antoine Beauchamp
# 
# Description
# -----------
# 

# Libraries ------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(optparse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--species',
              type = 'character',
              help = paste("The species whose data to process.")),
  make_option('--matrix',
              type = 'character',
              help = paste("Path to CSV file containing voxel- or sample-wise",
                           "expression matrix to label.")),
  make_option('--tree',
              type = 'character',
              help = paste("Path to .RData file containing neuroanatomical tree")),
  make_option('--treelabels',
              type = 'character'),
  make_option('--nlabels',
              type = 'numeric'),
  make_option('--greymatter',
              type = 'character',
              default = 'true',
              help = paste("Option to use only grey matter measurements.",
                           "[default %default]")),
  make_option('--labels',
              type = 'character'),
  make_option('--mask',
              type = 'character'),
  make_option('--outdir',
              type = 'character',
              default = 'data/',
              help = paste("[default %default]")),
  make_option('--verbose',
              type = 'character',
              default = 'true',
              help = paste("Verbosity option. [default %default]"))
)


# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset('--file=') %>% 
  str_remove('--file=') %>% 
  dirname()

path_tree_tools <- file.path(working_dir, 
                             script_dir, 
                             'functions', 
                             'tree_tools.R')
source(path_tree_tools)

path_processing_tools <- file.path(working_dir,
                                   script_dir,
                                   'functions',
                                   'processing_tools.R')
source(path_processing_tools)


#' Label voxels/samples with neuroanatomical regions
#'
#' @param measurements 
#' @param tree 
#' @param treefield 
#'
#' @return
labelRegions <- function(measurements, tree, treefield){
  
  #Get a list of all measurements and their corresponding ROIs
  listMeasurements <- tree$Get(treefield, filterFun = isLeaf)
  vecMeasurements <- unlist(listMeasurements)
  
  #Unlisting to vector will break the names. Fix them
  structNames <- names(listMeasurements)
  measurementNames <- names(vecMeasurements)
  
  for (struct in structNames){
    
    structRgx <- struct %>% 
      str_replace("\\(", "\\\\(") %>% 
      str_replace("\\)", "\\\\)")
    
    structRgx <- str_c("^", structRgx, "[0-9]+", "$")
    
    indStruct <- str_which(measurementNames, structRgx)
    
    measurementNames[indStruct] <- struct
  }
  
  #Proper names  
  names(vecMeasurements) <- measurementNames

  indMatchStructs <- match(measurements, vecMeasurements)
  
  measurementStructs <- names(vecMeasurements)[indMatchStructs]
  
  return(measurementStructs)
}


# Init -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

if (!(args[['greymatter']] %in% c('true', 'false'))) {
  stop(paste("--greymatter must be one of [true, false].",
             "Got:", args[['greymatter']]))
}

if (!(args[['verbose']] %in% c('true', 'false'))) {
  stop(paste("--verbose must be one of [true, false].",
             "Got:", args[['verbose']]))
}

verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)
gm <- ifelse(args[['greymatter']] == 'true', TRUE, FALSE)
species <- args[['species']]
fileMat <- args[['matrix']]
fileTree <- args[['tree']]
fileTreeLabels <- args[['treelabels']]
nlabels <- args[['nlabels']]

# species <- 'mouse'
# gm <- TRUE
# nlabels <- 67
# fileMat <- 'data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv'
# # fileMat <- "data/HumanExpressionMatrix_samples_pipeline_v1_homologs.csv"
# fileTree <- 'data/mouse/MouseExpressionTree_DSURQE.RData'
# # fileTree <- "data/human/HumanExpressionTree.RData"
# fileTreeLabels <- 'data/TreeLabels.RData'

if (!(species %in% c('mouse', 'human'))) {
  stop(paste("--species must be one of [mouse, human].",
             "Got:", species))
}

if (is.null(fileMat)){
  stop("--matrix empty with no default.")
}

if (is.null(fileTree)){
  stop("--tree empty with no default.")
}

if (is.null(fileTreeLabels)){
  stop("--treelabels empty with no default.")
}

if (is.null(nlabels)){
  stop("--nlabels empty with no default.")
}

if (!file.exists(fileMat)) {
  stop("Matrix file ", fileMat, " not found.")
}

if (!file.exists(fileTree)) {
  stop("Tree file ", fileTree, " not found.")
}

if (!file.exists(fileTreeLabels)) {
  stop("Tree label file ", fileTreeLabels, " not found.")
}

if (species == 'mouse') {
  
  labels <- args[['labels']]
  mask <- args[['mask']]
  
  # labels <- 'data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc'
  # mask <- 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc'
  
  if (is.null(labels)){
    stop("--labels empty with no default. Required when --species is mouse.")
  }
  
  if (is.null(mask)){
    stop("--mask empty with no default. Required when --species is mouse.")
  }
  
  if (!file.exists(labels)) {
    stop("Mouse label file ", labels, " not found.")
  }
  
  if (!file.exists(mask)) {
    stop("Mouse mask file ", mask, " not found.")
  }
  
}

if (verbose) {message("Labelling data from file: ", fileMat)}


# Importing and processing ---------------------------------------------------

if (verbose) {message("Importing the expression matrix...")}

#Import expression matrices
dfExpr <- suppressMessages(data.table::fread(fileMat,
                                             header = TRUE)) %>% 
  as_tibble()

#Extract genes and remove from df
genes <- dfExpr$Gene

dfExpr <- dfExpr %>% select(-Gene)

#Clean up mouse column names
if (species == 'mouse') {
  colnames(dfExpr) <- str_c("V", colnames(dfExpr))
}


# Build the data tree --------------------------------------------------------

if (verbose) {message("Importing the data tree...")}
  
load(fileTree)

if (species == 'mouse') {

  #Rename the tree
  tree <- Clone(treeMouseExpr)
  rm(treeMouseExpr)
  
  #Load atlas labels and mask
  label_vol <- mincGetVolume(labels)
  mask_vol <- mincGetVolume(mask)

  #Mask the atlas
  labels_masked <- label_vol[mask_vol == 1]
  
  treefield <- 'voxels'
  
  #Assign voxels to leaf nodes on the tree
  tree$Do(function(node){
    if(isLeaf(node)){
      node$voxels <- colnames(dfExpr)[labels_masked == node$label]
    }
  })
  
  #Aggregate voxel names up the tree
  tree$Do(function(node){
    node$voxels <- unlist(Aggregate(node, treefield, c))
  })

  #Remove white matter and ventricles
  if (gm) {
    cutAtNodes <- c('fiber tracts', 'ventricular systems')
    pruneAnatTree(tree, nodes = cutAtNodes, method = 'AtNode')
    tree <- FindNode(tree, 'Basic cell groups and regions')
  }
  
} else if (species == 'human') {
  
  #Rename
  tree <- Clone(treeHumanExpr)
  rm(treeHumanExpr)
  
  treefield <- 'samples'

  #Remove white matter and ventricles
  if (gm) {
    cutAtNodes <- c('white matter', 'sulci & spaces')
    pruneAnatTree(tree, cutAtNodes, method = 'AtNode')
    tree <- FindNode(tree, 'gray matter')
  }
  
} else {
  stop()
}
  
#Filter expression data for measurements that are in the pruned tree
tree_measurements <- tree[[treefield]]
dfExpr <- dfExpr[, (colnames(dfExpr) %in% tree_measurements)]


# Assign regional labels -----------------------------------------------------

if (verbose) {message("Assigning labels to the expression matrix...")}

#Transpose data frames
dfExpr <- dfExpr %>% as.matrix() %>% t()

#Assign gene names as columns. 
#Note: Necessary to do it this way if there are duplicate genes
colnames(dfExpr) <- genes

#Convert back to data frame
#Note: This will create new names if there are duplicated genes. 
#This is fine and desired.
dfExpr <- dfExpr %>% 
  as_tibble(rownames = "Measurement", .name_repair = "unique") %>% 
  column_to_rownames("Measurement")

#Load tree labels
load(fileTreeLabels)
if (species == 'mouse') {
  listLabels <- listLabelsMouse
} else if (species == 'human') {
  listLabels <- listLabelsHuman
} else {
  stop()
}

#Prune tree to given level of aggregation
region <- paste0('Region', nlabels)

if (region %in% names(listLabels)) {
  treePruned <- Clone(tree)
  pruneAnatTree(treePruned, nodes = listLabels[[region]], method = 'BelowNode')
  
  #Get the region labels for every voxel
  dfExpr[,'Region'] <- labelRegions(measurements = rownames(dfExpr), 
                                    tree = treePruned,
                                    treefield = treefield)
  
  #Remove identifiers from rownames
  rownames(dfExpr) <- NULL
  
} else {
  stop(paste("Label set", region, "not found in --treelabels file", 
             fileTreeLabels, "for --species", species))
}


# Write to file --------------------------------------------------------------

outfile <- fileMat %>% 
  basename() %>% 
  str_remove('.csv') %>% 
  str_c('_labelled') %>% 
  str_c('_', str_to_lower(region)) %>% 
  str_c('.csv')

outdir <- args[['outdir']]
outfile <- file.path(outdir, outfile)

if (verbose) {message(paste("Writing to file:", outfile, "..."))}

data.table::fwrite(dfExpr, file = outfile)
