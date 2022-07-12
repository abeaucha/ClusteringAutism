# ----------------------------------------------------------------------------
# generate_structure_masks.R
# Author: Antoine Beauchamp
# Created: July 12th, 2022
#
# Script to create AHBA microarray samples masked to a set of ROIs
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(MRIcrotome))
suppressPackageStartupMessages(library(data.tree))


# Functions ------------------------------------------------------------------

source('../../functions/tree_tools.R')

label_human_data <- function(x, tree, samples) {
  
  tree <- Clone(tree)
  
  tree$Do(function(node){
    node$struct <- rep(node$name, length(node$samples))
  }, traversal = 'post-order')
  
  tree_structs <- unlist(tree$Get('struct', filterFun = isLeaf))
  names(tree_structs) <- NULL
  
  tree_samples <- unlist(tree$Get('samples', filterFun = isLeaf))
  names(tree_samples) <- NULL
  
  sample_defs <- tibble(Structure = tree_structs,
                        sample_id = tree_samples)
  
  x <- x %>% 
    mutate(sample_id = samples) %>% 
    left_join(sample_defs, by = 'sample_id') %>% 
    select(-sample_id) %>% 
    rename(Region = Structure)
  
  return(x)
  
}


# Main -----------------------------------------------------------------------

#Human tree
human_treefile <- '../../data/human/HumanExpressionTree.RData'

#Tree labels
treelabels <- '../../data/TreeLabelsReordered.RData'

#Trees
load(human_treefile)
tree_human <- Clone(treeHumanExpr)
rm(treeHumanExpr)

#Tree labels
load(treelabels)

#Pruned human tree to desired level
human_regions <- c(listLabelsHumanReordered$Region88, 
                   'white matter', 
                   'sulci & spaces')
tree_human_pruned <- Clone(tree_human)
pruneAnatTree(tree_human_pruned, 
              nodes = human_regions, 
              method = 'BelowNode')

structs <- c('pons',
             'amygdala',
             'caudate nucleus',
             'thalamus',
             'dentate gyrus',
             'precentral gyrus',
             'myelencephalon')


# Create masks in MNI space --------------------------------------------------

samplefile <- 'data/AHBA_microarray_coordinates_mni.csv'
defsfile <- 'data/AHBA_microarray_coordinates_mni_defs.csv'
maskfile <- 'data/AHBA_microarray_mask_mni.mnc'
labelfile <- 'data/AHBA_microarray_labels_mni.mnc'

sample_coordinates <- suppressMessages(read_csv(samplefile))
defs <- suppressMessages(read_csv(defsfile))
sample_mask <- mincGetVolume(maskfile)
sample_labels <- mincGetVolume(labelfile)

defs_labelled <- label_human_data(x = defs,
                                  tree = tree_human_pruned,
                                  samples = defs$sample_id)

for (i in 1:length(structs)){
  
  struct <- structs[i]
  
  struct_labels <- defs_labelled %>% 
    filter(Region == struct) %>% 
    pull(label)
  
  struct_mask <- sample_mask
  struct_mask[!(sample_labels %in% struct_labels)] <- 0
  
  struct <- str_replace_all(structs[i], ' ', '_')
  
  outfile <- str_replace(maskfile, '.mnc', str_c('_', struct, '.mnc'))
  
  mincWriteVolume(buffer = struct_mask,
                  output.filename = outfile,
                  clobber = TRUE) 
}


# Create masks in study space ------------------------------------------------

samplefile <- 'data/AHBA_microarray_coordinates_mni_to_study_using_inverses.csv'
defsfile <- 'data/AHBA_microarray_coordinates_mni_defs.csv'
maskfile <- 'data/AHBA_microarray_mask_study_using_inverses.mnc'
labelfile <- 'data/AHBA_microarray_labels_study_using_inverses.mnc'

sample_coordinates <- suppressMessages(read_csv(samplefile))
defs <- suppressMessages(read_csv(defsfile))
sample_mask <- mincGetVolume(maskfile)
sample_labels <- mincGetVolume(labelfile)

defs_labelled <- label_human_data(x = defs,
                                  tree = tree_human_pruned,
                                  samples = defs$sample_id)

for (i in 1:length(structs)){
  
  struct <- structs[i]
  
  struct_labels <- defs_labelled %>% 
    filter(Region == struct) %>% 
    pull(label)
  
  struct_mask <- sample_mask
  struct_mask[!(sample_labels %in% struct_labels)] <- 0
  
  struct <- str_replace_all(structs[i], ' ', '_')
  
  outfile <- str_replace(maskfile, '.mnc', str_c('_', struct, '.mnc'))
  
  mincWriteVolume(buffer = struct_mask,
                  output.filename = outfile,
                  clobber = TRUE) 
}
