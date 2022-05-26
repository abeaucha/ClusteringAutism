# ----------------------------------------------------------------------------
# calculate_effect_sizes.R
# Author: Antoine Beauchamp
# Created: May 9th, 2022
#
# Calculate voxel-wise effect sizes for human patients
#
# Description
# -----------
#


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(doParallel))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--demographics',
              type = 'character',
              help = "Path to CSV file containing demographics data."),
  make_option('--imgdir',
              type = 'character',
              help = ""),
  make_option('--maskfile',
              type = 'character',
              help = ""),
  make_option('--outdir',
              type = 'character',
              help = ""),
  make_option('--ncontrols',
              type = 'numeric',
              help = ""),
  make_option('--parallel',
              type = 'character',
              help = ""),
  make_option('--nproc',
              type = 'numeric',
              help = "")
) 


# Functions ------------------------------------------------------------------


#' Get controls for a study participant
#'
#' @param participant 
#' @param controls 
#' @param ncontrols 
#' @param seed 
#'
#' @return
get_control_files <- function(participant, controls, imgfiles,
                              ncontrols = 10, seed = NULL){
  
  if (!is.null(seed)){set.seed(seed)}
  
  participant_age <- participant[, 'Age'][[1]]
  participant_sex <- participant[, 'Sex'][[1]]
  
  control_ids <- controls %>%
    filter(Sex == participant_sex) %>% 
    mutate(AgeDiff = abs(Age - participant_age)) %>% 
    top_n(n = -1*ncontrols, wt = AgeDiff) %>% 
    sample_n(size = ncontrols, replace = FALSE) %>% 
    pull(Extract_ID) %>% 
    str_remove('.mnc')
  
  control_files <- character(length(control_ids))
  for (j in 1:length(control_ids)){
    control_files[j] <- str_subset(imgfiles, control_ids[j])
  }
  
  return(control_files)
  
}


#' Compute the effect size for one participant
#'
#' @param participant 
#' @param controls 
#' @param mask 
#'
#' @return
compute_effect_size <- function(participant, controls, imgfiles, mask) {
  
  control_mean <- mincMean(filenames = controls)
  control_mean <- control_mean[,1]
  
  control_sd <- mincSd(filenames = controls)
  control_sd <- control_sd[,1]
  
  participant_id <- participant[, 'Extract_ID'][[1]] %>%
    str_remove('.mnc')
  
  participant_file <- str_subset(imgfiles, participant_id)
  
  participant_vol <- mincGetVolume(participant_file)
  
  participant_z <- (participant_vol - control_mean)/control_sd
  
  mask <- mincGetVolume(mask)
  participant_z[mask == 0] <- 0
  participant_z[is.infinite(participant_z)] <- 0
  participant_z[is.na(participant_z)] <- 0
  
  class(participant_z) <- class(participant_vol)
  attributes(participant_z) <- attributes(participant_vol)
  
  return(participant_z)
}


#' Title
#'
#' @param participant 
#' @param controls 
#' @param ncontrols 
#' @param mask 
#' @param outdir 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
executor <- function(participant, controls, ncontrols = 10, 
                     mask, imgdir, outdir, seed = NULL) {
  
  imgfiles <- list.files(imgdir, full.names = TRUE) %>% 
    str_subset('.mnc')
  
  control_files <- get_control_files(participant = participant,
                                     controls = controls,
                                     imgfiles = imgfiles,
                                     ncontrols = ncontrols,
                                     seed = seed)
  effect_size <- compute_effect_size(participant = participant,
                                     controls = control_files,
                                     imgfiles = imgfiles,
                                     mask = mask)

  outfile <- attributes(effect_size)[['filename']] %>%
    basename() %>%
    str_replace('.mnc', str_c('_effectsize_ncontrols', ncontrols, '.mnc'))
  outfile <- file.path(outdir, outfile)

  mincWriteVolume(effect_size,
                  output.filename = outfile,
                  clobber = TRUE)
  
  return(1)
} 


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
demofile <- args[['demographics']]
imgdir <- args[['imgdir']]
maskfile <- args[['maskfile']]
outdir <- args[['outdir']]
ncontrols <- args[['ncontrols']]
inparallel <- ifelse(args[['parallel']] == 'true', TRUE, FALSE)

# demofile <- 'data/human/registration/DBM_input_demo_passedqc.csv'
# imgdir <- 'data/human/registration/jacobians/absolute/smooth/'
# maskfile <- 'data/human/registration/reference_files/mask.mnc'
# outdir <- 'data/human/effect_sizes/absolute/'
# ncontrols <- 10
# inparallel <- TRUE

#Import demographics data
demographics <- data.table::fread(demofile, header = TRUE) %>% 
  as_tibble()

#Remove entries with no DX
demographics <- demographics %>% 
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex))

#Extract control participants
controls <- demographics %>% 
  filter(DX == 'Control')

#Extract participants
participants <- demographics %>% 
  filter(DX != 'Control')

#Option to run in parallel
if (inparallel) {
  nproc <- args[['nproc']]
  cl <- makeCluster(nproc)
  registerDoParallel(cl)
  tmp <- foreach(i = 1:nrow(participants), 
                 .packages = c('tidyverse', 'RMINC')) %dopar% {
                   executor(participant = participants[i,],
                            controls = controls,
                            ncontrols = ncontrols,
                            mask = maskfile,
                            imgdir = imgdir,
                            outdir = outdir,
                            seed = i)
                 } 
  stopCluster(cl)
} else {
  tmp <- foreach(i = 1:nrow(participants), 
                 .packages = c('tidyverse', 'RMINC')) %do% {
                   executor(participant = participants[i,],
                            controls = controls,
                            ncontrols = ncontrols,
                            mask = maskfile,
                            imgdir = imgdir,
                            outdir = outdir,
                            seed = i)
                 }
}