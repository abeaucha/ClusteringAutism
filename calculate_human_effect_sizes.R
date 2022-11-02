# ----------------------------------------------------------------------------
# calculate_human_effect_sizes.R
# Author: Antoine Beauchamp
# Created: May 9th, 2022
#
# Calculate voxel-wise effect sizes for human patients.
#
# Description
# -----------
# This script computes effect sizes for human study participants. Effect sizes
# are calculated on a voxel-wise basis by taking z-scores with respect to 
# study controls. Controls are matched to participants based on sex, site, 
# scanner and minimal absolute age difference. Effect size images are saved as 
# MINC files.


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(tcltk))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option('--demographics',
              type = 'character',
              help = "Path to CSV file containing human demographics data."),
  make_option('--imgdir',
              type = 'character',
              help = paste("Path to directory containing images to use",
                           "to compute effect sizes.")),
  make_option('--maskfile',
              type = 'character',
              help = "Path to mask file for the images."),
  make_option('--outdir',
              type = 'character',
              help = paste("Path to directory in which to save the effect",
                           "size images.")),
  make_option('--ncontrols',
              type = 'numeric',
              default = 10,
              help = paste("Number of controls to use when computing effect",
                           "sizes. [default %default]")),
  make_option('--threshold',
              type = 'numeric',
              default = 10,
              help = paste("[default %default]")),
  make_option('--parallel',
              type = 'character',
              default = 'false',
              help = "Option to run in parallel. [default %default]"),
  make_option('--nproc',
              type = 'numeric',
              help = paste("Number of processors to use in parallel.",
                           "Ignored if --parallel is false."))
) 


# Functions ------------------------------------------------------------------

#' Get propensity-matched controls for study participant
#'
#' @param participant (data.frame) A data.frame row containing demographic 
#' information for the participant of interest. 
#' @param controls (data.frame) Demographic information for the controls.
#' @param imgfiles (character vector) A set of paths to MINC image files that 
#' include (but are not limited to) all control images.
#' @param ncontrols (numeric scalar) Number of controls to use for propensity
#' matching
#' @param threshold (numeric scalar) 
#' @param seed (numeric scalar) Random seed to use for control sampling.
#'
#' @return (character vector) The paths to the images of the propensity-
#' matched controls.
get_control_files <- function(participant, controls, imgfiles,
                              ncontrols = 10, threshold = 10, seed = NULL){
  
  if (!is.null(seed)){set.seed(seed)}
  
  participant_age <- pull(participant, 'Age')
  participant_sex <- pull(participant, 'Sex')
  participant_site <- pull(participant, 'Site')
  participant_scanner <- pull(participant, 'Scanner')
  
  controls_tmp <- controls %>% 
    filter(Sex == participant_sex,
           Site == participant_site,
           Scanner == participant_scanner)
  
  ncontrols_test <- nrow(controls_tmp)
  
  if (ncontrols_test < threshold) {
    
    controls_tmp <- controls %>% 
      filter(Sex == participant_sex,
             Site == participant_site)
    
    ncontrols_test <- nrow(controls_tmp)
    
    if (ncontrols_test < threshold) {
      
      controls_tmp <- controls %>% 
        filter(Sex == participant_sex)
      
    }
  }
  
  controls <- controls_tmp
  
  control_ids <- controls %>%
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


#' Compute the effect size for a participant
#'
#' @param participant (data.frame) A data.frame row containing demographic 
#' information for the participant of interest. 
#' @param controls (character vector) Paths to the MINC images of the 
#' propensity-matched controls.
#' @param imgfiles (character vector) Paths to image MINC files including 
#' (but not limited to) the participant of interest.
#' @param mask (character scalar) Path to mask image MINC file. 
#'
#' @return (mincSingleDim) A vector of voxel-wise effect size values.
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


#' Calculate a participant's effect size
#'
#' @param participant (data.frame) A data.frame row containing demographic 
#' information for the participant of interest.
#' @param controls (data.frame) Demographic information for the controls.
#' @param ncontrols (numeric scalar) Number of controls to use for propensity
#' matching
#' @param threshold (numeric scalar) 
#' @param mask (character scalar) Path to mask image MINC file. 
#' @param imgdir (character vector) Paths to image MINC files.
#' @param outdir (character scalar) Path to directory in which to save effect 
#' size MINC image.
#' @param seed (numeric scalar) Random seed to use for control sampling.
#'
#' @return (numeric scalar) Value of 1 if function succeeds.
executor <- function(participant, controls, ncontrols = 10, 
                     threshold = 10, mask, imgdir, outdir, seed = NULL) {
  
  imgfiles <- list.files(imgdir, full.names = TRUE) %>% 
    str_subset('.mnc')
  
  control_files <- get_control_files(participant = participant,
                                     controls = controls,
                                     imgfiles = imgfiles,
                                     ncontrols = ncontrols,
                                     threshold = threshold,
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
threshold <- args[['threshold']]
inparallel <- ifelse(args[['parallel']] == 'true', TRUE, FALSE)

#Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Import demographics data
demographics <- data.table::fread(demofile, header = TRUE) %>% 
  as_tibble()

#Filter for POND and SickKids 
demographics <- demographics %>% 
  filter(Dataset %in% c("POND", "SickKids"))

#Remove entries with missing diagnosis, age, or sex
demographics <- demographics %>% 
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex),
         !is.na(Site),
         !is.na(Scanner))

#Extract control demographics 
controls <- demographics %>% 
  filter(DX == 'Control')

#Extract participant demographics
participants <- demographics %>% 
  filter(DX != 'Control')

#Compute effect sizes
pb <- txtProgressBar(max = nrow(participants), style = 3)
progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
if (inparallel) {
  nproc <- args[['nproc']]
  cl <- makeSOCKcluster(nproc)
  registerDoSNOW(cl)
  opts <- list(progress=progress)
  tmp <- foreach(i = 1:nrow(participants), 
                 .packages = c('tidyverse', 'RMINC'),
                 .options.snow=opts) %dopar% {
                   executor(participant = participants[i,],
                            controls = controls,
                            ncontrols = ncontrols,
                            threshold = threshold,
                            mask = maskfile,
                            imgdir = imgdir,
                            outdir = outdir,
                            seed = i)
                 } 
  close(pb)
  stopCluster(cl)
} else {
  tmp <- foreach(i = 1:nrow(participants), 
                 .packages = c('tidyverse', 'RMINC')) %do% {
                   progress(n = i)
                   executor(participant = participants[i,],
                            controls = controls,
                            ncontrols = ncontrols,
                            threshold = threshold,
                            mask = maskfile,
                            imgdir = imgdir,
                            outdir = outdir,
                            seed = i)
                 }
  close(pb)
}
