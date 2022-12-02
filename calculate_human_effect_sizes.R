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
  make_option("--demographics",
              type = "character",
              help = "Path to CSV file containing human demographics data."),
  make_option("--imgdir",
              type = "character",
              help = paste("Path to directory containing images to use",
                           "to compute effect sizes.")),
  make_option("--maskfile",
              type = "character",
              help = "Path to mask file for the images."),
  make_option("--outdir",
              type = "character",
              help = paste("Path to directory in which to save the effect",
                           "size images.")),
  make_option("--ncontrols",
              type = 'numeric',
              default = 10,
              help = paste("Number of controls to use when computing effect",
                           "sizes. [default %default]")),
  make_option("--threshold",
              type = "numeric",
              default = 10,
              help = paste("[default %default]")),
  make_option("--dataset",
              type = "numeric",
              default = 1,
              help = paste("A numeric flag indicating which data to use.", 
                           "Set to 1 for all data. Set to 2 for POND and",
                           "SickKids. Set to 3 for POND only. [default %default]")),
  make_option("--parallel",
              type = "character",
              default = 'false',
              help = "Option to run in parallel. [default %default]"),
  make_option("--nproc",
              type = "numeric",
              help = paste("Number of processors to use in parallel.",
                           "Ignored if --parallel is false."))
) 


# Functions ------------------------------------------------------------------

#' Get propensity-matched controls for study participant
#'
#' @param participant (data.frame) A data.frame row containing demographic 
#' information for the participant of interest. 
#' @param controls (data.frame) Demographic information for the controls.
#' @param ncontrols (numeric scalar) Number of controls to use for propensity
#' matching
#' @param threshold (numeric scalar) 
#' @param seed (numeric scalar) Random seed to use for control sampling.
#'
#' @return (data.frame) Demographic information for the propensity-matched 
#' controls.
get_matched_controls <- function(participant, controls, ncontrols = 10, 
                                 threshold = 10, seed = NULL){
  
  if (!is.null(seed)){set.seed(seed)}
  
  participant_age <- pull(participant, "Age")
  participant_sex <- pull(participant, "Sex")
  participant_site <- pull(participant, "Site")
  participant_scanner <- pull(participant, "Scanner")
  
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
  
  controls <- controls_tmp %>%
    mutate(AgeDiff = abs(Age - participant_age)) %>% 
    top_n(n = -1*ncontrols, wt = AgeDiff) %>% 
    sample_n(size = ncontrols, replace = FALSE) 
  
  return(controls)
  
}


#' Compute the effect size for a participant.
#'
#' @param participant (character scalar) Path to the MINC volume of 
#' interest.
#' @param controls (character vector) Paths to the MINC volumes of the 
#' controls.
#' @param mask (character scalar) Path to mask image MINC file. 
#'
#' @return (mincSingleDim) A vector of voxel-wise effect size values.
compute_effect_size <- function(participant, controls, mask) {
  
  sink(file = nullfile(), type = "output")
  control_mean <- mincMean(filenames = controls)
  control_mean <- control_mean[,1]
  
  control_sd <- mincSd(filenames = controls)
  control_sd <- control_sd[,1]
  sink(file = NULL)
  
  participant_vol <- mincGetVolume(participant)
  
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
#' matching.
#' @param threshold (numeric scalar) 
#' @param mask (character scalar) Path to mask volume. 
#' @param imgdir (character scalar) Path to the directory containing the 
#' image files.
#' @param outdir (character scalar) Path to directory in which to save effect 
#' sizes.
#' @param seed (numeric scalar) Random seed to use for control sampling.
#'
#' @return (numeric scalar) Value of 1 if function succeeds.
executor <- function(participant, controls, ncontrols = 10, 
                     threshold = 10, mask, imgdir, outdir,
                     dataset = NULL, seed = NULL) {
  
  #Get propensity-matched controls
  controls_matched <- get_matched_controls(participant = participant,
                                           controls = controls,
                                           ncontrols = ncontrols,
                                           threshold = threshold,
                                           seed = seed)
  
  #Extract image files
  imgfiles <- str_subset(list.files(imgdir, full.names = TRUE), ".mnc")
  
  #Image files for the controls
  control_ids = str_remove(controls_matched[["Extract_ID"]], ".mnc")
  control_files <- character(length(control_ids))
  for (j in 1:length(control_ids)){
    control_files[j] <- str_subset(imgfiles, control_ids[j])
  }
  
  #Image file for the participant
  participant_id <- str_remove(participant[[1, "Extract_ID"]], ".mnc")
  participant_file <- str_subset(imgfiles, participant_id)
  
  #Compute the effect size volume
  effect_size <- compute_effect_size(participant = participant_file,
                                     controls = control_files,
                                     mask = mask)
  
  #Get the resolution
  res <- max(minc.separation.sizes(participant_file))
  
  #Output file
  outfile <- participant_id %>% 
    basename() %>% 
    str_c("ES", 
          "res", res,
          "data", dataset, 
          "nc", ncontrols, 
          "thresh", threshold, sep = "_") %>% 
    str_c(".mnc")
  outfile <- file.path(outdir, outfile)
  mincWriteVolume(effect_size,
                  output.filename = outfile,
                  clobber = TRUE)

  return(1)
} 


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
demographics <- args[["demographics"]]
imgdir <- args[["imgdir"]]
maskfile <- args[["maskfile"]]
outdir <- args[["outdir"]]
ncontrols <- args[["ncontrols"]]
threshold <- args[["threshold"]]
dataset <- args[["dataset"]]
inparallel <- ifelse(args[["parallel"]] == "true", TRUE, FALSE)

#Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Import demographics data
demographics <- as_tibble(data.table::fread(demographics, header = TRUE))

#Filter data if desired
if (is.numeric(dataset)) {
  if (dataset != 1) {
    if (dataset == 2){
      demographics <- demographics %>% 
        filter(Dataset %in% c("POND", "SickKids"))
    } else if (dataset == 3) {
      demographics <- demographics %>% 
        filter(Dataset == "POND")
    } else {
      stop(str_c("Argument --dataset must be one of {1, 2, 3}: ", dataset))
    }
  }
} else {
  stop(str_c("Argument --dataset must be one of {1, 2, 3}: ", dataset))
}

#Remove entries with missing diagnosis, age, or sex
demographics <- demographics %>% 
  filter(!is.na(DX),
         !is.na(Age),
         !is.na(Sex),
         !is.na(Site),
         !is.na(Scanner))

#Extract control demographics 
controls <- demographics %>% 
  filter(DX == "Control")

#Extract participant demographics
participants <- demographics %>% 
  filter(DX != "Control")

#Compute effect sizes
pb <- txtProgressBar(max = nrow(participants), style = 3)
progress <- function(n) {setTxtProgressBar(pb = pb, value = n)}
if (inparallel) {
  nproc <- args[["nproc"]]
  cl <- makeSOCKcluster(nproc)
  registerDoSNOW(cl)
  opts <- list(progress=progress)
  tmp <- foreach(i = 1:nrow(participants), 
                 .packages = c("tidyverse", "RMINC"),
                 .options.snow=opts) %dopar% {
                   executor(participant = participants[i,],
                            controls = controls,
                            ncontrols = ncontrols,
                            threshold = threshold,
                            mask = maskfile,
                            imgdir = imgdir,
                            outdir = outdir,
                            dataset = dataset,
                            seed = i)
                 } 
  close(pb)
  stopCluster(cl)
} else {
  tmp <- foreach(i = 1:nrow(participants), 
                 .packages = c("tidyverse", "RMINC")) %do% {
                   progress(n = i)
                   executor(participant = participants[i,],
                            controls = controls,
                            ncontrols = ncontrols,
                            threshold = threshold,
                            mask = maskfile,
                            imgdir = imgdir,
                            outdir = outdir,
                            dataset = dataset,
                            seed = i)
                 }
  close(pb)
}
