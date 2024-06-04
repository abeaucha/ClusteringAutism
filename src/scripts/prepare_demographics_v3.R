# ----------------------------------------------------------------------------
# clean_v3_demographics.R
# Author: Antoine Beauchamp
#
# Clean human demographics information for v3 registration run.
#
# Description
# -----------


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Environment variables ------------------------------------------------------

PROJECTPATH <- Sys.getenv("PROJECTPATH")


# Main -----------------------------------------------------------------------

## Directories ---------------------------------------------------------------

# Registration base directories
registration_dir_v2 <- file.path(PROJECTPATH, "data/human/registration/v2/")
registration_dir_v3 <- file.path(PROJECTPATH, "data/human/registration/v3/")

# Jacobian image directories
imgdir_v3 <- file.path(registration_dir_v3, "jacobians")

# Demographics directories
demographics_dir_v2 <- file.path(registration_dir_v2, "subject_info/")
demographics_dir_v3 <- file.path(registration_dir_v3, "subject_info/")

# File suffix for the Jacobian images
file_suffix_v2 <- ".denoise.extracted_fwhm_2mm.mnc"
file_suffix_v3 <- ".denoise_fwhm_2mm.nii.gz"

# Images files and IDs
jacobians <- "absolute"
imgfiles <- imgdir_v3 %>% 
  file.path(jacobians, "smooth") %>% 
  list.files(pattern = "*.nii.gz")

# Create a data frame with image files and IDs
imgs <- tibble(file = imgfiles) %>% 
  mutate(Scan_ID = file %>% 
           str_remove(file_suffix_v3))


## POND + SickKids data ------------------------------------------------------

# Import v2 demographics
demographics_v2 <- "demographics.csv"
demographics_v2 <- file.path(demographics_dir_v2, demographics_v2)
demographics_v2 <- read_csv(demographics_v2, show_col_types = FALSE)

# Create a copy of v2 demographics for v3
demographics_v3 <- demographics_v2 

# Update image file names for v3
demographics_v3 <- demographics_v3 %>% 
  mutate(file = file %>% 
           str_remove(file_suffix_v2) %>% 
           str_c(file_suffix_v3)) 

# Rename POND subject IDs for v3
demographics_v3 <- demographics_v3 %>% 
  mutate(Subject_ID = str_replace(Subject_ID, "POND_", "sub-")) 

# Select desired columns
demographics_v3 <- demographics_v3 %>% 
  select(Subject_ID, Scan_ID, Dataset,
         Age, Sex, DX, 
         Site, Scanner,
         file)  


## HBN data ------------------------------------------------------------------

# HBN demographics directory
hbn_dir <- file.path(demographics_dir_v3, "HBN")

# Files containing basic demographics information (e.g. age, sex)
hbn_basic_files <- c("MRI_CBIC_participants.csv", 
                     "MRI_CUNY_participants.csv",
                     "MRI_RU_participants.csv")
hbn_basic_files <- file.path(hbn_dir, hbn_basic_files)
names(hbn_basic_files) <- c("CBIC", "CUNY", "RU")

# Import basic demographics data
hbn_basic <- hbn_basic_files %>% 
  map(read_csv, show_col_types = FALSE) %>% 
  map(function(x){
    x %>%
      select(Subject_ID = participant_id,
             Sex, Age)
  }) %>% 
  bind_rows(.id = "Site") %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male"))

# These CSVs seem to have duplicated entries. 
# Remove those.
hbn_basic <- distinct(hbn_basic)

# Create a data frame with image files and IDs
hbn_imgs_all <- imgs %>% 
  filter(str_starts(file, "sub-N")) %>% 
  mutate(Subject_ID = Scan_ID %>% 
           str_remove("sub-") %>% 
           str_remove("_.*"))

# Some participants have multiple scans
# Pull the IDs for those participants
hbn_participants_multiple_imgs <- hbn_imgs_all %>% 
  group_by(Subject_ID) %>% 
  count() %>% 
  ungroup() %>% 
  filter(n > 1) %>% 
  pull(Subject_ID)

# Extract imgs for participants with a single scan
hbn_imgs_unique <- hbn_imgs_all %>% 
  filter(!(Subject_ID %in% hbn_participants_multiple_imgs))

# Extract imgs for participants with multiple scans
hbn_imgs_multiple <- hbn_imgs_all %>% 
  filter(Subject_ID %in% hbn_participants_multiple_imgs)

# For participants with multiple scans, select Run 1 scans
hbn_imgs_multiple <- hbn_imgs_multiple %>% 
  filter(str_detect(Scan_ID, "run-01"))

# Data frame of unique scans
hbn_imgs <- bind_rows(hbn_imgs_unique, hbn_imgs_multiple)

# Join image files to basic demographics info
hbn_imgs <- hbn_imgs %>%
  left_join(hbn_basic, by = "Subject_ID")

# Import HBN diagnostic information
hbn_dx <- "9994_ConsensusDx_20220728.csv"
hbn_dx <- file.path(hbn_dir, "assessment_data", hbn_dx)
hbn_dx <- read_csv(hbn_dx, show_col_types = FALSE)
hbn_dx <- hbn_dx[2:nrow(hbn_dx),]

# Disorder groups to include
hbn_dx_include <- c("Neurodevelopmental Disorders",
                    "Obsessive Compulsive and Related Disorders",
                    "No Diagnosis Given")

# Neurodevelopmental disorders to exclude
hbn_dx_exclude <- c("Language Disorder",
                    "Specific Learning Disorder with Impairment in Mathematics",
                    "Specific Learning Disorder with Impairment in Reading",
                    "Specific Learning Disorder with Impairment in Written Expression",
                    "Speech Sound Disorder",
                    "Social (Pragmatic) Communication Disorder",
                    "Other Specified Neurodevelopmental Disorder",
                    "Unspecified Neurodevelopmental Disorder",
                    "Excoriation (Skin-Picking) Disorder")

# Select DX columns and rows
hbn_dx <- hbn_dx %>% 
  select(Subject_ID = EID, 
         DX_01_Cat,
         DX_01) %>% 
  filter(DX_01_Cat %in% hbn_dx_include,
         !(DX_01 %in% hbn_dx_exclude),
         !is.na(DX_01)) %>% 
  arrange(DX_01_Cat, DX_01)

# The diagnostic data frame seems to have duplicated entries.
# Remove those.
hbn_dx <- distinct(hbn_dx)

# Export full set of DX information
outfile <- "HBN_DX_categories.csv"
outfile <- file.path(hbn_dir, outfile)
write_csv(x = hbn_dx, file = outfile)

# Combine some DX using new labels
hbn_dx_labels <- tibble(DX_01 = c("ADHD-Combined Type", 
                                  "ADHD-Hyperactive/Impulsive Type",
                                  "ADHD-Inattentive Type",
                                  "Autism Spectrum Disorder",
                                  "Intellectual Disability-Mild",
                                  "Intellectual Disability-Moderate",
                                  "Intellectual Disability-Severe",
                                  "Other Specified Attention-Deficit/Hyperactivity Disorder", 
                                  "Other Specified Tic Disorder", 
                                  "Persistent (Chronic) Motor or Vocal Tic Disorder",
                                  "Provisional Tic Disorder",
                                  "Tourettes Disorder",
                                  "Unspecified Attention-Deficit/Hyperactivity Disorder",
                                  "Obsessive-Compulsive Disorder",
                                  "No Diagnosis Given"),
                        DX = c("ADHD", 
                               "ADHD",
                               "ADHD",
                               "ASD",
                               "ID",
                               "ID",
                               "ID",
                               "ADHD",
                               "Tourette Syndrome",
                               "Tourette Syndrome",
                               "Tourette Syndrome",
                               "Tourette Syndrome",
                               "ADHD", 
                               "OCD",
                               "Control"))

# Join new labels to DX data frame
hbn_dx <- hbn_dx %>% 
  left_join(hbn_dx_labels, by = "DX_01") %>% 
  select(Subject_ID, DX)

# Join image and diagnostic data frames
hbn <- left_join(hbn_imgs,
                 hbn_dx,
                 by = "Subject_ID") %>% 
  mutate(Dataset = "HBN",
         Scanner = "")


## Combine data --------------------------------------------------------------

# Combine demographics data frames
demographics <- bind_rows(demographics_v3, hbn)

# Replace file NIFTY extension with MINC
demographics <- demographics %>% 
  mutate(file = str_replace(file, ".nii.gz", ".mnc"))

# Remove participants with missing demographics
demographics <- demographics %>% 
  filter(!is.na(DX), 
         !is.na(Sex),
         !is.na(Age),
         !is.na(Site),
         !is.na(Scanner))

# These participants have messed up Jacobian images
ids_to_remove <- c("d8_0034_02",
                   "d8_0040_04",
                   "d8_0041_03",
                   "d8_0044_04")

# Remove broken images
demographics <- demographics %>% 
  filter(!(Subject_ID %in% ids_to_remove))

# Arrange columns and rows
demographics <- demographics %>% 
  select(Subject_ID, Scan_ID, Dataset,
         Age, Sex, DX, 
         Site, Scanner,
         file) %>% 
  arrange(Dataset, Subject_ID)

# Export
outfile <- "demographics.csv"
outfile <- file.path(demographics_dir_v3, outfile)
write_csv(x = demographics, file = outfile)
