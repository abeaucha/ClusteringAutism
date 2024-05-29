# ----------------------------------------------------------------------------
# prepare_HBN_demographics.R
# Author: Antoine Beauchamp
#


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Main -----------------------------------------------------------------------

## Directories ---------------------------------------------------------------

# Registration base directories
registration_dir_v3 <- "data/human/registration/v3/"

# Jacobian image directories
imgdir_v3 <- file.path(registration_dir_v3, "jacobians")

# Demographics directories
demographics_dir_v3 <- file.path(registration_dir_v3, "subject_info/")

# File suffix for the Jacobian images
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

# dx_group <- "NDD"
# dx_group <- "Anxiety"
# dx_group <- "AnxietyDepression"
# dx_group <- "Depression"
dx_group <- "Learning"
# dx_group <- "Language"
if (dx_group == "NDD") {
  
  #Disorder groups to include
  hbn_dx_include <- c("Neurodevelopmental Disorders",
                      "Obsessive Compulsive and Related Disorders",
                      "No Diagnosis Given")
  
  #Neurodevelopmental disorders to exclude
  hbn_dx_exclude <- c("Language Disorder",
                      "Specific Learning Disorder with Impairment in Mathematics",
                      "Specific Learning Disorder with Impairment in Reading",
                      "Specific Learning Disorder with Impairment in Written Expression",
                      "Speech Sound Disorder",
                      "Social (Pragmatic) Communication Disorder",
                      "Other Specified Neurodevelopmental Disorder",
                      "Unspecified Neurodevelopmental Disorder",
                      "Excoriation (Skin-Picking) Disorder")
  
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
  
} else if (dx_group == "Anxiety") {
  
  hbn_dx_include <- c("Anxiety Disorders",
                      "No Diagnosis Given")
  hbn_dx_exclude <- c("")
  
  hbn_dx_labels <- hbn_dx %>% 
    filter(DX_01_Cat %in% hbn_dx_include) %>% 
    select(DX_01) %>% 
    distinct() %>% 
    mutate(DX = ifelse(DX_01 == "No Diagnosis Given",
                       "Control", DX_01)) 
  
} else if (dx_group == "AnxietyDepression") {
  
  hbn_dx_include <- c("Anxiety Disorders",
                      "Depressive Disorders",
                      "No Diagnosis Given")
  hbn_dx_exclude <- c("")
  
  hbn_dx_labels <- hbn_dx %>% 
    filter(DX_01_Cat %in% hbn_dx_include) %>% 
    select(DX_01) %>% 
    distinct() %>% 
    mutate(DX = ifelse(DX_01 == "No Diagnosis Given",
                       "Control", DX_01)) 
  
} else if (dx_group == "Depression") {
  
  hbn_dx_include <- c("Depressive Disorders",
                      "No Diagnosis Given")
  hbn_dx_exclude <- c("")
  
  hbn_dx_labels <- hbn_dx %>% 
    filter(DX_01_Cat %in% hbn_dx_include) %>% 
    select(DX_01) %>% 
    distinct() %>% 
    mutate(DX = ifelse(DX_01 == "No Diagnosis Given",
                       "Control", DX_01)) 
  
} else if (dx_group == "Learning") {
  
  hbn_dx_include <- c("Specific Learning Disorder with Impairment in Mathematics",
                      "Specific Learning Disorder with Impairment in Reading",
                      "Specific Learning Disorder with Impairment in Written Expression",
                      "No Diagnosis Given")
  hbn_dx_exclude <- c("")
  
  hbn_dx_labels <- hbn_dx %>% 
    filter(DX_01 %in% hbn_dx_include) %>% 
    select(DX_01) %>% 
    distinct() %>% 
    mutate(DX = ifelse(DX_01 == "No Diagnosis Given",
                       "Control", DX_01)) 
  
} else if (dx_group == "Language") {
  
  hbn_dx_include <- c("Language Disorder",
                      "Speech Sound Disorder",
                      "Social (Pragmatic) Communication Disorder",
                      "No Diagnosis Given")
  hbn_dx_exclude <- c("")
  
  hbn_dx_labels <- hbn_dx %>% 
    filter(DX_01 %in% hbn_dx_include) %>% 
    select(DX_01) %>% 
    distinct() %>% 
    mutate(DX = ifelse(DX_01 == "No Diagnosis Given",
                       "Control", DX_01)) 
  
} else {
  stop()
}

# Select DX columns and rows
# hbn_dx <- hbn_dx %>% 
#   select(Subject_ID = EID, 
#          DX_01_Cat,
#          DX_01) %>% 
#   filter(DX_01_Cat %in% hbn_dx_include,
#          !(DX_01 %in% hbn_dx_exclude),
#          !is.na(DX_01)) %>% 
#   arrange(DX_01_Cat, DX_01)

hbn_dx <- hbn_dx %>% 
  select(Subject_ID = EID, 
         DX_01_Cat,
         DX_01) %>% 
  filter(DX_01 %in% hbn_dx_include,
         !(DX_01 %in% hbn_dx_exclude),
         !is.na(DX_01)) %>% 
  arrange(DX_01_Cat, DX_01)


# The diagnostic data frame seems to have duplicated entries.
# Remove those.
hbn_dx <- distinct(hbn_dx)

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

# Replace file NIFTY extension with MINC
hbn <- hbn %>% 
  mutate(file = str_replace(file, ".nii.gz", ".mnc"))

# Remove participants with missing demographics
hbn <- hbn %>% 
  filter(!is.na(DX), 
         !is.na(Sex),
         !is.na(Age),
         !is.na(Site),
         !is.na(Scanner))

ids_to_remove <- c("sub-NDAREK255DEE_acq-HCP_T1w",
                   "d8_0033_01",
                   "d8_0034_02",
                   "d8_0035_01",
                   "d8_0036_01",
                   "d8_0037_01",
                   "d8_0038_01",
                   "d8_0039_01",
                   "d8_0040_04",
                   "d8_0041_03",
                   "d8_0042_01",
                   "d8_0043_01",
                   "d8_0044_04")

hbn <- hbn %>% 
  filter(!(Subject_ID %in% ids_to_remove))

# Arrange columns and rows
hbn <- hbn %>% 
  select(Subject_ID, Scan_ID, Dataset,
         Age, Sex, DX, 
         Site, Scanner,
         file) %>% 
  arrange(Dataset, Subject_ID)

# Export
outfile <- paste("demographics", "HBN", dx_group, sep = "_")
outfile <- paste0(outfile, ".csv")
outfile <- file.path(demographics_dir_v3, outfile)
write_csv(x = hbn, file = outfile)
