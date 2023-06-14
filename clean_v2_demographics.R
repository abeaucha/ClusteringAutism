# ----------------------------------------------------------------------------
# clean_v2_demographics.R
# Author: Antoine Beauchamp
#
# Clean human demographics information for v2 registration run.
#
# Description
# -----------
# This script cleans and consolidates the demographic information for the
# version 2 registration run. POND and SickKids demographics and scan
# information were stored in separate data files. This program consolidates
# the various input files into a single demographics file "demographics.csv"

# Packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))

# Demographics directory
demographics_dir <- "data/human/registration/v2/subject_info/"

# File suffix for the Jacobian images
file_suffix <- ".denoise.extracted_fwhm_2mm.mnc"

# Clean up SickKids demographics
sickkids <- "Sickkids_dbm_input_demo_n500.csv"
sickkids <- file.path(demographics_dir, sickkids)
sickkids <- read_csv(sickkids, show_col_types = FALSE) %>% 
  select(Subject_ID, Scan_ID,
         Rating_Motion, Rating_IN3, QC_Status,
         Site, Scanner,
         Age, DX, Sex) %>%
  mutate(Dataset = "SickKids",
         file = str_c(Subject_ID, file_suffix))

# Clean up POND demographics
# Get registration information
pond_dbm <- "POND_release2_motion_in3_rescored_baseline_norepeats_pass_dbm_input_n791.csv"
pond_dbm <- file.path(demographics_dir, pond_dbm)
pond_dbm <- read_csv(pond_dbm, show_col_types = FALSE) %>% 
  select(Subject_ID, Scan_ID,
         Rating_Motion, Rating_IN3, QC_Status) %>%
  mutate(Dataset = "POND")

# Get demographics information
pond_metadata <- "/projects/POND/POND_BIDS/workflows/pond-metadata-20230111.xlsx"
pond_metadata <- read_excel(pond_metadata) %>% 
  mutate(Subject_ID = str_c("POND_", Subject)) %>% 
  select(Subject_ID,
         Sex = NSI_SEX,
         DX = PRIMARY_DIAGNOSIS) %>% 
  mutate(DX = ifelse(DX == "Typically Developing", "Control", DX),
         DX = ifelse(DX == "Unaffected Sibling", "Control", DX))

# Join registration info with demographics
pond <- pond_dbm %>% 
  left_join(pond_metadata, by = "Subject_ID")

# POND site and scanner information
pond_sites <- tibble(Scanner = c("SickKids_Trio", "SickKids_Prisma", "Queens_Trio", "Queens_Prisma", "Bloorview_Prisma"),
                     Site = c("Toronto", "Toronto", "Queens", "Queens", "Bloorview"))

# Get POND scanner information
pond_neuro <- "pond-neuroanatomy20230111.csv"
pond_neuro <- file.path(demographics_dir, pond_neuro)
pond_neuro <- read_csv(pond_neuro, show_col_types = FALSE) %>% 
  select(Subject_ID = subject, Scan_ID = scan,
         Scanner = scanned_on, Age = age_at_scan) %>%
  mutate(Subject_ID = str_replace(Subject_ID, "sub-", "POND_"),
         Scan_ID = str_c(Scan_ID, ".mnc")) %>% 
  left_join(pond_sites, by = "Scanner")

# Combine POND demographics
pond <- pond %>% 
  left_join(pond_neuro, by = c("Subject_ID", "Scan_ID")) %>% 
  mutate(file = str_replace(Scan_ID, ".mnc", file_suffix))

# Combined demographics
demographics <- bind_rows(pond, sickkids)

# Remove NAs
demographics <- demographics %>%
  filter(!is.na(DX), 
         !is.na(Sex),
         !is.na(Age),
         !is.na(Site),
         !is.na(Scanner))

# Patients to remove based on the QC README
to_remove <- c("sub-0880102",
               "sub-0881247",
               "sub-0881263",
               "sub-0881263",
               "sub-1050027",
               "sub-1050452",
               "sub-0881317",
               "sub-1050959",
               "sub-1050959")

# Remove patients
demographics <- demographics %>% 
  mutate(ID_temp = str_remove(Scan_ID, "_ses.*.mnc")) %>% 
  filter(!(ID_temp %in% to_remove)) %>% 
  select(-ID_temp)

# Export
outfile <- "demographics.csv"
outfile <- file.path(demographics_dir, outfile)
write_csv(x = demographics, file = outfile)
