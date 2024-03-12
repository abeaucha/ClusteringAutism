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


# Main -----------------------------------------------------------------------


## Directories ---------------------------------------------------------------

# Registration base directories
registration_dir_v2 <- "data/human/registration/v2/"
registration_dir_v3 <- "data/human/registration/v3/"

# Jacobian image directories
imgdir_v2 <- file.path(registration_dir_v2, "jacobians")
imgdir_v3 <- file.path(registration_dir_v3, "jacobians")

# Demographics directories
demographics_dir_v2 <- file.path(registration_dir_v2, "subject_info/")
demographics_dir_v3 <- file.path(registration_dir_v3, "subject_info/")

# File suffix for the Jacobian images
file_suffix <- ".denoise_fwhm_2mm.nii.gz"

# Images files and IDs
jacobians <- "absolute"
imgfiles <- imgdir_v3 %>% 
  file.path(jacobians, "smooth") %>% 
  list.files(pattern = "*.nii.gz")

# Create a data frame with image files and IDs
df_imgs <- tibble(file = imgfiles) %>% 
  mutate(Scan_ID = file %>% 
           str_remove(file_suffix))



## SickKids data -------------------------------------------------------------

demographics_v2 <- demographics_dir_v2 %>% 
  file.path("demographics.csv") %>% 
  read_csv(show_col_types = FALSE)

# Clean up SickKids demographics
sickkids <- "Sickkids_dbm_input_demo_n500.csv"
sickkids <- file.path(demographics_dir_v3, "SickKids", sickkids)
sickkids <- read_csv(sickkids, show_col_types = FALSE) %>% 
  select(Subject_ID, Scan_ID,
         Rating_Motion, Rating_IN3, QC_Status,
         Site, Scanner,
         Age, DX, Sex) %>%
  mutate(Dataset = "SickKids",
         file = str_c(Subject_ID, file_suffix))


## POND data -----------------------------------------------------------------

pond_metadata <- "pond-metadata-20240110.csv"
pond_metadata <- file.path(demographics_dir_v3, "POND", pond_metadata)
pond_metadata <- read_csv(pond_metadata) %>% 
  distinct() %>% 
  select(Subject_ID = SUBJECT,
         Sex = NSI_SEX,
         DX = PRIMARY_DIAGNOSIS) %>% 
  mutate(Subject_ID = as.character(Subject_ID),
         Subject_ID = ifelse(str_starts(Subject_ID, "88"), 
                             str_c("0", Subject_ID),
                             Subject_ID),
         Subject_ID = str_c("sub", Subject_ID, sep = "-"), 
         DX = ifelse(DX == "Typically Developing", "Control", DX),
         DX = ifelse(DX == "Unaffected Sibling", "Control", DX))

# POND site and scanner information
pond_sites <- tibble(Scanner = c("SickKids_Trio", "SickKids_Prisma", 
                                 "Queens_Trio", "Queens_Prisma", 
                                 "Bloorview_Prisma"),
                     Site = c("Toronto", "Toronto", 
                              "Queens", "Queens", 
                              "Bloorview"))

pond_neuro <- "pond-neuroanatomy-20240212.csv"
pond_neuro <- file.path(demographics_dir_v3, "POND", pond_neuro)
pond_neuro <- read_csv(pond_neuro) %>% 
  select(Subject_ID = subject, Scan_ID = scan,
         Scanner = scanned_on, Age = age_at_scan) %>% 
  left_join(pond_sites, by = "Scanner")

pond <- pond_neuro %>% 
  inner_join(df_imgs, by = "Scan_ID") %>% 
  left_join(pond_metadata, by = "Subject_ID")


# %>% 
filter(str_detect(Subject_ID, "^N")) %>% 
  mutate(Subject_ID = str_remove(Subject_ID, "_.*"))

# Join image files to basic demographics info
hbn_imgs <- hbn_imgs %>% 
  left_join(hbn_basic, by = "Subject_ID") 



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
             Release = release_number,
             Sex, Age)
  }) %>% 
  bind_rows(.id = "Site") %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male"))


# Create a data frame with image files and IDs
hbn_imgs <- tibble(file = imgfiles) %>% 
  mutate(Subject_ID = file %>% 
           str_remove(file_suffix) %>% 
           str_remove("sub-")) %>% 
  filter(str_detect(Subject_ID, "^N")) %>% 
  mutate(Subject_ID = str_remove(Subject_ID, "_.*"))

# Join image files to basic demographics info
hbn_imgs <- hbn_imgs %>% 
  left_join(hbn_basic, by = "Subject_ID") 





# 
# 
# hbn_dx <- "9994_ConsensusDx_20220728.csv"
# hbn_dx <- file.path(demographics_dir, "HBN", "assessment_data", hbn_dx)
# hbn_dx <- read_csv(hbn_dx)
# 
# hbn_dx_mv <- "9994_ConsensusDx_20220728_MV.csv"
# hbn_dx_mv <- file.path(demographics_dir, "HBN", "assessment_data", hbn_dx_mv)
# hbn_dx_mv <- read_csv(hbn_dx_mv)
# 
# 
# hbn_dx %>% 
#   select(DX_01_Cat, DX_01_Sub, DX_01) %>% 
#   distinct() %>% 
#   arrange(DX_01_Cat) %>% 
#   View()
# 
# hbn_dx %>% 
#   filter(DX_01_Cat == "Neurodevelopmental Disorders") %>% 
#   nrow()
# 
# hbn_dx %>% 
#   filter(DX_02_Cat == "Neurodevelopmental Disorders") %>% 
#   nrow()
# 
# hbn_dx %>% 
#   filter(DX_03_Cat == "Neurodevelopmental Disorders") %>% 
#   nrow()
# 
# hbn_dx %>% 
#   filter(DX_01_Cat == "Neurodevelopmental Disorders",
#          DX_02_Cat == "Neurodevelopmental Disorders") %>% 
#   select(DX_01_Cat, DX_01_Sub, DX_01, DX_02_Cat, DX_02_Sub, DX_02) %>% 
#   distinct() %>% 
#   View()
# 
# 
# hbn_dx %>% 
#   filter(is.na(DX_01_Cat))
