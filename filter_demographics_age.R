
suppressPackageStartupMessages(library(tidyverse))

demographics <- "data/human/registration/v3/subject_info/demographics.csv"
demographics <- read_csv(demographics, show_col_types = FALSE)

age_min <- 3
age_max <- 21

demographics_filt <-demographics %>%
  filter(Age >= age_min, Age <= age_max)

outfile <- "demographics_v3.1_age_filter.csv"
outfile <- file.path("data/human/registration/v3/subject_info/", outfile)
write_csv(x = demographics_filt, file = outfile)