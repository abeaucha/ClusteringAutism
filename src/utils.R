suppressPackageStartupMessages(library(tidyverse))


fetch_metadata <- function(metadata, ...) {
  kwargs <- list(...)
  
  if (!file.exists(metadata)) {
    stop(paste("Metadata file not found:", metadata))
  }
  
  df_metadata <- read_csv(metadata, show_col_types = FALSE)
  
  if (length(kwargs) == 0) {
    message("No parameters provided. Returning all metadata.")
    return(df_metadata)
  } else {
    df_params = as_tibble(kwargs)
  }
  
  df_match <- inner_join(df_metadata, df_params, by = colnames(df_params))
  nmatch <- nrow(df_match)
  if (nmatch == 0) {
    return(NULL)
  } else {
    return(df_match)
  }
}