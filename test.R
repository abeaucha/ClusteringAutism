rm(list = ls()); gc()
library(tidyverse)
library(parallel)

source("src/processing.R")

input <- "data/test/human/derivatives/v2/620/effect_sizes/resolution_3.0/absolute/effect_sizes.csv"

df <- as_tibble(data.table::fread(input, header = TRUE)) 

outfiles <- df[["file"]]

x <- df %>% 
  select(-file) %>% 
  as.matrix() %>% 
  t()

rm(list = "df"); gc()

outfiles <- file.path("test", outfiles)
margin <- 2
mask <- "data/human/registration/v2/reference_files/mask_0.8mm.mnc"
maskvol <- mincGetVolume(mask)
sum(maskvol)

#Check output files
if (is.null(outfiles)) {
  stop("Specify output files.")
}

#Check mask
if (is.null(mask)) {
  stop("Specify mask file.")
}

#Check that x is a matrix
if (!is.matrix(x)) {
  stop("x must be a matrix.")  
}

#Split matrix into array along margin
x <- asplit(x, MARGIN = margin)

#Check that number of output files matches the number of images
if (length(x) != length(outfiles)) {
  stop("Number of entries in x along margin ", margin, 
       " must be equal to the number of entries in outfiles")  
}

nproc <- 16
out <- parallel::mcmapply(vector_to_image, x, outfiles,
                          MoreArgs = list(mask = mask),
                          SIMPLIFY = TRUE, mc.cores = nproc)
cl <- makeCluster(nproc)
clusterEvalQ(cl, library(RMINC))
out <- clusterMap(cl = cl, fun = vector_to_image, x, outfiles, MoreArgs = list(mask = mask))
stopCluster(cl)

system("rm test/*")


nc <- c(2, 4, 8, 16)
df_sock <- tibble(nc = nc,
                  delta = c(0.5, 0.97, 1.91, 3.86),
                  type = "sock")

df_fork <- tibble(nc = nc, 
                  delta = c(1.21, 2.13, 4.61, 6.83),
                  type = "fork")

df_mem <- bind_rows(df_sock, df_fork)

ggplot(df_mem, aes(x = nc, y = delta, col = type)) + 
  geom_line() + 
  geom_point()
