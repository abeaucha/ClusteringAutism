suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))

nvoxels <- 10000
maskfile <- "data/human/registration/v2/reference_files/mask_0.8mm.mnc"
mask <- mincGetVolume(maskfile)

ind_mask <- mask > 0.5

ind_nvoxels <- which(ind_mask)[1:nvoxels]

mask_partial <- integer(length(mask))
mask_partial[ind_nvoxels] <- 1
attributes(mask_partial) <- attributes(mask)

outfile <- basename(maskfile) %>%
  str_remove(".mnc") %>%
  str_c("partial", "nvox", nvoxels, sep = "_") %>%
  str_c(".mnc")
outfile <- file.path(dirname(maskfile), outfile)

mincWriteVolume(
  buffer = mask_partial,
  output.filename = outfile,
  like.filename = maskfile,
  clobber = TRUE
)
