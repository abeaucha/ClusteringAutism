suppressPackageStartupMessages(library(RMINC))

imgdir <- "data/human/registration/reference_files/"
maskfile <- "mask_0.8mm.mnc"
maskfile <- file.path(imgdir, maskfile)

mask <- mincGetVolume(maskfile)
mask <- floor(mask)
mask_ind <- which(mask == 1)
mask_size <- length(mask_ind)

frac <- 0.75
mask_frac <- ceiling(frac*mask_size)
print(mask_frac)

mask_partial <- numeric(length(mask))
mask_partial[mask_ind[1:mask_frac]] <- 1
attributes(mask_partial) <- attributes(mask)

outfile <- paste0("mask_0.8mm_", frac, ".mnc")
outfile <- file.path(imgdir, outfile)

mincWriteVolume(buffer = mask_partial,
                output.filename = outfile,
                like.filename = maskfile,
                clobber = TRUE)
