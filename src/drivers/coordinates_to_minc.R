#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# coordinates_to_minc.R
# Author: Antoine Beauchamp
# Created: June 23rd, 2022
#
# Convert a table of world coordinates to a MINC image.
#
# Description
# -----------
# This script takes an input table of world coordinates and uses them to build
# a MINC image. The imaging space used to interpret the world coordinates is
# specified via a template MINC image.
#
# The output image can either be a mask image indicating the location of each
# row in the coordinates table, or a label image with a unique integer value
# for each row in the table.

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option(
    "--coordinates",
    type = "character",
    help = paste("Path to the file (.csv) containing the world", "coordinates.")
  ),
  make_option(
    "--template",
    type = "character",
    help = "Path to the template image (.mnc)"
  ),
  make_option(
    "--outfile",
    type = "character",
    help = paste(
      "Path to the file (.mnc) in which to save",
      "the output image."
    )
  ),
  make_option(
    "--type",
    type = "character",
    default = "labels",
    help = paste(
      "One of {mask, labels} indicating the type of",
      "image to return. [default %default]"
    )
  ),
  make_option(
    "--verbose",
    type = "character",
    default = "true",
    help = paste("Verbosity option. [default %default]")
  )
)


# Functions ------------------------------------------------------------------

#' Convert a set of world coordinates to voxel coordinates
#'
#' @param coords (character scalar) Path to the file (.csv) containing
#' the world coordinates.
#' @param template (character scalar) Path to the template image (.mnc)
#'
#' @return (matrix) Voxel coordinates
world_to_voxel <- function(coords, template) {
  coords <- suppressMessages(read_csv(coords)) %>%
    select(label, x, y, z) %>%
    column_to_rownames("label") %>%
    as.matrix()

  coords_voxel <- mincConvertWorldMatrix(
    world_matrix = t(coords),
    file = template,
    nearest_voxel = TRUE
  )
  coords_voxel <- t(coords_voxel)
  colnames(coords_voxel) <- c("x", "y", "z")
  rownames(coords_voxel) <- rownames(coords)

  return(coords_voxel)
}


#' Convert voxel coordinates to MINC image
#'
#' @param coords (matrix) Voxel coordinates
#' @param template (character scalar) Path to the template image (.mnc)
#' @param type (character scalar) Type of image to return
#'
#' @return (mincSingleDim) The coordinates image
coords_to_img <- function(coords, template, type = "labels") {
  template <- mincArray(mincGetVolume(template))

  img <- array(data = 0, dim = dim(template))
  for (s in 1:nrow(coords)) {
    i <- coords[s, 3]
    j <- coords[s, 2]
    k <- coords[s, 1]
    if (type == "mask") {
      img[i, j, k] <- 1
    } else if (type == "labels") {
      img[i, j, k] <- as.integer(rownames(coords)[s])
    } else {
      stop(paste("Argument `type` must be one of {'mask', 'labels'}:", type))
    }
  }
  attributes(img) <- attributes(template)

  return(img)
}


# Main -----------------------------------------------------------------------

# Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
coords <- args[["coordinates"]]
template <- args[["template"]]
outfile <- args[["outfile"]]
type <- args[["type"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

# Create output directory if needed
outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Convert world coordinates to voxel coordiantes
if (verbose) {
  message("Converting world coordinates to voxel...")
}
voxel_coordinates <- world_to_voxel(coords = coords, template = template)

# Construct image from voxel coordinates
if (verbose) {
  message("Creating image from coordinates...")
}
coords_img <- coords_to_img(
  coords = voxel_coordinates,
  template = template,
  type = type
)

# Export image
if (verbose) {
  message("Writing image to file...")
}
mincWriteVolume(buffer = coords_img, output.filename = outfile, clobber = TRUE)
