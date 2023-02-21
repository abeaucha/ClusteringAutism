# ----------------------------------------------------------------------------
# template.R
# Author: Antoine Beauchamp
# Created:
#
# Brief description
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--arg",
              type = "character",
              help = paste("Help message")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbosity option. [default %default]"))
) 


# Functions ------------------------------------------------------------------

#' Function description
#'
#' @param x (class) Parameter description.
#'
#' @return (class) Return description.
template <- function(x){
  return()
}


# Main -----------------------------------------------------------------------

#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))
arg <- args[["arg"]]
verbose <- ifelse(args[["verbose"]] == "true", TRUE, FALSE)

if (verbose) {message("")}

