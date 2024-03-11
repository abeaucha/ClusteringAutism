#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# preprocess_human_data.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2023

"""
Preprocess human data files.

Description
-----------
"""


# Packages --------------------------------------------------------------------

import argparse
import os
import utils


# Command line arguments ------------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--imgdir',
        type = str,
        help = ("Path to base directory containing Jacobians files. "
                "Expects sub-directories 'absolute' and 'relative', "
                "each containing a sub-directory named 'smooth', "
                "which holds the images.")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use. Executed in parallel when > 1."
    )

    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Verbosity.'
    )

    args = vars(parser.parse_args())

    return args


# Modules ---------------------------------------------------------------------

def main(imgdir, nproc = 1, verbose = True):
    """
    Execute the image pre-processing pipeline.

    Parameters
    ----------
    imgdir: str
        Path to base directory containing Jacobians files.
        Expects sub-directories 'absolute' and 'relative',
        each containing a subdirectory named 'smooth',
        which holds the images.
    nproc: int, default 1
        Number of processors to use. Executed in parallel when > 1.
    verbose: bool, default True
        Verbosity.

    Returns
    -------
    None
    """

    # Check existence of image directory
    if imgdir is None:
        raise Exception("Argument --imgdir must be specified")
    else:
        imgdir = os.path.join(imgdir, '')
        if not os.path.exists(imgdir):
            OSError("Image directory not found: {}".format(imgdir))

    # Iterate over jacobians
    jacobians = ['absolute', 'relative']
    for jac in jacobians:

        if verbose:
            print("Processing {} Jacobians files...".format(jac))

        # Get input images
        input_dir = os.path.join(imgdir, jac, 'smooth', '')
        if not os.path.exists(input_dir):
            raise OSError("Input directory not found: {}".format(input_dir))
        input_files = [os.path.join(input_dir, infile)
                       for infile in os.listdir(input_dir)
                       if '.gz' in infile]

        # Output directory
        output_dir = os.path.join(imgdir, jac, 'smooth_minc', '')

        # Unzip compressed files
        if verbose:
            print("\tExtracting compressed images...")
        unzipped_files = utils.gunzip_files(infiles = input_files,
                                            keep = True,
                                            nproc = nproc)

        # Convert images from NIFTY to MINC
        if verbose:
            print("\tConverting NIFTY files to MINC...")
        imgfiles = utils.convert_images(infiles = unzipped_files,
                                        input_format = 'nifty',
                                        output_format = 'minc',
                                        outdir = output_dir,
                                        keep = True,
                                        nproc = nproc)

    print("Pipeline complete.")

    return


# Execution -------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    main(**args)
