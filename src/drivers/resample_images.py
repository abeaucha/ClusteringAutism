#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# resample_images.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2022

"""
Resample a set of MINC images using the autocrop command line tool.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import utils


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--imgdir',
        type = str,
        help = ("Path to the directory containing the images (.mnc) to "
                "resample.")
    )

    parser.add_argument(
        '--outdir',
        type = str,
        help = ("Path to the directory in which to export the resampled "
                "images.")
    )

    parser.add_argument(
        '--isostep',
        type = str,
        help = ("Desired isotropic resolution of the resampled "
                "images (mm).")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = ("Number of processors to use. "
                "Executed in parallel when > 1.")
    )

    return vars(parser.parse_args())


# Main -----------------------------------------------------------------------

if __name__ == '__main__':

    # Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outdir = args['outdir']
    isostep = args['isostep']
    nproc = args['nproc']

    # Check isostep
    if isostep is None:
        raise Exception("Argument --isostep must be specified")

    # Ensure proper paths
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')

    # Create outdir if needed
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Get images in dir
    imgfiles = [os.path.join(imgdir, file)
                for file in os.listdir(imgdir)
                if '.mnc' in file]

    # Resample images
    outfiles = utils.resample_images(infiles = imgfiles,
                                     isostep = isostep,
                                     outdir = outdir,
                                     nproc = nproc)
