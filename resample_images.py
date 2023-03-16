#!.venv/bin/python3
# ----------------------------------------------------------------------------
# resample_images.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2022

"""
Resample a set of images using the autocrop command line tool.
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
        help = "Path to directory containing images to resample."
    )

    parser.add_argument(
        '--outdir',
        type = str,
        help = "Path to directory in which to save resampled images."
    )

    parser.add_argument(
        '--isostep',
        type = str,
        help = "Resolution of voxels in resampled images (mm)."
    )

    parser.add_argument(
        '--parallel',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = "Option to run in parallel."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = None,
        help = ("Number of processors to use in parallel. "
                "Ignored if --parallel set to false.")
    )

    return vars(parser.parse_args())


# Main -----------------------------------------------------------------------
if __name__ == '__main__':

    # Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outdir = args['outdir']
    isostep = args['isostep']
    parallel = True if args['parallel'] == 'true' else False
    nproc = args['nproc']

    # Parallel checks
    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified "
                            "when --parallel true")

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

    outfiles = utils.resample_images(infiles = imgfiles,
                                     isostep = isostep,
                                     outdir = outdir,
                                     parallel = parallel,
                                     nproc = nproc)
