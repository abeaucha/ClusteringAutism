#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# clean_mouse_pipeline_outputs.py
# Author: Antoine Beauchamp
# Created: April 30th, 2025

"""
Description
-----------
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

    # General arguments ---------------------------------------------------
    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/mouse/derivatives/',
        help = "Path to output directory."
    )

    parser.add_argument(
        '--resolution',
        type = int,
        default = 200,
        help = "Resolution (um) of images to link."
    )

    parser.add_argument(
        '--method',
        type = str,
        default = 'mean',
        help = "Method used to aggregate cluster maps."
    )

    parser.add_argument(
        '--transform',
        type = str,
        default = 'data/mouse/registration/transforms/MICe_scanbase.xfm',
        help = "Path to transform file (.xfm)."
    )

    parser.add_argument(
        '--transform-like',
        type = str,
        default = 'data/mouse/atlas/average_template_200um.mnc',
        help = "Path to transform likefile (.mnc)."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use."
    )

    args = vars(parser.parse_args())

    return args


# Modules --------------------------------------------------------------------

def main(pipeline_dir, resolution = 200,
         method = 'mean', transform = None,
         transform_like = None, nproc = 1):
    """
    Fetch mouse clustering outputs from Jacob's directory.

    Arguments
    ---------
    pipeline_dir: str
        Path to output directory.
    resolution: float
        Resolution (um) of images to link.
    method: str
        Method used to aggregate cluster maps.
    transform: str
        Path to transform file (.xfm) to use.
    transform_like: str
        Path to transform likefile (.mnc).
    nproc: int
        Number of processors to use in parallel.

    Returns
    -------
    None
    """

    # Create pipeline directories --------------------------------------------

    # Check existence of input directory
    if not os.path.exists(pipeline_dir):
        raise OSError("Pipeline directory not found: {}".format(pipeline_dir))

    # Path to cluster maps directory
    resolution_mm = float(resolution) / 1000

    # Path to clustering directory
    cluster_dir = os.path.join(pipeline_dir, 'clusters', 'resolution_{}'.format(resolution_mm), '')
    centroid_dir = os.path.join(pipeline_dir, 'centroids', 'resolution_{}'.format(resolution_mm), '')


    # Rename cluster files ---------------------------------------------------

    # Rename clusters CSV
    os.rename(os.path.join(cluster_dir, "Clusters.csv"),
              os.path.join(cluster_dir, "clusters.csv"))

    # Rename affinity matrix file
    os.rename(os.path.join(cluster_dir, "WMatrix.RData"),
              os.path.join(cluster_dir, "affinity.RData"))


    # Organize centroid images ------------------------------------------------

    # Iterate over jacobians
    jacobians = ('absolute', 'relative')
    for jtype in jacobians:

        print("Cleaning {} centroid images...".format(jtype))

        # Output centroid directory
        centroid_dir_j = os.path.join(centroid_dir, jtype, '')
        if not os.path.exists(centroid_dir_j):
            os.makedirs(centroid_dir_j)

        # Get input files
        jtype_short = jtype[:3]
        input_files = os.listdir(centroid_dir)
        input_files = [os.path.join(centroid_dir, file) for file in input_files]
        input_files = [file for file in input_files if '.mnc' in file]
        input_files = [file for file in input_files if str(resolution) in file]
        input_files = [file for file in input_files if jtype_short in file]
        input_files = [file for file in input_files if method in file]
        if len(input_files) == 0:
            raise OSError("Input files not found for specified parameters.")

        # Transform images if needed
        if transform is None:
            for file in input_files:
                os.rename(file, os.path.join(centroid_dir_j, os.path.basename(file)))
            outfiles = [os.path.join(centroid_dir_j, os.path.basename(file))
                        for file in input_files]
        else:
            if transform_like is None:
                raise Exception("Argument transform_like must be specified "
                                "when using transform.")
            outfiles = utils.transform_images(infiles = input_files,
                                              outdir = centroid_dir_j,
                                              like = transform_like,
                                              transform = transform,
                                              nproc = nproc)

        # Update file names
        for file in outfiles:
            file_split = os.path.basename(file).split('_')
            nk = file_split[3]
            k = file_split[1]
            file_new = 'centroid_nk_{}_k_{}.mnc'.format(nk, k)
            file_new = os.path.join(centroid_dir_j, file_new)
            os.rename(file, file_new)

    return


# Execution -------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    if args['transform'] == 'none':
        args['transform'] = None
        args['transform_like'] = None
    main(**args)
