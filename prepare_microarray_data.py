#!.venv/bin/python3
# ----------------------------------------------------------------------------
# prepare_microarray_data.py
# Author: Antoine Beauchamp

"""
Prepare the AHBA microarray sample coordinates for the similarity pipelines.

Description
-----------
This script transforms the AHBA microarray coordinates from the
MNI ICBM NLIN SYM 09c input space to the consensus average space of the
study. It also creates label and mask image in the study space for the
microarray samples.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import transcriptomic
from pyminc.volumes.factory import volumeFromFile


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/human/expression/v3/',
        help = "Path to the expression data directory."
    )

    parser.add_argument(
        '--transforms',
        nargs = '*',
        type = str,
        help = ("List of paths to the transforms from the registration ",
                "consensus average space to the MNI ICBM 152 NLIN 09c space.")
    )

    parser.add_argument(
        '--template',
        type = str,
        help = "Path to the consensus average template image (.mnc)."
    )

    parser.add_argument(
        '--metadata',
        type = str,
        default = 'data/human/expression/SampleInformation_pipeline_abagen.csv',
        help = ("Name of the metadata file (.csv) for the AHBA microarray ",
                "samples.")
    )

    parser.add_argument(
        '--annotations',
        type = str,
        default = 'data/human/expression/AHBA_microarray_sample_annotations.csv'
    )

    return vars(parser.parse_args())


# Modules --------------------------------------------------------------------


# # Directories
# expr_dir = 'data/human/expression'
# registration_dir = 'data/human/registration/'
# # version = 'v1'
# version = 'v2'
#
# # Expression inputs
# metadata = 'SampleInformation_pipeline_abagen.csv'
# annotations = 'AHBA_microarray_sample_annotations.csv'
#
# # Registration inputs
# template = 'model_0.8mm.mnc'
# # transforms = ['average_to_icbm_nlin_sym_09c0GenericAffine.mat',
# #               'average_to_icbm_nlin_sym_09c1Warp.nii.gz']
# transforms = ['to_target_0GenericAffine.mat',
#               'to_target_1Warp.nii']

def main(
    pipeline_dir, transforms, template,
    metadata = 'data/human/expression/SampleInformation_pipeline_abagen.csv',
    annotations = 'data/human/expression/AHBA_microarray_sample_annotations.csv'
):

    # Fetch microarray coordinates and transform to study space
    print("Fetching microarray coordinates...")
    coords = transcriptomic.prepare_microarray_coordinates(
        metadata = metadata,
        annotations = annotations,
        transforms = tuple(transforms)
    )

    # Rename coordinates file
    coords_new = coords.replace('.csv', '_{}.csv'.format(version))
    os.rename(coords, coords_new)
    coords = coords_new

    # Get template resolution
    vol = volumeFromFile(template)
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Create label image from microarray samples
    print("Creating label image...")
    labels = 'AHBA_microarray_labels_study_{}_{}mm.mnc'.format(version,
                                                               resolution)
    labels = os.path.join(expr_dir, labels)
    labels = transcriptomic.coordinates_to_minc(coordinates = coords,
                                                template = template,
                                                outfile = labels,
                                                type = 'labels')

    # Create mask image from microarray samples
    print("Creating mask image...")
    mask = 'AHBA_microarray_mask_study_{}_{}mm.mnc'.format(version, resolution)
    mask = os.path.join(expr_dir, mask)
    mask = transcriptomic.coordinates_to_minc(coordinates = coords,
                                              template = template,
                                              outfile = mask,
                                              type = 'mask')


# Main -----------------------------------------------------------------------
if __name__ == '__main__':

    args = parse_args()
    main()

