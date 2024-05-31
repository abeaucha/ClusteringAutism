#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# prepare_microarray_data.py
# Author: Antoine Beauchamp

"""
Prepare the AHBA microarray sample coordinates for the similarity pipelines.

Description
-----------
This script transforms the AHBA microarray sample coordinates from the
MNI ICBM NLIN SYM 09c input space to the consensus average space of the
imaging study. It also creates label and mask images in the study space
for the microarray samples.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import transcriptomic
import utils
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
        help = ("Path to the directory in which to export the microarray "
                "coordinates data.")
    )

    parser.add_argument(
        '--transforms',
        nargs = '*',
        type = str,
        help = ("List of paths to the transforms from the MNI ICBM 152 NLIN "
                "09c space to the registration consensus average space.")
    )

    parser.add_argument(
        '--template',
        type = str,
        default = 'data/human/registration/v3/reference_files/model_0.8mm.mnc',
        help = "Path to the consensus average template image (.mnc)."
    )

    parser.add_argument(
        '--metadata',
        type = str,
        default = 'data/human/expression/SampleInformation_pipeline_abagen.csv',
        help = ("Path to the file (.csv) containing the AHBA microarray "
                "sample metadata.")
    )

    parser.add_argument(
        '--annotations',
        type = str,
        default = 'data/human/expression/AHBA_microarray_sample_annotations.csv',
        help = ("Path to the file (.csv) containing AHBA microarray sample "
                "annotations indicating whether to keep or discard the "
                "sample.")
    )

    return vars(parser.parse_args())


# Modules --------------------------------------------------------------------

def main(
    pipeline_dir, transforms, template,
    metadata = 'data/human/expression/SampleInformation_pipeline_abagen.csv',
    annotations = 'data/human/expression/AHBA_microarray_sample_annotations.csv'
):
    """

    Parameters
    ----------
    pipeline_dir: str
        Path to the directory in which to export the microarray
        coordinates data.
    transforms: tuple of str
        Transforms from the MNI ICBM 152 NLIN 09c space to the
        registration consensus average space, to be passed to
        antsApplyTransformsToPoints.
    template: str
        Path to the consensus average template image (.mnc).
    metadata: str
        Path to file (.csv) containing the AHBA microarray sample
        metadata.
    annotations: str
        Path to the file (.csv) containing AHBA microarray sample
        annotations indicating whether to keep or discard the sample.

    Returns
    -------
    None
    """

    # Fetch the AHBA microarray coordinates and transform them
    # to the imaging space
    print("Fetching microarray coordinates...", flush = True)
    coords = transcriptomic.prepare_microarray_coordinates(
        outdir = pipeline_dir,
        metadata = metadata,
        annotations = annotations,
        transforms = transforms
    )

    # Fetch the image resolution
    vol = volumeFromFile(template)
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Create a label image from microarray samples
    print("Creating label image...", flush = True)
    labels = 'AHBA_microarray_labels_study_{}mm.mnc'.format(resolution)
    labels = os.path.join(pipeline_dir, labels)
    labels = utils.coordinates_to_minc(coordinates = coords,
                                       template = template,
                                       outfile = labels,
                                       type = 'labels')

    # Create a mask image from microarray samples
    print("Creating mask image...", flush = True)
    mask = 'AHBA_microarray_mask_study_{}mm.mnc'.format(resolution)
    mask = os.path.join(pipeline_dir, mask)
    mask = utils.coordinates_to_minc(coordinates = coords,
                                     template = template,
                                     outfile = mask,
                                     type = 'mask')


# Main -----------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    args['transforms'] = tuple(args['transforms'])
    main(**args)
