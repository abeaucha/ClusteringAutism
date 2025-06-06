#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# compute_image_similarity.py
# Author: Antoine Beauchamp
# Created: April 2nd, 2024

"""

Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import pandas as pd
import sys
import utils
from transcriptomic import transcriptomic_similarity


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--output-file',
        type = str,
        default = 'similarity.csv',
        help = ("Path to the file (.csv) in which to export the similarity "
                "values.")
    )

    parser.add_argument(
        '--input-file',
        type = str,
        help = ("Path to the file (.csv) containing image pairs for ",
                "which to evaluate the similarity.")
    )

    parser.add_argument(
        '--species',
        nargs = 2,
        type = str,
        help = "Strings indicating which species are being compared."
    )

    parser.add_argument(
        '--expr',
        nargs = 2,
        type = str,
        help = ("Paths to gene expression directories for the species "
                "being compared.")
    )

    parser.add_argument(
        '--masks',
        nargs = 2,
        type = str,
        help = ("Paths to the mask files (.mnc).")
    )

    parser.add_argument(
        '--microarray-coords',
        type = str,
        default = 'data/human/expression/AHBA_microarray_coordinates_study_v3.csv',
        help = ("Path to file (.csv) containing the world coordinates of "
                "the AHBA microarray samples.")
    )

    parser.add_argument(
        '--gene-space',
        type = str,
        default = 'average-latent-space',
        choices = ['average-latent-space', 'latent-space', 'homologous-genes'],
        help = ("The gene expression space to use to evaluate the similarity "
                "between clusters.")
    )

    parser.add_argument(
        '--n-latent-spaces',
        type = int,
        default = 100,
        help = ("Number of latent spaces to include when --gene-space is "
                "'average-latent-space'. Ignored otherwise.")
    )

    parser.add_argument(
        '--latent-space-id',
        type = int,
        default = 1,
        help = ("Latent space to use when --gene-space is 'latent-space'. "
                "Ignored otherwise.")
    )

    parser.add_argument(
        '--metric',
        type = str,
        default = 'correlation',
        help = ("The metric used to evaluate the similarity between cluster "
                "gene expression signature.")
    )

    parser.add_argument(
        '--signed',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to compute positive and negative similarity "
                "separately before averaging for a final value.")
    )

    parser.add_argument(
        '--threshold',
        type = str,
        default = 'top_n',
        choices = ['top_n', 'intensity', 'none'],
        help = ("Method used to threshold mouse and human cluster centroid "
                "images prior to constructing gene expression signatures.")
    )

    parser.add_argument(
        '--threshold-value',
        type = float,
        default = 0.2,
        help = ("Value used to threshold mouse and human images. Ignored if "
                "--threshold is 'none'")
    )

    parser.add_argument(
        '--threshold-symmetric',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to apply the threshold symmetrically to positive and "
                "negative image values.")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use."
    )

    return vars(parser.parse_args())


# Execution ------------------------------------------------------------------
if __name__ == '__main__':

    # Parse command line args
    kwargs = parse_args()
    kwargs['signed'] = True if kwargs['signed'] == 'true' else False
    kwargs['threshold_symmetric'] = (True if kwargs['threshold_symmetric'] == 'true'
                                     else False)

    # Import image pairs
    if kwargs['input_file'] is None: raise ValueError("Must provide input file")
    df_imgs = pd.read_csv(kwargs['input_file'])
    imgs = df_imgs.values.tolist()
    imgs = [tuple(x) for x in imgs]
    del kwargs['input_file']

    # Output file
    output_file = kwargs['output_file']
    del kwargs['output_file']

    # Extract kwargs for module
    kwargs['imgs'] = imgs

    # Compute pairwise similarity between cluster centroids
    # with utils.catchtime() as t:
    results = transcriptomic_similarity(**kwargs)
    results.to_csv(output_file, index = False)
    # print(f'Time: {t():.3f} seconds')
