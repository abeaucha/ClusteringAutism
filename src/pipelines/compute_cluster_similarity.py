#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# compute_cluster_similarity.py
# Author: Antoine Beauchamp
# Created: March 27th, 2024

"""
Pipeline to evaluate pairwise cluster similarity.

Description
-----------
This pipeline evaluates the pairwise gene expression similarity of mouse and
human clusters on the basis of the cluster centroid images. The centroid images
are obtained by specifying the mouse and human processing pipeline directories
as well as parameter set IDs in the command line arguments. Additional paths
must be specified to the mouse and human gene expression data directories.

Parameters governing the evaluation of the similarity computations are mappable.
"""

# Packages -------------------------------------------------------------------

import argparse


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/cross_species/v3/',
        help = "Path to the pipeline output directory."
    )

    parser.add_argument(
        '--species',
        nargs = 2,
        type = str,
        help = ("Strings indicating which species are being compared.")
    )

    parser.add_argument(
        '--input-dirs',
        nargs = 2,
        type = str,
        help = ("Paths to the processing pipeline directories containing "
                "images to compare.")
    )

    parser.add_argument(
        '--param-ids',
        nargs = 2,
        type = str,
        help = ("IDs specifying the processing pipeline parameter sets "
                "to use.")
    )

    parser.add_argument(
        '--expr-dirs',
        nargs = 2,
        type = str,
        help = ("Paths to gene expression directories for the species "
                "being compared.")
    )

    parser.add_argument(
        '--masks',
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
        choices = ['none', 'top_n', 'intensity'],
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
        '--jacobians',
        nargs = '*',
        type = str,
        default = ['absolute', 'relative'],
        help = "Jacobians to use."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use in parallel."
    )

    return vars(parser.parse_args())


# Modules --------------------------------------------------------------------


def main(pipeline_dir, species, input_dirs, param_ids, expr_dirs, masks,
         microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_study_v3.csv',
         gene_space = 'average-latent-space',
         n_latent_spaces = 100, latent_space_id = 1,
         metric = 'correlation', signed = True,
         threshold = 'top_n', threshold_value = 0.2,
         threshold_symmetric = True,
         jacobians = ('absolute', 'relative'),
         nproc = 1):
    pass


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    print(args)
    # main(**args)
