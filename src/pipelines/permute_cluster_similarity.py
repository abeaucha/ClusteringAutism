#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# permute_cluster_similarity.py
# Author: Antoine Beauchamp
# Created: April 25th, 2024

"""

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import shutil
import sys
import utils
import numpy as np
import pandas as pd
from itertools import product
import tempfile as tf


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
        help = "Strings indicating which species are being compared."
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
        nargs = 2,
        type = str,
        help = "Paths to the mask image files (.mnc)."
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
        '--jacobians',
        nargs = '*',
        type = str,
        default = ['absolute', 'relative'],
        help = "Jacobians to use."
    )

    parser.add_argument(
        '--execution',
        type = str,
        default = 'local',
        choices = ['local', 'slurm'],
        help = ("Flag indicating whether the pipeline should be executed "
                "or using the Slurm scheduler on a HPC cluster.")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use."
    )

    parser.add_argument(
        '--registry-name',
        type = str,
        default = "compute_cluster_similarity_registry",
        help = "Name of the registry directory for batched jobs."
    )

    parser.add_argument(
        '--registry-cleanup',
        type = str,
        default = "true",
        help = "Option to clean up registry after completion of batched jobs."
    )

    parser.add_argument(
        '--slurm-njobs',
        type = int,
        help = "Number of jobs to deploy on Slurm."
    )

    parser.add_argument(
        '--slurm-mem',
        type = str,
        help = "Memory per CPU."
    )

    parser.add_argument(
        '--slurm-time',
        type = str,
        help = "Walltime in hh:mm:ss format for Slurm jobs."
    )

    return vars(parser.parse_args())