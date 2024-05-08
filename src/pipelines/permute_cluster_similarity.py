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
import random
import os
import sys
import utils
import numpy as np
import pandas as pd
from process_human_images import centroids

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
        '--param-id',
        type = str,
        help = ("IDs specifying the similarity pipeline parameter set "
                "to use.")
    )

    parser.add_argument(
        '--input-dirs',
        nargs = 2,
        type = str,
        help = ("Paths to the processing pipeline directories containing "
                "images to compare.")
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
        default = 'data/human/expression/v3/AHBA_microarray_coordinates_study.csv',
        help = ("Path to file (.csv) containing the world coordinates of "
                "the AHBA microarray samples.")
    )

    parser.add_argument(
        '--permutations-n',
        type = int,
        default = 100,
        help = "Number of permutations to run."
    )

    parser.add_argument(
        '--permutations-start',
        type = int,
        default = 1,
        help = "Permutation to start on."
    )

    parser.add_argument(
        '--keep-centroids',
        type = str,
        default = 'false',
        help = "Option to keep permuted centroids."
    )

    # parser.add_argument(
    #     '--nproc',
    #     type = int,
    #     default = 1,
    #     help = "Number of processors to use."
    # )
    #
    # parser.add_argument(
    #     '--registry-name',
    #     type = str,
    #     default = "compute_cluster_similarity_registry",
    #     help = "Name of the registry directory for batched jobs."
    # )
    #
    # parser.add_argument(
    #     '--registry-cleanup',
    #     type = str,
    #     default = "true",
    #     help = "Option to clean up registry after completion of batched jobs."
    # )
    #
    # parser.add_argument(
    #     '--slurm-njobs',
    #     type = int,
    #     help = "Number of jobs to deploy on Slurm."
    # )
    #
    # parser.add_argument(
    #     '--slurm-mem',
    #     type = str,
    #     help = "Memory per CPU."
    # )
    #
    # parser.add_argument(
    #     '--slurm-time',
    #     type = str,
    #     help = "Walltime in hh:mm:ss format for Slurm jobs."
    # )

    return vars(parser.parse_args())


# Modules --------------------------------------------------------------------


def initialize(**kwargs):
    # Ensure proper paths
    pipeline_dir = os.path.join(kwargs['pipeline_dir'], '')
    input_dirs = [os.path.join(path, '') for path in kwargs['input_dirs']]
    expr_dirs = [os.path.join(path, '') for path in kwargs['expr_dirs']]

    # Fetch similarity pipeline parameters
    param_id = kwargs['param_id']
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    if not os.path.exists(metadata):
        raise OSError("Input pipeline metadata file not found: {}"
                      .format(metadata))
    params = utils.fetch_params_metadata(metadata = metadata,
                                         id = param_id)
    params = {key:val[0] for key, val in params.to_dict(orient = 'list').items()}

    # Build paths to input directories
    input_ids = (params['input_1_id'], params['input_2_id'])

    inputs = dict(effect_sizes = [], clusters = [], centroids = [])
    for i, x in enumerate(zip(input_dirs, input_ids)):
        metadata_i = os.path.join(x[0], 'metadata.csv')
        params_i = utils.fetch_params_metadata(metadata_i, id = x[1])

        res_i = params_i['resolution'][0]
        cluster_res_i = params_i['cluster_resolution'][0]

        res_i_str = 'resolution_{}'.format(res_i)
        cluster_res_i_str = 'resolution_{}'.format(cluster_res_i)

        input_dir_i = os.path.join(x[0], x[1], '')
        es_dir_i = os.path.join(input_dir_i, 'effect_sizes', res_i_str)
        cluster_dir_i = os.path.join(input_dir_i, 'clusters', cluster_res_i_str)
        centroid_dir_i = os.path.join(input_dir_i, 'centroids', res_i_str)

        inputs['effect_sizes'].append(es_dir_i)
        inputs['clusters'].append(cluster_dir_i)
        inputs['centroids'].append(centroid_dir_i)

    inputs = {key:tuple(val) for key, val in inputs.items()}

    # Pipeline directory
    pipeline_dir = os.path.join(pipeline_dir, param_id, 'permutations', '')
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)

    # Define pipeline sub-directories
    paths = dict(
        clusters = os.path.join(pipeline_dir, 'clusters', ''),
        centroids = os.path.join(pipeline_dir, 'centroids', ''),
        similarity = os.path.join(pipeline_dir, 'similarity', '')
    )

    # Create pipeline sub-directories
    for path in paths.values():
        if not os.path.exists(path):
            os.makedirs(path)

    return inputs, paths


def permute_cluster_labels(clusters, outdir, n = 100, start = 1,
                           min_per_k = None):
    outdir = os.path.join(outdir, '')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    df_clusters = pd.read_csv(clusters)
    cols = df_clusters.drop('ID', axis = 1).columns

    min_per_k = 0 if min_per_k is None else min_per_k

    outfiles = []
    for p in range(start, start + n):

        df_permute = df_clusters.copy()
        for col in cols:
            k = df_clusters[col].to_numpy()
            k_vals, k_freq = np.unique(k, return_counts = True)
            k_lt_min = k_vals[k_freq < min_per_k]
            lt_min = np.isin(k, k_lt_min)

            random.seed(p)
            np.random.seed(p)
            k[~lt_min] = np.random.choice(k[~lt_min],
                                          size = len(k[~lt_min]),
                                          replace = False)

            k[lt_min] = -1
            df_permute[col] = k

        outfile = 'clusters_permutation_{}.csv'.format(p)
        outfile = os.path.join(outdir, outfile)
        df_permute.to_csv(outfile, index = False)
        outfiles.append(outfile)

    return outfiles


def main(pipeline_dir, param_id, input_dirs, expr_dirs, masks,
         microarray_coords = 'data/human/expression/v3/AHBA_microarray_coordinates_study.csv',
         permutations_n = 100, permutations_start = 1, keep_centroids = False):
    # Get local kwargs
    kwargs = locals().copy()

    # Initialize pipeline
    print("Initializing pipeline...", flush = True)
    inputs, paths = initialize(**kwargs)

    # Create cluster permutations
    print("Generating cluster permutations...", flush = True)
    clusters = os.path.join(inputs['clusters'][0], 'clusters.csv')
    permutations = permute_cluster_labels(clusters = clusters,
                                          outdir = paths['clusters'],  # Path to pipeline cluster dir
                                          n = permutations_n,
                                          start = permutations_start)

    permutations_end = permutations_start + permutations_n
    permutations_range = range(permutations_start, permutations_end)
    for p, f in zip(permutations_range, permutations):
        print("Permutation {} of {}".format(p, permutations_end-1))

        centroid_kwargs = dict(
            clusters = f,
            imgdir = inputs['effect_sizes'][0],
            outdir = paths['centroids'],
            mask = masks[0],
            method = ...,
            execution = ...,
            nproc = ...,
            registry_name = ...,
            registry_cleanup = ...,
            slurm_mem = ...,
            slurm_time = ...
        )
        centroid_outputs = centroids(**centroid_kwargs)

    # Iterate over permutations
    # For each permutation:
    ## 1. Regenerate centroids
    ## 2. Put together list of cluster pairs on diagonal +/- n off diagonal
    ## 3. Compute cluster pair similarity

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    args['input_dirs'] = tuple(args['input_dirs'])
    args['expr_dirs'] = tuple(args['expr_dirs'])
    args['masks'] = tuple(args['masks'])
    args['keep_centroids'] = (True if args['keep_centroids'] == 'true'
                              else False)
    main(**args)
