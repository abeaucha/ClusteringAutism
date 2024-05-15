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
from compute_cluster_similarity import generate_cluster_pairs


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
        default = ['data/human/derivatives/v3/', 'data/mouse/derivatives/v3/'],
        help = ("Paths to the processing pipeline directories containing "
                "images to compare.")
    )

    parser.add_argument(
        '--expr-dirs',
        nargs = 2,
        type = str,
        default = ['data/human/expression', 'data/mouse/expression'],
        help = ("Paths to gene expression directories for the species "
                "being compared.")
    )

    parser.add_argument(
        '--masks',
        nargs = 2,
        type = str,
        default = ['data/human/registration/v3/reference_files/mask_0.8mm.mnc',
                   'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc'],
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
        '--off-diagonal',
        type = int,
        default = 1,
        help = "Number of off diagonal elements to evaluate similarity"
    )

    parser.add_argument(
        '--keep-centroids',
        type = str,
        default = 'false',
        help = "Option to keep permuted centroids."
    )

    parser.add_argument(
        '--execution',
        type = str,
        default = 'local',
        choices = ['local', 'slurm'],
        help = ("Flag indicating whether the pipeline should be executed "
                "or using the Slurm scheduler on a HPC cluster.")
    )

    # parser.add_argument(
    #     '--nproc',
    #     type = int,
    #     default = 1,
    #     help = "Number of processors to use."
    # )

    parser.add_argument(
        '--registry-name',
        type = str,
        default = "permute_cluster_similarity_registry",
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
    params['latent_space_id'] = (str(1) if params['latent_space_id'] is np.nan
                                 else params['latent_space_id'])

    # Build paths to input directories
    input_ids = (params['input_1_id'], params['input_2_id'])

    inputs = dict(effect_sizes = [], clusters = [], centroids = [])
    for i, x in enumerate(zip(input_dirs, input_ids)):
        metadata_i = os.path.join(x[0], 'metadata.csv')
        params_i = utils.fetch_params_metadata(metadata_i, id = x[1])

        res_i = params_i['resolution'][0]
        cluster_res_i = params_i['cluster_resolution'][0]
        params['input_{}_centroid_method'.format(i + 1)] = params_i['centroid_method'][0]

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

    return inputs, paths, params


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


def subset_cluster_pairs(centroid_pairs, off_diagonal = 1):
    centroid_pairs_subset = []
    for pair in centroid_pairs:
        nk = [os.path.basename(x) for x in pair]
        nk = [os.path.splitext(x)[0] for x in nk]
        nk = [x.split('_') for x in nk]
        nk = [int(x[-3]) for x in nk]
        nk_diff = abs(nk[0] - nk[1])
        if nk_diff < off_diagonal + 1:
            centroid_pairs_subset.append(pair)
    return centroid_pairs_subset


def main(pipeline_dir, param_id, input_dirs, expr_dirs, masks,
         microarray_coords = 'data/human/expression/v3/AHBA_microarray_coordinates_study.csv',
         permutations_n = 100, permutations_start = 1,
         off_diagonal = 1, keep_centroids = False,
         execution = 'local',
         registry_name = 'permute_cluster_similarity_registry',
         registry_cleanup = True, slurm_njobs = None, slurm_mem = None,
         slurm_time = None):
    # Get local kwargs
    kwargs = locals().copy()

    # Initialize pipeline
    print("Initializing pipeline...", flush = True)
    inputs, paths, params = initialize(**kwargs)

    # Create cluster permutations
    print("Generating cluster permutations...", flush = True)
    clusters = os.path.join(inputs['clusters'][0], 'clusters.csv')
    permutations = permute_cluster_labels(clusters = clusters,
                                          outdir = paths['clusters'],  # Path to pipeline cluster dir
                                          n = permutations_n,
                                          start = permutations_start)

    # Driver script and kwargs
    driver = 'transcriptomic_similarity.py'
    driver_kwargs = {key:val for key, val in params.items()
                     if 'input' not in key}
    driver_kwargs.update(dict(
        species = (params['input_1_species'], params['input_2_species']),
        expr = kwargs['expr_dirs'],
        masks = kwargs['masks'],
        microarray_coords = kwargs['microarray_coords']
    ))
    del driver_kwargs['id']

    # Iterate over permutations
    permutations_end = permutations_start + permutations_n
    permutations_range = range(permutations_start, permutations_end)
    for p, f in zip(permutations_range, permutations):
        print("Permutation {} of {}".format(p, permutations_end - 1),
              flush = True)

        outdir = os.path.join(paths['centroids'], 'permutation_{}'.format(p), '')
        centroid_kwargs = dict(
            clusters = f,
            imgdir = inputs['effect_sizes'][0],
            outdir = outdir,
            mask = masks[0],
            method = params['input_1_centroid_method'],
            execution = 'slurm',
            nproc = 8,
            registry_cleanup = False,
            slurm_mem = '16G',
            slurm_time = 60
        )
        centroid_outputs = centroids(**centroid_kwargs)
        # centroid_outputs = dict(
        #     absolute = os.path.join(outdir, 'absolute', ''),
        #     relative = os.path.join(outdir, 'relative', '')
        # )

        # Generate all pairs of centroid images
        print("Generating centroid image pairs...", flush = True)
        centroid_dirs = [outdir, inputs['centroids'][1]]
        centroid_pairs = generate_cluster_pairs(centroid_dirs = centroid_dirs)

        # Identify subset of centroid pairs on and off the nk diagonal
        centroid_pairs = subset_cluster_pairs(centroid_pairs = centroid_pairs,
                                              off_diagonal = off_diagonal)
        centroid_pairs = pd.DataFrame(centroid_pairs)

        # Execution mode
        print("Evaluating cluster similarity...", flush = True)
        if execution == 'local':

            # Export cluster pairs
            outfile = 'centroid_pairs_permutation_{}.csv'.format(p)
            outfile = os.path.join(paths['similarity'], outfile)
            centroid_pairs.to_csv(outfile, index = False)
            input_file = outfile

            # Output file for permuted similarity
            output_file = 'similarity_permutation_{}.csv'.format(p)
            output_file = os.path.join(paths['similarity'], output_file)

            # Update driver kwargs
            driver_kwargs['input-file'] = input_file
            driver_kwargs['output-file'] = output_file
            driver_kwargs = {key.replace('_', '-'):val
                             for key, val in driver_kwargs.items()}

            # Execute driver
            utils.execute_local(script = driver, kwargs = driver_kwargs)

        elif execution == 'slurm':

            # Slurm job resources
            resources = dict(
                nodes = 1,
                mem = slurm_mem,
                time = slurm_time
            )

            # Create registry
            registry = utils.Registry(resources = resources,
                                      name = '{}_{}'.format(registry_name, p))

            # Create data batches
            registry.create_batches(x = centroid_pairs,
                                    nbatches = slurm_njobs,
                                    prefix = 'centroid_pairs_batch')

            # Create jobs for batches
            kwargs['nproc'] = 1
            registry.create_jobs(script = driver,
                                 kwargs = driver_kwargs)

            # Submit jobs
            out = registry.submit_jobs(wait = True, cleanup = registry_cleanup)

            # Export results
            print("Exporting results...", flush = True)
            output_file = 'similarity_permutation_{}.csv'.format(p)
            output_file = os.path.join(paths['similarity'], output_file)
            out.to_csv(output_file, index = False)

        else:
            raise ValueError("Argument `execution` must be one of "
                             "{'local', 'slurm'}")

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    args['input_dirs'] = tuple(args['input_dirs'])
    args['expr_dirs'] = tuple(args['expr_dirs'])
    args['masks'] = tuple(args['masks'])
    args['keep_centroids'] = (True if args['keep_centroids'] == 'true'
                              else False)
    args['registry_cleanup'] = (True if args['registry_cleanup'] == 'true'
                                else False)
    main(**args)
