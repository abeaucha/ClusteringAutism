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
import utils
import numpy as np
import pandas as pd
from compute_cluster_similarity import generate_cluster_pairs
from process_human_images import centroids
from shutil import rmtree
from transcriptomic import transcriptomic_similarity
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/cross_species/',
        help = ("Path to the directory in which to export pipeline "
                "outputs. A uniquely identified sub-directory will be "
                "created using the specified set of pipeline parameters.")
    )

    parser.add_argument(
        '--params-id',
        type = str,
        help = ("ID specifying the similarity pipeline parameter set "
                "to use.")
    )

    parser.add_argument(
        '--input-dirs',
        nargs = 2,
        type = str,
        default = ['data/human/derivatives/', 'data/mouse/derivatives/'],
        help = ("Paths to the processing pipeline directories containing "
                "the images used to build the permutations. Expects "
                "subdirectories 'effect_sizes', 'clusters', and 'centroids'.")
    )

    parser.add_argument(
        '--expr-dirs',
        nargs = 2,
        type = str,
        default = ['data/human/expression', 'data/mouse/expression'],
        help = ("Paths to the gene expression directories for the species "
                "being compared.")
    )

    parser.add_argument(
        '--masks',
        nargs = 2,
        type = str,
        default = ['data/human/registration/reference_files/mask_0.8mm.mnc',
                   'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc'],
        help = "Paths to the mask image files (.mnc)."
    )

    parser.add_argument(
        '--microarray-coords',
        type = str,
        default = 'data/human/expression/AHBA_microarray_coordinates_study.csv',
        help = ("Path to file (.csv) containing the world coordinates of "
                "the AHBA microarray samples in the human imaging study space.")
    )

    parser.add_argument(
        '--permutations-n',
        type = int,
        default = 100,
        help = "Number of permutations to evaluate."
    )

    parser.add_argument(
        '--permutations-start',
        type = int,
        default = 1,
        help = "Starting permutation ID."
    )

    parser.add_argument(
        '--permutations-ids',
        type = int,
        nargs = '*',
        help = ("Optional list of specific permutation IDs to evaluate. "
                "If passed, this overrides --permutations-n and "
                "--permutations-start.")
    )

    parser.add_argument(
        '--off-diagonal',
        type = int,
        default = 1,
        help = ("Number of off-diagonal cluster solution combinations for "
                "which to evaluate the similarity permutations.")
    )

    parser.add_argument(
        '--keep-centroids',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = "Option to keep the permuted centroid images."
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
    """
    Initialize the similarity permutation pipeline.

    Parameters
    ----------
    kwargs: dict
        All arguments passed to the main() function.

    Returns
    -------
    inputs: dict of tuple
        Dictionary with keys 'effect_sizes', 'clusters', and 'centroids'
        containing tuples of paths to the respective input directories.
    paths: dict of str
        Dictionary with keys 'clusters', 'centroids', and 'similarity'
        containing paths to the pipeline output directories.
    params: dict
        Dictionary containing pipeline parameters.
    """

    # Ensure proper paths
    pipeline_dir = os.path.join(kwargs['pipeline_dir'], '')
    input_dirs = [os.path.join(path, '') for path in kwargs['input_dirs']]
    expr_dirs = [os.path.join(path, '') for path in kwargs['expr_dirs']]

    # Fetch the similarity pipeline parameters
    params_id = kwargs['params_id']
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    if not os.path.exists(metadata):
        raise OSError("Input pipeline metadata file not found: {}"
                      .format(metadata))
    params = utils.fetch_params_metadata(metadata = metadata,
                                         id = params_id)
    params = {key:val[0] for key, val in params.to_dict(orient = 'list').items()}

    # If latent_space_id is NaN, set to 1 (subsequently unused)
    params['latent_space_id'] = (str(1) if params['latent_space_id'] is np.nan
                                 else params['latent_space_id'])

    # Build paths to the input directories
    input_ids = (params['input_1_id'], params['input_2_id'])
    inputs = dict(effect_sizes = [], clusters = [], centroids = [])
    for i, x in enumerate(zip(input_dirs, input_ids)):

        # Import input i params
        metadata_i = os.path.join(x[0], 'metadata.csv')
        params_i = utils.fetch_params_metadata(metadata_i, id = x[1])

        # Extract resolutions and centroid method
        res_i = params_i['resolution'][0]
        cluster_res_i = params_i['cluster_resolution'][0]
        params['input_{}_centroid_method'.format(i + 1)] = params_i['centroid_method'][0]

        # Resolution strings for paths
        res_i_str = 'resolution_{}'.format(res_i)
        cluster_res_i_str = 'resolution_{}'.format(cluster_res_i)

        # Input paths
        input_dir_i = os.path.join(x[0], x[1], '')
        es_dir_i = os.path.join(input_dir_i, 'effect_sizes', res_i_str)
        cluster_dir_i = os.path.join(input_dir_i, 'clusters', cluster_res_i_str)
        centroid_dir_i = os.path.join(input_dir_i, 'centroids', res_i_str)

        # Append input paths to lists in dictionary
        inputs['effect_sizes'].append(es_dir_i)
        inputs['clusters'].append(cluster_dir_i)
        inputs['centroids'].append(centroid_dir_i)

    # Convert lists to tuple
    inputs = {key:tuple(val) for key, val in inputs.items()}

    # Create pipeline directory
    pipeline_dir = os.path.join(pipeline_dir, params_id, 'permutations', '')
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


def permute_cluster_labels(clusters, outdir, n = 100, start = 1, ids = None,
                           min_per_k = None):
    """
    Permute the cluster assignment labels.

    Parameters
    ----------
    clusters: str
        Path to the file (.csv) containing cluster assignments.
    outdir: str
        Path to the output directory.
    n: int, default 100
        Number of permutations to generate.
    start: int, default 1
        Starting permutation seed.
    ids: None or list of int
        Optional list of specific permutation IDs to evaluate.
        If passed, this overrides arguments `n` and `start`.
    min_per_k: None or int
        Minimum number of patients per cluster to include in
        permutation. Clusters below threshold are ignored.

    Returns
    -------
    outfiles: list of str
        Paths to the permuted cluster assignment files (.csv).
    """

    # Create output directory if needed
    outdir = os.path.join(outdir, '')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Import cluster assignments
    df_clusters = pd.read_csv(clusters)

    # Fetch cluster assignment columns
    cols = df_clusters.drop('ID', axis = 1).columns

    # Set minimum assigned min cluster to 0 if not specified
    min_per_k = 0 if min_per_k is None else min_per_k

    # Iterate over permutations
    outfiles = []
    if ids is None:
        ids = list(range(start, start + n))
    for p in ids:

        # Iterate over clusters
        df_permute = df_clusters.copy()
        for col in cols:

            # Fetch cluster labels
            k = df_clusters[col].to_numpy()

            # Identify clusters with fewer patients than the minimum specified
            k_vals, k_freq = np.unique(k, return_counts = True)
            k_lt_min = k_vals[k_freq < min_per_k]
            lt_min = np.isin(k, k_lt_min)

            # Shuffle cluster assignments
            random.seed(p)
            np.random.seed(p)
            k[~lt_min] = np.random.choice(k[~lt_min],
                                          size = len(k[~lt_min]),
                                          replace = False)

            # Assign new cluster labels to data frame
            k[lt_min] = -1
            df_permute[col] = k

        # Export permuted cluster labels
        outfile = 'clusters_permutation_{}.csv'.format(p)
        outfile = os.path.join(outdir, outfile)
        df_permute.to_csv(outfile, index = False)
        outfiles.append(outfile)

    return outfiles


def subset_cluster_pairs(centroid_pairs, off_diagonal = 1):
    """
    Subset cluster pairs for diagonal and off-diagonal cluster solutions.

    Parameters
    ----------
    centroid_pairs: list of list of str
        List containing paths to the centroid image pairs.
    off_diagonal: int, default 1
        Number of off-diagonal cluster solution combinations for which
        to evaluate the similarity permutations.

    Returns
    -------
    centroid_pairs_subset: list of list of str
        List containing paths to the subsetted centroid image pairs.
    """

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


@utils.timing
def main(pipeline_dir, params_id, input_dirs, expr_dirs, masks,
         microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_study.csv',
         permutations_n = 100, permutations_start = 1, permutations_ids = None,
         off_diagonal = 1, keep_centroids = False,
         execution = 'local', nproc = 1, slurm_mem = None, slurm_time = None):
    """
    Execute the similarity permutation pipeline.

    Parameters
    ----------
    pipeline_dir: str
        Path to the directory in which to export pipeline outputs.
        A uniquely identified subdirectory will be created using the
        specified set of pipeline parameters.
    params_id: str
        Parameter set ID for the similarity pipeline to permute.
    input_dirs: tuple of str
        Paths to the processing pipeline directories containing
        the images used to build the permutations. Expects
        subdirectories 'effect_sizes', 'clusters', and 'centroids'.
    expr_dirs: tuple of str
        Paths to the gene expression directories for the species
        being compared.
    masks: tuple of str
        Paths to the mask image files (.mnc).
    microarray_coords: str
        Path to file (.csv) containing the world coordinates of the
        AHBA microarray samples in the human imaging study space.
    permutations_n: int, default 100
        Number of cluster permutations to generate.
    permutations_start: int, default 1
        Starting permutation ID.
    permutations_ids: list of int or None
        Optional list of specific permutation IDs to evaluate.
        If passed, this overrides arguments `n` and `start`.
    off_diagonal: int, default 1
        Number of off-diagonal cluster solution combinations for which
        to evaluate the similarity permutations.
    keep_centroids: bool, default False
        Options to keep permuted centroid images.
    execution: {'local', 'slurm'}
        Flag indicating whether the pipeline should be executed locally or
        using the SLURM scheduler on an HPC cluster.
    nproc: int, default 1
        Number of jobs to deploy.
    slurm_mem: str, default None
        Memory per CPU on Slurm. Ignored when execution = 'local'.
    slurm_time: str, default = None
        Walltime in hh:mm:ss format for Slurm jobs.
        Ignored when execution = 'local'.

    Returns
    -------
    None
    """

    # Get local kwargs
    kwargs = locals().copy()

    # Initialize the pipeline
    print("Initializing pipeline...", flush = True)
    inputs, paths, params = initialize(**kwargs)

    # Generate cluster permutations
    print("Generating cluster permutations...", flush = True)
    clusters = os.path.join(inputs['clusters'][0], 'clusters.csv')
    permutations = permute_cluster_labels(clusters = clusters,
                                          outdir = paths['clusters'],
                                          n = permutations_n,
                                          start = permutations_start,
                                          ids = permutations_ids)

    # Driver script kwargs
    driver_kwargs = {key:val for key, val in params.items()
                     if 'input' not in key}
    driver_kwargs.update(dict(
        species = (params['input_1_species'], params['input_2_species']),
        expr = kwargs['expr_dirs'],
        masks = kwargs['masks'],
        microarray_coords = kwargs['microarray_coords']
    ))
    del driver_kwargs['id']
    driver_kwargs['n_latent_spaces'] = int(driver_kwargs['n_latent_spaces'])
    driver_kwargs['latent_space_id'] = int(driver_kwargs['latent_space_id'])
    driver_kwargs['signed'] = True if driver_kwargs['signed'] == 'true' else False
    driver_kwargs['threshold_value'] = float(driver_kwargs['threshold_value'])
    driver_kwargs['threshold_symmetric'] = True if driver_kwargs['threshold_symmetric'] == 'true' else False

    # Initialize Dask client for execution
    if execution == 'local':
        cluster = LocalCluster(n_workers = nproc,
                               threads_per_worker=1)
    elif execution == 'slurm':
        cluster = SLURMCluster(
            cores = 1,
            memory = slurm_mem,
            walltime = slurm_time
        )
        cluster.scale(jobs = nproc)
    else:
        raise ValueError("Argument `--execution` must be one of {'local', 'slurm'}")

    client = Client(cluster)

    # Add client to driver kwargs
    driver_kwargs['client'] = client

    # Iterate over permutations
    # If permutation IDs given, use those. Otherwise, iterate in sequence.
    if permutations_ids is None:
        permutations_end = permutations_start + permutations_n
        permutations_range = range(permutations_start, permutations_end)
        permutations_ids = list(permutations_range)
    for p, f in zip(permutations_ids, permutations):
        print("Permutation {}".format(p), flush = True)

        # Compute permuted cluster centroid images
        print("Generating permuted centroids...", flush = True)
        outdir = os.path.join(paths['centroids'], 'permutation_{}'.format(p), '')
        # centroid_kwargs = dict(
        #     clusters = f,
        #     imgdir = inputs['effect_sizes'][0],
        #     outdir = outdir,
        #     mask = masks[0],
        #     method = params['input_1_centroid_method'],
        #     execution = 'slurm',
        #     nproc = 8,
        #     registry_cleanup = True,
        #     registry_name = 'permutation_{}'.format(p),
        #     slurm_mem = '16G',
        #     slurm_time = 120
        # )
        centroid_kwargs = dict(
            clusters = f,
            imgdir = inputs['effect_sizes'][0],
            outdir = outdir,
            mask = masks[0],
            method = params['input_1_centroid_method'],
            execution = 'local',
            nproc = 8,
            registry_cleanup = True,
            registry_name = 'permutation_{}'.format(p)
        )
        centroid_outputs = centroids(**centroid_kwargs)

        # Generate all pairs of centroid images
        print("Generating centroid image pairs...", flush = True)
        centroid_dirs = [outdir, inputs['centroids'][1]]
        centroid_pairs = generate_cluster_pairs(centroid_dirs = centroid_dirs)

        # Identify subset of centroid pairs on and off the nk diagonal
        centroid_pairs = subset_cluster_pairs(centroid_pairs = centroid_pairs,
                                              off_diagonal = off_diagonal)

        # Add image pairs to driver kwargs
        driver_kwargs['imgs'] = [tuple(x) for x in centroid_pairs]

        # Evaluate the cluster similarity using specified execution
        print("Evaluating cluster similarity...", flush = True)
        results = transcriptomic_similarity(**driver_kwargs)

        # Output file for permuted similarity
        output_file = 'similarity_permutation_{}.csv'.format(p)
        output_file = os.path.join(paths['similarity'], output_file)
        results.to_csv(output_file, index = False)

        # Remove permuted centroids if specified
        if not keep_centroids: rmtree(centroid_dirs[0])

    # Close Dask client
    client.close()

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
