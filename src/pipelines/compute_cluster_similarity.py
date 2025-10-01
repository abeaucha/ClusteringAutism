#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# compute_cluster_similarity.py
# Author: Antoine Beauchamp
# Created: March 27th, 2024

"""
Evaluate the similarity between pairs of clusters.

Description
-----------
This pipeline evaluates the similarity between pairs of imaging clusters
based on the transcriptomic similarity of their centroid images.
"""


# Packages -------------------------------------------------------------------

import argparse
import os
import utils
import pandas as pd
from itertools import product
from transcriptomic import transcriptomic_similarity
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster
from time import sleep

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
        help = ("Optional parameter set ID for pipeline. "
                "ID will be randomly generated if not specified.")
    )

    parser.add_argument(
        '--species',
        nargs = 2,
        type = str,
        default = ['human', 'mouse'],
        help = ("List of strings indicating which species are being compared. "
                "Entries must be 'human' or 'mouse'.")
    )

    parser.add_argument(
        '--input-dirs',
        nargs = 2,
        type = str,
        default = ['data/human/derivatives/', 'data/mouse/derivatives/'],
        help = ("Paths to the processing pipeline directories containing "
                "the centroid images to compare. Expects a sub-directory "
                "'centroids'.")
    )

    parser.add_argument(
        '--input-params-ids',
        nargs = 2,
        type = str,
        help = ("List of integer IDs specifying the processing pipeline "
                "parameter sets to use.")
    )

    parser.add_argument(
        '--expr-dirs',
        nargs = 2,
        type = str,
        default = ['data/human/expression/', 'data/mouse/expression/'],
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
        '--gene-space',
        type = str,
        default = 'avg-mlp-latent-space',
        choices = ['avg-mlp-latent-space', 'mlp-latent-space', 'vae-latent-space', 'homologous-genes'],
        help = ("The gene expression space to use to evaluate the similarity "
                "between cluster centroid images.")
    )

    parser.add_argument(
        '--n-latent-spaces',
        type = int,
        default = 100,
        help = ("The number of latent spaces to include when --gene-space is "
                "'avg-mlp-latent-space'. Ignored otherwise.")
    )

    parser.add_argument(
        '--latent-space-id',
        type = int,
        default = 1,
        help = ("The ID of the latent space to use when --gene-space is "
                "'mlp-latent-space'. Ignored otherwise.")
    )

    parser.add_argument(
        '--metric',
        type = str,
        default = 'correlation',
        help = ("The metric used to evaluate the similarity between cluster "
                "centroid gene expression signature.")
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
        help = ("The method used to threshold the cluster centroid images "
                "prior to constructing the gene expression signatures.")
    )

    parser.add_argument(
        '--threshold-value',
        type = float,
        default = 0.2,
        help = ("The value used to threshold the centroid images. Ignored if "
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
        help = "The class of Jacobian images to use."
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
    Initialize the similarity pipeline.

    Parameters
    ----------
    kwargs: dict
        All arguments passed to the main() function.

    Returns
    -------
    paths: dict
        Dictionary containing the pipeline path and the paths to the
        input centroid directories.
    """

    # Ensure proper paths
    pipeline_dir = os.path.join(kwargs['pipeline_dir'], '')
    input_dirs = [os.path.join(path, '') for path in kwargs['input_dirs']]

    # Convert bool to str for multilingual consistency
    for key, val in kwargs.items():
        if type(val) is bool:
            kwargs[key] = 'true' if val else 'false'

    # Fetch the input pipeline parameters
    params = dict()
    param_ids = kwargs['input_params_ids']
    metadata = [os.path.join(path, 'metadata.csv') for path in input_dirs]
    for i in range(len(metadata)):
        if not os.path.exists(metadata[i]):
            raise OSError("Input pipeline metadata file not found: {}"
                          .format(metadata[i]))
        params_i = utils.fetch_params_metadata(metadata[i],
                                               id = param_ids[i])
        params_i = params_i[['id', 'dataset', 'resolution']]
        params_i = params_i.to_dict(orient = 'list')
        params_i['species'] = [kwargs['species'][i]]
        params_i = {'_'.join(['input', str(i + 1), key]):val[0]
                    for key, val in params_i.items()}
        params.update(params_i)

    # Update the parameter sets with the current pipeline parameters
    params.update(
        dict(gene_space = kwargs['gene_space'],
             n_latent_spaces = kwargs['n_latent_spaces'],
             latent_space_id = kwargs['latent_space_id'],
             metric = kwargs['metric'],
             signed = kwargs['signed'],
             threshold = kwargs['threshold'],
             threshold_value = kwargs['threshold_value'],
             threshold_symmetric = kwargs['threshold_symmetric'])
    )

    # Create the pipeline directory
    if kwargs['params_id'] is None:
        params_id = utils.random_id(3)
    else:
        params_id = kwargs['params_id']
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir,
                                           params_id = params_id)
    pipeline_dir = os.path.join(pipeline_dir, 'similarity', '')
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)

    # Extract the input pipeline image resolutions
    resolutions = list([params['input_1_resolution'],
                        params['input_2_resolution']])

    # Extract the input pipeline centroid directories
    centroid_dirs = [os.path.join(input_dirs[i], param_ids[i],
                                  'centroids',
                                  'resolution_{}'.format(resolutions[i]), '')
                     for i in range(len(input_dirs))]

    # Dictionary containing the pipeline paths
    paths = dict(
        pipeline = pipeline_dir,
        centroids = centroid_dirs
    )

    return paths


def generate_cluster_pairs(centroid_dirs, jacobians = ('absolute', 'relative')):
    """
    Generate pairs of cluster centroid images.

    Parameters
    ----------
    centroid_dirs: list of str
        List of paths to the directories containing the centroid images.
        Expects subdirectories named 'absolute', 'relative', or both.
    jacobians: tuple of str
        Strings indicating which Jacobian images to use.

    Returns
    -------
    centroid_pairs: list of tuple of str
        List of tuples containing the pairs of centroid images.
    """

    # Iterate over Jacobians
    for j, jac in enumerate(jacobians):

        # Update input paths with jacobians
        centroid_dirs_j = [os.path.join(path, jac, '')
                           for path in centroid_dirs]

        # Get input centroid image files
        centroids_j = [os.listdir(path) for path in centroid_dirs_j]

        # Prepend directory path to centroid image files
        centroids_j = [[os.path.join(centroid_dirs_j[i], file)
                        for file in centroids_j[i]]
                       for i in range(len(centroids_j))]

        # Expand centroid combinations for current Jacobians
        cluster_pairs_j = [list(pair) for pair in
                           list(product(centroids_j[0], centroids_j[1]))]

        # Concatenate Jacobian image pairs
        if j == 0:
            cluster_pairs = cluster_pairs_j
        else:
            cluster_pairs = cluster_pairs + cluster_pairs_j

    return cluster_pairs


def main(pipeline_dir, species, input_dirs, input_params_ids, expr_dirs, masks,
         params_id = None,
         microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_study.csv',
         gene_space = 'avg-mlp-latent-space',
         n_latent_spaces = 100, latent_space_id = 1,
         metric = 'correlation', signed = True,
         threshold = 'top_n', threshold_value = 0.2,
         threshold_symmetric = True,
         jacobians = ('absolute', 'relative'),
         execution = 'local', nproc = 1,
         slurm_mem = None, slurm_time = None):

    """
    Execute the cluster similarity pipeline.

    Parameters
    ----------
    pipeline_dir: str
        Path to the directory in which to export pipeline outputs.
        A uniquely identified subdirectory will be created using the
        specified set of pipeline parameters.
    species: tuple of str
        Strings indicating which species are being compared.
    input_dirs: tuple of str
        Paths to the processing pipeline directories containing the
        centroid images to compare. Expects a subdirectory 'centroids'.
    input_params_ids: tuple of str
        List of integer IDs specifying the processing pipeline
        parameter sets to use.
    expr_dirs: tuple of str
        Paths to the gene expression directories for the species
        being compared.
    masks: tuple of str
        Paths to the mask image files (.mnc).
    params_id: str
        Optional parameter set ID for pipeline. ID will be randomly
        generated if not specified.
    microarray_coords: str
        Path to file (.csv) containing the world coordinates of the
        AHBA microarray samples in the human imaging study space.
    gene_space: {'avg-mlp-latent-space', 'mlp-latent-space', 'vae-latent-space', 'homologous-genes'}
        The gene expression space to use to evaluate the similarity
        between cluster centroid images.
    n_latent_spaces: int, default 100
        The number of latent spaces to include when
        gene_space = 'avg-mlp-latent-space'. Ignored otherwise.
    latent_space_id: int, default 1
        The ID of the latent space to use when
        gene_space = 'mlp-latent_space'. Ignored otherwise.
    metric: str, default 'correlation'
        The metric used to evaluate the similarity between cluster
        centroid gene expression signature.
    signed: bool, default True
        Option to compute positive and negative similarity separately
        before averaging for a final value.
    threshold: {'top_n', 'intensity', None}
        The method used to threshold the cluster centroid images prior
        to constructing the gene expression signatures.
    threshold_value: float, default 0.2
        The value used to threshold the centroid images. Ignored if
        threshold = None.
    threshold_symmetric: bool, default True
        Option to apply the threshold symmetrically to positive and
        negative image values.
    jacobians: str or tuple of str, default ('absolute', 'relative')
        The class of Jacobian images to use.
    execution: {'local', 'slurm'}
        Flag indicating whether the pipeline should be executed or using
        the Slurm scheduler on an HPC cluster.
    nproc: int, default 1
        Number of processors to use.
    slurm_mem: str, default None
        Memory per CPU on Slurm. Ignored when execution = 'local'.
    slurm_time: str, default None
        Walltime in hh:mm:ss format for Slurm jobs.
        Ignored when execution = 'local'.

    Returns
    -------
    None
    """

    # Adapt gene space parameters to the selected gene space
    if gene_space == 'avg-mlp-latent-space':
        latent_space_id = None
    elif gene_space == 'mlp-latent-space':
        n_latent_spaces = None
    elif gene_space == 'vae-latent-space':
        latent_space_id = None
        n_latent_spaces = None
    elif gene_space == 'homologous-genes':
        latent_space_id = None
        n_latent_spaces = None
    else:
        raise ValueError("Argument `gene_space` must be one of "
                         "{'avg-mlp-latent-space', 'mlp-latent-space', "
                         "'vae-latent-space', 'homologous-genes'}.")

    # Adapt thresholding parameters if no thresholding specified
    if threshold is None:
        threshold_value = None
        threshold_symmetric = None

    # Get local kwargs
    kwargs = locals().copy()

    # Initialize pipeline directory tree
    print("Initializing pipeline...", flush = True)
    paths = initialize(**kwargs)

    # Generate pairs of centroid images
    print("Generating centroid image pairs...", flush = True)
    cluster_pairs = generate_cluster_pairs(centroid_dirs = paths['centroids'],
                                           jacobians = jacobians)

    # Export centroid pairs
    output_file = os.path.join(paths['pipeline'], 'centroid_pairs.csv')
    pd.DataFrame(cluster_pairs).to_csv(output_file, index = False)

    # Set imgs arg in driver kwargs
    kwargs['imgs'] = [tuple(x) for x in cluster_pairs[:1000]]

    # Initialize Dask client for execution
    if kwargs['execution'] == 'local':
        cluster = LocalCluster(n_workers = kwargs['nproc'],
                               threads_per_worker = 1)
    elif kwargs['execution'] == 'slurm':
        cluster = SLURMCluster(
            cores = 1,
            memory = kwargs['slurm_mem'],
            walltime = kwargs['slurm_time']
        )
        cluster.scale(jobs = kwargs['nproc'])
    else:
        raise ValueError("Argument `--execution` must be one of {'local', 'slurm'}")

    client = Client(cluster)

    print("Dask dashboard at: {}".format(client.dashboard_link))

    # Add client to driver kwargs
    kwargs['client'] = client

    # Clean up driver kwargs
    kwargs['expr'] = kwargs.pop('expr_dirs')
    del kwargs['pipeline_dir']
    del kwargs['params_id']
    del kwargs['input_dirs']
    del kwargs['input_params_ids']
    del kwargs['jacobians']
    del kwargs['execution']
    del kwargs['nproc']
    del kwargs['slurm_mem']
    del kwargs['slurm_time']

    # Compute pairwise similarity between cluster centroids
    print("Evaluating transcriptomic similarity of clusters...", flush = True)
    results = transcriptomic_similarity(**kwargs)

    # Close Dask client
    print("Closing client", flush = True)
    client.close()
    print("Waiting 10 seconds...", flush = True)
    sleep(10)
    print("Closing cluster", flush = True)
    cluster.close()

    # Export similarity values
    print("Exporting similarity values...", flush = True)
    output_file = os.path.join(paths['pipeline'], 'similarity.csv')
    results.to_csv(output_file, index = False)

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    args['input_dirs'] = tuple(args['input_dirs'])
    if args['input_params_ids'] is None:
        raise ValueError("Argument `--input-params-ids` is required.")
    else:
        args['input_params_ids'] = tuple(args['input_params_ids'])
    args['species'] = tuple(args['species'])
    args['expr_dirs'] = tuple(args['expr_dirs'])
    args['masks'] = tuple(args['masks'])
    args['jacobians'] = tuple(args['jacobians'])
    args['signed'] = True if args['signed'] == 'true' else False
    args['threshold'] = (None if args['threshold'] == 'none'
                         else args['threshold'])
    args['threshold_symmetric'] = (True if args['threshold_symmetric'] == 'true'
                                   else False)
    #with utils.catchtime() as t:
    main(**args)
    #print(f'Time: {t():.3f} seconds')
