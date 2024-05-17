#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# compute_cluster_similarity.py
# Author: Antoine Beauchamp
# Created: March 27th, 2024

"""
Pipeline to evaluate pairwise cluster similarity.

Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import os
import sys
import utils
import pandas as pd
from itertools import product



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
        default = 'data/human/expression/v3/AHBA_microarray_coordinates_study.csv',
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
    paths: list of str
        List of paths to the pipeline sub-directories
    """

    # Ensure proper paths
    pipeline_dir = os.path.join(kwargs['pipeline_dir'], '')
    input_dirs = [os.path.join(path, '') for path in kwargs['input_dirs']]
    expr_dirs = [os.path.join(path, '') for path in kwargs['expr_dirs']]

    # Convert bool to str for multi-language consistency
    for key, val in kwargs.items():
        if type(val) is bool:
            kwargs[key] = 'true' if val else 'false'

    # Fetch pipeline parameters
    params = dict()
    param_ids = kwargs['param_ids']
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

    # Update parameter sets with pipeline parameters
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

    # Create pipeline directory for parameter set
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir)
    pipeline_dir = os.path.join(pipeline_dir, 'similarity', '')
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)

    # Image resolutions
    resolutions = list([params['input_1_resolution'],
                        params['input_2_resolution']])

    # Centroid directories
    centroid_dirs = [os.path.join(input_dirs[i], param_ids[i],
                                  'centroids',
                                  'resolution_{}'.format(resolutions[i]), '')
                     for i in range(len(input_dirs))]

    # Dictionary containing pipeline paths
    paths = dict(
        pipeline = pipeline_dir,
        centroids = centroid_dirs
    )

    return paths


def generate_cluster_pairs(centroid_dirs, jacobians = ('absolute', 'relative')):
    """
    Generate pairs of cluster centroid images

    Parameters
    ----------
    centroid_dirs: list of str
        List of paths to the directories containing the centroid images.
        Expects sub-directories named 'absolute', 'relative', or both.
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


def main(pipeline_dir, species, input_dirs, param_ids, expr_dirs, masks,
         microarray_coords = 'data/human/expression/v3/AHBA_microarray_coordinates_study.csv',
         gene_space = 'average-latent-space',
         n_latent_spaces = 100, latent_space_id = 1,
         metric = 'correlation', signed = True,
         threshold = 'top_n', threshold_value = 0.2,
         threshold_symmetric = True,
         jacobians = ('absolute', 'relative'),
         execution = 'local', nproc = 1,
         registry_name = 'compute_cluster_similarity_registry',
         registry_cleanup = True, slurm_njobs = None, slurm_mem = None,
         slurm_time = None):
    # Adapt gene space parameters to selected gene space
    if gene_space == 'average-latent-space':
        latent_space_id = None
    elif gene_space == 'latent-space':
        n_latent_spaces = None
    elif gene_space == 'homologous-genes':
        latent_space_id = None
        n_latent_spaces = None
    else:
        raise ValueError

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
    cluster_pairs = pd.DataFrame(cluster_pairs)

    # Driver script
    driver = 'transcriptomic_similarity.py'

    # Update kwargs for driver
    kwargs['expr'] = kwargs.pop('expr_dirs')
    kwargs['signed'] = 'true' if signed else 'false'
    kwargs['threshold_symmetric'] = 'true' if threshold_symmetric else 'false'
    del kwargs['pipeline_dir']
    del kwargs['input_dirs']
    del kwargs['param_ids']
    del kwargs['jacobians']
    del kwargs['execution']
    del kwargs['registry_name']
    del kwargs['registry_cleanup']
    del kwargs['slurm_njobs']
    del kwargs['slurm_mem']
    del kwargs['slurm_time']

    # Execution mode
    print("Evaluating cluster similarity...", flush = True)
    if execution == 'local':

        # Export cluster pairs
        outfile = os.path.join(paths['pipeline'], 'centroid_pairs.csv')
        cluster_pairs.to_csv(outfile, index = False)

        # Update kwargs for driver
        kwargs['input-file'] = outfile
        kwargs['output-file'] = os.path.join(paths['pipeline'], 'similarity.csv')
        kwargs = {key.replace('_', '-'):val for key, val in kwargs.items()}

        # Execute driver
        utils.execute_local(script = driver, kwargs = kwargs)


    elif execution == 'slurm':

        # Slurm job resources
        resources = dict(
            nodes = 1,
            mem = slurm_mem,
            time = slurm_time
        )

        # Create registry
        registry = utils.Registry(resources = resources,
                                  name = registry_name)

        # Create data batches
        registry.create_batches(x = cluster_pairs,
                                nbatches = slurm_njobs,
                                prefix = 'centroid_pairs_batch')

        # Create jobs for batches
        kwargs['nproc'] = 1
        registry.create_jobs(script = driver,
                             kwargs = kwargs)

        # Submit jobs
        out = registry.submit_jobs(wait = True, cleanup = registry_cleanup)

        # Export results
        print("Exporting results...", flush = True)
        output_file = os.path.join(paths['pipeline'], 'similarity.csv')
        out.to_csv(output_file, index = False)

    else:
        raise ValueError("Argument `execution` must be one of "
                         "{'local', 'slurm'}")

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    args['input_dirs'] = tuple(args['input_dirs'])
    args['param_ids'] = tuple(args['param_ids'])
    args['species'] = tuple(args['species'])
    args['expr_dirs'] = tuple(args['expr_dirs'])
    args['masks'] = tuple(args['masks'])
    args['jacobians'] = tuple(args['jacobians'])
    args['signed'] = True if args['signed'] == 'true' else False
    args['threshold'] = (None if args['threshold'] == 'none'
                         else args['threshold'])
    args['threshold_symmetric'] = (True if args['threshold_symmetric'] == 'true'
                                   else False)
    args['registry_cleanup'] = (True if args['registry_cleanup'] == 'true'
                                else False)
    with utils.catchtime() as t:
        main(**args)
    print(f'Time: {t():.3f} seconds')
