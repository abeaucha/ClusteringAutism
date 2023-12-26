#!/usr/bin/env python3


# Packages -------------------------------------------------------------------

import argparse
import os
import sys
import utils
import pandas as pd
from glob import glob
from pyminc.volumes.factory import volumeFromFile


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    # General arguments ---------------------------------------------------
    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/human/derivatives/v2/',
        help = ("Path to the directory in which to store pipeline outputs. "
                "A sub-directory will be created based on the datasets "
                "specified using --datasets. Directories for the various "
                "pipeline stages will be created within this sub-directory.")
    )

    parser.add_argument(
        '--input-dir',
        type = str,
        default = 'data/human/registration/v2/jacobians_resampled/resolution_0.8/',
        help = ("Path to the directory containing Jacobian images. The program "
                "will look for a sub-directory 'resolution_##' using the value "
                "passed to --resolution.")
    )

    parser.add_argument(
        '--demographics',
        type = str,
        default = 'data/human/registration/v2/subject_info/demographics.csv',
        help = "Path to file (.csv) containing demographics information."
    )

    parser.add_argument(
        '--mask',
        type = str,
        default = 'data/human/registration/v2/reference_files/mask_0.8mm.mnc',
        help = "Path to mask file (.mnc) associated with the study images."
    )

    parser.add_argument(
        '--datasets',
        nargs = '*',
        type = str,
        default = ['POND', 'SickKids'],
        help = ("List of strings indicating which datasets to include in "
                "processing.")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        help = "Number of processors to use in parallel."
    )

    # Effect size arguments ---------------------------------------------------
    parser.add_argument(
        '--es-method',
        type = str,
        default = 'normative-growth',
        choices = ['normative-growth', 'propensity-matching'],
        help = "Method to use to compute effect size images."
    )

    parser.add_argument(
        '--es-group',
        type = str,
        default = 'patients',
        choices = ['patients', 'controls', 'all'],
        help = "Group of participants for which to compute effect sizes."
    )

    parser.add_argument(
        '--es-nbatches',
        type = int,
        default = 1,
        help = "Number of batches to use to process effect sizes."
    )

    parser.add_argument(
        '--es-df',
        type = int,
        default = 3,
        help = ("The number of degrees of freedom for the natural splines used "
                "in normative growth modelling. Ignored if --es-method is "
                "'propensity-matching'.")
    )

    parser.add_argument(
        '--es-batch',
        type = str,
        default = 'Site-Scanner',
        help = ("Batch variables to use for normalization prior to normative "
                "growth modelling. Variables must be found in --demographics. "
                "Multiple batch variables must be specified in a single string, "
                "separated by a hyphen (e.g. 'Site-Scanner').")
    )

    parser.add_argument(
        '--es-ncontrols',
        nargs = '*',
        type = int,
        default = 10,
        help = ("The number of controls to use for propensity-matching. "
                "Ignored if --es-method is 'normative-growth'.")
    )

    parser.add_argument(
        '--es-matrix-file',
        type = str,
        default = 'effect_sizes.csv',
        help = ("The basename of the file (.csv) in which to store voxel-wise "
                "effect size matrices.")
    )

    # Clustering arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-resolution',
        type = float,
        default = 3.0,
        help = ""
    )

    parser.add_argument(
        '--cluster-nk-max',
        type = int,
        default = 10,
        help = ("The maximum number of clusters to identify when clustering. "
                "The program will obtain cluster solutions from 2 up to the "
                "value provided.")
    )

    parser.add_argument(
        '--cluster-metric',
        type = str,
        default = 'correlation',
        help = "The distance metric to use in similarity network fusion."
    )

    parser.add_argument(
        '--cluster-K',
        type = int,
        default = 10,
        help = ("The number of nearest-neighbours to use in similarity network "
                "fusion.")
    )

    parser.add_argument(
        '--cluster-sigma',
        type = float,
        default = 0.5,
        help = "The variance for the local model in similarity network fusion."
    )

    parser.add_argument(
        '--cluster-t',
        type = int,
        default = 20,
        help = ("The number of iterations for the diffusion process in "
                "similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-file',
        type = str,
        default = 'clusters.csv',
        help = ("The basename of the file (.csv) in which to store cluster "
                "assignments.")
    )

    parser.add_argument(
        '--cluster-affinity-file',
        type = str,
        default = 'affinity.csv',
        help = (
            "The basename of the file (.csv) in which to store the affinity "
            "matrix produced by similarity network fusion.")
    )

    # Cluster maps arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-map-method',
        type = str,
        default = 'mean',
        help = "The method to use to compute cluster centroid images."
    )

    args = vars(parser.parse_args())

    return args


# Modules ----------------------------------------------------------------------

def initialize(**kwargs):
    # Effect size calculation method parameters
    if kwargs['es_method'] == 'normative-growth':
        kwargs['es_ncontrols'] = None
    elif kwargs['es_method'] == 'propensity-matching':
        kwargs['es_df'] = None
        kwargs['es_batch'] = None
    else:
        raise ValueError

    # Image resolution
    vol = volumeFromFile(kwargs['mask'])
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Pipeline parameters
    params = dict(
        dataset = '-'.join(kwargs['datasets']),
        resolution = resolution,
        es_method = kwargs['es_method'],
        es_group = kwargs['es_group'],
        es_df = kwargs['es_df'],
        es_batch = (None if kwargs['es_batch'] is None
                    else '-'.join(kwargs['es_batch'])),
        es_ncontrols = kwargs['es_ncontrols'],
        cluster_resolution = kwargs['cluster_resolution'],
        cluster_nk_max = kwargs['cluster_nk_max'],
        cluster_metric = kwargs['cluster_metric'],
        cluster_K = kwargs['cluster_K'],
        cluster_sigma = kwargs['cluster_sigma'],
        cluster_t = kwargs['cluster_t'],
        cluster_map_method = kwargs['cluster_map_method']
    )

    pipeline_dir = kwargs['pipeline_dir']
    input_dir = kwargs['input_dir']

    # Create pipeline directory
    params_id = utils.random_id(3)
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir,
                                           params_id = params_id)
    params_id = utils.fetch_params_id(metadata = metadata,
                                      params = params)

    # Directories for pipeline stages
    imgdir = os.path.join(pipeline_dir, 'jacobians', '')
    es_dir = os.path.join(pipeline_dir, 'effect_sizes',
                          'resolution_{}'.format(resolution), '')
    cluster_dir = os.path.join(pipeline_dir, 'clusters',
                               'resolution_{}'.format(
                                   kwargs['cluster_resolution']),
                               '')
    cluster_map_dir = os.path.join(pipeline_dir, 'cluster_maps',
                                   'resolution_{}'.format(resolution), '')

    # Check existence of input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: ".format(input_dir))

    # Create pipeline sub-directories
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)
    if not os.path.exists(es_dir):
        os.makedirs(es_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    if not os.path.exists(cluster_map_dir):
        os.makedirs(cluster_map_dir)

    # Filter for data sets ----------------------------------------------------

    demographics = kwargs['demographics']
    datasets = kwargs['datasets']

    # Import demographics
    df_demographics = pd.read_csv(demographics)

    # Filter individuals for data subset
    df_demographics = (df_demographics
                       .loc[df_demographics['Dataset'].isin(datasets)]
                       .copy())

    # Write out demographics subset to subset directory
    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)

    # Create symlinks to Jacobian images
    print("Creating symlinks to Jacobian images...")
    jacobians = ['absolute', 'relative']
    for j, jac in enumerate(jacobians):

        # Create symlinks to Jacobian images
        input_files = glob(os.path.join(input_dir, jac, '') + '*.mnc')
        if len(input_files) == 0:
            raise OSError("No input files in directory: ".format(input_dir))
        input_files_in_dataset = [[f for f in input_files if g in f][0]
                                  for g in df_demographics['file'].to_list()]
        imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                     dst = os.path.join(imgdir, jac, ''))

    # Dictionary containing pipeline paths
    paths = dict(
        pipeline = pipeline_dir,
        jacobians = imgdir,
        effect_sizes = es_dir,
        clusters = cluster_dir,
        centroids = cluster_map_dir
    )

    return paths


def effect_sizes(imgdir, demographics, mask, outdir,
                 method = 'normative-growth', group = 'patients',
                 nbatches = 1, df = 3, batch = 'Site-Scanner',
                 ncontrols = None, matrix_file = 'effect_sizes.csv',
                 execution = 'local', nproc = 1, **kwargs):

    jacobians = ['absolute', 'relative']
    if execution == 'local':
        for j, jac in enumerate(jacobians):
            utils.execute_local(script = 'compute_effect_sizes.R',
                                args = ...)
    elif execution == 'slurm':
        utils.execute_slurm(script = 'compute_effect_sizes.R',
                            args = ...,
                            slurm_args = ...)
    else:
        raise ValueError

    return


def clustering():
    if env == 'local':
        utils.execute_local(script = 'generate_clusters.R',
                            args = ...)
    elif env == 'slurm':
        utils.execute_slurm(script = 'generate_clusters.R',
                            args = ...,
                            slurm_args = ...)
    else:
        raise ValueError

    return


def centroids():
    if env == 'local':
        utils.execute_local(script = 'create_cluster_centroids.R',
                            args = ...)
    elif env == 'slurm':
        utils.execute_slurm(script = 'create_cluster_centroids.R',
                            args = ...,
                            slurm_args = ...)
    else:
        raise ValueError
    return


def main(pipeline_dir, input_dir, demographics, mask,
         datasets = ('POND', 'SickKids'), nproc = 1,
         es_method = 'normative-growth', es_group = 'patients',
         es_nbatches = 1, es_df = 3,
         es_batch = ('Site', 'Scanner'), es_ncontrols = 10,
         es_matrix_file = 'effect_sizes.csv',
         cluster_resolution = 3.0,
         cluster_nk_max = 10, cluster_metric = 'correlation',
         cluster_K = 10, cluster_sigma = 0.5, cluster_t = 20,
         cluster_file = 'clusters.csv',
         cluster_affinity_file = 'affinity.csv',
         cluster_map_method = 'mean'):
    # Get dictionary of function kwargs
    kwargs = locals().copy()

    # Initialize pipeline directory
    # Get paths to pipeline
    paths = initialize(**kwargs)

    print(paths)
    sys.exit()

    # Compute effect sizes
    # Arguments for effect sizes
    es_kwargs = {key.replace('es_', ''): val
                 for key, val in kwargs.items() if 'es_' in key}
    es_kwargs.update(
        dict(imgdir = paths['jacobians'],
             demographics = demographics,
             mask = mask,
             outdir = paths['effect_sizes'])
    )
    effect_sizes(**es_kwargs)

    slurm_args = dict(
        job_name = 'process_human_data_v2_0.8mm',
        nodes = 1,
        cpus_per_task = 8,
        mem = '64G',
        time = '72:00:00',
        chdir = os.getcwd(),
        output = 'logs/process_human_data_v2_0.8mm_%j.out'
    )

    # clustering()
    # centroids()


if __name__ == '__main__':
    args = parse_args()
    args['datasets'] = tuple(args['datasets'])
    main(**args)
