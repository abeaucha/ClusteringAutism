#!.venv/bin/python3
# ----------------------------------------------------------------------------
# pl_process_initialize.py
# Author: Antoine Beauchamp
# Created: November 29th, 2023

"""
Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import os
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

    # parser.add_argument(
    #     '--nproc',
    #     type = int,
    #     help = "Number of processors to use in parallel."
    # )

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
        default = 4,
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
        help = ("The basename of the file (.csv) in which to store the affinity "
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


if __name__ == '__main__':

    args = parse_args()

    # Effect size calculation method parameters
    if args['es_method'] == 'normative-growth':
        args['es_ncontrols'] = None
    elif args['es_method'] == 'propensity-matching':
        args['es_df'] = None
        args['es_batch'] = None
    else:
        raise ValueError

    # Image resolution
    vol = volumeFromFile(args['mask'])
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Pipeline parameters
    params = dict(
        dataset = '-'.join(args['datasets']),
        resolution = resolution,
        es_method = args['es_method'],
        es_group = args['es_group'],
        es_df = args['es_df'],
        es_batch = args['es_batch'],
        es_ncontrols = args['es_ncontrols'],
        cluster_resolution = args['cluster_resolution'],
        cluster_nk_max = args['cluster_nk_max'],
        cluster_metric = args['cluster_metric'],
        cluster_K = args['cluster_K'],
        cluster_sigma = args['cluster_sigma'],
        cluster_t = args['cluster_t'],
        cluster_map_method = args['cluster_map_method']
    )

    pipeline_dir = args['pipeline_dir']
    input_dir = args['input_dir']

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
                               'resolution_{}'.format(args['cluster_resolution']), '')
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

    demographics = args['demographics']
    datasets = args['datasets']

    # Import demographics
    df_demographics = pd.read_csv(demographics)

    # Filter individuals for data subset
    df_demographics = (df_demographics
                       .loc[df_demographics['Dataset'].isin(datasets)]
                       .copy())

    # Write out demographics subset to subset directory
    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)

    # Iterate over jacobians to compute effect sizes
    print("Creating symlinks to Jacobian images...")
    jacobians = ['absolute', 'relative']
    for j, jac in enumerate(jacobians):

        # Create symlinks to Jacobian images ----------------------------------
        input_files = glob(os.path.join(input_dir, jac, '') + '*.mnc')
        if len(input_files) == 0:
            raise OSError("No input files in directory: ".format(input_dir))
        input_files_in_dataset = [[f for f in input_files if g in f][0]
                                  for g in df_demographics['file'].to_list()]
        imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                     dst = os.path.join(imgdir, jac, ''))
