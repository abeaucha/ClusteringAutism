#!.venv/bin/python3
# ----------------------------------------------------------------------------
# process_human_data.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2023

"""
Pipeline to process and cluster human neuroimaging data. 

Description
-----------
The pipeline executes three main stages:
1. Compute effect size images using absolute and relative Jacobian images
2. Cluster individuals used similarity network fusion and spectral clustering
3. Create cluster centroid images using absolute and relative effect size images

The pipeline is mappable over multiple parameter sets.
"""

# Packages -------------------------------------------------------------------

import argparse
from pipelines import process_human_data
from itertools import product
from functools import reduce


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # General arguments ---------------------------------------------------
    parser.add_argument(
        '--pipeline-dir',
        type=str,
        default='data/human/derivatives/',
        help=("Path to the directory in which to store pipeline outputs. "
              "A sub-directory will be created based on the datasets specified "
              "using --datasets. Directories for the various pipeline stages "
              "will be created within this sub-directory.")
    )

    parser.add_argument(
        '--input-dir',
        type=str,
        default='data/human/registration/v1/jacobians_resampled/',
        help=("Path to the directory containing Jacobian images. The program "
              "will look for a sub-directory 'resolution_##' using the value "
              "passed to --resolution.")
    )

    parser.add_argument(
        '--resolution',
        type=float,
        default=3.0,
        help=("Resolution (mm) of Jacobian images to use.")
    )

    parser.add_argument(
        '--demographics',
        type=str,
        default='data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv',
        help="Path to file (.csv) containing demographics information."
    )

    parser.add_argument(
        '--mask',
        type=str,
        default='data/human/registration/v1/reference_files/mask_3.0mm.mnc',
        help="Path to mask file (.mnc) associated with the study images."
    )

    parser.add_argument(
        '--datasets',
        nargs='*',
        type=str,
        default=['POND', 'SickKids'],
        help=("List of strings indicating which datasets to include in processing.")
    )

    parser.add_argument(
        '--nproc',
        type=int,
        help="Number of processors to use in parallel."
    )

    # Effect size arguments ---------------------------------------------------
    parser.add_argument(
        '--es-method',
        nargs=1,
        type=str,
        default=['normative-growth'],
        choices=['normative-growth', 'propensity-matching'],
        help="Method to use to compute effect size images."
    )

    parser.add_argument(
        '--es-nbatches',
        type = int,
        default = 1,
        help = "Number of batches to use to process effect sizes."
    )

    parser.add_argument(
        '--es-df',
        nargs='*',
        type=int,
        default=[3],
        help=("The number of degrees of freedom for the natural splines used "
              "in normative growth modelling. Ignored if --es-method is "
             "'propensity-matching'.")
    )

    parser.add_argument(
        '--es-batch',
        nargs='*',
        type=str,
        default=['Site-Scanner'],
        help=("Batch variables to use for normalization prior to normative "
              "growth modelling. Variables must be found in --demographics. "
              "Multiple batch variables must be specified in a single string, "
              "separated by a hyphen (e.g. 'Site-Scanner').")
    )

    parser.add_argument(
        '--es-ncontrols',
        nargs='*',
        type=int,
        default=[10],
        help=("The number of controls to use for propensity-matching. "
             "Ignored if --es-method is 'normative-growth'.")
    )

    parser.add_argument(
        '--es-matrix-file',
        type=str,
        default='effect_sizes.csv',
        help=("The basename of the file (.csv) in which to store voxel-wise "
              "effect size matrices.")
    )

    # Clustering arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-nk-max',
        nargs=1,
        type=int,
        default=[10],
        help=("The maximum number of clusters to identify when clustering. "
             "The program will obtain cluster solutions from 2 up to the "
             "value provided.")
    )

    parser.add_argument(
        '--cluster-metric',
        nargs='*',
        type=str,
        default=['correlation'],
        help="The distance metric to use in similarity network fusion."
    )

    parser.add_argument(
        '--cluster-K',
        nargs='*',
        type=int,
        default=[10],
        help=("The number of nearest-neighbours to use in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-sigma',
        nargs='*',
        type=float,
        default=[0.5],
        help=("The variance for the local model in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-t',
        nargs='*',
        type=int,
        default=[20],
        help=("The number of iterations for the diffusion process in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-file',
        type=str,
        default='clusters.csv',
        help=("The basename of the file (.csv) in which to store cluster assignments.")
    )

    parser.add_argument(
        '--cluster-affinity-file',
        type=str,
        default='affinity.csv',
        help=("The basename of the file (.csv) in which to store the affinity matrix "
             "produced by similarity network fusion.")
    )

    # Cluster maps arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-map-method',
        nargs='*',
        type=str,
        default=['mean'],
        help=("The method to use to compute cluster centroid images.")
    )

    args = vars(parser.parse_args())

    return args


# Pipeline --------------------------------------------------------------------
if __name__ == '__main__':

    args = parse_args()
    args['datasets'] = tuple(args['datasets'])
    params = {key: val for key, val in args.items()
              if 'es_' in key or 'cluster_' in key}
    del params['es_matrix_file']
    del params['es_nbatches']
    del params['cluster_file']
    del params['cluster_affinity_file']
    param_keys = list(params.keys())
    for param_vals in product(*params.values()):
        param_set = dict(list(zip(param_keys, param_vals)))
        param_set['es_batch'] = tuple(param_set['es_batch'].split('-'))
        param_msg = [': '.join([str(key), str(val)])
                     for key, val in param_set.items()]
        param_msg = reduce(lambda x, y: x + '\n\t' + y, param_msg)
        param_msg = "Running parameter set:\n\t{}\n".format(param_msg)
        print(param_msg)
        args.update(param_set)
        process_human_data(**args)
