#!.venv/bin/python3
# ----------------------------------------------------------------------------
# process_human_data.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2023

"""
Pipeline to process human neuroimaging data.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import sys
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
        help=("Path to pipeline directory.")
    )

    parser.add_argument(
        '--input-dir',
        type=str,
        default='data/human/registration/jacobians_resampled/',
        help=("Path to Jacobians directory.")
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
        default='data/human/registration/DBM_input_demo_passedqc_wfile.csv',
        help=("Path to CSV file containing demographics information.")
    )

    parser.add_argument(
        '--mask',
        type=str,
        default='data/human/registration/reference_files/mask_3.0mm.mnc',
        help=("Path to image mask")
    )

    parser.add_argument(
        '--datasets',
        nargs='*',
        type=str,
        default=['POND', 'SickKids'],
        help=("List of strings indicating datasets to include in processing.")
    )

    parser.add_argument(
        '--parallel',
        type=str,
        default='true',
        choices=['true', 'false'],
        help=("Option to run in parallel.")
    )

    parser.add_argument(
        '--nproc',
        type=int,
        help=("Number of processors to use in parallel.")
    )

    # parser.add_argument(
    #     '--verbose',
    #     type=str,
    #     default='true',
    #     choices=['true', 'false'],
    #     help='Verbosity.'
    # )

    # Effect size arguments ---------------------------------------------------
    parser.add_argument(
        '--es-method',
        nargs=1,
        type=str,
        default=['normative-growth'],
        choices=['normative-growth', 'propensity-matching'],
        help=("Method to use to compute effect sizes.")
    )

    parser.add_argument(
        '--es-df',
        nargs='*',
        type=int,
        default=[3],
        help=("Degrees of freedom hyperparameter in normative growth modelling.")
    )

    parser.add_argument(
        '--es-combat',
        nargs='*',
        type=str,
        default=['true'],
        help=("Option to run ComBat normalization prior to normative growth modelling.")
    )

    parser.add_argument(
        '--es-combat-batch',
        nargs='*',
        type=str,
        default=['Site-Scanner'],
        help=("Batch variables to use in ComBat normalization.")
    )

    parser.add_argument(
        '--es-ncontrols',
        nargs='*',
        type=int,
        default=[10],
        help=("Number of controls to use for propensity-matching.")
    )

    parser.add_argument(
        '--es-matrix-file',
        type=str,
        default='effect_sizes.csv',
        help=("Name of CSV file in which to store voxel-wise effect size matrices.")
    )

    # Clustering arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-nk-max',
        nargs=1,
        type=int,
        default=[10],
        help=("Maximum number of clusters to identify in cluster solutions.")
    )

    parser.add_argument(
        '--cluster-metric',
        nargs='*',
        type=str,
        default=['correlation'],
        help=("Distance metric to use in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-K',
        nargs='*',
        type=int,
        default=[10],
        help=("Number of nearest-neighbours to use in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-sigma',
        nargs='*',
        type=float,
        default=[0.5],
        help=("Variance for the local model in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-t',
        nargs='*',
        type=int,
        default=[20],
        help=("Number of iterations for the diffusion process in similarity network fusion.")
    )

    parser.add_argument(
        '--cluster-file',
        type=str,
        default='clusters.csv',
        help=("Name of the CSV file in which to store cluster assignments.")
    )

    parser.add_argument(
        '--cluster-affinity-file',
        type=str,
        default='affinity.csv',
        help=("Name of the file in which to store similarity network fusion affinity matrix.")
    )

    # Cluster maps arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-map-method',
        nargs='*',
        type=str,
        default=['mean'],
        help=("Method to use to aggregate images within each cluster.")
    )

    args = vars(parser.parse_args())

    return args


# Pipeline --------------------------------------------------------------------
if __name__ == '__main__':

    args = parse_args()
    args['parallel'] = True if args['parallel'] == 'true' else False
    # args['verbose'] = True if args['verbose'] == 'true' else False
    params = {key: val for key, val in args.items()
              if 'es_' in key or 'cluster_' in key}
    del params['es_matrix_file']
    del params['cluster_file']
    del params['cluster_affinity_file']
    param_keys = list(params.keys())
    for param_vals in product(*params.values()):
        param_set = dict(list(zip(param_keys, param_vals)))
        param_set['es_combat'] = (True if param_set['es_combat'] == 'true'
                                  else False)
        param_set['es_combat_batch'] = param_set['es_combat_batch'].split('-')
        param_msg = [': '.join([str(key), str(val)])
                     for key, val in param_set.items()]
        param_msg = reduce(lambda x, y: x + '\n\t' + y, param_msg)
        param_msg = "Running parameter set:\n\t{}\n".format(param_msg)
        print(param_msg)
        args.update(param_set)
        process_human_data(**args)
