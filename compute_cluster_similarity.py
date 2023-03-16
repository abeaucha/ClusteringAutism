#!.venv/bin/python3
# ----------------------------------------------------------------------------
# cluster_similarity.py
# Author: Antoine Beauchamp
# Created: March 15th, 2023

"""
Pipeline to evaluate pairwise cluster similarity.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
from pipelines import compute_cluster_similarity
from itertools import product
from functools import reduce


# Command line arguments -----------------------------------------------------
def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/cross_species/similarity/',
        help = "Path to pipeline output directory."
    )

    parser.add_argument(
        '--mouse-cluster-dir',
        type = str,
        help = "Path to directory containing mouse cluster maps."
    )

    parser.add_argument(
        '--human-cluster-dir',
        type = str,
        help = "Path to directory containing human cluster maps."
    )

    parser.add_argument(
        '--mouse-expr-dir',
        type = str,
        default = 'data/mouse/expression/',
        help = "Path to mouse expression directory."
    )

    parser.add_argument(
        '--human-expr-dir',
        type = str,
        default = 'data/human/expression/',
        help = "Path to human expression directory."
    )

    parser.add_argument(
        '--mouse-mask',
        type = str,
        default = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
        help = ("Path to mouse mask (.mnc) used to construct the gene "
                "expression matrix.")
    )

    parser.add_argument(
        '--human-mask',
        type = str,
        default = 'data/human/registration/reference_files/mask_3.0mm.mnc',
        help = "Path to human mask (.mnc)."
    )

    parser.add_argument(
        '--mouse-resolution',
        type = float,
        default = 0.2,
        help = "Resolution (mm) of mouse images."
    )

    parser.add_argument(
        '--human-resolution',
        type = float,
        default = 3.0,
        help = "Resolution (mm) of human images."
    )

    parser.add_argument(
        '--human-microarray-coords',
        type = str,
        default = 'data/human/expression/AHBA_microarray_coordinates_studyspace.csv',
        help = ("Path to file (.csv) containing the world coordinates of "
                "the AHBA microarray samples.")
    )

    parser.add_argument(
        '--gene-space',
        type = str,
        default = 'average-latent-space',
        choices = ['average-latent-space', 'latent-space', 'homologous-genes'],
        help = "Gene expression space to use when computing similarity."
    )

    parser.add_argument(
        '--n-latent-spaces',
        nargs = '*',
        type = int,
        default = [100],
        help = ("Number of latent spaces to include when --gene-space is "
                "'average-latent-space'. Ignored otherwise.")
    )

    parser.add_argument(
        '--latent-space-id',
        nargs = '*',
        type = int,
        default = [1],
        help = ("Latent space to use when --gene-space is 'latent-space'. "
                "Ignored otherwise.")
    )

    parser.add_argument(
        '--metric',
        nargs = '*',
        type = str,
        default = ['correlation'],
        help = "Similarity metric."
    )

    parser.add_argument(
        '--signed',
        nargs = '*',
        type = str,
        default = ['true'],
        help = ("Option to compute positive and negative similarity "
                "separately before averaging for a final value.")
    )

    parser.add_argument(
        '--threshold',
        type = str,
        default = 'top_n',
        choices = ['none', 'top_n', 'intensity'],
        help = "Method used to threshold mouse and human images."
    )

    parser.add_argument(
        '--threshold-value',
        nargs = '*',
        type = float,
        default = [0.2],
        help = "Value used to threshold mouse and human images."
    )

    parser.add_argument(
        '--threshold-symmetric',
        nargs = '*',
        type = str,
        default = ['true'],
        help = ("Option to apply threshold symmetrically to positive and "
                "negative values.")
    )

    parser.add_argument(
        '--parallel',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = "Option to run in parallel."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        help = "Number of processors to use in parallel."
    )

    return vars(parser.parse_args())


# Pipeline --------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()

    args['parallel'] = True if args['parallel'] == 'true' else False

    if args['gene_space'] == 'average-latent-space':
        args['latent_space_id'] = None
    elif args['gene_space'] == 'latent-space':
        args['n_latent_spaces'] = None
    else:
        args['latent_space_id'] = None
        args['n_latent_spaces'] = None

    if args['threshold'] == 'none':
        args['threshold'] = None
        args['threshold_value'] = None
        args['threshold_symmetric'] = None

    params = {key: val for key, val in args.items() if (type(val) is list)}
    param_keys = list(params.keys())
    for param_vals in product(*params.values()):
        param_set = dict(list(zip(param_keys, param_vals)))
        param_set['signed'] = (True if param_set['signed'] == 'true'
                               else 'false')
        param_set['threshold_symmetric'] = (
            True if param_set['threshold_symmetric'] == 'true' else 'false'
        )
        param_msg = [': '.join([str(key), str(val)])
                     for key, val in param_set.items()]
        param_msg = reduce(lambda x, y: x + '\n\t' + y, param_msg)
        param_msg = "Running parameter set:\n\t{}\n".format(param_msg)
        print(param_msg)
        args.update(param_set)
        compute_cluster_similarity(**args)
