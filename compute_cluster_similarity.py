#!.venv/bin/python3
# ----------------------------------------------------------------------------
# compute_cluster_similarity.py
# Author: Antoine Beauchamp
# Created: March 15th, 2023

"""
Pipeline to evaluate pairwise cluster similarity.

Description
-----------
This pipeline evaluates the pairwise gene expression similarity of mouse and
human clusters on the basis of the cluster centroid images. The centroid images
are obtained by specifying the mouse and human processing pipeline directories
as well as parameter set IDs in the command line arguments. Additional paths
must be specified to the mouse and human gene expression data directories.

Parameters governing the evaluation of the similarity computations are mappable.
"""

# Packages -------------------------------------------------------------------

import argparse
from functools import reduce
from itertools import product
from pipelines import compute_cluster_similarity


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/cross_species/v2/',
        help = "Path to the pipeline output directory."
    )
    
    parser.add_argument(
        '--human-pipeline-dir',
        type = str,
        default = 'data/human/derivatives/v2/',
        help = "Path to the human processing pipeline directory."
    )

    parser.add_argument(
        '--mouse-pipeline-dir',
        type = str,
        default = 'data/mouse/derivatives/v2/',
        help = "Path to the mouse processing pipeline directory."
    )
    
    parser.add_argument(
        '--human-params-id',
        type = str,
        help = ("ID specifying the human processing pipeline parameter "
                "set to use.")
    )
    
    parser.add_argument(
        '--mouse-params-id',
        type = str,
        help = ("ID specifying the mouse processing pipeline parameter "
                "set to use.")
    )
    
    parser.add_argument(
        '--human-expr-dir',
        type = str,
        default = 'data/human/expression/',
        help = "Path to the human gene expression directory."
    )

    parser.add_argument(
        '--mouse-expr-dir',
        type = str,
        default = 'data/mouse/expression/',
        help = "Path to mouse gene expression directory."
    )

    parser.add_argument(
        '--human-mask',
        type = str,
        default = 'data/human/registration/v2/reference_files/mask_0.8mm.mnc',
        help = "Path to the human mask (.mnc)."
    )
    
    parser.add_argument(
        '--mouse-mask',
        type = str,
        default = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
        help = ("Path to the mouse mask (.mnc) used to construct the mouse "
                "gene expression matrix.")
    )

    parser.add_argument(
        '--human-microarray-coords',
        type = str,
        default = 'data/human/expression/AHBA_microarray_coordinates_study_v2.csv',
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
        help = ("The metric used to evaluate the similarity between cluster "
                "gene expression signature.")
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
        help = ("Method used to threshold mouse and human cluster centroid "
                "images prior to constructing gene expression signatures.")
    )

    parser.add_argument(
        '--threshold-value',
        nargs = '*',
        type = float,
        default = [0.2],
        help = ("Value used to threshold mouse and human images. Ignored if "
                "--threshold is 'none'")
    )

    parser.add_argument(
        '--threshold-symmetric',
        nargs = '*',
        type = str,
        default = ['true'],
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
    args['jacobians'] = tuple(args['jacobians'])
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
