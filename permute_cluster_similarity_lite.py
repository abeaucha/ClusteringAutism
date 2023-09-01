#!.venv/bin/python3
# ----------------------------------------------------------------------------
# permute_cluster_similarity_lite.py
# Author: Antoine Beauchamp
# Created: September 1st, 2023

"""
Pipeline to compute permuted cluster similarity values for specific human and
mouse cluster solutions.

Description
-----------
This pipeline iteratively permutes the cluster assignments for the human
participants and evaluates the pairwise mouse-human cluster similarity values
using the permuted human clusters. The number of permutations can be specified
as well as the permutation to begin with. Each permutation has a fixed random
seed corresponding to the permutation integer. The pipeline is mappable over
absolute and relative Jacobians.

Note that cluster similarity parameters are inherited from the output of
a mouse-human cluster similarity pipeline. Thus, the similarity pipeline but
be completed before running permutations.
"""

# Packages -------------------------------------------------------------------

import argparse
from pipelines import permute_cluster_similarity_lite


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
        '--params-id',
        type = str,
        help = ("ID specifying the mouse-human cluster similarity pipeline "
                "parameter set to use.")
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
        '--human-nk',
        type = int,
        default = 2,
        help = "Human cluster solution."
    )

    parser.add_argument(
        '--mouse-nk',
        type = int,
        default = 2,
        help = "Mouse cluster solution."
    )

    parser.add_argument(
        '--npermutations',
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
        '--jacobians',
        nargs = '*',
        type = str,
        default = ['absolute', 'relative'],
        help = "Jacobians to use"
    )

    parser.add_argument(
        '--keep-cluster-maps',
        type = str,
        default = 'false',
        help = "Option to keep permuted cluster maps."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use. Executes in parallel if > 1."
    )

    return vars(parser.parse_args())


# Main -----------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()
    args['jacobians'] = tuple(args['jacobians'])
    args['keep_cluster_maps'] = (True if args['keep_cluster_maps'] == 'true'
                                 else False)
    permute_cluster_similarity_lite(**args)
