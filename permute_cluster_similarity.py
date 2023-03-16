# ----------------------------------------------------------------------------
# permute_cluster_similarity.py
# Author: Antoine Beauchamp
# Created: March 16th, 2023

"""
Brief description.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
from pipelines import permute_cluster_similarity


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--human-pipeline-dir',
        type = str,
        default = 'data/human/derivatives/',
        help = "Path to human pipeline directory."
    )

    parser.add_argument(
        '--mouse-pipeline-dir',
        type = str,
        default = 'data/mouse/derivatives/',
        help = "Path to mouse pipeline directory."
    )

    parser.add_argument(
        '--human-expr-dir',
        type = str,
        default = 'data/human/expression/',
        help = "Path to human expression directory."
    )

    parser.add_argument(
        '--mouse-expr-dir',
        type = str,
        default = 'data/mouse/expression/',
        help = "Path to mouse expression directory."
    )

    parser.add_argument(
        '--human-mask',
        type = str,
        default = 'data/human/registration/reference_files/mask_3.0mm.mnc',
        help = "Path to human mask (.mnc)."
    )

    parser.add_argument(
        '--mouse-mask',
        type = str,
        default = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
        help = ("Path to mouse mask (.mnc) used to construct the gene "
                "expression matrix.")
    )

    parser.add_argument(
        '--human-microarray-coords',
        type = str,
        default = 'data/human/expression/AHBA_microarray_coordinates_studyspace.csv',
        help = ("Path to file (.csv) containing the world coordinates of "
                "the AHBA microarray samples.")
    )

    parser.add_argument(
        '--output-dir',
        type = str,
        default = 'data/cross_species/',
        help = "Path to pipeline output directory."
    )

    parser.add_argument(
        '--human-cluster-params-id',
        type = str,
        help = "Parameter set ID for human clusters."
    )

    parser.add_argument(
        '--human-dataset',
        type = str,
        default = 'POND_SickKids',
        help = "Human dataset to use."
    )
    
    parser.add_argument(
        '--mouse-dataset',
        type = str,
        default = 'Models_135',
        help = "Mouse dataset to use."
    )

    parser.add_argument(
        '--npermutations',
        type = int,
        default = 100,
        help = "Number of permutations to run."
    )

    parser.add_argument(
        '--cluster-map-method',
        type = str,
        default = 'mean',
        choices = ['mean', 'median'],
        help = "Method to use to aggregate images within each cluster."
    )

    parser.add_argument(
        '--sim-gene-space',
        type = str,
        default = 'average-latent-space',
        choices = ['average-latent-space', 'latent-space', 'homologous-genes'],
        help = "Gene expression space to use when computing similarity."
    )

    parser.add_argument(
        '--sim-n-latent-spaces',
        type = int,
        default = 100,
        help = ("Number of latent spaces to include when --gene-space is "
                "'average-latent-space'. Ignored otherwise.")
    )

    parser.add_argument(
        '--sim-latent-space-id',
        type = int,
        default = 1,
        help = ("Latent space to use when --gene-space is 'latent-space'. "
                "Ignored otherwise.")
    )

    parser.add_argument(
        '--sim-metric',
        type = str,
        default = 'correlation',
        help = "Similarity metric."
    )

    parser.add_argument(
        '--sim-signed',
        type = str,
        default = 'true',
        help = ("Option to compute positive and negative similarity "
                "separately before averaging for a final value.")
    )

    parser.add_argument(
        '--sim-threshold',
        type = str,
        default = 'top_n',
        choices = ['none', 'top_n', 'intensity'],
        help = "Method used to threshold mouse and human images."
    )

    parser.add_argument(
        '--sim-threshold-value',
        type = float,
        default = 0.2,
        help = "Value used to threshold mouse and human images."
    )

    parser.add_argument(
        '--sim-threshold-symmetric',
        type = str,
        default = 'true',
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


# Main -----------------------------------------------------------------------

if __name__ == '__main__':

    args = parse_args()

    args['parallel'] = True if args['parallel'] == 'true' else False
    args['sim_signed'] = True if args['sim_signed'] == 'true' else False

    if args['sim_gene_space'] == 'average-latent-space':
        args['sim_latent_space_id'] = None
    elif args['sim_gene_space'] == 'latent-space':
        args['sim_n_latent_spaces'] = None
    else:
        args['sim_latent_space_id'] = None
        args['sim_n_latent_spaces'] = None

    if args['sim_threshold'] == 'none':
        args['sim_threshold'] = None
        args['sim_threshold_value'] = None
        args['sim_threshold_symmetric'] = None
    else:
        args['sim_threshold_symmetric'] = (
            True if args['sim_threshold_symmetric'] == 'true' else False
        )

    permute_cluster_similarity(**args)
