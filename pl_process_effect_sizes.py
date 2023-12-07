#!.venv/bin/python3
# ----------------------------------------------------------------------------
# process_effect_sizes.py
# Author: Antoine Beauchamp
# Created: November 29th, 2023

"""
Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import os
import processing
import utils


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
        '--params-id',
        type = str,
        help = ("ID specifying the human processing pipeline parameter "
                "set to use.")
    )

    parser.add_argument(
        '--mask',
        type = str,
        default = 'data/human/registration/v2/reference_files/mask_0.8mm.mnc',
        help = "Path to mask file (.mnc) associated with the study images."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        help = "Number of processors to use in parallel."
    )

    parser.add_argument(
        '--jacobians',
        type = str,
        choices = ['absolute', 'relative'],
        help = "Jacobians to use."
    )

    parser.add_argument(
        '--nbatches',
        type = int,
        default = 1,
        help = "Number of batches to use to process effect sizes."
    )

    parser.add_argument(
        '--matrix-file',
        type = str,
        default = 'effect_sizes.csv',
        help = ("The basename of the file (.csv) in which to store voxel-wise "
                "effect size matrices.")
    )

    args = vars(parser.parse_args())

    return args


if __name__ == '__main__':

    # Parse command line argument
    args = parse_args()

    # Extract args
    pipeline_dir = args['pipeline_dir']
    params_id = args['params_id']
    mask = args['mask']
    nproc = args['nproc']
    jacobians = args['jacobians']
    nbatches = args['nbatches']
    matrix_file = args['matrix_file']

    if jacobians is None:
        raise ValueError("jacobians must be one of {'absolute', 'relative'}")

    # Get pipeline metadata
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    params = utils.fetch_params_metadata(metadata, id = args['params_id'])
    if params is None:
        raise Exception("Pipeline not initialized: {}. "
                        "Run pl_process_initialize.py with desired parameters."
                        .format(params_id))
    if params.shape[0] > 1:
        raise Exception()

    # Check existence of pipeline directory
    pipeline_dir = os.path.join(pipeline_dir, params_id)
    if not os.path.exists(pipeline_dir):
        raise Exception("Pipeline not initialized: {}. "
                        "Run pl_process_initialize.py with desired parameters."
                        .format(params_id))

    # Extract pipeline parameters
    resolution = params['resolution'][0]
    cluster_resolution = params['cluster_resolution'][0]
    method = params['es_method'][0]
    group = params['es_group'][0]
    df = params['es_df'][0]
    batch = params['es_batch'][0]
    ncontrols = params['es_ncontrols'][0]

    # Input directory
    imgdir = os.path.join(pipeline_dir, 'jacobians', jacobians, '')

    # Effect size directory
    es_dir = os.path.join(pipeline_dir, 'effect_sizes',
                          'resolution_{}'.format(resolution))

    # Demographics file
    demographics = os.path.join(pipeline_dir, 'demographics.csv')

    # Arguments for effect size module
    kwargs = dict(imgdir = imgdir,
                  demographics = demographics,
                  mask = mask,
                  outdir = os.path.join(es_dir, jacobians, ''),
                  nproc = nproc,
                  method = method,
                  group = group)

    # Update arguments as needed
    if method == 'normative-growth':
        kwargs.update(
            dict(df = df,
                 batch = batch.split('-'),
                 nbatches = nbatches)
        )
    else:
        kwargs.update(
            dict(ncontrols = ncontrols)
        )

    # Calculate effect sizes
    es_files = processing.calculate_human_effect_sizes(**kwargs)

    # Resample effect size images to clustering resolution if needed
    if resolution != cluster_resolution:
        print("Downsampling effect sizes to {}mm...".format(
            cluster_resolution))
        es_dir_downsampled = es_dir.replace(
            'resolution_{}'.format(resolution),
            'resolution_3.0')
        es_files_downsampled = utils.resample_images(
            infiles = es_files,
            outdir = os.path.join(es_dir_downsampled, jacobians, ''),
            isostep = cluster_resolution,
            parallel = True if nproc > 1 else False,
            nproc = nproc
        )
    else:
        es_dir_downsampled = es_dir
        es_files_downsampled = es_files

    # Build effect size matrix
    print("Building effect size matrix...")
    df_es = processing.build_voxel_matrix(imgfiles = es_files_downsampled,
                                          mask = mask, file_col = True,
                                          sort = True, parallel = True,
                                          nproc = nproc)
    df_es['file'] = [os.path.basename(file) for file in df_es['file']]
    df_es.to_csv(os.path.join(es_dir_downsampled, jacobians, matrix_file),
                 index = False)
