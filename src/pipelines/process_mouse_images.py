#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# process_mouse_images.py
# Author: Antoine Beauchamp
# Created: July 16th, 2025

"""
Execute the mouse image processing pipeline.

Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import os
import utils
import pandas as pd
from shutil import rmtree
from pyminc.volumes.factory import volumeFromFile

# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # General arguments ------------------------------------------------------
    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/mouse/derivatives/',
        help = ("Path to the directory in which to export pipeline "
                "outputs. A uniquely identified sub-directory will be "
                "created using the specified set of pipeline parameters.")
    )

    parser.add_argument(
        '--params-id',
        type = str,
        help = ("Optional ID for pipeline parameter set.")
    )

    parser.add_argument(
        '--input-dir',
        type = str,
        default = 'data/mouse/registration/',
        help = ("Path to the directory containing input Jacobian "
                "MINC images.")
    )

    parser.add_argument(
        '--models',
        type = str,
        default = 'data/mouse/registration/Names_Paper.csv',
        help = ("Path to the file (.csv) containing mouse model "
                "information.")
    )

    parser.add_argument(
        '--mask',
        type = str,
        default = 'data/mouse/registration/reference_files/scanbase_second_level-nlin-3_mask_200um.mnc',
        help = "Path to the mask file (.mnc)."
    )

    # Clustering arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-resolution',
        type = float,
        default = 0.2,
        help = ("The isotropic resolution (mm) at which to run the clustering "
                "stage. Effect size images will be resampled to this "
                "resolution if it is different from the resolution of the "
                "input images.")
    )

    parser.add_argument(
        '--cluster-nk-max',
        type = int,
        default = 10,
        help = ("The maximum number of clusters to identify when clustering. "
                "Solutions will be obtained from 2 clusters to the value "
                "provided.")
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
        help = ("The basename of the file (.csv) in which to store the "
                "affinity matrix produced by similarity network fusion.")
    )

    # Centroid arguments -----------------------------------------------------
    parser.add_argument(
        '--centroid-method',
        type = str,
        default = 'mean',
        choices = ['mean', 'median'],
        help = "The method to use to compute cluster centroid images."
    )

    parser.add_argument(
        '--transform',
        type = str,
        default = 'data/mouse/registration/transforms/MICe_scanbase.xfm',
        help = "Path to transform file (.xfm)."
    )

    parser.add_argument(
        '--transform-like',
        type = str,
        default = 'data/mouse/atlas/average_template_200um.mnc',
        help = "Path to transform likefile (.mnc)."
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use."
    )

    return vars(parser.parse_args())



# Modules --------------------------------------------------------------------

def initialize(**kwargs):
    """
    Initialize the pipeline

    Parameters
    ----------
    kwargs: dict
        All arguments passed to the main() function.

    Returns
    -------
    paths: dict
        Dictionary containing the paths to the pipeline sub-directories.
    """

    # Fetch the image resolution
    vol = volumeFromFile(kwargs['mask'])
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Pipeline parameters
    params = dict(
        dataset = 'MICe',
        resolution = f'{resolution:.1f}',
        cluster_resolution = f'{kwargs["cluster_resolution"]:.1f}',
        cluster_nk_max = kwargs['cluster_nk_max'],
        cluster_metric = kwargs['cluster_metric'],
        cluster_K = kwargs['cluster_K'],
        cluster_sigma = kwargs['cluster_sigma'],
        cluster_t = kwargs['cluster_t'],
        centroid_method = kwargs['centroid_method']
    )

    # Extract the pipeline and input directories
    pipeline_dir = kwargs['pipeline_dir']
    input_dir = kwargs['input_dir']

    # Create pipeline directory based on the pipeline parameters
    print("Creating pipeline directories...", flush = True)
    if kwargs['params_id'] is None:
        params_id = utils.random_id(3)
    else:
        params_id = kwargs['params_id']
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir,
                                           params_id = params_id)
    params_id = utils.fetch_params_id(metadata = metadata,
                                      params = params)

    # Define sub-directories for pipeline stages
    es_dir = os.path.join(pipeline_dir, 'effect_sizes', '')
    cluster_dir = os.path.join(pipeline_dir, 'clusters',
                               'resolution_{}'.format(
                                   kwargs['cluster_resolution']),
                               '')
    centroid_dir = os.path.join(pipeline_dir, 'centroids',
                                'resolution_{}'.format(resolution), '')

    # Check existence of the input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: {}".format(input_dir))

    # Create pipeline sub-directories
    if not os.path.exists(es_dir):
        os.makedirs(es_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    if not os.path.exists(centroid_dir):
        os.makedirs(centroid_dir)

    # Copy model names file into pipeline directory
    df_models = pd.read_csv(kwargs['models'], index_col = 0)
    df_models.columns = ['file', 'ID']
    df_models['file'] = df_models['file'] + '.mnc'
    df_models.to_csv(os.path.join(pipeline_dir, 'model_names.csv'),
                     index = False)

    # Dictionary containing pipeline paths
    paths = dict(
        pipeline = pipeline_dir,
        effect_sizes = es_dir,
        clusters = cluster_dir,
        centroids = centroid_dir,
    )

    return paths


def main(pipeline_dir, input_dir, models, mask, params_id = None,
         cluster_resolution = 0.2, cluster_nk_max = 10,
         cluster_metric = 'correlation', cluster_K = 10,
         cluster_sigma = 0.5, cluster_t = 20,
         cluster_file = 'clusters.csv',
         cluster_affinity_file = 'affinity.csv',
         centroid_method = 'mean',
         transform = None, transform_like = None, nproc = 1):

    # Initialize pipeline
    kwargs = locals().copy()
    paths = initialize(**kwargs)

    del kwargs['models']
    del kwargs['params_id']
    del kwargs['transform']
    del kwargs['transform_like']
    del kwargs['nproc']
    kwargs['pipeline_dir'] = paths['pipeline']
    kwargs['registration_dir'] = kwargs.pop('input_dir')
    kwargs = {key.replace('_', '-'): val for key, val in kwargs.items()}

    # Execute processing pipeline
    print("Executing processing pipeline...", flush = True)
    script = 'process_mouse_images.R'
    utils.execute_local(script = script, kwargs = kwargs)

    # Transform images
    if transform is not None:
        print("Transforming centroids...", flush = True)
        if transform_like is None:
            raise Exception("Argument transform_like must be specified "
                            "when using transform.")

        # Iterate over Jacobians
        jacobians = ('absolute', 'relative')
        for jtype in jacobians:

            # Temporary centroids directory
            centroid_dir_tmp_j = os.path.join(paths['pipeline'], 'centroids_tmp', jtype)
            if not os.path.exists(centroid_dir_tmp_j):
                os.makedirs(centroid_dir_tmp_j)

            # Get input centroid images
            centroid_dir_j = os.path.join(paths['centroids'], jtype)
            input_files = os.listdir(centroid_dir_j)
            input_files = [os.path.join(centroid_dir_j, file)
                           for file in input_files]

            # Transform centroid images
            output_files = utils.transform_images(infiles = input_files,
                                                  outdir = centroid_dir_tmp_j,
                                                  like = transform_like,
                                                  transform = transform,
                                                  nproc = nproc)

            for file in output_files:
                os.rename(file, os.path.join(centroid_dir_j, os.path.basename(file)))

        rmtree(os.path.join(paths['pipeline'], 'centroids_tmp'))

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    main(**args)