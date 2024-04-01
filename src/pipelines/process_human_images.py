#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# process_human_images.py
# Author: Antoine Beauchamp
# Created: March 11th, 2024

"""
Execute human image processing pipeline.

Description
-----------
"""


# Packages -------------------------------------------------------------------

import argparse
import os
import sys
import utils
import processing
import pandas as pd
from glob import glob
from pyminc.volumes.factory import volumeFromFile


# Command line arguments -----------------------------------------------------

def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    # General arguments ------------------------------------------------------
    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/human/derivatives/v2/',
        help = ("Path to the directory in which to store pipeline outputs. "
                "A sub-directory will be created based on the parameters "
                "used to run the pipeline.")
    )

    parser.add_argument(
        '--input-dir',
        type = str,
        default = 'data/human/registration/v2/jacobians_resampled/resolution_0.8/',
        help = ("Path to the directory containing Jacobian images. "
                "This directory must contain sub-directories named 'absolute' "
                "and 'relative' containing the Jacobian images.")
    )

    parser.add_argument(
        '--demographics',
        type = str,
        default = 'data/human/registration/v2/subject_info/demographics.csv',
        help = "Path to the file (.csv) containing demographics information."
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

    # Effect size arguments --------------------------------------------------
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
        '--es-df',
        type = int,
        default = 3,
        help = ("The number of degrees of freedom for the natural splines used "
                "in normative growth modelling. Ignored if --es-method is "
                "'propensity-matching'.")
    )

    parser.add_argument(
        '--es-batch',
        nargs = '*',
        type = str,
        default = ['Site', 'Scanner'],
        help = ("Batch variables to use for normalization prior to normative "
                "growth modelling. Variables must be found in --demographics.")
    )

    parser.add_argument(
        '--es-ncontrols',
        type = int,
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
        help = ("The resolution (mm) at which to run the clustering stage. "
                "Effect size images will be resampled to this resolution "
                "if it is different from the resolution of the input images.")
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

    # Execution arguments ----------------------------------------------------
    parser.add_argument(
        '--stages',
        nargs = '*',
        type = str,
        default = ['effect-sizes', 'clusters', 'centroids'],
        help = "List of strings indicating which pipeline stages to execute."
    )

    parser.add_argument(
        '--execution',
        type = str,
        default = 'local',
        choices = ['local', 'slurm'],
        help = ("Flag indicating whether the pipeline should be executed "
                "or using the Slurm scheduler on a HPC cluster.")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use in parallel."
    )

    parser.add_argument(
        '--registry-name',
        type = str,
        default = "process_human_images_registry",
        help = "Name of the registry directory for batched jobs."
    )

    parser.add_argument(
        '--registry-cleanup',
        type = str,
        default = "true",
        help = "Option to clean up registry after completion of batched jobs."
    )

    parser.add_argument(
        '--slurm-njobs',
        type = int,
        help = "Number of jobs to deploy on Slurm."
    )

    parser.add_argument(
        '--slurm-mem',
        type = str,
        help = "Memory per CPU."
    )

    parser.add_argument(
        '--slurm-time',
        type = int,
        help = "Walltime in minutes for Slurm jobs."
    )

    return vars(parser.parse_args())


# Modules --------------------------------------------------------------------

@utils.timing
def initialize(**kwargs):
    """
    Initialize the image processing pipeline.

    Parameters
    ----------
    kwargs: dict
        All arguments passed to the main() function.

    Returns
    -------
    paths: list of str
        List of paths to the pipeline sub-directories
    """
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
        resolution = f'{resolution:.1f}',
        es_method = kwargs['es_method'],
        es_group = kwargs['es_group'],
        es_df = kwargs['es_df'],
        es_batch = (None if kwargs['es_batch'] is None
                    else '-'.join(kwargs['es_batch'])),
        es_ncontrols = kwargs['es_ncontrols'],
        cluster_resolution = f'{kwargs["cluster_resolution"]:.1f}',
        cluster_nk_max = kwargs['cluster_nk_max'],
        cluster_metric = kwargs['cluster_metric'],
        cluster_K = kwargs['cluster_K'],
        cluster_sigma = kwargs['cluster_sigma'],
        cluster_t = kwargs['cluster_t'],
        centroid_method = kwargs['centroid_method']
    )

    pipeline_dir = kwargs['pipeline_dir']
    input_dir = kwargs['input_dir']

    # Create pipeline directory based on the pipeline parameters
    print("Creating pipeline directories...")
    params_id = utils.random_id(3)
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir,
                                           params_id = params_id)
    params_id = utils.fetch_params_id(metadata = metadata,
                                      params = params)

    # Define sub-directories for pipeline stages
    imgdir = os.path.join(pipeline_dir, 'jacobians', '')
    es_dir = os.path.join(pipeline_dir, 'effect_sizes',
                          'resolution_{}'.format(resolution), '')
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
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)
    if not os.path.exists(es_dir):
        os.makedirs(es_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    if not os.path.exists(centroid_dir):
        os.makedirs(centroid_dir)

    # Filter for data sets ---------------------------------------------------

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
            raise OSError("No input files in directory: {}".format(input_dir))
        input_files_in_dataset = [[f for f in input_files if g in f][0]
                                  for g in df_demographics['file'].to_list()]
        imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                     dst = os.path.join(imgdir, jac, ''))

    # Dictionary containing pipeline paths
    paths = dict(
        pipeline = pipeline_dir,
        demographics = demographics,
        jacobians = imgdir,
        effect_sizes = es_dir,
        clusters = cluster_dir,
        centroids = centroid_dir
    )

    return paths


@utils.timing
def effect_sizes(imgdir, demographics, mask, outdir,
                 method = 'normative-growth', group = 'patients',
                 df = 3, batch = ('Site', 'Scanner'), ncontrols = None,
                 matrix_file = 'effect_sizes.csv', matrix_resolution = 3.0,
                 execution = 'local', nproc = 1,
                 registry_name = None, registry_cleanup = True,
                 slurm_njobs = None, slurm_mem = None, slurm_time = None):
    """

    Parameters
    ----------
    imgdir: str
        Path to the directory containing the input Jacobian images. Must
        contain sub-directories 'absolute' and 'relative', which hold the
        image files (.mnc).
    demographics: str
        Path to the file (.csv) containing the demographics information
    mask: str
        Path to the mask file (.mnc)
    outdir: str
        Path to the directory in which to export the effect size sub-directories
        and images.
    method: {'normative-growth', 'propensity-matching'}
        Method used to compute the effect size images.
    group: {'patients', 'controls', 'all'}
        Group of participants for which to compute effect sizes.
    df: int, default 3
        Number of degrees of freedom to use when `method`='normative-growth'
    batch: str or tuple of str, default ('Site', 'Scanner')
        Batch variables to normalize against when `method` =
        'normative-growth'
    ncontrols: int, default 10
        Number of controls to use for propensity matching when `method` =
        'propensity-matching'
    matrix_file: str, default 'effect_sizes.csv'
        Basename of the file (.csv) in which to write the absolute and relative
        effect size matrices.
    matrix_resolution: float, default 3.0
        Resolution (mm) of the images used for computing the effect size
        matrices. Images will be resampled to this resolution if it is not
        equal to the resolution of the input images.
    execution: {'local', 'slurm'}
        Flag indicating whether to run locally or using Slurm.
    nproc: int, default 1
        Number of processors to use. Executed in parallel if > 1.
    registry_name: str, default None
        Name of the registry directory for batched jobs.
    registry_cleanup: bool, default True
        Option to clean up registry after completion of batched jobs.
    slurm_njobs: int, default None
        Number of jobs to deploy when execution = 'slurm'.
    slurm_mem: str, default None
        Memory per CPU when execution = 'slurm', e.g. '16G'.
    slurm_time: int, default = None
        Walltime (minutes) for Slurm jobs when execution = 'slurm'

    Returns
    -------
    out: dict of str
    """

    # Clean up arguments
    kwargs = locals().copy()
    kwargs = {key.replace('_', '-'):val for key, val in kwargs.items()}
    kwargs['batch'] = (None if kwargs['batch'] is None
                       else '-'.join(kwargs['batch']))
                       
    # Driver script
    script = 'compute_effect_sizes.R'

    # Dictionary to store outputs
    out = dict(
        absolute = dict(imgdir = '',
                        matrix = ''),
        relative = dict(imgdir = '',
                        matrix = '')
    )

    # Iterate over Jacobians
    for j in out.keys():

        print("Computing {} effect size images...".format(j))
        kwargs['imgdir'] = os.path.join(imgdir, j, '')
        kwargs['outdir'] = os.path.join(outdir, j, '')
        if registry_name is not None:
            registry_name_es = registry_name + "_es_" + j
            kwargs['registry-name'] = registry_name_es
        kwargs['registry-cleanup'] = "true" if registry_cleanup else "false"
        utils.execute_local(script = script, kwargs = kwargs)

        # Create the effect size matrix
        print("Building {} effect size matrix...".format(j))

        # Path to effect size images
        imgfiles = os.listdir(os.path.join(outdir, j, ''))
        imgfiles = [file for file in imgfiles if '.mnc' in file]
        imgfiles = [os.path.join(outdir, j, file) for file in imgfiles]

        # If the matrix resolution is specified, resample to the matrix res
        if matrix_resolution is None:
            outdir_f = outdir
            mask_f = mask
            imgfiles_f = imgfiles

        else:

            # Get the image resolution
            vol = volumeFromFile(mask)
            resolution = vol.getSeparations()
            if len(set(resolution)) == 1:
                resolution = resolution[0]
            else:
                raise Exception
            vol.closeVolume()

            # Resample if needed
            if (matrix_resolution != resolution):

                print("Resampling effect size images to {}mm..."
                      .format(matrix_resolution))

                # Directory for resampled images
                outdir_f = outdir.replace(
                    'resolution_{}'.format(resolution),
                    'resolution_{}'.format(matrix_resolution)
                )

                # Resample images
                imgfiles_f = utils.resample_images(
                    infiles = imgfiles,
                    outdir = os.path.join(outdir_f, j, ''),
                    isostep = matrix_resolution,
                    nproc = nproc
                )

                # Resample mask
                mask_f = utils.resample_image(
                    infile = mask,
                    isostep = matrix_resolution,
                    outdir = outdir_f,
                    suffix = f'_autocrop_{matrix_resolution:.1f}mm'
                )

            else:
                outdir_f = outdir
                mask_f = mask
                imgfiles_f = imgfiles

        # Build effect size matrix and export
        df_es = processing.build_voxel_matrix(imgfiles = imgfiles_f,
                                              mask = mask_f,
                                              file_col = True,
                                              sort = True,
                                              version = "v1",
                                              parallel = True,
                                              nproc = nproc)
        df_es['file'] = [os.path.basename(file) for file in df_es['file']]
        df_es.to_csv(os.path.join(outdir_f, j, matrix_file), index = False)

        # Add outputs to dictionary
        out[j]['imgdir'] = os.path.join(outdir, j, '')
        out[j]['matrix'] = os.path.join(outdir_f, j, matrix_file)

    return out


@utils.timing
def clustering(infiles, rownames = 'file', nk_max = 10,
               metric = 'correlation', K = 10, sigma = 0.5, t = 20,
               cluster_file = 'clusters.csv', affinity_file = 'affinity.csv'):
    """
    Identify clusters of patients.

    Parameters
    ----------
    infiles: list of str
        Paths to the files (.csv) containing the effect size matrices.
    rownames: str, default 'file'
        Column in the input files containing row names.
    nk_max: int, default 10
        Maximum number of clusters to identify. Solutions will be
        obtained for nk = 2 to nk = nk_max
    metric: str, default 'correlation'
        Similarity metric used to compute the SNF affinity matrices.
    K: int, default 10
        Number of nearest-neighbours used to compute the SNF affinity
        matrices.
    sigma: float, default 0.5
        Variance for the local model in the SNF affinity matrices.
    t: int, default 20
        Number of iterations for the diffusion process in SNF.
    cluster_file: str, default 'clusters.csv'
        Path to the file (.csv) in which to save the cluster assignments.
    affinity_file: str, default 'affinity.csv'
        Path to the file (.csv) in which to save the SNF affinity matrix.

    Returns
    -------
    cluster_file: str
        Path to the cluster file.
    """

    # Clean up arguments
    kwargs = locals().copy()
    kwargs = {key.replace('_', '-'):val for key, val in kwargs.items()}
    kwargs['file1'] = kwargs['infiles'][0]
    kwargs['file2'] = kwargs['infiles'][1]
    del kwargs['infiles']

    # Driver script
    script = 'generate_clusters.R'
    utils.execute_local(script = script, kwargs = kwargs)

    return cluster_file


@utils.timing
def centroids(clusters, imgdir, outdir, mask,
              method = 'mean', execution = 'local', nproc = 1,
              registry_name = None, registry_cleanup = True,
              slurm_mem = None, slurm_time = None):
    """
    Compute cluster centroid images.

    Parameters
    ----------
    clusters: str
        Path to the file (.csv) containing cluster assignments.
    imgdir: str
        Path to the directory containing the individual images (.mnc).
    outdir: str
        Path to the output directory.
    mask: str
        Path to a mask image (.mnc).
    method: {'mean', 'median'}
        Method used to compute the cluster centroids.
    execution: {'local', 'slurm'}
        Flag indicating whether to run locally or using Slurm
    nproc: int, default 1
        Number of processors to use. Executed in parallel if > 1.
    registry_name: str, default None
        Name of the registry directory for batched jobs.
    registry_cleanup: bool, default True
        Option to clean up registry after completion of batched jobs.
    slurm_mem: str, default None
        Memory per CPU when execution = 'slurm', e.g. '16G'.
    slurm_time: int, default = None
        Walltime (minutes) for Slurm jobs when execution = 'slurm'

    Returns
    -------
    out: str
        The directory containing the centroid image files.
    """

    # Script args
    kwargs = locals().copy()
    kwargs = {key.replace('_', '-'):val for key, val in kwargs.items()}
    kwargs['cluster-file'] = kwargs.pop('clusters')

    # Driver script
    script = 'compute_cluster_centroids.R'

    # Dictionary to store outputs
    out = dict(absolute = '', relative = '')

    # Iterate over Jacobians
    for j in out.keys():
        print("Computing {} cluster centroid images...".format(j))
        kwargs['imgdir'] = os.path.join(imgdir, j, '')
        kwargs['outdir'] = os.path.join(outdir, j, '')
        if registry_name is not None:
            registry_name_es = registry_name + "_centroids_" + j
            kwargs['registry-name'] = registry_name_es
        kwargs['registry-cleanup'] = "true" if registry_cleanup else "false"
        utils.execute_local(script = script, kwargs = kwargs)
        out[j] = os.path.join(outdir, j, '')

    return out


@utils.timing
def main(pipeline_dir, input_dir, demographics, mask,
         datasets = ('POND', 'SickKids'),
         es_method = 'normative-growth', es_group = 'patients',
         es_df = 3, es_batch = ('Site', 'Scanner'), es_ncontrols = 10,
         es_matrix_file = 'effect_sizes.csv',
         cluster_resolution = 3.0,
         cluster_nk_max = 10, cluster_metric = 'correlation',
         cluster_K = 10, cluster_sigma = 0.5, cluster_t = 20,
         cluster_file = 'clusters.csv',
         cluster_affinity_file = 'affinity.csv',
         centroid_method = 'mean',
         stages = ('effect-sizes', 'clusters', 'centroids'),
         execution = 'local', nproc = 1,
         registry_name = None, registry_cleanup = True,
         slurm_njobs = None, slurm_mem = None, slurm_time = None):
    """
    Execute the image processing pipeline.

    Parameters
    ----------
    pipeline_dir: str
        Path to the directory in which to store pipeline outputs.
    input_dir: str
        Path to the directory containing Jacobian images.
    demographics: str
        Path to the file (.csv) containing the demographics information.
    mask: str
        Path to the mask file (.mnc) for the registration.
    datasets: tuple of str
        Datasets to include in processing.
    es_method: {'normative-growth', 'propensity-matching'}
        Method used to compute the effect size images.
    es_group: {'patients', 'controls', 'all'}
        Group of participants for which to compute effect sizes.
    es_df: int, default 3
        Number of degrees of freedom to use when `es_method` =
        'normative-growth'
    es_batch: str or tuple of str
        Batch variables to normalize against when `es_method` =
        'normative-growth'
    es_ncontrols: int, default 10
        Number of controls to use for propensity matching when
        `es_method` = 'propensity-matching'
    es_matrix_file: str, default 'effect_sizes.csv'
        Basename of the file (.csv) in which to write the absolute and
        relative effect size matrices.
    cluster_resolution: float, default 3.0
        Resolution (mm) of effect size images at which to execute
        clustering.
    cluster_nk_max: int, default 10
        Maximum number of clusters to identify. Clustering solutions
        are identified from nk = 2 to nk = `cluster_nk_max`
    cluster_metric: str, default 'correlation'
        Distance metric used to generate affinity matrices during
        clustering.
    cluster_K: int, default 10
        Number of nearest-neighbours to consider when building affinity
        matrices during clustering.
    cluster_sigma: float, default 0.5
        Variance of the local model used to generate affinity matrices
        during clustering.
    cluster_t: int, default 20
        Number of iterations for the diffusion process in similarity
        network fusion during clustering.
    cluster_file: str, default 'clusters.csv'
        Basename of the file (.csv) in which to write the cluster
        assignments.
    cluster_affinity_file: str, default 'affinity.csv'
        Basename of the file (.csv) in which to write the fused
        affinity matrix.
    centroid_method: {'mean', 'median'}
        Method used to compute cluster centroid images.
    stages: tuple of str
        Strings indicating which pipeline stages to execute. Must be a
        combination of 'effect-sizes', 'clusters', and 'centroids'.
    execution: {'local', 'slurm'}
        Flag indicating whether the pipeline should be executed or
        using the Slurm scheduler on an HPC cluster.
    nproc: int, default 1
        Number of processors to use in parallel.
    registry_name: str, default None
        Name of the registry directory to use in batched jobs.
    registry_cleanup: bool, default True
        Option to clean up registry after completion of batched jobs.
    slurm_njobs: int, default None
        Number of jobs to deploy on Slurm
    slurm_mem: str, default None
        Memory per CPU for Slurm jobs.
    slurm_time: int, default = None
        Walltime (minutes) for Slurm jobs.

    Returns
    -------
    None
    """

    # Get dictionary of function kwargs
    kwargs = locals().copy()

    # Error if > 3 stages passed
    if len(stages) > 3:
        raise Exception("Got more than three stages: {}".format(stages))

    # Error if stages not in options list
    stages_opt = ('effect-sizes', 'clusters', 'centroids')
    for stage in stages:
        if stage not in stages_opt:
            raise ValueError("Pipeline stage '{}' not recognized."
                             .format(stage))

    # Generate stages flags
    stages = {stage:True if stage in stages else False
              for stage in stages_opt}

    # Initialize pipeline directory tree
    print("Initializing pipeline...")
    paths = initialize(**kwargs)
    
    # Compute effect size images
    if stages['effect-sizes']:
        print("Computing effect sizes...")
        es_kwargs = {key.replace('es_', ''):val
                     for key, val in kwargs.items() if 'es_' in key}
        es_kwargs.update(
            dict(imgdir = paths['jacobians'], demographics = paths['demographics'],
                 mask = mask, outdir = paths['effect_sizes'],
                matrix_resolution = cluster_resolution,
                execution = execution, nproc = nproc,
                registry_name = registry_name,
                registry_cleanup = registry_cleanup,
                slurm_njobs = slurm_njobs, slurm_mem = slurm_mem,
                slurm_time = slurm_time)
        )
        es_outputs = effect_sizes(**es_kwargs)

    # Generate clusters
    if stages['clusters']:
        print("Generating clusters...")
        cluster_kwargs = dict(
            infiles = [es_outputs[key]['matrix'] for key in es_outputs.keys()],
            nk_max = cluster_nk_max, metric = cluster_metric, K = cluster_K,
            sigma = cluster_sigma, t = cluster_t,
            cluster_file = os.path.join(paths['clusters'], cluster_file),
            affinity_file = os.path.join(paths['clusters'], cluster_affinity_file)
        )
        clusters = clustering(**cluster_kwargs)

    # Compute cluster centroid images
    if stages['centroids']:
        print("Generating cluster centroids...")
        centroid_kwargs = dict(
            clusters = clusters, imgdir = paths['effect_sizes'],
            outdir = paths['centroids'], mask = mask,
            method = centroid_method,
            execution = execution, nproc = nproc,
            registry_name = registry_name,
            registry_cleanup = registry_cleanup,
            slurm_mem = slurm_mem,
            slurm_time = slurm_time
        )
        centroid_outputs = centroids(**centroid_kwargs)

    print("Pipeline complete.")

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    args['datasets'] = tuple(args['datasets'])
    args['es_batch'] = tuple(args['es_batch'])
    args['stages'] = tuple(args['stages'])
    args['registry_cleanup'] = True if args['registry_cleanup'] == 'true' \
        else False
    main(**args)
