#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# process_human_images.py
# Author: Antoine Beauchamp
# Created: March 11th, 2024

"""
Execute the human image processing pipeline.

Description
-----------
The human image processing pipeline consists of three stages:
1. Compute effect size images for the desired participant group based
   on absolute and relative Jacobian images.
2. Generate clusters of participants based on absolute and relative
   Jacobian effect size images.
3. Generate cluster centroid images for all cluster solutions.

The pipeline can be deployed locally or on an HPC cluster that uses Slurm.
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
from shutil import rmtree
from process_human_images import initialize, effect_sizes, clustering, centroids


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
        default = 'data/human/derivatives/v3/',
        help = ("Path to the directory in which to export pipeline "
                "outputs. A uniquely identified sub-directory will be "
                "created using the specified set of pipeline parameters.")
    )

    parser.add_argument(
        '--input-dir',
        type = str,
        default = 'data/human/registration/v3/jacobians_resampled/resolution_0.8/',
        help = ("Path to the directory containing the input Jacobian "
                "MINC images. Expects sub-directories 'absolute' and "
                "'relative' containing the images.")
    )

    parser.add_argument(
        '--demographics',
        type = str,
        default = 'data/human/registration/v3/subject_info/demographics.csv',
        help = ("Path to the file (.csv) containing participant "
                "demographics information.")
    )

    parser.add_argument(
        '--mask',
        type = str,
        default = 'data/human/registration/v3/reference_files/mask_0.8mm.mnc',
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

    parser.add_argument(
        '--cv-n',
        type = int,
        help = "Number of cross-validation samples."
    )

    parser.add_argument(
        '--cv-start',
        type = int,
        default = 1,
        help = "Starting CV sample ID"
    )

    # Effect size arguments --------------------------------------------------

    parser.add_argument(
        '--es-method',
        type = str,
        default = 'normative-growth',
        choices = ['normative-growth', 'propensity-matching'],
        help = "The method to use to compute the effect size images."
    )

    parser.add_argument(
        '--es-group',
        type = str,
        default = 'controls',
        choices = ['patients', 'controls', 'all'],
        help = ("The group of participants for which to compute effect "
                "sizes.")
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
                "locally or using the Slurm scheduler on an HPC cluster.")
    )

    parser.add_argument(
        '--nproc',
        type = int,
        default = 1,
        help = "Number of processors to use when executing locally."
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
def main(pipeline_dir, input_dir, demographics, mask,
         datasets = ('POND', 'SickKids'), 
         cv_n = None, cv_start = 1,
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
    nproc: int, default 1
        Number of processors to use in parallel.
    execution: {'local', 'slurm'}
        Flag indicating whether the pipeline should be executed or
        using the Slurm scheduler on an HPC cluster.
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
    print("Initializing pipeline...", flush = True)
    paths = initialize(**kwargs)

    # Create cross-validation subdirectory
    cv_dir = os.path.join(paths['pipeline'], 'cross_validation')
    if not os.path.exists(cv_dir):
        os.makedirs(cv_dir)

    # Iterate over CV samples
    cv_start = kwargs['cv_start']
    cv_n = kwargs['cv_n']
    for sample in range(cv_start, cv_start + cv_n):

        print("Running cross-validation sample {} ...".format(sample))

        # Create CV pipeline output sub-directories
        cv_sample_dir = os.path.join(cv_dir, 'sample_{}'.format(sample), '')
        cv_sample_paths = dict(
            effect_sizes = paths['effect_sizes'].replace(paths['pipeline'], cv_sample_dir),
            clusters = paths['clusters'].replace(paths['pipeline'], cv_sample_dir),
            centroids = paths['centroids'].replace(paths['pipeline'], cv_sample_dir)
        )

        for path in cv_sample_paths.values():
            if not os.path.exists(path):
                os.makedirs(path)

        # Registry name for CV sample
        registry_name_sample = '{}_{}'.format(registry_name, sample)
    
        # Compute effect size images
        if stages['effect-sizes']:
            print("Computing effect sizes...", flush = True)
            es_kwargs = {key.replace('es_', ''):val
                        for key, val in kwargs.items() if 'es_' in key}
            es_kwargs.update(
                dict(imgdir = paths['jacobians'], demographics = paths['demographics'],
                    mask = mask, outdir = cv_sample_paths['effect_sizes'],
                    cv_seed = sample,
                    matrix_resolution = cluster_resolution,
                    execution = execution, nproc = nproc,
                    registry_name = registry_name_sample,
                    registry_cleanup = registry_cleanup,
                    slurm_njobs = slurm_njobs, slurm_mem = slurm_mem,
                    slurm_time = slurm_time)
            )
            es_outputs = effect_sizes(**es_kwargs)

        # Generate clusters
        if stages['clusters']:
            print("Generating clusters...", flush = True)
            cluster_kwargs = dict(
                infiles = [es_outputs[key]['matrix'] for key in es_outputs.keys()],
                nk_max = 2, metric = cluster_metric, K = cluster_K,
                sigma = cluster_sigma, t = cluster_t,
                cluster_file = os.path.join(cv_sample_paths['clusters'], cluster_file),
                affinity_file = os.path.join(cv_sample_paths['clusters'], cluster_affinity_file)
            )
            clusters = clustering(**cluster_kwargs)

        # Compute cluster centroid images
        if stages['centroids']:
            print("Generating cluster centroids...", flush = True)
            centroid_kwargs = dict(
                clusters = clusters, imgdir = cv_sample_paths['effect_sizes'],
                outdir = cv_sample_paths['centroids'], mask = mask,
                method = centroid_method,
                execution = execution, nproc = nproc,
                registry_name = registry_name_sample,
                registry_cleanup = registry_cleanup,
                slurm_mem = slurm_mem,
                slurm_time = slurm_time
            )
            centroid_outputs = centroids(**centroid_kwargs)

        rmtree(os.path.join(cv_sample_dir, 'effect_sizes'))

    print("Pipeline complete.", flush = True)

    return


# Execution ------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    if args['cv_n'] is None:
        raise ValueError("Specify the number of cross-validation samples.")
    args['datasets'] = tuple(args['datasets'])
    args['es_batch'] = tuple(args['es_batch'])
    args['stages'] = tuple(args['stages'])
    args['registry_cleanup'] = True if args['registry_cleanup'] == 'true' \
        else False
    main(**args)
