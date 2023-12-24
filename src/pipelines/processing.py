import os
import sys
import utils
import pandas as pd
from glob import glob
from pyminc.volumes.factory import volumeFromFile


def initialize(**kwargs):
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
        resolution = resolution,
        es_method = kwargs['es_method'],
        es_group = kwargs['es_group'],
        es_df = kwargs['es_df'],
        es_batch = (None if kwargs['es_batch'] is None
                    else '-'.join(kwargs['es_batch'])),
        es_ncontrols = kwargs['es_ncontrols'],
        cluster_resolution = kwargs['cluster_resolution'],
        cluster_nk_max = kwargs['cluster_nk_max'],
        cluster_metric = kwargs['cluster_metric'],
        cluster_K = kwargs['cluster_K'],
        cluster_sigma = kwargs['cluster_sigma'],
        cluster_t = kwargs['cluster_t'],
        cluster_map_method = kwargs['cluster_map_method']
    )

    pipeline_dir = kwargs['pipeline_dir']
    input_dir = kwargs['input_dir']

    # Create pipeline directory
    params_id = utils.random_id(3)
    metadata = os.path.join(pipeline_dir, 'metadata.csv')
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir,
                                           params_id = params_id)
    params_id = utils.fetch_params_id(metadata = metadata,
                                      params = params)

    # Directories for pipeline stages
    imgdir = os.path.join(pipeline_dir, 'jacobians', '')
    es_dir = os.path.join(pipeline_dir, 'effect_sizes',
                          'resolution_{}'.format(resolution), '')
    cluster_dir = os.path.join(pipeline_dir, 'clusters',
                               'resolution_{}'.format(
                                   kwargs['cluster_resolution']),
                               '')
    cluster_map_dir = os.path.join(pipeline_dir, 'cluster_maps',
                                   'resolution_{}'.format(resolution), '')

    # Check existence of input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: ".format(input_dir))

    # Create pipeline sub-directories
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)
    if not os.path.exists(es_dir):
        os.makedirs(es_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    if not os.path.exists(cluster_map_dir):
        os.makedirs(cluster_map_dir)

    # Filter for data sets ----------------------------------------------------

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
            raise OSError("No input files in directory: ".format(input_dir))
        input_files_in_dataset = [[f for f in input_files if g in f][0]
                                  for g in df_demographics['file'].to_list()]
        imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                     dst = os.path.join(imgdir, jac, ''))

    # Dictionary containing pipeline paths
    paths = dict(
        pipeline = pipeline_dir,
        jacobians = imgdir,
        effect_sizes = es_dir,
        clusters = cluster_dir,
        centroids = cluster_map_dir
    )

    return paths


def clustering():
    if env == 'local':
        utils.execute_local(script = 'generate_clusters.R',
                            args = ...)
    elif env == 'slurm':
        utils.execute_slurm(script = 'generate_clusters.R',
                            args = ...,
                            slurm_args = ...)
    else:
        raise ValueError

    return


def centroids():
    return


def main(pipeline_dir = 'data/human/derivatives/v2/',
         input_dir = 'data/human/registration/v2/jacobians_resampled/resolution_0.8/',
         demographics = 'data/human/registration/v2/subject_info/demographics.csv',
         mask = 'data/human/registration/v2/reference_files/mask_0.8mm.mnc',
         datasets = ('POND', 'SickKids'),
         parallel = True, nproc = None,
         es_method = 'normative-growth', es_group = 'patients',
         es_nbatches = 1, es_df = 3,
         es_batch = ('Site', 'Scanner'), es_ncontrols = 10,
         es_matrix_file = 'effect_sizes.csv',
         cluster_resolution = 3.0,
         cluster_nk_max = 10, cluster_metric = 'correlation',
         cluster_K = 10, cluster_sigma = 0.5, cluster_t = 20,
         cluster_file = 'clusters.csv',
         cluster_affinity_file = 'affinity.csv',
         cluster_map_method = 'mean'):
    # Get dictionary of function kwargs
    kwargs = locals().copy()

    # Initialize pipeline directory
    # Get paths to pipeline
    paths = initialize(**kwargs)

    # Compute effect sizes
    def effect_sizes(imgdir, demographics, mask, outdir,
                     method = 'normative-growth', env = 'local',
                     nproc = 1, **kwargs):

        # Create output dir if needed
        imgdir = os.path.join(imgdir, '')
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        jacobians = ['absolute', 'relative']
        if env == 'local':
            for j, jac in enumerate(jacobians):
                utils.execute_local(script = 'compute_effect_sizes.R',
                                    args = ...)
        elif env == 'slurm':
            utils.execute_slurm(script = 'compute_effect_sizes.R',
                                args = ...,
                                slurm_args = ...)
        else:
            raise ValueError

        return

    # Arguments for effect sizes
    es_kwargs = {key.replace('es_', ''): val
                 for key, val in kwargs.items() if 'es_' in key}
    es_kwargs.update(
        dict(imgdir = paths['jacobians'],
             demographics = demographics,
             mask = mask,
             outdir = paths['effect_sizes'])
    )
    effect_sizes(**es_kwargs)

    slurm_args = dict(
        job_name = 'process_human_data_v2_0.8mm',
        nodes = 1,
        cpus_per_task = 8,
        mem = '64G',
        time = '72:00:00',
        chdir = os.getcwd(),
        output = 'logs/process_human_data_v2_0.8mm_%j.out'
    )

    # clustering()
    # centroids()
