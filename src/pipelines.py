import os
import processing
import transcriptomic
import utils
import pandas as pd
from glob import glob
from itertools import product
from shutil import rmtree, copyfile
from pyminc.volumes.factory import volumeFromFile
import sys


def fetch_mouse_clustering_outputs(pipeline_dir, input_dir, resolution = 200,
                                   method = 'mean', transform = None,
                                   transform_like = None, parallel = False,
                                   nproc = None):
    """
    Fetch mouse clustering outputs from Jacob's directory.

    Arguments
    ---------
    pipeline_dir: str
        Path to output directory.
    input_dir:
        Path to input directory.
    resolution: float
        Resolution (um) of images to link.
    method: str
        Method used to aggregate cluster maps.
    transform: str
        Path to transform file (.xfm) to use.
    transform_like: str
        Path to transform likefile (.mnc).
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    None
    """

    # Create pipeline directories ---------------------------------------------

    # Path to cluster maps directory
    resolution_mm = float(resolution) / 1000

    params = dict(resolution = resolution_mm,
                  cluster_map_method = method)
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir)

    # Path to clustering directory
    cluster_dir = os.path.join(pipeline_dir, 'clusters', '')
    cluster_map_dir = os.path.join(pipeline_dir, 'cluster_maps', 'resolution_{}'.format(resolution_mm), '')

    # Check existence of input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: {}".format(input_dir))

    # Get cluster files -------------------------------------------------------

    # Make link to cluster file
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    #Clusters and affinity matrix
    cluster_file = os.path.join(input_dir, 'Clusters.csv')
    affinity_file = '/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/SNFMatrix_Paper/WMatrix.RData'

    #Copy files
    copyfile(cluster_file, os.path.join(cluster_dir, 'clusters.csv'))
    copyfile(affinity_file, os.path.join(cluster_dir, 'affinity.RData'))

    # Iterate over jacobians
    jacobians = ['absolute', 'relative']
    for jac in jacobians:

        print("Getting {} Jacobian cluster maps...".format(jac))

        # Cluster map directory
        cluster_map_dir_jac = os.path.join(cluster_map_dir, jac, '')
        if not os.path.exists(cluster_map_dir_jac):
            os.makedirs(cluster_map_dir_jac)

        # Get input files
        jac_short = jac[:3]
        input_files = os.listdir(input_dir)
        input_files = [os.path.join(input_dir, file) for file in input_files]
        input_files = [file for file in input_files if '.mnc' in file]
        input_files = [file for file in input_files if str(resolution) in file]
        input_files = [file for file in input_files if jac_short in file]
        input_files = [file for file in input_files if method in file]
        if len(input_files) == 0:
            raise OSError("Input files not found for specified parameters.")

        # Option to transform images
        if transform is None:
            outfiles = utils.mk_symlinks(src = input_files,
                                         dst = cluster_map_dir_jac)
        else:
            if transform_like is None:
                raise Exception("Argument transform_like must be specified "
                                "when using transform.")
            outfiles = utils.transform_images(infiles = input_files,
                                              outdir = cluster_map_dir_jac,
                                              like = transform_like,
                                              transform = transform,
                                              parallel = parallel,
                                              nproc = nproc)

        # Update file names
        for file in outfiles:
            file_split = os.path.basename(file).split('_')
            nk = file_split[3]
            k = file_split[1]
            file_new = 'cluster_map_nk_{}_k_{}.mnc'.format(nk, k)
            file_new = os.path.join(cluster_map_dir_jac, file_new)
            os.rename(file, file_new)

    return


def process_human_data(pipeline_dir = 'data/human/derivatives/v1/',
                       input_dir = 'data/human/registration/v1/jacobians_resampled/resolution_0.8/',
                       demographics = 'data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv',
                       mask = 'data/human/registration/v1/reference_files/mask_0.8mm.mnc',
                       datasets = ('POND', 'SickKids'),
                       parallel = True, nproc = None,
                       es_method = 'normative-growth', es_nbatches = 1,
                       es_df = 3, es_batch = ('Site', 'Scanner'),
                       es_ncontrols = 10, es_matrix_file = 'effect_sizes.csv',
                       cluster_resolution = 3.0,
                       cluster_nk_max = 10, cluster_metric = 'correlation',
                       cluster_K = 10, cluster_sigma = 0.5, cluster_t = 20,
                       cluster_file = 'clusters.csv',
                       cluster_affinity_file = 'affinity.csv',
                       cluster_map_method = 'mean'):
    """
    Docstring
    
    """

    # Convert parameter tuples to lists
    datasets = list(datasets)
    es_batch = list(es_batch)

    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the "
                             "number of processors to use in parallel.")

    # Effect size calculation method parameters
    if es_method == 'normative-growth':
        es_ncontrols = None
    elif es_method == 'propensity-matching':
        es_df = None
        es_batch = None
    else:
        raise ValueError

    # Image resolution
    vol = volumeFromFile(mask)
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Create pipeline directories ---------------------------------------------

    # Pipeline parameters
    params = dict(
        dataset = '-'.join(datasets),
        resolution = resolution,
        es_method = es_method,
        es_df = es_df,
        es_batch = (None if es_batch is None
                    else '-'.join(es_batch)),
        es_ncontrols = es_ncontrols,
        cluster_resolution = cluster_resolution,
        cluster_nk_max = cluster_nk_max,
        cluster_metric = cluster_metric,
        cluster_K = cluster_K,
        cluster_sigma = cluster_sigma,
        cluster_t = cluster_t,
        cluster_map_method = cluster_map_method
    )

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
    es_dir = os.path.join(pipeline_dir, 'effect_sizes', 'resolution_{}'.format(resolution), '')
    cluster_dir = os.path.join(pipeline_dir, 'clusters', 'resolution_{}'.format(cluster_resolution), '')
    cluster_map_dir = os.path.join(pipeline_dir, 'cluster_maps', 'resolution_{}'.format(resolution), '')

    # Check existence of input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: ".format(input_dir))

    # Filter for data sets ----------------------------------------------------

    # Import demographics
    df_demographics = pd.read_csv(demographics)

    # Filter individuals for data subset
    df_demographics = (df_demographics
                       .loc[df_demographics['Dataset'].isin(datasets)]
                       .copy())

    # Write out demographics subset to subset directory
    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)

    # Execute pipeline --------------------------------------------------------

    # Iterate over jacobians to compute effect sizes
    jacobians = ['absolute', 'relative']
    for j, jac in enumerate(jacobians):

        print("Processing {} Jacobian images...".format(jac))

        # Create symlinks to Jacobian images ----------------------------------
        print("Creating symlinks to Jacobian images...")
        input_files = glob(os.path.join(input_dir, jac, '') + '*.mnc')
        if len(input_files) == 0:
            raise OSError("No input files in directory: ".format(input_dir))
        input_files_in_dataset = [[f for f in input_files if g in f][0]
                                  for g in df_demographics['file'].to_list()]
        imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                     dst = os.path.join(imgdir, jac, ''))

        # Compute effect sizes ------------------------------------------------
        print("Computing effect size images...")
        es_kwargs = dict(imgdir = os.path.join(imgdir, jac, ''),
                         demographics = demographics,
                         mask = mask,
                         outdir = os.path.join(es_dir, jac, ''),
                         nproc = nproc,
                         method = es_method)

        if es_method == 'normative-growth':
            es_kwargs.update(
                dict(df = es_df,
                     batch = es_batch,
                     nbatches = es_nbatches)
            )
        else:
            es_kwargs.update(
                dict(ncontrols = es_ncontrols)
            )

        es_files = processing.calculate_human_effect_sizes(**es_kwargs)

        # Resample effect size images -----------------------------------------
        if resolution != cluster_resolution:
            print("Downsampling effect sizes to {}mm...".format(cluster_resolution))
            es_dir_downsampled = es_dir.replace(
                'resolution_{}'.format(resolution),
                'resolution_3.0')
            es_files_downsampled = utils.resample_images(
                infiles = es_files,
                outdir = os.path.join(es_dir_downsampled, jac, ''),
                isostep = cluster_resolution,
                parallel = parallel,
                nproc = nproc
            )
        else:
            es_dir_downsampled = es_dir
            es_files_downsampled = es_files

        print("Building effect size matrix...")
        df_es = processing.build_voxel_matrix(imgfiles = es_files_downsampled,
                                              mask = mask, file_col = True,
                                              sort = True, parallel = True,
                                              nproc = nproc)
        df_es['file'] = [os.path.basename(file) for file in df_es['file']]
        df_es.to_csv(os.path.join(es_dir_downsampled, jac, es_matrix_file),
                     index = False)

    # Cluster effect sizes ----------------------------------------------------
    print("Clustering absolute and relative effect size images...")
    cluster_kwargs = dict(
        infiles = [os.path.join(es_dir_downsampled, jac, es_matrix_file)
                   for jac in jacobians],
        rownames = 'file',
        nk_max = cluster_nk_max,
        metric = cluster_metric,
        K = cluster_K,
        sigma = cluster_sigma,
        t = cluster_t,
        cluster_file = os.path.join(cluster_dir, cluster_file),
        affinity_file = (os.path.join(cluster_dir, cluster_affinity_file)
                         if cluster_affinity_file is not None else None)
    )

    cluster_file = processing.cluster_human_data(**cluster_kwargs)

    # Create cluster maps -----------------------------------------------------
    for j, jac in enumerate(jacobians):
        print("Creating representative cluster maps for {} images..."
              .format(jac))
        cluster_map_kwargs = dict(
            clusters = cluster_file,
            imgdir = os.path.join(es_dir, jac, ''),
            outdir = os.path.join(cluster_map_dir, jac, ''),
            mask = mask,
            method = cluster_map_method,
            nproc = nproc
        )
        cluster_maps = processing.create_cluster_maps(**cluster_map_kwargs)

    return


def compute_cluster_similarity(human_pipeline_dir, mouse_pipeline_dir,
                               human_params_id, mouse_params_id,
                               pipeline_dir = 'data/cross_species/',
                               human_expr_dir = 'data/human/expression/',
                               mouse_expr_dir = 'data/mouse/expression/',
                               human_mask = 'data/human/registration/v2/reference_files/mask_0.8mm.mnc',
                               mouse_mask = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
                               human_microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_study_v2.csv',
                               gene_space = 'average-latent-space',
                               n_latent_spaces = 100, latent_space_id = 1,
                               metric = 'correlation', signed = True,
                               threshold = 'top_n', threshold_value = 0.2,
                               threshold_symmetric = True,
                               jacobians = ('absolute', 'relative'),
                               parallel = True, nproc = None):
    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified "
                            "when --parallel true")

    #Check jacobians type
    if type(jacobians) is tuple:
        jacobians = list(jacobians)
    elif type(jacobians) is str:
        jacobians = [jacobians]
    else:
        raise TypeError("`jacobians` must be a string or a tuple of strings.")

    # Ensure proper paths
    human_pipeline_dir = os.path.join(human_pipeline_dir, '')
    mouse_pipeline_dir = os.path.join(mouse_pipeline_dir, '')

    # Pipeline metadata
    human_metadata = os.path.join(human_pipeline_dir, 'metadata.csv')
    mouse_metadata = os.path.join(mouse_pipeline_dir, 'metadata.csv')

    # Fetch human pipeline parameters
    human_params = utils.fetch_params_metadata(human_metadata,
                                               id = human_params_id)
    human_params = human_params.to_dict(orient = 'list')
    human_params = {'_'.join(['human', key]): val[0] for key, val in
                    human_params.items()}

    # Fetch mouse pipeline parameters
    mouse_params = utils.fetch_params_metadata(mouse_metadata,
                                               id = mouse_params_id)
    mouse_params = mouse_params.to_dict(orient = 'list')
    mouse_params = {'_'.join(['mouse', key]): val[0] for key, val in
                    mouse_params.items()}

    # Combine parameter sets
    params = human_params.copy()
    params.update(mouse_params)
    params.update(dict(gene_space = gene_space,
                       n_latent_spaces = n_latent_spaces,
                       latent_space_id = latent_space_id,
                       metric = metric,
                       signed = signed,
                       threshold = threshold,
                       threshold_value = threshold_value,
                       threshold_symmetric = threshold_symmetric))

    # Create pipeline directory for parameter set
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir)
    pipeline_dir = os.path.join(pipeline_dir, 'similarity', '')
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)

    # Image resolutions
    human_resolution = params['human_resolution']
    mouse_resolution = params['mouse_resolution']

    # Cluster map directories
    human_cluster_dir = os.path.join(human_pipeline_dir, human_params_id, 'cluster_maps',
                                     'resolution_{}'.format(human_resolution),
                                     '')
    mouse_cluster_dir = os.path.join(mouse_pipeline_dir, mouse_params_id, 'cluster_maps',
                                     'resolution_{}'.format(mouse_resolution),
                                     '')

    # Combine mouse and human data into tuples
    expr = (human_expr_dir, mouse_expr_dir)
    masks = (human_mask, mouse_mask)

    # Iterate over Jacobians
    for j, jac in enumerate(jacobians):
        print(
            "Evaluating similarity of {} Jacobian cluster maps...".format(jac)
        )

        # Update input paths with jacobians
        human_cluster_dir_jac = os.path.join(human_cluster_dir, jac, '')
        mouse_cluster_dir_jac = os.path.join(mouse_cluster_dir, jac, '')

        # Get mouse and human cluster map files
        human_cluster_maps = os.listdir(human_cluster_dir_jac)
        mouse_cluster_maps = os.listdir(mouse_cluster_dir_jac)

        # Update cluster map files with directory paths
        human_cluster_maps = [os.path.join(human_cluster_dir_jac, file)
                              for file in human_cluster_maps]
        mouse_cluster_maps = [os.path.join(mouse_cluster_dir_jac, file)
                              for file in mouse_cluster_maps]

        # Expand mouse and human cluster map combinations
        cluster_pairs = list(product(human_cluster_maps, mouse_cluster_maps))

        # Compute pairwise similarity between cluster maps
        out = transcriptomic.transcriptomic_similarity(
            imgs = cluster_pairs,
            expr = expr,
            masks = masks,
            microarray_coords = human_microarray_coords,
            gene_space = gene_space,
            n_latent_spaces = n_latent_spaces,
            latent_space_id = latent_space_id,
            metric = metric,
            signed = signed,
            threshold = threshold,
            threshold_value = threshold_value,
            threshold_symmetric = threshold_symmetric,
            parallel = parallel,
            nproc = nproc
        )

        # Export similarity
        outfile = 'similarity_{}.csv'.format(jac)
        outfile = os.path.join(pipeline_dir, outfile)
        out.to_csv(outfile, index = False)

    return


def permute_cluster_similarity(pipeline_dir, params_id,
                               human_pipeline_dir = 'data/human/derivatives/v2/',
                               mouse_pipeline_dir = 'data/mouse/derivatives/v2/',
                               human_expr_dir = 'data/human/expression/',
                               mouse_expr_dir = 'data/mouse/expression/',
                               human_mask = 'data/human/registration/v2/reference_files/mask_0.8mm.mnc',
                               mouse_mask = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
                               human_microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_study_v2.csv',
                               npermutations = 100, permutations_start = 1,
                               jacobians = ('absolute', 'relative'),
                               keep_cluster_maps = False, nproc = 1):

    #Check jacobians type
    if type(jacobians) is tuple:
        jacobians = list(jacobians)
    elif type(jacobians) is str:
        jacobians = [jacobians]
    else:
        raise TypeError("`jacobians` must be a string or a tuple of strings.")

    #Get parameter set with specified ID
    params_metadata = os.path.join(pipeline_dir, 'metadata.csv')
    if not os.path.exists(params_metadata):
        raise OSError("File not found: {}".params_metadata)
    params = utils.fetch_params_metadata(params_metadata, id = str(params_id))

    #Pipeline directory
    pipeline_dir = os.path.join(pipeline_dir, params_id, 'permutations', '')
    if not os.path.exists(pipeline_dir):
        os.makedirs(pipeline_dir)

    # Output directory
    output_dir = os.path.join(pipeline_dir, 'similarity', '')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #Human and mouse parameter set IDs
    human_params_id = params['human_id'].values[0]
    mouse_params_id = params['mouse_id'].values[0]

    #Human and mouse pipeline resolutions
    human_resolution = params['human_resolution'].values[0]
    human_cluster_resolution = params['human_cluster_resolution'].values[0]
    mouse_resolution = params['mouse_resolution'].values[0]

    #Human directories
    human_pipeline_dir = os.path.join(human_pipeline_dir, human_params_id, '')
    human_es_dir = os.path.join(human_pipeline_dir, 'effect_sizes', 'resolution_{}'.format(human_resolution), '')
    human_cluster_dir = os.path.join(human_pipeline_dir, 'clusters', 'resolution_{}'.format(human_cluster_resolution), '')

    #Mouse directories
    mouse_pipeline_dir = os.path.join(mouse_pipeline_dir, mouse_params_id, '')
    mouse_cluster_map_dir = os.path.join(mouse_pipeline_dir, 'cluster_maps', 'resolution_{}'.format(mouse_resolution),)

    # Permute cluster labels
    human_cluster_file = os.path.join(human_cluster_dir, 'clusters.csv')
    human_cluster_perm_dir = os.path.join(pipeline_dir, 'clusters', '')
    human_perm_files = processing.permute_cluster_labels(
        cluster_file = human_cluster_file,
        outdir = human_cluster_perm_dir,
        npermutations = npermutations,
        start = permutations_start)

    #Human cluster centroid method
    cluster_map_method = params['human_cluster_map_method'].values[0]

    # Transcriptomic similarity args
    similarity_kwargs = params[['gene_space', 'n_latent_spaces',
                                'latent_space_id', 'metric', 'signed',
                                'threshold', 'threshold_value',
                                'threshold_symmetric']]
    similarity_kwargs = similarity_kwargs.to_dict(orient = 'list')
    similarity_kwargs = {key:val[0] for key, val in similarity_kwargs.items()}
    similarity_kwargs['n_latent_spaces'] = int(similarity_kwargs['n_latent_spaces'])
    similarity_kwargs['signed'] = bool(similarity_kwargs['signed'])
    similarity_kwargs['threshold_value'] = float(similarity_kwargs['threshold_value'])
    similarity_kwargs['threshold_symmetric'] = bool(similarity_kwargs['threshold_symmetric'])

    # Iterate over permutations
    p_range = range(permutations_start, permutations_start + npermutations)
    for p, pfile in zip(p_range, human_perm_files):
        print("Permutation {} of {}".format(p, p_range[len(p_range)-1]))

        # Iterate over jacobians
        for j in jacobians:

            print("\tProcessing {} Jacobians...".format(j))

            # Create permuted human cluster maps
            print("\t\tCreating cluster maps...")
            human_imgdir = os.path.join(pipeline_dir, 'cluster_maps',
                                        'permutation_{}'.format(p),
                                        'resolution_{}'.format(human_resolution),
                                        j, '')
            cluster_map_kwargs = dict(
                clusters = pfile,
                imgdir = os.path.join(human_es_dir, j, ''),
                outdir = human_imgdir,
                mask = human_mask,
                method = cluster_map_method,
                nproc = nproc
            )
            human_cluster_maps = processing.create_cluster_maps(
                **cluster_map_kwargs)

            # Get mouse cluster maps
            mouse_imgdir = os.path.join(mouse_cluster_map_dir, j, '')
            mouse_cluster_maps = os.listdir(mouse_imgdir)
            mouse_cluster_maps = [os.path.join(mouse_imgdir, img) for img in
                                  mouse_cluster_maps]

            # Identify cluster pairs to evaluate
            cluster_pairs_all = list(
                product(human_cluster_maps, mouse_cluster_maps)
            )
            cluster_pairs = []
            for pair in cluster_pairs_all:

                human = pair[0]
                human = os.path.basename(human).replace('.mnc', '').split('_')
                human_nk = int(human[-3])

                mouse = pair[1]
                mouse = os.path.basename(mouse).replace('.mnc', '').split('_')
                mouse_nk = int(mouse[-3])

                nk_diff = abs(human_nk - mouse_nk)
                if nk_diff < 2:
                    cluster_pairs.append(pair)

            # Evaluate cluster similarity
            print("\t\tComputing pairwise cluster similarity...")
            similarity_kwargs.update(dict(
                imgs = cluster_pairs,
                expr = (human_expr_dir, mouse_expr_dir),
                masks = (human_mask, mouse_mask),
                microarray_coords = human_microarray_coords,
                parallel = True if nproc > 1 else False,
                nproc = nproc
            ))
            sim = transcriptomic.transcriptomic_similarity(**similarity_kwargs)

            # Export permuted similarity
            outfile = 'similarity_permutation_{}_{}.csv'.format(p, j)
            outfile = os.path.join(output_dir, outfile)
            sim.to_csv(outfile, index = False)

        # Remove permuted cluster maps
        if not keep_cluster_maps:
            rmtree(os.path.join(pipeline_dir, 'cluster_maps',
                                'permutation_{}'.format(p)))
