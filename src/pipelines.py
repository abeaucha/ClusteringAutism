import os
import processing
import transcriptomic
import utils
import pandas as pd
from glob import glob
from itertools import product
from shutil import rmtree


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

    # Check existence of input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: {}".format(input_dir))

    # Path to clustering directory
    cluster_dir = os.path.join(pipeline_dir, 'clusters', '')

    # Make link to cluster file
    cluster_file = os.path.join(input_dir, 'Clusters.csv')
    utils.mk_symlinks(src = [cluster_file],
                      dst = cluster_dir)
    os.rename(os.path.join(cluster_dir, os.path.basename(cluster_file)),
              os.path.join(cluster_dir, 'clusters.csv'))

    # Path to cluster maps directory
    resolution_mm = float(resolution) / 1000
    cluster_map_dir = os.path.join(pipeline_dir, 'cluster_maps')
    cluster_map_dir = utils.mkdir_from_params(
        params = dict(cluster_map_method = method),
        outdir = cluster_map_dir
    )
    cluster_map_dir = os.path.join(cluster_map_dir,
                                   'resolution_{}'.format(resolution_mm),
                                   '')

    # Iterate over jacobians
    jacobians = ['absolute', 'relative']
    for jac in jacobians:

        print("Getting {} Jacobian cluster maps...".format(jac))

        # Cluster map directory
        cluster_map_dir_jac = os.path.join(cluster_map_dir, jac, '')

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


def process_human_data(pipeline_dir = 'data/human/derivatives/',
                       input_dir = 'data/human/registration/jacobians_resampled/',
                       resolution = 3.0,
                       demographics = 'data/human/registration/DBM_input_demo_passedqc.csv',
                       mask = 'data/human/registration/reference_files/mask_3.0mm.mnc',
                       datasets = ['POND', 'SickKids'], parallel = True,
                       nproc = None, es_method = 'normative-growth', es_df = 3,
                       es_combat = True, es_combat_batch = ['Site', 'Scanner'],
                       es_ncontrols = 10, es_matrix_file = 'effect_sizes.csv',
                       cluster_nk_max = 10, cluster_metric = 'correlation',
                       cluster_K = 10, cluster_sigma = 0.5, cluster_t = 20,
                       cluster_file = 'clusters.csv',
                       cluster_affinity_file = 'affinity.csv',
                       cluster_map_method = 'mean'):
    """
    Docstring
    
    """

    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified "
                            "when --parallel true")

    if es_method == 'normative-growth':
        if not es_combat:
            es_combat_batch = None
        es_ncontrols = None
    else:
        es_df = None
        es_combat = False
        es_combat_batch = None

    # Filter for data sets ----------------------------------------------------

    # Create output directory for specified datasets
    pipeline_dir = utils.mkdir_from_list(inlist = datasets,
                                         basedir = pipeline_dir)

    # Import demographics
    df_demographics = pd.read_csv(demographics)

    # Filter individuals for data subset
    df_demographics = (df_demographics
                       .loc[df_demographics['Dataset'].isin(datasets)]
                       .copy())

    # Write out demographics subset to subset directory
    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)

    # Create pipeline directories ---------------------------------------------

    # Paths to input directory
    input_dir = os.path.join(input_dir,
                             'resolution_{}'.format(resolution),
                             '')
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: ".format(input_dir))

    # Paths to pipeline image directory
    imgdir = os.path.join(pipeline_dir,
                          'jacobians',
                          'resolution_{}'.format(resolution),
                          '')

    # Pipeline parameters
    params = dict(
        effect_sizes = dict(es_method = es_method,
                            es_df = es_df,
                            es_combat = es_combat,
                            es_combat_batch = (None if es_combat_batch is None
                                               else '-'.join(es_combat_batch)),
                            es_ncontrols = es_ncontrols),
        clusters = dict(cluster_nk_max = cluster_nk_max,
                        cluster_metric = cluster_metric,
                        cluster_K = cluster_K,
                        cluster_sigma = cluster_sigma,
                        cluster_t = cluster_t),
        cluster_maps = dict(cluster_map_method = cluster_map_method)
    )

    # Create directories for pipeline stages based on parameters
    stage_params = {}
    stage_dirs = {}
    for key, val in params.items():
        if key == list(params.keys())[0]:
            params_id = utils.random_id(3)
        else:
            params_id = '-'.join([params_id, utils.random_id(3)])
        stage_params.update(val)
        stage_dir = os.path.join(pipeline_dir, key, '')
        metadata = os.path.join(stage_dir, 'metadata.csv')
        stage_dir = utils.mkdir_from_params(params = stage_params,
                                            outdir = stage_dir,
                                            params_id = params_id)
        stage_dir = os.path.join(stage_dir,
                                 'resolution_{}'.format(resolution),
                                 '')
        stage_dirs[key] = stage_dir
        params_id = utils.get_params_id(params = stage_params,
                                        metadata = metadata)

    # Extract directories for pipeline stages
    es_dir = stage_dirs['effect_sizes']
    cluster_dir = stage_dirs['clusters']
    cluster_map_dir = stage_dirs['cluster_maps']

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
                         parallel = parallel,
                         nproc = nproc,
                         method = es_method)

        if es_method == 'normative-growth':
            es_kwargs.update(
                dict(df = es_df,
                     combat = es_combat,
                     combat_batch = es_combat_batch)
            )
        else:
            es_kwargs.update(
                dict(ncontrols = es_ncontrols)
            )

        es_files = processing.calculate_human_effect_sizes(**es_kwargs)

        # Build effect size matrix --------------------------------------------
        print("Building effect size voxel matrix...")
        df_es = processing.build_voxel_matrix(infiles = es_files,
                                              mask = mask,
                                              sort = True,
                                              file_col = True,
                                              parallel = parallel,
                                              nproc = nproc)
        df_es['file'] = [os.path.basename(file) for file in df_es['file']]
        df_es.to_csv(os.path.join(es_dir, jac, es_matrix_file), index = False)

    # Cluster effect sizes ----------------------------------------------------
    print("Clustering absolute and relative effect size images...")
    cluster_kwargs = dict(
        infiles = [os.path.join(es_dir, jac, es_matrix_file)
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
            method = cluster_map_method
        )
        cluster_maps = processing.create_cluster_maps(**cluster_map_kwargs)

    print("Done.")

    return


def compute_cluster_similarity(mouse_cluster_dir, human_cluster_dir,
                               pipeline_dir = 'data/similarity/',
                               mouse_expr_dir = 'data/mouse/expression/',
                               human_expr_dir = 'data/human/expression/',
                               mouse_mask = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
                               human_mask = 'data/human/registration/reference_files/mask_3.0mm.mnc',
                               mouse_resolution = 0.2, human_resolution = 3.0,
                               human_microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_studyspace.csv',
                               gene_space = 'average-latent-space',
                               n_latent_spaces = 100, latent_space_id = 1,
                               metric = 'correlation', signed = True,
                               threshold = 'top_n', threshold_value = 0.2,
                               threshold_symmetric = True,
                               parallel = True, nproc = None):
    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified "
                            "when --parallel true")

    # Ensure proper paths
    human_cluster_dir = os.path.join(human_cluster_dir, '')
    mouse_cluster_dir = os.path.join(mouse_cluster_dir, '')

    # Get parameter IDs for mouse and human pipelines
    human_input_params_id = human_cluster_dir.split('/')[-2]
    mouse_input_params_id = mouse_cluster_dir.split('/')[-2]

    # Update paths with resolution information
    human_cluster_dir = os.path.join(human_cluster_dir,
                                     'resolution_{}'.format(human_resolution),
                                     '')
    mouse_cluster_dir = os.path.join(mouse_cluster_dir,
                                     'resolution_{}'.format(mouse_resolution),
                                     '')

    # Create output directory from parameters
    params = dict(
        human_input_id=human_input_params_id,
        human_resolution=human_resolution,
        mouse_input_id=mouse_input_params_id,
        mouse_resolution=mouse_resolution,
        gene_space=gene_space,
        n_latent_spaces=n_latent_spaces,
        metric=metric,
        signed=signed,
        threshold=threshold,
        threshold_value=threshold_value,
        threshold_symmetric=threshold_symmetric
    )
    pipeline_dir = utils.mkdir_from_params(params = params,
                                           outdir = pipeline_dir)

    # Combine mouse and human data into tuples
    expr = (mouse_expr_dir, human_expr_dir)
    masks = (mouse_mask, human_mask)

    # Iterate over Jacobians
    jacobians = ['absolute', 'relative']
    for j, jac in enumerate(jacobians):
        print(
            "Evaluating similarity of {} Jacobian cluster maps...".format(jac)
        )

        # Update input paths with jacobians
        human_cluster_dir_jac = os.path.join(human_cluster_dir, jac, '')
        mouse_cluster_dir_jac = os.path.join(mouse_cluster_dir, jac, '')

        # Get mouse and human cluster map files
        mouse_cluster_maps = os.listdir(mouse_cluster_dir_jac)
        human_cluster_maps = os.listdir(human_cluster_dir_jac)

        # Update cluster map files with directory paths
        mouse_cluster_maps = [os.path.join(mouse_cluster_dir_jac, file)
                              for file in mouse_cluster_maps]
        human_cluster_maps = [os.path.join(human_cluster_dir_jac, file)
                              for file in human_cluster_maps]

        # Expand mouse and human cluster map combinations
        cluster_pairs = list(product(mouse_cluster_maps, human_cluster_maps))

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


def permute_cluster_similarity(human_pipeline_dir = 'data/human/derivatives/',
                               mouse_pipeline_dir = 'data/mouse/derivatives/',
                               human_expr_dir = 'data/human/expression/',
                               mouse_expr_dir = 'data/mouse/expression/',
                               human_mask = 'data/human/registration/reference_files/mask_3.0mm.mnc',
                               mouse_mask = 'data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc',
                               human_microarray_coords = 'data/human/expression/AHBA_microarray_coordinates_studyspace.csv',
                               output_dir = 'data/cross_species/',
                               human_cluster_params_id = None,
                               human_dataset = 'POND_SickKids',
                               mouse_dataset = 'Models_135',
                               npermutations = 100,
                               cluster_map_method = 'mean',
                               sim_gene_space = 'average-latent-space',
                               sim_n_latent_spaces = 5,
                               sim_latent_space_id = 1,
                               sim_metric = 'correlation',
                               sim_signed = True,
                               sim_threshold = 'top_n',
                               sim_threshold_value = 0.2,
                               sim_threshold_symmetric = True,
                               parallel = True,
                               nproc = None):

    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified "
                            "when --parallel true")
    
    # Paths to human directories
    if human_cluster_params_id is None:
        raise Exception("Specify parameter set ID for human clusters.")
    human_pipeline_dir = os.path.join(human_pipeline_dir, human_dataset, '')
    human_es_params_id = human_cluster_params_id.split('-')[0]
    human_es_dir = os.path.join(human_pipeline_dir, 'effect_sizes', human_es_params_id, 'resolution_3.0', '')
    human_cluster_dir = os.path.join(human_pipeline_dir, 'clusters', human_cluster_params_id, 'resolution_3.0', '')

    # Paths to mouse directories
    mouse_pipeline_dir = os.path.join(mouse_pipeline_dir, mouse_dataset, '')
    mouse_cluster_map_dir = os.path.join(mouse_pipeline_dir, 'cluster_maps', '')
    mouse_metadata = os.path.join(mouse_cluster_map_dir, 'metadata.csv')
    df_mouse_metadata = utils.fetch_metadata(metadata = mouse_metadata,
                                             cluster_map_method = cluster_map_method)
    mouse_cluster_map_params_id = df_mouse_metadata['id'].values[0]
    mouse_cluster_map_dir = os.path.join(mouse_cluster_map_dir, mouse_cluster_map_params_id, 'resolution_0.2', '')

    # Output directory
    datasets = '-'.join([mouse_dataset, human_dataset])
    output_dir = os.path.join(output_dir, datasets, 'permutations', 'similarity', '')
    params = dict(
        human_input_id = human_cluster_params_id,
        cluster_map_method = cluster_map_method,
        gene_space = sim_gene_space,
        n_latent_spaces = sim_n_latent_spaces,
        metric = sim_metric,
        signed = sim_signed,
        threshold = sim_threshold,
        threshold_value = sim_threshold_value,
        threshold_symmetric = sim_threshold_symmetric
    )
    output_dir = utils.mkdir_from_params(params = params, outdir = output_dir)

    # Permute cluster labels
    human_perm_dir = os.path.join(human_pipeline_dir, 'permutations', '')
    human_cluster_file = os.path.join(human_cluster_dir, 'clusters.csv')
    human_cluster_perm_dir = os.path.join(human_perm_dir, 'clusters',
                                          human_cluster_params_id,
                                          'resolution_3.0',
                                          '')
    human_perm_files = processing.permute_cluster_labels(cluster_file = human_cluster_file,
                                                         outdir = human_cluster_perm_dir,
                                                         npermutations = npermutations)

    # Iterate over permutations
    for p, pfile in enumerate(human_perm_files):

        print("Permutation {} of {}".format(p+1, len(human_perm_files)))

        perm = os.path.basename(pfile)
        perm = os.path.splitext(perm)[0]
        perm = perm.replace('clusters_', '')

        #Iterate over jacobians
        jacobians = ['absolute', 'relative']
        for j in jacobians:

            print("\tProcessing {} Jacobians...".format(j))

            print("\t\tCreating cluster maps...")
            human_imgdir = os.path.join(human_es_dir, j, '')
            human_cluster_map_dir = os.path.join(human_perm_dir, 'cluster_maps', 'resolution_3.0', j, '')
            cluster_map_kwargs = dict(
                clusters = pfile,
                imgdir = human_imgdir,
                outdir = human_cluster_map_dir,
                mask = human_mask,
                method = cluster_map_method
            )
            human_cluster_maps = processing.create_cluster_maps(**cluster_map_kwargs)

            #TODO: Potentially need to resample cluster maps to 1.0mm here

            mouse_imgdir = os.path.join(mouse_cluster_map_dir, j, '')
            mouse_cluster_maps = os.listdir(mouse_imgdir)
            mouse_cluster_maps = [os.path.join(mouse_imgdir, img) for img in
                                  mouse_cluster_maps]

            cluster_pairs_all = list(product(mouse_cluster_maps, human_cluster_maps))
            cluster_pairs = []
            for pair in cluster_pairs_all:

                m = pair[0]
                m = os.path.basename(m).replace('.mnc', '').split('_')
                m_nk = int(m[-3])

                h = pair[1]
                h = os.path.basename(h).replace('.mnc', '').split('_')
                h_nk = int(h[-3])

                if m_nk == h_nk:
                    m_k = int(m[-1])
                    h_k = int(h[-1])
                    if abs(m_k - h_k) < 2:
                        cluster_pairs.append(pair)

            expr = (mouse_expr_dir, human_expr_dir)
            masks =(mouse_mask, human_mask)

            print("\t\tComputing pairwise cluster similarity...")

            sim = transcriptomic.transcriptomic_similarity(
                imgs = cluster_pairs,
                expr = expr,
                masks = masks,
                microarray_coords = human_microarray_coords,
                gene_space = sim_gene_space,
                n_latent_spaces = sim_n_latent_spaces,
                latent_space_id = sim_latent_space_id,
                metric = sim_metric,
                signed = sim_signed,
                threshold = sim_threshold,
                threshold_value = sim_threshold_value,
                threshold_symmetric = sim_threshold_symmetric,
                parallel = parallel,
                nproc = nproc
            )

            outfile = 'similarity_{}_{}.csv'.format(perm, j)
            outfile = os.path.join(output_dir, outfile)
            sim.to_csv(outfile, index = False)

    rmtree(os.path.join(human_perm_dir, 'cluster_maps'))