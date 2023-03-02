import os
import pandas as pd
from src import utils, processing
from glob import glob
from re import sub

def process_human_data(pipeline_dir = 'data/human/derivatives/', 
                       input_dir = 'data/human/registration/jacobians_resampled/', 
                       resolution = 3.0, 
                       demographics = 'data/human/registration/DBM_input_demo_passedqc.csv', 
                       mask = 'data/human/registration/reference_files/mask_3.0mm.mnc', 
                       datasets = ['POND', 'SickKids'], parallel = True, 
                       nproc = None, verbose = True, 
                       es_method = 'normative-growth', es_df = 3, 
                       es_combat = True, es_combat_batch = ['Site', 'Scanner'], 
                       es_ncontrols = 10, es_matrix_file = 'effect_sizes.csv', 
                      cluster_nk_max = 10, cluster_metric = 'correlation',
                      cluster_K = 10, cluster_sigma = 0.5, cluster_t = 20,
                      cluster_file = 'clusters.csv', cluster_affinity_file = None,
                      cluster_map_method = 'mean'):
    
    """
    Docstring
    
    """
    
    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified when --parallel true")
            
    if es_method == 'normative-growth':
        if not es_combat:
            es_combat_batch = None
        es_ncontrols = None
    else:
        es_df = None
        es_combat = False
        es_combat_batch = None

        
    # Filter for data sets ---------------------------------------------------------
        
    #Create output directory for specified datasets
    pipeline_dir = utils.mkdir_from_list(inlist = datasets,
                                         basedir = pipeline_dir)

    #Import demographics
    df_demographics = pd.read_csv(demographics)

    #Filter individuals for data subset
    df_demographics = (df_demographics
                       .loc[df_demographics['Dataset'].isin(datasets)]
                       .copy())

    #Write out demographics subset to subset directory
    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)

    
    # Create pipeline directories ---------------------------------------------
    
    #Paths to input directory
    input_dir = os.path.join(input_dir, 'resolution_{}'.format(resolution), '')
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: ".format(input_dir))

    #Paths to pipeline image directory
    imgdir = os.path.join(pipeline_dir, 'jacobians', 'resolution_{}'.format(resolution),'')

    #Pipeline parameters
    params = {'effect_sizes': {'es_method':es_method,
                               'es_df':es_df,
                               'es_combat':es_combat,
                               'es_combat_batch':(None if es_combat_batch is None 
                                                  else '-'.join(es_combat_batch)),
                               'es_ncontrols':es_ncontrols},
              'clusters': {'cluster_nk_max':cluster_nk_max,
                           'cluster_metric':cluster_metric,
                           'cluster_K':cluster_K,
                           'cluster_sigma':cluster_sigma,
                           'cluster_t':cluster_t},
              'cluster_maps':{'cluster_map_method':cluster_map_method}}

    #Create directories for pipeline stages based on parameters
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
        stage_dir = os.path.join(stage_dir, 'resolution_{}'.format(resolution), '')
        stage_dirs[key] = stage_dir
        params_id = utils.get_params_id(params = stage_params, metadata = metadata)

    #Extract directories for pipeline stages
    es_dir = stage_dirs['effect_sizes']
    cluster_dir = stage_dirs['clusters']
    cluster_map_dir = stage_dirs['cluster_maps']
    
    # Execute pipeline -------------------------------------------------------------
    
    #Iterate over jacobians
    jacobians = ['absolute', 'relative']
    for j, jac in enumerate(jacobians):

        print("Processing {} Jacobian images...".format(jac))

        # Create symlinks to Jacobian images -----------------------------------
        print("Creating symlinks to Jacobian images...")
        input_files = glob(os.path.join(input_dir, jac, '')+'*.mnc')
        if len(input_files) == 0:
            raise OSError("No input files in directory: ".format(input_dir))
        input_files_in_dataset = [[f for f in input_files if g in f][0] 
                                  for g in df_demographics['file'].to_list()]
        imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                     dst = os.path.join(imgdir, jac, ''))
        
        
        # Compute effect sizes ----------------------------------------------------------
        print("Computing effect size images...")
        es_kwargs = {'imgdir':os.path.join(imgdir, jac, ''),
                     'demographics':demographics,
                     'mask':mask,
                     'outdir':os.path.join(es_dir, jac, ''),
                     'parallel':parallel,
                     'nproc':nproc,
                     'method':es_method}
        if es_method == 'normative-growth':
            es_kwargs.update({'df':es_df, 
                              'combat':es_combat, 
                              'combat_batch':es_combat_batch})
        else:
            es_kwargs.update({'ncontrols':es_ncontrols})
        es_files = processing.calculate_human_effect_sizes(**es_kwargs)
        

        # Build effect size matrix ---------------------------------------------------
        print("Building effect size voxel matrix...")
        df_es = processing.build_voxel_matrix(infiles = es_files,
                                              mask = mask,
                                              sort = True,
                                              file_col = True,
                                              parallel = parallel,
                                              nproc = nproc)
        df_es['file'] = [os.path.basename(file) for file in df_es['file']]
        df_es.to_csv(os.path.join(es_dir, jac, es_matrix_file), index = False)
    
    
    # Cluster effect sizes ----------------------------------------------------------------
    print("Clustering absolute and relative effect size images...")
    cluster_kwargs = {'infiles':[os.path.join(es_dir, jac, es_matrix_file) 
                                 for jac in jacobians],
                      'rownames':'file',
                      'nk_max':cluster_nk_max,
                      'metric':cluster_metric,
                      'K':cluster_K,
                      'sigma':cluster_sigma,
                      't':cluster_t,
                      'cluster_file':os.path.join(cluster_dir, cluster_file),
                      'affinity_file':(os.path.join(cluster_dir, cluster_affinity_file) 
                                       if cluster_affinity_file is not None else None)}
    cluster_file = processing.cluster_human_data(**cluster_kwargs)

    
    # Create cluster maps -------------------------------------------------------------------------
    for j, jac in enumerate(jacobians):
        print("Creating representative cluster maps for {} images...".format(jac))
        cluster_map_kwargs = {'clusters':cluster_file,
                             'imgdir':os.path.join(es_dir, jac, ''),
                             'outdir':os.path.join(cluster_map_dir, jac, ''),
                             'mask':mask,
                             'method':cluster_map_method}
        cluster_maps = processing.create_cluster_maps(**cluster_map_kwargs)
    
    print("Done.")
    
    return