# ----------------------------------------------------------------------------
# process_human_data.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2023

"""
Pipeline to process human neuroimaging data.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import pandas as pd
from src import utils, processing
from glob import glob
from re import sub


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
        default = 'data/human/derivatives/',
        help = ("Path to pipeline directory.")
    )
    
    parser.add_argument(
        '--input-dir',
        type = str,
        default = 'data/human/registration/jacobians_resampled/',
        help = ("Path to Jacobians directory.")
    )
    
    parser.add_argument(
        '--resolution',
        type = float,
        default = 3.0,
        help = ("Resolution (mm) of Jacobian images to use.")
    )
    
    parser.add_argument(
        '--demographics',
        type = str,
        default = 'data/human/registration/DBM_input_demo_passedqc.csv',
        help = ("Path to CSV file containing demographics information.")
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        default = 'data/human/registration/reference_files/mask_3.0mm.mnc',
        help = ("Path to image mask")
    )
    
    parser.add_argument(
        '--datasets',
        nargs = '*',
        type = str,
        default = ['POND', 'SickKids'],
        help = ("List of strings indicating datasets to include in processing.")
    )
    
    parser.add_argument(
        '--parallel',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to run in parallel.")
    )
    
    parser.add_argument(
        '--nproc',
        type = int,
        help = ("Number of processors to use in parallel.")
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Verbosity.'
    )
    
    # Effect size arguments ---------------------------------------------------
    parser.add_argument(
        '--es-method',
        type = str,
        default = 'normative-growth',
        choices = ['normative-growth', 'propensity-matching'],
        help = ("Method to use to compute effect sizes.")
    )
        
    parser.add_argument(
        '--es-df',
        type = int,
        default = 3,
        help = ("Degrees of freedom hyperparameter in normative growth modelling.")
    )
        
    parser.add_argument(
        '--es-combat',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to run ComBat normalization prior to normative growth modelling.")
    )
        
    parser.add_argument(
        '--es-combat-batch',
        nargs = '*',
        type = str,
        default = ['Site', 'Scanner'],
        help = ("Batch variables to use in ComBat normalization.")
    )
        
    parser.add_argument(
        '--es-ncontrols',
        type = int,
        default = 10,
        help = ("Number of controls to use for propensity-matching.")
    )
        
    parser.add_argument(
        '--es-matrix-file',
        type = str,
        default = 'effect_sizes.csv',
        help = ("Name of CSV file in which to store voxel-wise effect size matrices.")
    )
        
        
    # Clustering arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-nk-max',
        type = int,
        default = 10,
        help = ("Maximum number of clusters to identify in cluster solutions.")
    )
        
    parser.add_argument(
        '--cluster-metric',
        type = str,
        default = 'correlation',
        choices = ['correlation'],
        help = ("Distance metric to use in similarity network fusion.")
    )
        
    parser.add_argument(
        '--cluster-K',
        type = int,
        default = 10,
        help = ("Number of nearest-neighbours to use in similarity network fusion.")
    )
        
    parser.add_argument(
        '--cluster-sigma',
        type = float,
        default = 0.5,
        help = ("Variance for the local model in similarity network fusion.")
    )
        
    parser.add_argument(
        '--cluster-t',
        type = int,
        default = 20,
        help = ("Number of iterations for the diffusion process in similarity network fusion.")
    )
        
    parser.add_argument(
        '--cluster-file',
        type = str,
        default = 'clusters.csv',
        help = ("Name of the CSV file in which to store cluster assignments.")
    )
        
    parser.add_argument(
        '--cluster-affinity-file',
        type = str,
        help = ("Name of the file in which to store similarity network fusion affinity matrix.")
    )
        
    # Cluster maps arguments ---------------------------------------------------
    parser.add_argument(
        '--cluster-map-method',
        type = str,
        default = 'mean',
        choices = ['mean', 'median'],
        help = ("Method to use to aggregate images within each cluster.")
    )
    
    args = vars(parser.parse_args())
    
    return args


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    
    #General parameters
    pipeline_dir = args['pipeline_dir']
    input_dir = args['input_dir']
    resolution = args['resolution']
    demographics = args['demographics']
    datasets = args['datasets']
    mask = args['mask']
    parallel = True if args['parallel'] == 'true' else False
    nproc = args['nproc']
    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified when --parallel true")
    
    #Effect size parameters..
    es_method = args['es_method']
    if es_method == 'normative-growth':
        es_df = args['es_df']
        es_combat = True if args['es_combat'] == 'true' else False
        if es_combat:
            es_combat_batch = args['es_combat_batch']
        else:
            es_combat_batch = None
        es_ncontrols = None
    else:
        es_df = None
        es_combat = False
        es_combat_batch = None
        es_ncontrols = args['es_ncontrols']
    es_matrix_file = args['es_matrix_file']

    #Clustering parameters
    cluster_nk_max = args['cluster_nk_max']
    cluster_metric = args['cluster_metric']
    cluster_K = args['cluster_K']
    cluster_sigma = args['cluster_sigma']
    cluster_t = args['cluster_t']
    cluster_file = args['cluster_file']
    cluster_affinity_file = args['cluster_affinity_file']
        
    #Cluster map parameters
    cluster_map_method = args['cluster_map_method']
    
    
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

    #Create a new column for patient image files
    df_demographics['file'] = (df_demographics['Extract_ID']
                               .str.replace('.mnc', 
                                            '_fwhm_4vox.mnc', 
                                            regex = True))

    #Write out demographics subset to subset directory
    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)
    
    
    # Create pipeline directories ---------------------------------------------
    
    #Paths to input directory
    input_dir = os.path.join(input_dir, 'resolution_{}'.format(resolution), '')
    
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
    
    
    return
    
if __name__=='__main__':
    main()
