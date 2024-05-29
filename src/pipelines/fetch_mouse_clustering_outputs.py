#!.venv/bin/python3
# ----------------------------------------------------------------------------
# get_mouse_clustering_data.py 
# Author: Antoine Beauchamp
# Created: March 7th, 2023

"""
Pipeline to create links to mouse clustering data.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
from itertools import product
from functools import reduce


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    # General arguments ---------------------------------------------------
    parser.add_argument(
        '--input-dir',
        type = str,
        default = '/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters_Paper/',
        help = "Path to input directory."
    )
    
    parser.add_argument(
        '--pipeline-dir',
        type = str,
        default = 'data/mouse/derivatives/v2/',
        help = "Path to output directory."
    )
    
    parser.add_argument(
        '--resolution',
        type = int,
        default = 200,
        help = "Resolution (um) of images to link."
    )
    
    parser.add_argument(
        '--method',
        nargs = '*',
        type = str,
        default = ['mean'],
        help = "Method used to aggregate cluster maps."
    )
    
    parser.add_argument(
        '--transform',
        type = str,
        default = 'data/mouse/registration/MICe_scanbase.xfm',
        help = "Path to transform file (.xfm)."
    )
    
    parser.add_argument(
        '--transform-like',
        type = str,
        default = 'data/mouse/atlas/average_template_200um.mnc',
        help = "Path to transform likefile (.mnc)."
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
    
    args = vars(parser.parse_args())
    
    return args


# Modules ----------------------

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
    cluster_map_dir = os.path.join(pipeline_dir, 'cluster_maps',
                                   'resolution_{}'.format(resolution_mm), '')

    # Check existence of input directory
    if not os.path.exists(input_dir):
        raise OSError("Input directory not found: {}".format(input_dir))

    # Get cluster files -------------------------------------------------------

    # Make link to cluster file
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    # Clusters and affinity matrix
    cluster_file = os.path.join(input_dir, 'Clusters.csv')
    affinity_file = '/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/SNFMatrix_Paper/WMatrix.RData'

    # Copy files
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


# Pipeline --------------------------------------------------------------------

if __name__ == '__main__':
    
    args = parse_args()
    args['parallel'] = True if args['parallel'] == 'true' else False
    if args['transform'] == 'none':
        args['transform'] = None
        args['transform_like'] = None
    params = {key:val for key, val in args.items() 
              if type(val) is list}
    param_keys = [key for key in params.keys()]
    for param_vals in product(*params.values()):
        param_set = dict(list(zip(param_keys, param_vals)))
        param_msg = [': '.join([str(key), str(val)])
                     for key, val in param_set.items()]
        param_msg = reduce(lambda x, y: x+'\n\t'+y, param_msg)
        param_msg = "Running parameter set:\n\t{}\n".format(param_msg)
        print(param_msg)
        args.update(param_set)
        fetch_mouse_clustering_outputs(**args)
    
