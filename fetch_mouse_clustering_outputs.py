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

import sys
import argparse
from pipelines import fetch_mouse_clustering_outputs
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
        help = ("Path to input directory.")
    )
    
    parser.add_argument(
        '--output-dir',
        type = str,
        default = 'data/mouse/derivatives/Models_135/',
        help = ("Path to output directory.")
    )
    
    parser.add_argument(
        '--resolution',
        type = int,
        default = 200,
        help = ("Resolution (um) of images to link.")
    )
    
    parser.add_argument(
        '--method',
        nargs = '*',
        type = str,
        default = ['mean'],
        help = ("Method used to aggregate cluster maps.")
    )
    
    parser.add_argument(
        '--transform',
        type = str,
        default = 'data/mouse/registration/MICe_scanbase.xfm',
        help = ("Path to transform file (.xfm).")
    )
    
    parser.add_argument(
        '--like',
        type = str,
        default = 'data/mouse/atlas/average_template_200um.mnc',
        help = ("Path to transform likefile (.mnc).")
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
    
    args = vars(parser.parse_args())
    
    return args


# Pipeline --------------------------------------------------------------------
if __name__ == '__main__':
    
    args = parse_args()
    args['parallel'] = bool(args['parallel'])
    if args['transform'] == 'none':
        args['transform'] = None
        args['like'] = None
    params = {key:val for key, val in args.items() 
              if type(val) is list}
    param_keys = [key for key in params.keys()]
    for param_vals in product(*params.values()):
        param_set = dict(list(zip(param_keys, param_vals)))
        param_msg = [': '.join([str(key), str(val)]) for key, val in param_set.items()]
        param_msg = reduce(lambda x, y: x+'\n\t'+y, param_msg)
        param_msg = "Running parameter set:\n\t{}\n".format(param_msg)
        print(param_msg)
        args.update(param_set)
        fetch_mouse_clustering_outputs(**args)
    
