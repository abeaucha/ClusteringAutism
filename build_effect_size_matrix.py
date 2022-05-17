# ----------------------------------------------------------------------------
# build_effect_size_matrix.py 
# Author: Antoine Beauchamp
# Created: May 16th, 2022


# Packages -------------------------------------------------------------------

import argparse
import os
import numpy                as np
import pandas               as pd
import multiprocessing      as mp
from pyminc.volumes.factory import volumeFromFile
from glob                   import glob
from tqdm                   import tqdm
from functools              import partial

# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--datadir',
        type = str,
        help = ("Directory containing effect size MINC images.")
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        help = ("Mask file."))
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = ("Directory in which to save effect size CSV file.")
    )
        
    parser.add_argument(
        '--outfile',
        type = str,
        default = 'effect_sizes.csv',
        help = ("Name of CSV file in which to export effect sizes.")
    )
    
    parser.add_argument(
        '--parallel',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = "Option to run in parallel."
    )
    
    parser.add_argument(
        '--nproc',
        type = int,
        default = mp.cpu_count(),
        help = ("Number of processors to use in parallel. "
                "Ignored if --parallel set to false.")
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Verbosity.'
    )
    
    args = vars(parser.parse_args())
    
    return args


# Functions ------------------------------------------------------------------



# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    datadir = args['datadir']
    outdir = args['outdir']
    outfile = args['outfile']
    mask = args['mask']
    parallel = True if args['parallel'] == 'true' else False
#     verbose = True if args['verbose'] == 'true' else False

    

    
if __name__=='__main__':
    main()