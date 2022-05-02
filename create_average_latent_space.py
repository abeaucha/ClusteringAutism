# ----------------------------------------------------------------------------
# create_average_latent_space.py 
# Author: Antoine Beauchamp
# Created: May 2nd, 2022

"""
Create an average gene expression latent space.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import numpy as np
import pandas as pd
from glob import glob
from datatable import fread
from re import sub

# Functions ------------------------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/MLP_outcomes/',
        help = ("Directory containing latent space data sets.")
    )
    
    parser.add_argument(
        '--fileglob',
        type = str,
        help = ("Glob to get the latent space CSV files.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        help = ("Output directory.")
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


def average_latent_spaces(files):
    
    """ 
    """
    
    latent_space = fread(files[0], header = True).to_numpy()
    for file in files[1:]:
        
        latent_space_tmp = fread(file, header = True).to_numpy()
        
        latent_space = np.mean(np.array([latent_space, 
                                         latent_space_tmp]), 
                               axis = 0)
        
    return latent_space
    

# Main -----------------------------------------------------------------------

def main():
    
    args = parse_args()
    datadir = args['datadir']
    outdir = args['outdir']
    file_glob = args['fileglob']
    verbose = True if args['verbose'] == 'true' else False
    
    if outdir is None:
        outdir = datadir
    
    datadir = os.path.join(datadir, '')
    
    files = glob(datadir+file_glob)
    
    if verbose:
        print("Computing average latent space...")
    
    latent_space = average_latent_spaces(files = files)

    outfile = sub(r'[*]', '_average', file_glob)
    outfile = os.path.join(outdir, outfile)
    
    if verbose:
        print("Writing to file {} ...".format(outfile))
    
    pd.DataFrame(latent_space).to_csv(outfile)
    
    return


if __name__=='__main__':
    main()
