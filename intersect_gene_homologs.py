# ----------------------------------------------------------------------------
# intersect_gene_homologs.py 
# Author: Antoine Beauchamp
# Created: April 4th, 2022

"""
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import numpy  as np
import pandas as pd
from re import sub

# Functions ------------------------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--mouse',
        type = str,
        help = ("Path to csv file containing mouse voxel-wise expression "
                "matrix.")
    )
    
    parser.add_argument(
        '--human',
        type = str,
        help = ("Path to csv file containing human sample-wise expression "
                "matrix.")
    )
    
    parser.add_argument(
        '--homologs',
        type = str,
        help = ("Path to csv file containing mouse-human homologous genes")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = ("Output directory")
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


def intersect_gene_homologs(mouse, human, homologs):
    
    """ 
    """
    
    ind_mouse = homologs['Mouse'].isin(mouse.index)
    ind_human = homologs['Human'].isin(human.index)
    homologs = homologs.loc[ind_mouse & ind_human]
    
    mouse = mouse.loc[mouse.index.isin(homologs['Mouse'])]
    human = human.loc[human.index.isin(homologs['Human'])]
    
    mouse = mouse.sort_index()
    homologs = homologs.sort_values(by = 'Mouse')
    mouse.index = homologs['Human']
    mouse.index.name = 'Gene'
    
    mouse = mouse.sort_index()
    human = human.sort_index()
    
    return mouse, human


# Main -----------------------------------------------------------------------

def main():
    
    #Load command line arguments
    args = parse_args()
    mousefile = args['mouse']
    humanfile = args['human']
    homologs = args['homologs']
    outdir = args['outdir']
    verbose = True if args['verbose'] == 'true' else False
    
    if verbose:
        print("Importing files...")
    
    #Import data frames
    mouse = pd.read_csv(mousefile, index_col = 'Gene')
    human = pd.read_csv(humanfile, index_col = 'Gene')
    homologs = pd.read_csv(homologs)
    
    if verbose:
        print("Intersecting gene homologs...")
    
    #Intersect gene homologs
    mouse, human = intersect_gene_homologs(mouse, human, homologs)
    
    #Create output file names
    mousefile_out = os.path.basename(mousefile)
    humanfile_out = os.path.basename(humanfile)
    
    mousefile_out = sub('.csv', '_homologs.csv', mousefile_out)
    humanfile_out = sub('.csv', '_homologs.csv', humanfile_out)
    
    mousefile_out = os.path.join(outdir, mousefile_out)
    humanfile_out = os.path.join(outdir, humanfile_out)
    
    if verbose:
        print("Writing to file...")
    
    #Write to file
    mouse.to_csv(mousefile_out)
    human.to_csv(humanfile_out)
    
    return

if __name__=='__main__':
    main()
