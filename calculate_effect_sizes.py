# ----------------------------------------------------------------------------
# calculate_effect_sizes.py 
# Author: Antoine Beauchamp
# Created: May 9th, 2022

"""


Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import numpy as np
import pandas as pd
from datatable import fread


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--demographics',
        type = str,
        help = ("Path to CSV file containing demographics data.")
    )
    
    parser.add_argument(
        '--imgdir',
        type = str,
        help = ()
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        help = ()
    )
    
    args = vars(parser.parse_args())
    
    return args


# Functions ------------------------------------------------------------------


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    demographics_file = args['demographics']
    imgdir = args['imgdir']
    outdir = args['outdir']
    
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    
    demographics = (fread(demographics_file, header = True)
                    .to_pandas())
    
    
    demographics.dropna(subset = ['DX'],
                        inplace = True)
    
    controls = demographics.loc[demographics['DX'] == 'Control']
    participants = demographics.loc[demographics['DX'] != 'Control']

    print(participants.columns)
    
    i = 0
    
    print(participants.loc[0, 'Extract_ID'])
    
    quit()
    
    
    return

if __name__=='__main__':
    main()