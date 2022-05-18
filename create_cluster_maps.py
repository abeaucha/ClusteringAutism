# ----------------------------------------------------------------------------
# create_cluster_maps.py 
# Author: Antoine Beauchamp
# Created: May 18th, 2022

"""

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
        '--clusterfile',
        type = str,
        help = ("")
    )
    
        
    args = vars(parser.parse_args())
    
    return args


# Functions ------------------------------------------------------------------


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    clusterfile = args['clusterfile']
    
    df_clusters = fread(clusterfile, header = True).to_pandas()
    print(df_clusters.head())
    quit()
    
    
if __name__=='__main__':
    main()