# ----------------------------------------------------------------------------
# create_cluster_signatures.py 
# Author: Antoine Beauchamp
# Created: May 3rd, 2022

"""


Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import sys
import contextlib
import numpy                as np
import pandas               as pd
import multiprocessing      as mp
from glob                   import glob
from datatable              import fread
from pyminc.volumes.factory import volumeFromFile
from functools              import partial
from itertools              import product
from tqdm                   import tqdm
from io                     import StringIO


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--clusterdir',
        type = str,
        help = ("Directory containing cluster masks.")
    )
    
    parser.add_argument(
        '--exprdir',
        type = str,
        help = ("Directory containing expression data sets.")
    )
    
    parser.add_argument(
        '--exprglob',
        type = str,
        help = ("Glob to get expression CSV files.")
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        help = ("Mask file used to generated the expression data sets."))
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = ("Directory in which to save cluster signatures file.")
    )
        
    parser.add_argument(
        '--outfile',
        type = str,
        default = 'cluster_signatures.csv',
        help = ("Name of CSV file in which to export cluster signatures.")
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

@contextlib.contextmanager
def silence():
    save_stdout = sys.stdout
    sys.stdout = StringIO()
    yield
    sys.stdout = save_stdout
    

def import_cluster_mask(cluster, mask):
    
    """
    """
        
    mask_vol = volumeFromFile(mask)
    mask_array = np.array(mask_vol.data.flatten())
    mask_vol.closeVolume()
    
    cluster_vol = volumeFromFile(cluster)
    cluster_array = np.array(cluster_vol.data.flatten())
    cluster_vol.closeVolume()
    
    cluster_masked = cluster_array[mask_array == 1]
    
    return cluster_masked


def aggregate_expression(cluster, expression):

    """
    """
        
    return np.array(expression
                    .loc[cluster == 1]
                    .mean())


def create_cluster_signature(infiles, mask):
    
    """
    """
        
    cluster = import_cluster_mask(cluster = infiles[0],
                                  mask = mask)
    
    with silence():
        expression = (fread(infiles[1], 
                            header = True)
                      .to_pandas())
    
    signature = aggregate_expression(cluster = cluster,
                                     expression = expression)
    
    return signature


def build_signatures_table(clusters, expression, mask, 
                           parallel = False, nproc = None):
    
    """
    """
        
    infiles = list(product(clusters, expression))
    
    create_cluster_signature_partial = partial(create_cluster_signature, 
                                               mask = mask)
    
    if parallel:
        
        if nproc is None:
            nproc = mp.cpu_count()
            
        pool = mp.Pool(nproc)
        
        signatures = []
        for signature in tqdm(pool.imap(create_cluster_signature_partial,
                                        infiles),
                              total = len(infiles)):
            signatures.append(signature)
            
        pool.close()
        pool.join()
        
    else:
    
        signatures = list(map(create_cluster_signature_partial,
                               tqdm(infiles)))
        
    df_signatures = pd.DataFrame(np.asarray(signatures))
    df_signatures['clusterfile'] = [os.path.basename(i[0]) for i in infiles]
    df_signatures['exprfile'] = [os.path.basename(i[1]) for i in infiles]
    
    return df_signatures


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    cluster_dir = args['clusterdir']
    expr_dir = args['exprdir']
    expr_glob = args['exprglob']
    parallel = True if args['parallel'] == 'true' else False
    verbose = True if args['verbose'] == 'true' else False
    
    #Ensure directories are proper paths
    cluster_dir = os.path.join(cluster_dir, '')
    expr_dir = os.path.join(expr_dir, '')
    
    #Input files
    expr_files = glob(expr_dir+expr_glob)
    cluster_files = glob(cluster_dir+'*.mnc')
    mask_file = args['mask']
        
    #Build signatures data frame        
    if verbose:
        print("Creating cluster signatures...")
    
    df_signatures = build_signatures_table(clusters = cluster_files,
                                           expression = expr_files,
                                           mask = mask_file,
                                           parallel = parallel,
                                           nproc = args['nproc'])

    #Write to file
    if verbose:
        print("Writing to file...")
        
    outdir = args['outdir']
    outfile = args['outfile']
    outfile = os.path.join(outdir, outfile)
        
    df_signatures.to_csv(outfile, index = False)
    
    return
    
    
if __name__=='__main__':
    main()
