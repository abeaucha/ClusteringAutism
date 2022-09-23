# ----------------------------------------------------------------------------
# mouse_cluster_signatures.py 
# Author: Antoine Beauchamp
# Created: May 3rd, 2022

"""
Create gene expression signatures for mouse imaging clusters.

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
        '--cluster-dir',
        type = str,
        help = ("Directory containing cluster mask images.")
    )
    
    parser.add_argument(
        '--expr-dir',
        type = str,
        help = ("Directory containing gene expression .csv files.")
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        help = ("Mask file used to generate the expression data sets."))
    
    parser.add_argument(
        '--sign',
        type = str,
        choices = ['positive', 'negative'],
        help = ()
    )
    
    parser.add_argument(
        '--outfile',
        type = str,
        default = 'cluster_signatures.csv',
        help = ("Path to the .csv file in which to export cluster ",
                "signatures.")
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
    

def import_image(img, mask):
    
    """
    Import and mask a MINC image.
    
    Arguments
    ---------
    img: str
        Path to the MINC file to import.
    mask: str
        Path to the mask MINC file. 
    
    Returns
    -------
    img_masked: numpy.ndarray
        A 1-dimensional NumPy array containing the masked image voxels.
    """
    
    mask_vol = volumeFromFile(mask)
    mask_array = np.array(mask_vol.data.flatten())
    mask_vol.closeVolume()
    
    img_vol = volumeFromFile(img)
    img_array = np.array(img_vol.data.flatten())
    img_vol.closeVolume()
    
    img_masked = img_array[mask_array == 1]
    
    return img_masked


def aggregate_expression(cluster, expression):

    """
    Aggregate expression values in a cluster.
    
    Arguments
    ---------
    cluster: numpy.ndarray
        A NumPy array of dimension 1 containing a voxel-wise binary 
        mask with values {0,1}. 
    expression: pandas.core.frame.DataFrame
        A DataFrame with rows corresponding to voxels and columns 
        corresponding to gene expression variables.
    
    Returns
    -------
    numpy.ndarray
        A NumPy array of dimension 1 containing the aggregated 
        expression values.
    """
        
    return np.array(expression
                    .loc[cluster == 1]
                    .mean())


def create_cluster_signature(infiles, mask, sign = None):
    
    """
    Create an expression signature for a cluster. 
    
    Arguments
    ---------
    infiles: tuple of str
        Tuple with two elements containing 1. the path to the cluster 
        mask image and 2. the path to the gene expression .csv file.
    mask: str
        Path to mask image file used to create the gene expression 
        files.
    sign: str
        One of {'positive', 'negative'} indicating whether to use
        only negative or positive mask values. All values are used
        if None. 
    
    Returns
    -------
    signature: numpy.ndarray
        A NumPy array of dimension 1 containing the expression 
        signature. 
    """
        
    if all([sign != val for val in [None, 'positive', 'negative']]):
        raise ValueError("sign must be one of {None, 'positive', 'negative'}")

    cluster = import_image(img = infiles[0],
                           mask = mask)

    cluster = np.int8(np.floor(cluster))

    if sign == 'positive':
        cluster[cluster < 0] = 0
    elif sign == 'negative':
        cluster[cluster > 0] = 0

    cluster = np.abs(cluster)

    with silence():
        expression = (fread(infiles[1], 
                            header = True)
                      .to_pandas())

    signature = aggregate_expression(cluster = cluster,
                                     expression = expression)
    
    return signature


def build_signatures_table(clusters, expression, mask, sign = None,
                           parallel = False, nproc = None):
    
    """
    Create a DataFrame containing cluster expression signatures.
    
    Arguments
    ---------
    clusters: list of str
        Paths to cluster mask image files.
    expression: list of str
        Paths to gene expression .csv files.
    mask: str
        Path to mask image file used to create the gene expression 
        files.
    sign: str
        One of {'positive', 'negative'} indicating whether to use
        only negative or positive mask values. All values are used
        if None.
    parallel: bool (optional)
        Run in parallel.
    nproc: int (optional)
        Number of processors to use in parallel.
    
    Returns
    -------
    df_signatures: pandas.core.frame.DataFrame
        A DataFrame containing expression signatures for each cluster.
    """
        
    infiles = list(product(clusters, expression))
    
    create_cluster_signature_partial = partial(create_cluster_signature, 
                                               mask = mask,
                                               sign = sign)
    
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
    df_signatures['cluster_file'] = [os.path.basename(i[0]) for i in infiles]
    df_signatures['expr_file'] = [os.path.basename(i[1]) for i in infiles]
    df_signatures['mask_sign'] = sign
    
    return df_signatures

# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    cluster_dir = args['cluster_dir']
    expr_dir = args['expr_dir']
    mask_file = args['mask']
    sign = args['sign']
    outfile = args['outfile']
    parallel = True if args['parallel'] == 'true' else False
    verbose = True if args['verbose'] == 'true' else False
    
    #Create outdir if needed
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #Ensure directories are proper paths
    cluster_dir = os.path.join(cluster_dir, '')
    expr_dir = os.path.join(expr_dir, '')
    
    #Input files
    expr_files = glob(expr_dir+'*.csv')
    cluster_files = glob(cluster_dir+'*.mnc')
    
    #Build signatures data frame        
    if verbose:
        print("Creating cluster signatures...")
    
    df_signatures = build_signatures_table(clusters = cluster_files,
                                           expression = expr_files,
                                           mask = mask_file,
                                           sign = sign,
                                           parallel = parallel,
                                           nproc = args['nproc'])

    #Write to file
    if verbose:
        print("Writing to file...")
        
    outfile = args['outfile']
    df_signatures.to_csv(outfile, index = False)
    
    return
    
if __name__=='__main__':
    main()
