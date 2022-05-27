# ----------------------------------------------------------------------------
# create_cluster_mask.py 
# Author: Antoine Beauchamp
# Created: May 18th, 2022

"""
Create a mask indicating which voxels are in a cluster.

Description
-----------
This script creates a set of masks to indicate which voxels belong to 
a cluster. This done by thresholding a representative cluster map image 
to retain only those voxels that pass the threshold.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import numpy                as np
from pyminc.volumes.factory import volumeLikeFile, volumeFromFile
from glob                   import glob
from re                     import sub
from functools              import partial
from tqdm                   import tqdm


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--imgdir',
        type = str,
        help = ("Path to directory containing representative cluster map "
                "images.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        help = ("Path to directory in which to save cluster masks.")
    )
    
    parser.add_argument(
        '--threshold',
        type = float,
        default = 0.5,
        help = ("Threshold to determine which voxels are in the mask.")
    )
    
    parser.add_argument(
        '--symmetric',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to apply threshold symmetrically to negative and ",
                "positive values.")
    )
    
    args = vars(parser.parse_args())
    
    return args


# Functions ------------------------------------------------------------------

def create_cluster_mask(infile, outfile, threshold = 0.2, symmetric = True, 
                        comparison = '>'):
    
    """
    Create a cluster mask.
    
    Arguments
    ---------
    infile: str
        Path to the MINC file containing the cluster map to use to 
        create the mask.
    outfile: str
        Path to the MINC file in which to save the cluster mask.
    threshold: float
        Threshold to determine which voxels are in the mask.
    symmetric: bool
        Option to apply threshold symmetrically to negative and 
        positive values.
    comparison: str
        Symbol indicating which logical comparison to use to determine 
        cluster voxels.
    
    Returns
    -------
    None
    """
    
    cluster_mask = volumeLikeFile(likeFilename = infile,
                                 outputFilename = outfile,
                                 labels = True)
    
    cluster_map = volumeFromFile(infile)
    
    if symmetric:
        if comparison == '>':
            ind = np.abs(cluster_map.data) > threshold
        else:
            ind = np.abs(cluster_map.data) < threshold
    else:
        if comparison == '>':
            ind = cluster_map.data > threshold
        else:
            ind = cluster_map.data < threshold
            
    cluster_mask.data[ind] = 1

    cluster_mask.writeFile()
    cluster_mask.closeVolume()
    
    cluster_map.closeVolume()
    
    return


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outdir = args['outdir']
    threshold = args['threshold']
    symmetric = True if args['symmetric'] == 'true' else False
    
    #Ensure proper paths
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    
    #Create outdir if needed
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #Get input images
    infiles = glob(imgdir+'*.mnc')
    
    #Create output file names
    outfiles = [os.path.basename(file) for file in infiles]
    outfiles = [sub(r'.mnc', '', file) for file in outfiles]
    outfiles = [file+'_mask_threshold{}.mnc'.format(threshold) 
                for file in outfiles]
    outfiles = [os.path.join(outdir, file) for file in outfiles]
    
    #Partial function for iteration
    create_cluster_mask_partial = partial(create_cluster_mask,
                                          threshold = threshold,
                                          symmetric = symmetric)
    
    #Iterate over images
    list(map(create_cluster_mask_partial, tqdm(infiles), outfiles))
    
    return
    
if __name__=='__main__':
    main()
