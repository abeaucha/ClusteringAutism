# ----------------------------------------------------------------------------
# build_effect_size_matrix.py 
# Author: Antoine Beauchamp
# Created: May 16th, 2022

"""
Build a matrix containing voxel-wise effect sizes.
"""


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
        '--imgdir',
        type = str,
        help = ("Directory containing effect size MINC images.")
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        help = ("Path to mask file."))
    
    parser.add_argument(
        '--outfile',
        type = str,
        default = 'effect_sizes.csv',
        help = ("Path to .csv file in which to export effect size matrix.")
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


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outfile = args['outfile']
    mask = args['mask']
    parallel = True if args['parallel'] == 'true' else False
    verbose = True if args['verbose'] == 'true' else False

    #Ensure proper paths
    imgdir = os.path.join(imgdir, '')
    
    #Get MINC images in dir
    imgfiles = glob(imgdir+'*.mnc')
    
    #Create partial function for mapping
    import_image_partial = partial(import_image,
                                   mask = mask)
    
    if verbose:
        print("Importing images...")
    
    #Import images
    if parallel:
        if nproc is None:
            nproc = mp.cpu_count()
            
        pool = mp.Pool(nproc)
        
        imgs = []
        for img in tqdm(pool.imap(import_image_partial, imgfiles),
                        total = len(imgfiles)):
            imgs.append(img)
            
        pool.close()
        pool.join()
        
    else:
        
        imgs = list(map(import_image_partial, tqdm(imgfiles)))

    if verbose:
        print("Building matrix...")
        
    #Convert to data frame
    df_imgs = pd.DataFrame(np.asarray(imgs))
    df_imgs['file'] = imgfiles
    
    if verbose:
        print("Writing to file...")

    #Save matrix
    df_imgs.sort_values(by = 'file').to_csv(outfile, index = False)

    return
    
if __name__=='__main__':
    main()