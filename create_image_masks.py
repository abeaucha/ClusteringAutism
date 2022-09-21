# ----------------------------------------------------------------------------
# create_image_masks.py 
# Author: Antoine Beauchamp
# Created: May 18th, 2022

"""
Create masks for a set of MINC images.

Description
-----------
This script creates a mask for each image in a set using one of two 
thresholding methods: 

1. Using an intensity threshold 
2. Using the top number or fraction of voxels. 
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import sys
import numpy                as np
import pandas               as pd
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
        help = ("Path to directory containing images to mask.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        help = ("Path to directory in which to masks.")
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        help = ("Path to optional mask to apply before thresholding.")
    )
        
    parser.add_argument(
        '--method',
        type = str,
        default = 'intensity',
        choices = ['intensity', 'top_n'],
        help = "Method to use for thresholding."
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
    
    parser.add_argument(
        '--comparison',
        type = str,
        default = 'gt',
        choices = ['gt', 'lt'],
        help = ("How to apply the intensity threshold. Ignored if method ",
                "is 'top_n'.")
    )
    
    args = vars(parser.parse_args())
    
    return args


# Functions ------------------------------------------------------------------

def threshold_intensity(data, threshold, symmetric = True, comparison = 'gt'):
    
    """
    Create a mask using an intensity threshold.
    
    Arguments
    --------
    data: numpy.ndarray
        Data to threshold
    threshold: float
        Intensity value at which to threshold.  
    symmetric: bool
        Option to apply threshold symmetrically to negative and 
        positive values.
    comparison: str
        One of 'gt' or 'lt' indicating whether the mask should 
        return values greater than or lesser than the threshold.
    
    Returns
    -------
    out: numpy.ndarray
        Mask corresponding to thresholded data.
    """
    
    if symmetric:
        if comparison == 'gt':
            ind = np.abs(data) > threshold
        else:
            ind = np.abs(data) < threshold
    else:
        if comparison == 'gt':
            ind = data > threshold
        else:
            ind = data < threshold
    out = np.zeros_like(data)
    out[ind] = 1
    return out
        

def threshold_top_n(data, n, symmetric = True, tolerance = 1e-5):
    
    """
    Create a mask using the top n values.
    
    Arguments
    ---------
    data: numpy.ndarray
        Data to threshold
    n: int or float
        If |n| > 1, the top n values will be used. 
        If 0 < |n| < 1, the top fraction of values will be used.
        If n < 0, the top negative values will be used. 
    symmetric: bool
        Option to apply the threshold symmetrically.
    tolerance: float
        Absolute values below this threshold will be excluded.
    
    Returns
    -------
    out: numpy.ndarray
        Mask corresponding to thresholded data.
    """

    #Volume indices
    ni, nj, nk = data.shape
    indices =[[i,j,k] for i in range(ni) 
              for j in range(nj) 
              for k in range(nk)]
    indices = np.asarray(indices)

    #Volume values
    values = np.array(data.flatten())
    values = values.reshape(len(values), 1)

    #Combine indices and values
    df = pd.DataFrame(np.concatenate((indices, values), axis = 1),
                columns = ['i', 'j', 'k', 'value'])

    #Symmetric option
    if symmetric:    
        df['value'] = np.abs(df['value'])

    #Tolerance filter
    ind_tolerance = np.abs(df['value']) > tolerance

    #Filter top/bottom values
    if n > 0:
        ind_sign = df['value'] > 0
        df_top_n = df.loc[ind_tolerance & ind_sign]
        if n < 1:
            n = int(np.floor(n*df_top_n.shape[0]))
        df_top_n = df_top_n.nlargest(n = n, columns = 'value')
    elif n < 0:
        ind_sign = df['value'] < 0
        df_top_n = df.loc[ind_tolerance & ind_sign]
        n = abs(n)
        if n < 1:
            n = int(np.floor(n*df_top_n.shape[0]))
        df_top_n = df_top_n.nsmallest(n = n, columns = 'value')
    else:
        raise ValueError("n cannot be 0")

    #Extract top indices
    top_indices = df_top_n.loc[:,['i', 'j', 'k']].to_numpy()

    #Create mask
    n = top_indices.shape[0]
    out = np.zeros_like(data)
    for row in range(n):
        i = int(top_indices[row,0])
        j = int(top_indices[row,1])
        k = int(top_indices[row,2])
        out[i,j,k] = 1
        
    return out
        

def create_image_mask(infile, outfile, maskfile = None, method = 'intensity',
                      threshold = 0.5, symmetric = True, comparison = 'gt'):
    
    """
    Create a threshold mask based on an input image. 
    
    Arguments
    ---------
    infile: str
        Path to the MINC file containing the image to mask.
    outfile: str
        Path to the MINC file in which to save the mask.
    maskfile: str
        Path to a MINC file containing a mask to apply before
        thresholding.
    method: str
        One of {'intensity', 'top_n'} indicating which method 
        to use for thresholding.
    threshold: float
        Threshold to determine which voxels are in the mask. 
        If method = 'intensity', this is the intensity value at
        which to threshold. 
        If method = 'top_n', this is the number of fraction of top
        voxels to use.
    symmetric: bool
        Option to apply threshold symmetrically to negative and 
        positive values.
    comparison: str
        String indicating how to apply the threshold if 
        method = 'intensity'. This is ignored if method = 'top_n'.
    
    Returns
    -------
    None
    """
    
    img_vol = volumeFromFile(infile)
    img_data = np.array(img_vol.data)
    img_vol.closeVolume()
    
    if maskfile is not None:
        mask_vol = volumeFromFile(maskfile)
        mask_data = np.array(mask_vol.data)
        mask_vol.closeVolume()
        img_data[mask_data == 0] = 0
    
    if method == 'intensity':
        img_mask_data = threshold_intensity(data = img_data,
                                            threshold = threshold,
                                            symmetric = symmetric,
                                            comparison = comparison)
    elif method == 'top_n':
        img_mask_data = threshold_top_n(data = img_data,
                                        n = threshold, 
                                        symmetric = symmetric)
    else:
        raise ValueError("method must be one of {'intensity', 'top_n'}")
        
    img_mask_vol = volumeLikeFile(likeFilename = infile,
                                  outputFilename = outfile,
                                  labels = True)
        
    img_mask_vol.data = img_mask_data
    img_mask_vol.writeFile()
    img_mask_vol.closeVolume()
    
    return


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outdir = args['outdir']
    mask = args['mask']
    method = args['method']
    threshold = args['threshold']
    symmetric = True if args['symmetric'] == 'true' else False
    comparison = args['comparison']
    
    #Cast threshold to integer if number of voxels used
    if (method == 'top_n') & (abs(threshold) > 1):
        threshold = int(threshold)
    
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
    outfiles = [file+'_mask' for file in outfiles]
    if method == 'intensity':
        outfiles = [file+'_intensity' for file in outfiles]
    elif method == 'top_n':
        outfiles = [file+'_topn' for file in outfiles]
    if symmetric:
        outfiles = [file+'_symmetric' for file in outfiles]
    else:
        if threshold > 0:
            outfiles = [file+'_positive' for file in outfiles]
        elif threshold < 0:
            outfiles = [file+'_negative' for file in outfiles]
    outfiles = [file+'_threshold{}.mnc'.format(abs(threshold)) 
                for file in outfiles]
    outfiles = [os.path.join(outdir, file) for file in outfiles]
    
    #Partial function for iteration
    create_image_mask_partial = partial(create_image_mask,
                                        maskfile = mask,
                                        method = method,
                                        threshold = threshold,
                                        symmetric = symmetric,
                                        comparison = comparison)
    
    #Iterate over images
    list(map(create_image_mask_partial, tqdm(infiles), outfiles))
    
    return
    
if __name__=='__main__':
    main()
