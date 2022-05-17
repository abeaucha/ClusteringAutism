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
        '--imgdir',
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

def import_image(img, mask):
    
    """
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
    outdir = args['outdir']
    outfile = args['outfile']
    mask = args['mask']
    parallel = True if args['parallel'] == 'true' else False
#     verbose = True if args['verbose'] == 'true' else False

    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    
    imgfiles = glob(imgdir+'*.mnc')
    
    import_image_partial = partial(import_image,
                                   mask = mask)
    
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

    df_imgs = pd.DataFrame(np.asarray(imgs))
    df_imgs['imagefile'] = [os.path.basename(i) for i in imgfiles]

    outfile = os.path.join(outdir, outfile)
    df_imgs.to_csv(outfile)

    return
    
    
if __name__=='__main__':
    main()