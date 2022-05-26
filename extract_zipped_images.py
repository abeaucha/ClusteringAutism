# ----------------------------------------------------------------------------
# extract_jacobian_images.py 
# Author: Antoine Beauchamp
# Created: May 9th, 2022

"""
Extract images from files compressed using GNU zip.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import sys
import multiprocessing      as mp
from glob                   import glob
from re                     import sub
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
        help = ("Path to directory containing compressed images.")
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
    
    args = vars(parser.parse_args())
    
    return args


# Functions ------------------------------------------------------------------

def extract_minc(gzfile):
    
    """
    Extract MINC image from compressed file
    
    Arguments
    ---------
    gzfile: str
        Path to image file to extract.
        
    Returns
    -------
    None
    """
    
    os.system('gunzip -f -k {}'.format(gzfile))
    gzfile = sub(r'.gz', '', gzfile)
    os.system('nii2mnc -quiet -clobber {}'.format(gzfile))
    os.system('rm {}'.format(gzfile))
    
    return


# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    parallel = True if args['parallel'] == 'true' else False
    nproc = args['nproc']
    
    #Get paths to compressed files
    imgdir = os.path.join(imgdir, '')
    files = glob(imgdir+'*.nii.gz')
    
    #Extract MINC files in parallel or sequentially
    if parallel:
        
        pool = mp.Pool(nproc)
        tmp_list = []
        for tmp in tqdm(pool.imap(extract_minc, files), total = len(files)):
            tmp_list.append(tmp)
            
        pool.close()
        pool.join()
        
    else:
        
        tmp_list = list(map(extract_minc, tqdm(files)))
    
    return


if __name__=='__main__':
    main()
