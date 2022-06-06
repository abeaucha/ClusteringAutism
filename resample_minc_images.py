# ----------------------------------------------------------------------------
# resample_minc_images.py 
# Author: Antoine Beauchamp
# Created: May 16th, 2022

"""
Resample a set of MINC images using the autocrop command line tool.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import multiprocessing      as mp
from glob                   import glob
from re                     import sub
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
        help = ("Path to directory containing images to resample.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        help = ("Path to directory in which to save resampled images.")
    )
    
    parser.add_argument(
        '--isostep',
        type = str,
        help = ("Resolution of voxels in resampled images (mm).")
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

def resample_image(infile, outdir, isostep):
    
    """
    Resample a MINC image
    
    Arguments
    ---------
    infile: str
        Path to the MINC file to resample.
    outdir: str
        Path to the directory in which to save the resampled image.
    isostep: float
        Isotropic dimension of voxels in resampled image (millimeters).
    
    
    Returns
    -------
    None
    """

    outfile = sub(r'.mnc', '_resampled_{}.mnc'.format(isostep), infile)
    outfile = os.path.basename(outfile)
    outfile = os.path.join(outdir, outfile)
    cmd_autocrop = ('autocrop -quiet -clobber -isostep {} {} {}'
                    .format(isostep, infile, outfile))
    os.system(cmd_autocrop)

    return
    

# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outdir = args['outdir']
    isostep = args['isostep']
    parallel = True if args['parallel'] == 'true' else False
    nproc = args['nproc']
    
    #Ensure proper paths
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    
    #Create outdir if needed
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #Get MINC images in dir
    files = glob(imgdir+'*.mnc')
    
    #Create partial function for mapping
    resample_partial = partial(resample_image, 
                               outdir = outdir,
                               isostep = isostep)

    if parallel:
        
        pool = mp.Pool(nproc)
        tmp_list = []
        for tmp in tqdm(pool.imap(resample_partial, files), 
                        total = len(files)):
            tmp_list.append(tmp)
            
        pool.close()
        pool.join()
        
    else:
        
        tmp_list = list(map(resample_partial, tqdm(files)))
    
    return

if __name__=='__main__':
    main()
