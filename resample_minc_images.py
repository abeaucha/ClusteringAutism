# ----------------------------------------------------------------------------
# resample_minc_images.py 
# Author: Antoine Beauchamp
# Created: May 16th, 2022

"""
Downsample a set of MINC images. 
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
        help = ("Path to directory containing images to downsample.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        help = ("Path to directory in which to save downsampled images.")
    )
    
    parser.add_argument(
        '--isostep',
        type = str,
        help = ("Resolution of voxels in downsampled images (mm).")
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

def downsample_image(infile, outdir, isostep):
    
    """
    """
    
    outfile = sub(r'.mnc', '_downsampled_{}.mnc'.format(isostep), infile)
    outfile = os.path.basename(outfile)
    outfile = os.path.join(outdir, outfile)
    cmd_autocrop = ('autocrop -quiet -clobber -isostep {} {} {}'
                    .format(isostep, infile, outfile))
    os.system(cmd_autocrop)

    return
    

# Main ------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    imgdir = args['imgdir']
    outdir = args['outdir']
    isostep = args['isostep']
    parallel = True if args['parallel'] == 'true' else False
    nproc = args['nproc']
    
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    
    files = glob(imgdir+'*.mnc')
    
    downsample_partial = partial(downsample_image, 
                                 outdir = outdir,
                                 isostep = isostep)

    if parallel:
        
        pool = mp.Pool(nproc)
        tmp_list = []
        for tmp in tqdm(pool.imap(downsample_partial, files), 
                        total = len(files)):
            tmp_list.append(tmp)
            
        pool.close()
        pool.join()
        
    else:
        
        tmp_list = list(map(downsample_partial, tqdm(files)))
    
    return


if __name__=='__main__':
    main()