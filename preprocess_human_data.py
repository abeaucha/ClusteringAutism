# ----------------------------------------------------------------------------
# preprocess_human_data.py 
# Author: Antoine Beauchamp
# Created: February 22nd, 2023

"""
Preprocess human data files.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os
from src import utils


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--basedir',
        type = str,
        help = ("Path to base directory containing Jacobians files.")
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

def template(x):
    
    """
    Template
    
    Arguments
    ---------
    x: dtype
        Argument description.
        
    Returns
    -------
    None
    """
    
    return
    

# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    basedir = args['basedir']
    parallel = True if args['parallel'] == 'true' else False
    nproc = args['nproc']
    verbose = True if args['verbose'] == 'true' else False
  
    #Base directory checks
    if basedir is None:
        raise Exception("Argument --basedir must be specified")
    else:
        basedir = os.path.join(basedir, '')
        if not os.path.exists(basedir):
            OSError("Base directory not found: {}".format(basedir))

    #Parallel checks
    if parallel:
        if nproc is None:
            raise Exception("Argument --nproc must be specified when --parallel true")

    #Iterate over jacobians
    jacobians = ['absolute', 'relative']
    for jac in jacobians:
    
        if verbose:
            print("Processing {} Jacobians files...".format(jac))
    
        #Input files
        input_dir = os.path.join(basedir, jac, 'smooth', '')
        if not os.path.exists(input_dir):
            raise OSError("Input directory not found: {}".format(input_dir))
        input_files = [os.path.join(input_dir, infile) 
                       for infile in os.listdir(input_dir) 
                       if '.gz' in infile] 
    
        #Output directory
        output_dir = os.path.join(basedir, jac, 'smooth_minc', '')
    
        if verbose:
            print("\tExtracting compressed images...")
        unzipped_files = utils.gunzip_files(infiles = input_files,
                                            keep = True,
                                            parallel = parallel,
                                            nproc = nproc)

        if verbose:
            print("\tConverting NIFTY files to MINC...")
        imgfiles = utils.convert_images(infiles = unzipped_files,
                                        input_format = 'nifty',
                                        output_format = 'minc',
                                        outdir = output_dir,
                                        keep = True,
                                        parallel = parallel,
                                        nproc = nproc)

    return
    
if __name__=='__main__':
    main()
