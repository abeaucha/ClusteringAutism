# ----------------------------------------------------------------------------
# template.py 
# Author: Antoine Beauchamp
# Created: 

"""
Brief description.

Description
-----------
"""

# Packages -------------------------------------------------------------------

import argparse
import os


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--arg',
        type = str,
        help = ("Help message")
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
    
#     #Parse command line arguments
#     args = parse_args()
#     arg = args['arg1']
#     verbose = True if args['verbose'] == 'true' else False
    
#     if verbose:
#         print("")

    input_dir =     

    
    return
    
if __name__=='__main__':
    main()
