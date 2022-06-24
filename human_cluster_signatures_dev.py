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
from mouse_cluster_signatures import import_image
from glob import glob
from datatable import fread
import numpy as np
from pyminc.volumes.factory import volumeFromFile


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
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
    verbose = True if args['verbose'] == 'true' else False

    mask_file = 'data/human/registration/reference_files/mask_1.0mm.mnc'
    
    cluster_dir = 'data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_0.1/'
    cluster_dir = os.path.join(cluster_dir, '')
    
    cluster_files = glob(cluster_dir+'*.mnc')
    
    microarray_labels = 'data/human/expression/AHBA_microarray_samples_mask_studyspace_1.0mm.mnc'
    
    labels_vol = volumeFromFile(microarray_labels)
    labels = np.array(labels_vol.data.flatten())
    print(np.sum(labels))
    quit()
    
    microarray_labels = import_image(img = microarray_labels, mask = mask_file)
    
    microarray_defs = 'data/human/expression/AHBA_microarray_samples_defs.csv'
    microarray_defs = (fread(microarray_defs, 
                            header = True)
                      .to_pandas())
    print(microarray_defs.shape)
    print(np.sum(microarray_labels != 0))
    quit()
    
    cluster = import_image(img = cluster_files[0], mask = mask_file)
    
    print(len(cluster))
    print(len(microarray_labels))
    
    microarray_in_cluster = microarray_labels[cluster == 1]
    microarray_in_cluster = microarray_in_cluster[microarray_in_cluster != 0]

    print(len(microarray_in_cluster))
    
    print(np.unique(microarray_in_cluster))
    
#     expr_file = 'data/human/expression/input_space/HumanExpressionMatrix_samples_pipeline_v1_homologs_scaled.csv'
#     expression = (fread(expr_file, 
#                             header = True)
#                       .to_pandas())
    
    
    
    if verbose:
        print("")
    
    return
    
if __name__=='__main__':
    main()
