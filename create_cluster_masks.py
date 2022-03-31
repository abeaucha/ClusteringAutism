import argparse
import os
import numpy as np
from pyminc.volumes.factory import volumeLikeFile, volumeFromFile
from glob import glob
from re import sub
from functools import partial
from tqdm import tqdm


def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/mouse/cluster_maps/',
        help = ("")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/mouse/cluster_masks/',
        help = ("")
    )
    
    parser.add_argument(
        '--threshold',
        type = float,
        default = 0.5,
        help = ("")
    )
    
    parser.add_argument(
        '--symmetric',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("")
    )
    
    args = vars(parser.parse_args())
    
    return args


def create_cluster_mask(infile, outfile, threshold = 0.2, symmetric = True, comparison = '>'):
    
    """
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


def main():
    
    args = parse_args()
    datadir = args['datadir']
    outdir = args['outdir']
    threshold = args['threshold']
    symmetric = True if args['symmetric'] == 'true' else False
    
    datadir = os.path.join(datadir, '')
    outdir = os.path.join(outdir, '')
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    infiles = glob(datadir+'*.mnc')
    
    outfiles = [os.path.basename(file) for file in infiles]
    outfiles = [sub(r'.mnc', '', file) for file in outfiles]
    outfiles = [file+'_mask_threshold{}.mnc'.format(threshold) for file in outfiles]
    outfiles = [os.path.join(outdir, file) for file in outfiles]
    
    create_cluster_mask_partial = partial(create_cluster_mask,
                                          threshold = threshold,
                                          symmetric = symmetric)
    
    list(map(create_cluster_mask_partial, tqdm(infiles), outfiles))
    
    
    return
    
    
if __name__=='__main__':
    main()