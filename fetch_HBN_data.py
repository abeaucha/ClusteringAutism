#!.venv/bin/python3
# ----------------------------------------------------------------------------
# fetch_HBN_data.py
# Author: Antoine Beauchamp

"""
Fetch Healthy Brain Network T1w images from Bilal's directory
"""


# Packages -------------------------------------------------------------------
import multiprocessing as mp
import os
from shutil import copyfile
from tqdm import tqdm


# Functions ------------------------------------------------------------------
def copier(paths): copyfile(src = paths[0], dst = paths[1])


# Main -----------------------------------------------------------------------
if __name__ == '__main__':

    # Directories
    input_dir = '/projects/bilal/HBN/'
    output_dir = 'data/human/HBN/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # HBN subject directories
    subjects = [d for d in os.listdir(input_dir) if 'sub' in d]

    # Collate input file paths
    input_files = []
    for i, s in enumerate(subjects):
        imgdir = os.path.join(input_dir, s, 'anat')
        if os.path.exists(imgdir):
            imgfiles = [f for f in os.listdir(imgdir) if 'T1w' in f]
            imgfiles = [f for f in imgfiles if 'nii.gz' in f]
            imgfiles = [f for f in imgfiles if 'VNav' not in f]
            if len(imgfiles) > 0:
                for file in imgfiles:
                    input_files.append(os.path.join(imgdir, file))

    # Output file paths
    output_files = [os.path.join(output_dir, os.path.basename(f))
                    for f in input_files]

    # Zip input and outputs
    paths = list(zip(input_files, output_files))

    # Copy images
    pool = mp.Pool(8)
    results = []
    for result in tqdm(pool.imap(copier, paths), total = len(paths)):
        results.append(result)
    pool.close()
    pool.join()