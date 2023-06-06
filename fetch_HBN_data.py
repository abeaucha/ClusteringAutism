#!.venv/bin/python3

# Packages
import os
from shutil import copyfile

# Directories
input_dir = '/projects/bilal/HBN/'
output_dir = 'data/human/HBN/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Input sub-directories
input_dirs = [os.path.join(input_dir, d, 'anat', '')
              for d in os.listdir(input_dir) if 'sub' in d]

# Copy T1w images
for d in input_dirs:
    if os.path.exists(d):
        files = [f for f in os.listdir(d) if 'T1w' in f]
        files = [f for f in files if 'nii.gz' in f]
        if len(files) != 0:
            for f in files:
                src = os.path.join(d, f)
                dst = os.path.join(output_dir, f)
                copyfile(src = src, dst = dst)

