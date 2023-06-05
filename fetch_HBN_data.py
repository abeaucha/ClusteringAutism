import os

input_dir = '/projects/bilal/HBN/'
output_dir = 'data/human/HBN/'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

input_subdirs = [f for f in os.listdir(input_dir) if 'sub' in f]

input_subdirs[0]