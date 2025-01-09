#!/usr/bin/env python3

import os
import utils
from itertools import product
import pandas as pd

def generate_cluster_pairs(centroid_dirs, jacobians = ('absolute', 'relative')):
    """
    Generate pairs of cluster centroid images.

    Parameters
    ----------
    centroid_dirs: list of str
        List of paths to the directories containing the centroid images.
        Expects sub-directories named 'absolute', 'relative', or both.
    jacobians: tuple of str
        Strings indicating which Jacobian images to use.

    Returns
    -------
    centroid_pairs: list of tuple of str
        List of tuples containing the pairs of centroid images.
    """

    # Iterate over Jacobians
    for j, jac in enumerate(jacobians):

        # Update input paths with jacobians
        centroid_dirs_j = [os.path.join(path, jac, '')
                           for path in centroid_dirs]

        # Get input centroid image files
        # centroids_j = [os.listdir(path) for path in centroid_dirs_j]
        centroid_imgs = ['centroid_nk_2_k_1.mnc', 'centroid_nk_2_k_2.mnc']
        centroids_j = [centroid_imgs, centroid_imgs]

        # Prepend directory path to centroid image files
        centroids_j = [[os.path.join(centroid_dirs_j[i], file)
                        for file in centroids_j[i]]
                       for i in range(len(centroids_j))]

        # Expand centroid combinations for current Jacobians
        cluster_pairs_j = [list(pair) for pair in
                           list(product(centroids_j[0], centroids_j[1]))]

        # Concatenate Jacobian image pairs
        if j == 0:
            cluster_pairs = cluster_pairs_j
        else:
            cluster_pairs = cluster_pairs + cluster_pairs_j

    return cluster_pairs




centroid_dir = 'data/human/derivatives/v3/700/centroids/resolution_0.8/'
sample_dir = 'data/human/derivatives/v3/916/cross_validation/'

output_dir = 'data/cross_species/v3/control_cv/similarity/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# nsamples = len(os.listdir(sample_dir))
nsamples = 50

for i in range(1, nsamples+1):
    sample_centroid_dir = 'data/human/derivatives/v3/916/cross_validation/sample_{}/centroids/resolution_0.8'.format(i)
    cluster_pairs_i = generate_cluster_pairs(centroid_dirs = [centroid_dir, sample_centroid_dir])
    if i == 1:
        cluster_pairs = cluster_pairs_i
    else:
        cluster_pairs = cluster_pairs + cluster_pairs_i

cluster_pairs = pd.DataFrame(cluster_pairs)

# Driver script
driver = 'transcriptomic_similarity.py'

kwargs = dict(
    species = ('human', 'human'),
    expr = ('data/human/expression', 'data/human/expression'),
    masks = ('data/human/registration/v3/reference_files/mask_0.8mm.mnc',
             'data/human/registration/v3/reference_files/mask_0.8mm.mnc'),
    microarray_coords = 'data/human/expression/v3/AHBA_microarray_coordinates_study.csv',
    gene_space = 'average-latent-space',
    n_latent_spaces = 50,
    metric = 'correlation',
    signed = 'true',
    threshold = 'top_n',
    threshold_value = 0.2,
    threshold_symmetric = 'true'
)

# outfile = os.path.join(output_dir, 'centroid_pairs.csv')
# cluster_pairs.to_csv(outfile, index = False)
#
# kwargs['input-file'] = outfile
# kwargs['output-file'] = os.path.join(output_dir, 'similarity.csv')
# kwargs = {key.replace('_', '-'):val for key, val in kwargs.items()}
#
# # Execute the driver
# utils.execute_local(script = driver, kwargs = kwargs)

resources = dict(
    nodes = 1,
    mem = slurm_mem,
    time = slurm_time
)

# Create the registry
registry = utils.Registry(resources = resources,
                          name = 'control_samples_similarity_registry')

# Create data batches
registry.create_batches(x = cluster_pairs,
                        nbatches = slurm_njobs,
                        prefix = 'centroid_pairs_batch')

# Create jobs for batches
kwargs['nproc'] = 1
registry.create_jobs(script = driver,
                     kwargs = kwargs)

# Submit the jobs
out = registry.submit_jobs(wait = True, cleanup = True)

# Export the results
print("Exporting results...", flush = True)
output_file = os.path.join(output_dir, 'similarity.csv')
out.to_csv(output_file, index = False)