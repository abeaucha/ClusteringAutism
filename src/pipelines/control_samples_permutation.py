#!/usr/bin/env python3

print("Packages...")
import os
import pandas as pd
from permute_cluster_similarity import permute_cluster_labels
from process_human_images import centroids
from control_samples_similarity import generate_cluster_pairs
import sys

permutations_start = 1
permutations_n = 2
nsamples = 50

inputs_es_dir = 'data/human/derivatives/v3/700/effect_sizes/resolution_0.8/'
inputs_cluster_dir = 'data/human/derivatives/v3/700/clusters/resolution_3.0/'

inputs_sample_dir = 'data/human/derivatives/v3/916/cross_validation/'

print("Initializing...")

pipeline_dir = 'data/cross_species/v3/control_cv/permutations/'
if not os.path.exists(pipeline_dir):
    os.makedirs(pipeline_dir)

# Define pipeline sub-directories
paths = dict(
    clusters = os.path.join(pipeline_dir, 'clusters', ''),
    centroids = os.path.join(pipeline_dir, 'centroids', ''),
    similarity = os.path.join(pipeline_dir, 'similarity', '')
)

# Create pipeline sub-directories
for path in paths.values():
    if not os.path.exists(path):
        os.makedirs(path)


print("Permuting labels...")
clusters = os.path.join(inputs_cluster_dir, 'clusters.csv')
permutations = permute_cluster_labels(clusters = clusters,
                                      outdir = paths['clusters'],
                                      n = permutations_n,
                                      start = permutations_start)

for file in permutations:
    pd.read_csv(file).iloc[:, :2].to_csv(file, index = False)

sys.exit()

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

permutations_end = permutations_start + permutations_n
permutations_range = range(permutations_start, permutations_end)
permutations_ids = list(permutations_range)

for p, f in zip(permutations_ids, permutations):

    registry_name_p = 'control_samples_permutation_registry_{}'.format(p)

    print("Generating permuted centroids...", flush = True)
    outdir = os.path.join(paths['centroids'], 'permutation_{}'.format(p), '')
    centroid_kwargs = dict(
        clusters = f,
        imgdir = inputs_es_dir,
        outdir = outdir,
        mask = 'data/human/registration/v3/reference_files/mask_0.8mm.mnc',
        method = 'mean',
        execution = 'slurm',
        nproc = 8,
        registry_cleanup = True,
        registry_name = registry_name_p,
        slurm_mem = '16G',
        slurm_time = 120
    )
    centroid_outputs = centroids(**centroid_kwargs)

    for i in range(1, nsamples+1):
        sample_centroid_dir = os.path.join(inputs_sample_dir, 'sample_{}'.format(i), 'centroids', 'resolution_0.8')
        cluster_pairs_i = generate_cluster_pairs(centroid_dirs = [outdir, sample_centroid_dir])
        if i == 1:
            cluster_pairs = cluster_pairs_i
        else:
            cluster_pairs = cluster_pairs + cluster_pairs_i

    cluster_pairs = pd.DataFrame(cluster_pairs)

    # Evaluate the cluster similarity using specified execution
    print("Evaluating cluster similarity...", flush = True)

    resources = dict(
        nodes = 1,
        mem = '16G',
        time = '8:00:00'
    )

    # Create the registry
    registry = utils.Registry(resources = resources,
                              name = registry_name_p)

    # Create data batches
    registry.create_batches(x = cluster_pairs,
                            nbatches = 300,
                            prefix = 'centroid_pairs_batch')

    # Create jobs for batches
    kwargs['nproc'] = 1
    registry.create_jobs(script = driver,
                         kwargs = kwargs)

    # Submit the jobs
    out = registry.submit_jobs(wait = True, cleanup = True)

    # Export results
    print("Exporting results...", flush = True)
    output_file = 'similarity_permutation_{}.csv'.format(p)
    output_file = os.path.join(paths['similarity'], output_file)
    out.to_csv(output_file, index = False)

    print()