#!/bin/bash
#SBATCH --job-name=compute_cluster_similarity
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=192
#SBATCH --time=1:00:00
#SBATCH --chdir=/scratch/abeaucha/ClusteringAutism/main/

# Activate virtual environment
source activate_venv.sh

# Execute pipeline
compute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--species human mouse \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--input-params-ids 013 107 \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--gene-space vae-latent-space \
--jacobians absolute relative \
--execution local \
--nproc 384

deactivate