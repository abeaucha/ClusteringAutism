#!/bin/bash
#SBATCH --job-name=permute_cluster_similarity
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=192
#SBATCH --time=24:00:00
#SBATCH --chdir=/scratch/abeaucha/ClusteringAutism/main/

# Activate virtual environment
source activate_venv.sh

# Execute pipeline
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--params-id 852 \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 500 \
--off-diagonal 1 \
--execution local \
--nproc 576

# --permutations-ids 67 71 \
