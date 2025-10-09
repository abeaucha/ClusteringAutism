#!/bin/bash
#SBATCH --job-name=permute_cluster_similarity
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=192
#SBATCH --time=24:00:00
#SBATCH --chdir=/scratch/abeaucha/ClusteringAutism/main/
#SBATCH --output=logs/permute_similarity_PONDSK_MICe_%j.out

# Activate virtual environment
source activate_venv.sh

# Execute pipeline
# Took about 15 hrs for 100 permutations on Trillium
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--params-id 852 \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--permutations-start 401 \
--permutations-n 100 \
--off-diagonal 1 \
--execution local \
--nproc 384

# --permutations-ids 67 71 \
