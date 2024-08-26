#!/bin/bash
#SBATCH --job-name=permute_cluster_similarity
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/permute_cluster_similarity_v3_HBN_MICe_%j.out
#SBATCH --dependency=afterok:12690724

# Activate virtual environment
source activate_venv_hpc.sh

# Pipeline registry directory
REGISTRY="permutations_registry_${SLURM_JOB_ID}"

# Execute pipeline
# Maximum of 300 jobs allowed on HPF (??)
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--param-id 861 \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--permutations-start 201 \
--permutations-n 10 \
--off-diagonal 1 \
--execution slurm \
--registry-name $REGISTRY \
--registry-cleanup true \
--slurm-njobs 300 \
--slurm-mem 16G \
--slurm-time 8:00:00

# --permutations-ids 95 98 \
