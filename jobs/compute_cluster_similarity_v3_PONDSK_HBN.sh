#!/bin/bash
#SBATCH --job-name=compute_cluster_similarity
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/compute_cluster_similarity_v3_POND_HBN_%j.out
#SBATCH --qos=abeauchamp_q
#SBATCH --dependency=afterok:10305196

# Activate virtual environment
source activate_venv_hpc.sh

# Pipeline registry directory
REGISTRY="compute_cluster_similarity_registry_${SLURM_JOB_ID}"

# Execute pipeline
# Maximum of 300 jobs allowed on HPF (??)
compute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--species human human \
--input-dirs data/human/derivatives/v3/ data/human/derivatives/v3/ \
--param-ids 700 013 \
--expr-dirs data/human/expression data/human/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--gene-space average-latent-space \
--n-latent-spaces 50 \
--jacobians absolute relative \
--execution slurm \
--registry-name $REGISTRY \
--registry-cleanup false \
--slurm-njobs 300 \
--slurm-mem 16G \
--slurm-time 8:00:00
