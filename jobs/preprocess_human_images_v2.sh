#!/bin/bash
#SBATCH --job-name=preprocess_human_images_v2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/preprocess_human_images_v2_%j.out

# This pipeline ran in ___ minutes with --time=__:__:__
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Execute pipeline
preprocess_human_images.py \
--imgdir data/human/registration/v2/jacobians/ \
--nproc $SLURM_CPUS_PER_TASK

deactivate
