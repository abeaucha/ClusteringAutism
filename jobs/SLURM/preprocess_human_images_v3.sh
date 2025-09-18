#!/bin/bash
#SBATCH --job-name=preprocess_human_images_v3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/preprocess_human_images_v3_%j.out

# This pipeline ran in about 200 minutes with --time=12:00:00
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Execute pipeline
preprocess_human_images.py \
--imgdir data/human/registration/v3/jacobians/ \
--nproc $SLURM_CPUS_PER_TASK

deactivate
