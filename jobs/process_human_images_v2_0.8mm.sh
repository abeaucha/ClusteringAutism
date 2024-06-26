#!/bin/bash
#SBATCH --job-name=process_human_images_0.8mm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_images_v2_0.8mm_%j.out

# This pipeline ran in 148 minutes with --time=48:00:00.
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Execute pipeline
process_human_images.py \
--pipeline-dir data/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v2/subject_info/demographics.csv \
--mask data/human/registration/v2/reference_files/mask_0.8mm.mnc \
--execution slurm \
--nproc 8 \
--slurm-njobs 300 \
--slurm-time 60 \
--slurm-mem 16G

deactivate
