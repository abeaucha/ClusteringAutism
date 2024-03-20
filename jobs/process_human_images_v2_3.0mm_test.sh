#!/bin/bash
#SBATCH --job-name=process_human_images_3.0mm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_images_v2_3.0mm_%j.out

# This pipeline ran in 20 minutes with --time=24:00:00.
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Execute pipeline
process_human_images.py \
--pipeline-dir data/test/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_3.0/ \
--demographics data/test/human/demographics_v2_test.csv \
--mask data/human/registration/v2/reference_files/mask_3.0mm.mnc \
--execution slurm \
--nproc 8 \
--slurm-njobs 50 \
--slurm-time 30 \
--slurm-mem 8G

deactivate