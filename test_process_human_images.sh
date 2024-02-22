#!/bin/bash
#SBATCH --job-name=test_process_human_images
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_images_%j.out

#source activate_venv_hpc.sh

#process_human_images.py \
#--pipeline-dir data/test/human/derivatives/v2/ \
#--input-dir data/human/registration/v2/jacobians_resampled/resolution_3.0/ \
#--demographics data/human/registration/v2/subject_info/demographics.csv \
#--mask data/human/registration/v2/reference_files/mask_3.0mm.mnc \
#--es-nbatches 1 \
#--execution slurm \
#--nproc 8 \
#--slurm-njobs 100 \
#--slurm-time 30 \
#--slurm-mem 16G

process_human_images.py \
--pipeline-dir data/test/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v2/subject_info/demographics.csv \
--mask data/human/registration/v2/reference_files/mask_0.8mm.mnc \
--es-nbatches 1 \
--execution slurm \
--nproc 8 \
--slurm-njobs 200 \
--slurm-time 30 \
--slurm-mem 16G

# 128GB seems to work. Runs pretty quickly with 200 jobs.
# Exporting images works for 3.0mm but broke for 0.8mm