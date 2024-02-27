#!/bin/bash
#SBATCH --job-name=test_process_human_images
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_images_%j.out

#source activate_venv_hpc.sh

# 3.0mm executed locally
#process_human_images.py \
#--pipeline-dir data/test/human/derivatives/v2/ \
#--input-dir data/human/registration/v2/jacobians_resampled/resolution_3.0/ \
#--demographics data/human/registration/v2/subject_info/demographics.csv \
#--mask data/human/registration/v2/reference_files/mask_3.0mm.mnc \
#--es-nbatches 1 \
#--execution local \
#--nproc 8

# 3.0mm executed on Slurm
process_human_images.py \
--pipeline-dir data/test/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_3.0/ \
--demographics data/human/registration/v2/subject_info/demographics.csv \
--mask data/human/registration/v2/reference_files/mask_3.0mm.mnc \
--es-nbatches 1 \
--execution slurm \
--nproc 8 \
--slurm-njobs 50 \
--slurm-time 30 \
--slurm-mem 8G

# 0.8mm executed on Slurm
#process_human_images.py \
#--pipeline-dir data/test/human/derivatives/v2/ \
#--input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
#--demographics data/human/registration/v2/subject_info/demographics.csv \
#--mask data/human/registration/v2/reference_files/mask_0.8mm.mnc \
#--es-nbatches 1 \
#--execution slurm \
#--nproc 8 \
#--slurm-njobs 300 \
#--slurm-time 30 \
#--slurm-mem 16G

# 0.8mm executed on Slurm with partial mask
#process_human_images.py \
#--pipeline-dir data/test/human/derivatives/v2/ \
#--input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
#--demographics data/human/registration/v2/subject_info/demographics.csv \
#--mask data/human/registration/v2/reference_files/mask_0.8mm_test.mnc \
#--es-nbatches 1 \
#--execution slurm \
#--nproc 8 \
#--slurm-njobs 100 \
#--slurm-time 30 \
#--slurm-mem 16G
