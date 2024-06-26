#!/bin/bash
#SBATCH --job-name=process_human_images_0.8mm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_images_v3_0.8mm_HBN_%j.out

# This pipeline ran in ___ minutes with --time=12:00:00.
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Pipeline registry directory
REGISTRY="process_human_images_registry_${SLURM_JOB_ID}"

# What node is executing?
echo $SLURM_NODELIST

# Execute pipeline
process_human_images.py \
--pipeline-dir data/human/derivatives/v3/ \
--input-dir data/human/registration/v3/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v3/subject_info/demographics.csv \
--mask data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--datasets HBN \
--es-method normative-growth \
--es-batch Site \
--es-group patients \
--es-df 3 \
--cluster-resolution 3.0 \
--execution slurm \
--nproc 8 \
--registry-name $REGISTRY \
--registry-cleanup false \
--slurm-njobs 300 \
--slurm-time 60 \
--slurm-mem 16G

deactivate
