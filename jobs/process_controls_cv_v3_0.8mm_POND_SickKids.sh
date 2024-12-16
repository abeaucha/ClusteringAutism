#!/bin/bash
#SBATCH --job-name=controls_cv
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_controls_cv_v3_0.8mm_POND_SK_%j.out
##SBATCH --dependency=afterok:

# This pipeline ran in ___ minutes with --time=12:00:00.
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Pipeline registry directory
REGISTRY="process_controls_cv_registry_${SLURM_JOB_ID}"

# What node is executing?
echo $SLURM_NODELIST

# Execute pipeline
process_control_images_resampled.py \
--pipeline-dir data/human/derivatives/v3/ \
--input-dir data/human/registration/v3/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v3/subject_info/demographics.csv \
--mask data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--datasets POND SickKids \
--cv-n 20 \
--cv-start 23 \
--es-method normative-growth \
--es-group controls \
--es-df 3 \
--cluster-resolution 3.0 \
--execution slurm \
--nproc 8 \
--registry-name $REGISTRY \
--registry-cleanup true \
--slurm-njobs 300 \
--slurm-time 120 \
--slurm-mem 16G

deactivate
