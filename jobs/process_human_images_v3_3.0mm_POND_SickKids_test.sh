#!/bin/bash
#SBATCH --job-name=process_human_images_3.0mm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_images_v3_3.0mm_POND_SK_%j.out

# This pipeline ran in 30 minutes with --time=2:00:00.
# Should be able to run it with a shorter walltime.

# Activate virtual environment
source activate_venv_hpc.sh

# Pipeline registry directory
REGISTRY="process_human_images_registry_${SLURM_JOB_ID}"

# Execute pipeline
process_human_images.py \
--pipeline-dir data/test/human/derivatives/v3/ \
--input-dir data/human/registration/v3/jacobians_resampled/resolution_3.0/ \
--demographics data/test/human/demographics_v3_test.csv \
--mask data/human/registration/v3/reference_files/mask_3.0mm.mnc \
--datasets POND SickKids \
--es-method normative-growth \
--es-group patients \
--es-df 3 \
--cluster-resolution 3.0 \
--execution slurm \
--nproc $SLURM_CPUS_PER_TASK \
--slurm-njobs 50 \
--slurm-time 30 \
--slurm-mem 8G

deactivate

#--registry-name $REGISTRY \
#--registry-cleanup true \