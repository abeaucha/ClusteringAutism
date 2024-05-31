#!/bin/bash
#SBATCH --job-name=prepare_microarray_data
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --chdir=/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main
#SBATCH --output=logs/prepare_microarray_data_v3_%j.out

# Activate virtual environment
source activate_venv_hpc.sh

# Execute pipeline
prepare_microarray_data.py \
  --pipeline-dir data/human/expression/v3/ \
  --transforms \
  data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_1InverseWarp.nii.gz \
  [data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_0GenericAffine.mat,1] \
  --template data/human/registration/v3/reference_files/model_0.8mm.mnc

deactivate
