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
source activate_venv.sh

# Execute pipeline
prepare_microarray_data.py \
  --pipeline-dir data/human/expression/v3/ \
  --transforms \
  data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_1InverseWarp.nii.gz \
  [data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_0GenericAffine.mat,1] \
  --template data/human/registration/v3/reference_files/model_0.8mm.mnc

deactivate
