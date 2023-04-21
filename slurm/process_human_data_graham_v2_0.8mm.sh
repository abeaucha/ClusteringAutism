#!/bin/bash
#SBATCH --job-name=process_human_data_v2_0.8mm
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_data_v2_0.8mm_%j.out

source activate_venv.sh

#Takes about 10 hours
ti=$(date)
echo "Start time: $ti"
python3 process_human_data.py \
  --pipeline-dir data/human/derivatives/v2/ \
  --input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
  --demographics data/human/registration/v2/subject_info/demographics.csv \
  --mask data/human/registration/v2/reference_files/mask_0.8mm.mnc \
  --datasets POND SickKids \
  --nproc $SLURM_CPUS_PER_TASK \
  --es-method normative-growth \
  --es-nbatches 8 \
  --es-df 3 \
  --es-batch Site-Scanner \
  --cluster-resolution 3.0 \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --cluster-map-method mean
tf=$(date)
echo "End time: $tf"
