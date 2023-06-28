#!/bin/bash
#SBATCH --job-name=process_human_data_v2_0.8mm_test
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_data_v2_0.8mm_test_%j.out

source activate_venv.sh

#Takes about 8 hours
ti=$(date)
echo "Start time: $ti"
python3 process_human_data.py \
  --pipeline-dir data/human/derivatives/v2_test/ \
  --input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
  --demographics data/human/registration/v2/subject_info/demographics_filter_v1.csv \
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