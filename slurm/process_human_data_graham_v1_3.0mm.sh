#!/bin/bash
#SBATCH --job-name=process_human_data_v1_3.0mm
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_data_v1_3.0mm_%j.out

source activate_venv.sh

ti=$(date)
echo "Start time: $ti"
python3 process_human_data.py \
  --pipeline-dir data/human/derivatives/v1/ \
  --input-dir data/human/registration/v1/jacobians_resampled/resolution_3.0/ \
  --demographics data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv \
  --mask data/human/registration/v1/reference_files/mask_3.0mm.mnc \
  --datasets POND SickKids \
  --nproc $SLURM_CPUS_PER_TASK \
  --es-method normative-growth \
  --es-nbatches 1 \
  --es-df 3 \
  --es-batch Site-Scanner \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --cluster-map-method mean
tf=$(date)
echo "End time: $tf"
