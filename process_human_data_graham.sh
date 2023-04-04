#!/bin/bash
#SBATCH --job-name=process_human_data
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_data_%j.out

source activate_venv.sh

python3 process_human_data.py \
  --pipeline-dir /scratch/abeaucha/data/human/derivatives/ \
  --input-dir data/human/registration/jacobians_resampled/ \
  --resolution 3.0 \
  --demographics data/human/registration/DBM_input_demo_passedqc_wfile.csv \
  --mask data/human/registration/reference_files/mask_3.0mm.mnc \
  --datasets POND SickKids \
  --nproc 8 \
  --es-method normative-growth \
  --es-df 3 \
  --es-batch Site-Scanner \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --cluster-map-method mean
