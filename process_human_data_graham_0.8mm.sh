#!/bin/bash
#SBATCH --job-name=process_human_data
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_data_0.8mm_%j.out

source activate_venv.sh

outdir=/scratch/abeaucha/data/human/derivatives/v1/

python3 process_human_data.py \
  --pipeline-dir $outdir \
  --input-dir data/human/registration/v1/jacobians_resampled/ \
  --resolution 0.8 \
  --demographics data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv \
  --mask data/human/registration/v1/reference_files/mask_0.8mm.mnc \
  --datasets POND SickKids \
  --nproc $SLURM_CPUS_PER_TASK \
  --es-method normative-growth \
  --es-nbatches 4 \
  --es-df 3 \
  --es-batch Site-Scanner \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --cluster-map-method mean
