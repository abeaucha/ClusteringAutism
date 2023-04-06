#!/bin/bash
#SBATCH --job-name=process_human_data
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --chdir=/project/def-jlerch/abeaucha/Paper_ClusteringAutism/main
#SBATCH --output=logs/process_human_data_%j.out

source activate_venv.sh

outdir=data/human/derivatives/v1

python3 process_human_data.py \
  --pipeline-dir $SLURM_TMPDIR \
  --input-dir data/human/registration/v1/jacobians_resampled/ \
  --resolution 3.0 \
  --demographics data/human/registration/v1/DBM_input_demo_passedqc_wfile.csv \
  --mask data/human/registration/v1/reference_files/mask_3.0mm.mnc \
  --datasets POND SickKids \
  --nproc $SLURM_CPUS_PER_TASK \
  --es-method normative-growth \
  --es-df 3 \
  --es-batch Site-Scanner \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --cluster-map-method mean

cp -r ${SLURM_TMPDIR}/* ${outdir}
