#!/usr/bin/env bash

export QBATCH_SYSTEM=slurm

cd /hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main

source activate_hpfenv.sh

python3 process_human_data.py \
  --pipeline-dir data/human/derivatives/ \
  --input-dir data/human/registration/jacobians_resampled/ \
  --resolution 0.8 \
  --demographics data/human/registration/DBM_input_demo_passedqc_wfile.csv \
  --mask data/human/registration/reference_files/mask_0.8mm.mnc \
  --datasets POND SickKids \
  --parallel true \
  --nproc 8 \
  --es-method normative-growth \
  --es-df 3 \
  --es-combat true \
  --es-combat-batch Site-Scanner \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --cluster-map-method mean
