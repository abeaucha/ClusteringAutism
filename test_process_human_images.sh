#!/bin/bash

#source activate_venv.sh

#process_human_images.py \
#--pipeline-dir data/test/human/derivatives/v2/ \
#--input-dir data/human/registration/v2/jacobians_resampled/resolution_3.0/ \
#--demographics data/human/registration/v2/subject_info/demographics.csv \
#--mask data/human/registration/v2/reference_files/mask_3.0mm.mnc \
#--es-nbatches 1 \
#--execution slurm \
#--nproc 8 \
#--slurm-njobs 32 \
#--slurm-time 30 \
#--slurm-mem 8G

process_human_images.py \
--pipeline-dir data/test/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v2/subject_info/demographics.csv \
--mask data/human/registration/v2/reference_files/mask_0.8mm.mnc \
--es-nbatches 1 \
--execution slurm \
--nproc 8 \
--slurm-njobs 100 \
--slurm-time 30 \
--slurm-mem 16G