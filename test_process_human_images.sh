#!/bin/bash

#source activate_venv.sh

process_human_images.py \
--pipeline-dir data/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_3.0/ \
--demographics data/human/registration/v2/subject_info/demographics.csv \
--mask data/human/registration/v2/reference_files/mask_3.0mm.mnc \
--es-nbatches 2 \
--nproc 4