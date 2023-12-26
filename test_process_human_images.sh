#!/bin/bash

source activate_venv.sh

process_human_images.py \
--pipeline-dir data/human/derivatives/v2/ \
--input-dir data/human/registration/v2/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v2/subject_info/demographics.csv \
--mask data/human/registration/v2/reference_files/mask_0.8mm.mnc