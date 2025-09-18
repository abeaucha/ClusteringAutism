#!/bin/bash

process_human_images.py \
--pipeline-dir data/human/derivatives/v3/ \
--input-dir data/human/registration/v3/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v3/subject_info/demographics.csv \
--mask data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--datasets HBN \
--es-method normative-growth \
--es-group patients \
--es-df 3 \
--es-batch Site \
--cluster-resolution 3.0 \
--nproc 8
