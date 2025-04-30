#!/bin/bash

# Execute pipeline
process_human_images.py \
--params-id 100 \
--pipeline-dir data/human/derivatives/ \
--input-dir data/human/registration/jacobians/ \
--demographics data/human/registration/subject_info/demographics.csv \
--mask data/human/registration/reference_files/mask_3.0mm.mnc \
--datasets POND SickKids \
--es-method normative-growth \
--es-batch Site Scanner \
--es-group patients \
--es-df 3 \
--cluster-resolution 3.0 \
--execution local \
--nproc 8
