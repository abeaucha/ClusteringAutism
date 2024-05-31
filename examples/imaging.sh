#!/bin/bash
# An example shell script of programs used to process the human imaging
# data.

# Activate virtual environment
source activate_venv.sh

# Preprocess the human Jacobian images
preprocess_human_images.py \
--imgdir data/human/registration/v3/jacobians/ \
--nproc 8

# Resample the absolute Jacobian images
resample_images.py \
--imgdir data/human/registration/v3/jacobians/absolute/smooth_minc/ \
--outdir data/human/registration/v3/jacobians_resampled/resolution_0.8/absolute/ \
--isostep 0.8 \
--nproc 8

# Resample the relative Jacobian images
resample_images.py \
--imgdir data/human/registration/v3/jacobians/relative/smooth_minc/ \
--outdir data/human/registration/v3/jacobians_resampled/resolution_0.8/relative/ \
--isostep 0.8 \
--nproc 8

# Process the human images
process_human_images.py \
--pipeline-dir data/human/derivatives/v3/ \
--input-dir data/human/registration/v3/jacobians_resampled/resolution_0.8/ \
--demographics data/human/registration/v3/subject_info/demographics.csv \
--mask data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--datasets POND SickKids \
--es-method normative-growth \
--es-group patients \
--es-df 3 \
--cluster-resolution 3.0 \
--execution local \
--nproc 8

deactivate