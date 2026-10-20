#!/bin/bash

# Resample absolute images
resample_images.py \
--imgdir $PROJECTPATH/data/human/registration/v3/jacobians/absolute/smooth_minc/ \
--outdir $PROJECTPATH/data/human/registration/v3/jacobians_resampled/resolution_0.8/absolute/ \
--isostep 0.8 \
--nproc 8

# Resample relative images
resample_images.py \
 --imgdir $PROJECTPATH/data/human/registration/v3/jacobians/relative/smooth_minc/ \
 --outdir $PROJECTPATH/data/human/registration/v3/jacobians_resampled/resolution_0.8/relative/ \
 --isostep 0.8 \
 --nproc 8


