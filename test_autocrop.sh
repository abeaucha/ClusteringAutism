#!/bin/bash

autocrop -quiet -clobber -isostep 3.0 \
/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/human/registration/v2/reference_files/mask_0.8mm.mnc \
mask_test.mnc

#python3 resample_images.py \
#--imgdir test_input \
#--outdir test_output \
#--isostep 3.0