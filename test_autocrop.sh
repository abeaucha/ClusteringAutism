#!/bin/bash

autocrop -quiet -clobber -isostep 3.0 \
test_input/mask_0.8mm.mnc \
test_output/mask_test.mnc

#python3 resample_images.py \
#--imgdir test_input \
#--outdir test_output \
#--isostep 3.0