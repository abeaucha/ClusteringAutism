#!/bin/bash
#Transform AHBA microarray coordinates from MNI space to a study space.

input=data/human/expression/AHBA_microarray_coordinates_mni.csv
output=data/human/expression/AHBA_microarray_coordinates_studyspace_v2.csv
affine=data/human/registration/v2/average_to_MNI/to_target_0GenericAffine.mat
warp=data/human/registration/v2/average_to_MNI/to_target_1Warp.nii

source activate_venv.sh

antsApplyTransformsToPoints \
  -d 3 \
  -i $input \
  -o $output \
  -t $warp \
  -t $affine

deactivate



