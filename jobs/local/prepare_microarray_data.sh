#!/bin/bash

prepare_microarray_data.py \
  --pipeline-dir data/human/expression/v3/ \
  --transforms \
  data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_1InverseWarp.nii.gz \
  [data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_0GenericAffine.mat,1] \
  --template data/human/registration/v3/reference_files/model_0.8mm.mnc \
  --metadata data/human/expression/SampleInformation_pipeline_abagen.csv \
  --annotations data/human/expression/AHBA_microarray_sample_annotations.csv

