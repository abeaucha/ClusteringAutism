#!/bin/bash
# An example shell script of programs used to process the human imaging
# data.

# Activate virtual environment
source activate_venv.sh

# Prepare the AHBA microarray data
prepare_microarray_data.py \
  --pipeline-dir data/human/expression/v3/ \
  --transforms \
  data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_1InverseWarp.nii.gz \
  [data/human/registration/v3/transforms/MNI_09c_0.5mm_to_MNI_09b_hires_0GenericAffine.mat,1] \
  --template data/human/registration/v3/reference_files/model_0.8mm.mnc

# Evaluate the similarity between mouse and human imaging clusters
compute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--species human mouse \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--param-ids 700 107 \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--gene-space average-latent-space \
--n-latent-spaces 100 \
--jacobians absolute relative \
--execution local \
--nproc 16

