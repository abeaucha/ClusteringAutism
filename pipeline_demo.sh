#!/bin/bash

# Increase limit on file descriptors
ulimit -Sn 100000

# Human image processing pipeline
echo "Processing human images..."
process_human_images.py \
--params-id 001 \
--pipeline-dir data/human/derivatives/ \
--input-dir data/human/registration/jacobians/ \
--demographics data/human/registration/subject_info/demographics.csv \
--mask data/human/registration/reference_files/mask_3.0mm.mnc \
--datasets POND SickKids \
--es-method normative-growth \
--es-batch Site Scanner \
--cluster-nk-max 2 \
--execution local \
--nproc 8

# Mouse image processing pipeline
echo "Processing mouse images..."
process_mouse_images.R

# Clean mouse pipeline outputs
echo "Organizing mouse pipeline outputs..."
clean_mouse_pipeline_outputs.py \
  --pipeline-dir data/mouse/derivatives/001 \
  --resolution 200 \
  --method mean \
  --transform data/mouse/registration/transforms/MICe_scanbase.xfm \
  --transform-like data/mouse/atlas/average_template_200um.mnc \
  --nproc 4

# Mouse-human cluster similarity pipeline
echo "Evaluating similarity of mouse and human clusters..."
compute_cluster_similarity.py \
  --pipeline-dir data/cross_species/ \
  --params-id 001 \
  --species human mouse \
  --input-dirs data/human/derivatives/ data/mouse/derivatives/ \
  --input-params-ids 001 001 \
  --expr-dirs data/human/expression data/mouse/expression \
  --masks data/human/registration/reference_files/mask_3.0mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
  --microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
  --gene-space average-latent-space \
  --n-latent-spaces 5 \
  --jacobians absolute relative \
  --execution local \
  --nproc 2