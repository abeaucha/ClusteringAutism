#!/bin/bash

permute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--param-id 861 \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 50 \
--off-diagonal 1 \
--nproc 16

# --permutations-ids 95 98 \
