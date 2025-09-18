#!/bin/bash

permute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--param-id 779 \
--input-dirs data/human/derivatives/v3/ data/human/derivatives/v3/ \
--expr-dirs data/human/expression data/human/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 50 \
--off-diagonal 1 \
--nproc 16

# --permutations-ids 71 95 \
