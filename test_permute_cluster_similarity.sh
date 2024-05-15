#!/usr/bin/env bash

permute_cluster_similarity.py \
--pipeline-dir data/test/cross_species/v3/ \
--param-id 984 \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--permutations-n 2 \
--permutations-start 1 \
--off-diagonal 2 \
--execution slurm \
--registry-name permutation_test \
--registry-cleanup false \
--slurm-njobs 2 \
--slurm-mem 16G \
--slurm-time 8:00:00

