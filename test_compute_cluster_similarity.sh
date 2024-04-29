#!/usr/bin/env bash

compute_cluster_similarity.py \
--pipeline-dir data/test/cross_species/v3/ \
--species human mouse \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--param-ids 700 107 \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--gene-space average-latent-space \
--n-latent-spaces 10 \
--jacobians absolute relative \
--execution slurm \
--registry-cleanup false \
--slurm-njobs 10 \
--slurm-mem 16G \
--slurm-time 2:00:00
