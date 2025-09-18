#!/bin/bash

compute_cluster_similarity.py \
--pipeline-dir data/cross_species/v3/ \
--species human human \
--input-dirs data/human/derivatives/v3/ data/human/derivatives/v3/ \
--param-ids 700 013 \
--expr-dirs data/human/expression data/human/expression \
--masks data/human/registration/v3/reference_files/mask_0.8mm.mnc data/human/registration/v3/reference_files/mask_0.8mm.mnc \
--microarray-coords data/human/expression/v3/AHBA_microarray_coordinates_study.csv \
--gene-space avg-mlp-latent-space \
--n-latent-spaces 50 \
--jacobians absolute relative \
--nproc 8