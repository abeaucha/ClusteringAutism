#!/usr/bin/env bash

compute_cluster_similarity.py \
--pipeline-dir data/test/cross_species/v3/ \
--species human mouse \
--input-dirs data/human/derivatives/v3/ data/mouse/derivatives/v3/ \
--param-ids 547 107 \
--expr-dirs data/human/expression/ data/mouse/expression \
--masks data/human/registration/v3/reference_files/mask_3.0mm.mnc \
data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc