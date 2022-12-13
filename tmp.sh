#!/bin/bash

Rscript calculate_human_effect_sizes_tmp.R \
	--demographics data/human/registration/DBM_input_demo_passedqc.csv \
	--imgdir data/human/registration/jacobians/absolute/smooth/ \
	--outdir data/human/effect_sizes/absolute/ \
	--maskfile data/human/registration/reference_files/mask.mnc \
	--ncontrols 10 \
	--threshold 5 \
	--dataset 1 \
	--parallel true \
	--nproc 4

Rscript anatomical_similarity_matrix.R \
    --mouse-dir data/mouse/expression/latent_space/ \
    --mouse-labels data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc \
    --mouse-mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
    --mouse-tree data/mouse/expression/MouseExpressionTree_DSURQE.RData \
    --human-dir data/human/expression/latent_space/ \
    --human-tree data/human/expression/HumanExpressionTree.RData \
    --human-metadata data/human/expression/SampleInformation_pipeline_abagen.csv \
    --tree-labels data/TreeLabelsReordered.RData \
    --metric correlation \
    --outfile data/similarity/anatomical_similarity.csv \
    --parallel true \
    --nproc 5
