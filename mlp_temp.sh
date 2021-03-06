#!/bin/bash

source activate_venv

python3 train_multilayer_perceptron.py \
	--outdir data/MLP_outcomes/ \
	--training data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_GM_labelled_67_scaled.csv \
	--mousetransform data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv \
	--humantransform data/HumanExpressionMatrix_samples_pipeline_v1_homologs_scaled.csv \
	--nunits 200 \
	--L2 1e-6 \
	--nepochs 5 \
	--learningrate 1e-5 \
	--seed 1 \
	--confusionmatrix true \
	--saveparams true \
	--paramsheader true \
	--verbose true

cat data/MLP_outcomes/MLP_params.csv

deactivate
