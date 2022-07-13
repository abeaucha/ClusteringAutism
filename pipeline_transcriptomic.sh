#!/bin/bash

source activate_venv.sh

#Get transcriptomic files from transcriptomic similarity repo
# source get_transcriptomic_files.sh

#Filter gene sets for homologous genes
echo "Filtering homologous genes..."
python3 intersect_gene_homologs.py \
	--mouse data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed.csv \
	--human data/human/expression/HumanExpressionMatrix_samples_pipeline_v1.csv \
	--homologs data/MouseHumanGeneHomologs.csv \
	--outdir data/ \
	--verbose true
    
#Move homologous expression matrices to respective directories
mv data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv data/mouse/expression/
mv data/HumanExpressionMatrix_samples_pipeline_v1_homologs.csv data/human/expression/

#Label mouse matrix with 67 grey matter labels
echo "Labelling expression matrices..."
Rscript label_expression_matrix.R \
	--species mouse \
	--matrix data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv \
	--tree data/mouse/expression/MouseExpressionTree_DSURQE.RData \
	--treelabels data/TreeLabels.RData \
	--nlabels 67 \
	--labels data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc \
	--mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--greymatter true \
    --outdir data/mouse/expression/ \
	--verbose true

echo "Processing expression matrices..."

#Normalize labelled mouse matrix
Rscript process_expression_matrix.R \
	--infile data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_GM_labelled_67.csv \
	--transpose false \
	--scale true \
	--aggregate false \
	--outdir data/mouse/expression/ \
	--verbose true

#Normalize unlabelled mouse matrix with all voxels
Rscript process_expression_matrix.R \
	--infile data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv \
	--transpose true \
	--scale true \
	--aggregate false \
	--outdir data/mouse/expression/ \
	--verbose true

#Normalize unlabelled human matrix with all samples
Rscript process_expression_matrix.R \
	--infile data/human/expression/HumanExpressionMatrix_samples_pipeline_v1_homologs.csv \
	--transpose true \
	--scale true \
	--aggregate false \
	--outdir data/human/expression/ \
	--verbose true

#Create symlinks to input space matrices
ln -s $(pwd)/data/mouse/expression/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv $(pwd)/data/mouse/expression/input_space/

ln -s $(pwd)/data/human/expression/HumanExpressionMatrix_samples_pipeline_v1_homologs_scaled.csv $(pwd)/data/human/expression/input_space/

#Generate 500 latent spaces    
source generate_latent_spaces.sh

deactivate