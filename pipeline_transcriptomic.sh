#!/bin/bash

source activate_venv.sh

python3 build_voxel_matrix.py \
	--datadir /projects/abeauchamp/Projects/MouseHumanMapping/Paper_TranscriptomicSimilarity/AMBA/data/expression/ \
	--dataset coronal \
	--outdir data/mouse/ \
	--mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--log2 true \
	--groupexp true \
	--threshold 0.2 \
	--impute true \
	--parallel true \
	--nproc 4

#Filter gene sets for homologous genes
echo "Filtering homologous genes..."
python3 intersect_gene_homologs.py \
	--mouse data/mouse/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed.csv \
	--human data/human/HumanExpressionMatrix_samples_pipeline_v1.csv \
	--homologs data/MouseHumanGeneHomologs.csv \
	--outdir data/ \
	--verbose true

#Label mouse matrix with 67 grey matter labels
echo "Labelling expression matrices..."
Rscript label_expression_matrix.R \
	--species mouse \
	--matrix data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv \
	--tree data/mouse/MouseExpressionTree_DSURQE.RData \
	--treelabels data/TreeLabels.RData \
	--nlabels 67 \
	--labels data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc \
	--mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--greymatter true \
	--verbose true

echo "Processing expression matrices..."

#Normalize labelled mouse matrix
Rscript process_expression_matrix.R \
	--infile data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_GM_labelled_67.csv \
	--transpose false \
	--scale true \
	--aggregate false \
	--outdir data/ \
	--verbose true

#Normalize unlabelled mouse matrix with all voxels
Rscript process_expression_matrix.R \
	--infile data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv \
	--transpose true \
	--scale true \
	--aggregate false \
	--outdir data/ \
	--verbose true

#Normalize unlabelled human matrix with all samples
Rscript process_expression_matrix.R \
	--infile data/HumanExpressionMatrix_samples_pipeline_v1_homologs.csv \
	--transpose true \
	--scale true \
	--aggregate false \
	--outdir data/ \
	--verbose true

#Generate 500 latent spaces    
source generate_latent_spaces.sh

deactivate