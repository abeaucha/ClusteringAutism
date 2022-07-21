#!/bin/bash

#Activate virtual environment
source activate_venv.sh

#Get transcriptomic files from transcriptomic similarity repo
echo "Obtaining transcriptomic files..."
source get_transcriptomic_files.sh

#Transcriptomic directories
mouse_dir=data/mouse/expression/
human_dir=data/human/expression/

#Filter gene sets for homologous genes
echo "Filtering homologous genes..."
python3 intersect_gene_homologs.py \
	--mouse ${mouse_dir}MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed.csv \
	--human ${human_dir}HumanExpressionMatrix_samples_pipeline_abagen.csv \
	--homologs data/MouseHumanGeneHomologs.csv \
	--outdir data/ \
	--verbose true
    
#Move homologous expression matrices to respective directories
mv data/MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv ${mouse_dir}
mv data/HumanExpressionMatrix_samples_pipeline_abagen_homologs.csv ${human_dir}

#Label mouse matrix with 67 grey matter labels
echo "Labelling expression matrices..."
Rscript label_expression_matrix.R \
	--species mouse \
	--matrix ${mouse_dir}MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv \
	--tree ${mouse_dir}MouseExpressionTree_DSURQE.RData \
	--treelabels data/TreeLabels.RData \
	--nlabels 67 \
	--labels data/mouse/atlas/DSURQE_CCFv3_labels_200um.mnc \
	--mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
	--greymatter true \
	--outdir ${mouse_dir} \
	--verbose true

echo "Processing expression matrices..."

#Normalize labelled mouse matrix
Rscript process_expression_matrix.R \
	--infile ${mouse_dir}MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_GM_labelled_67.csv \
	--transpose false \
	--scale true \
	--aggregate false \
	--outdir ${mouse_dir} \
	--verbose true

#Normalize unlabelled mouse matrix with all voxels
Rscript process_expression_matrix.R \
	--infile ${mouse_dir}MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs.csv \
	--transpose true \
	--scale true \
	--aggregate false \
	--outdir ${mouse_dir}/input_space/ \
	--verbose true

#Normalize unlabelled human matrix with all samples
Rscript process_expression_matrix.R \
	--infile ${human_dir}HumanExpressionMatrix_samples_pipeline_abagen_homologs.csv \
	--transpose true \
	--scale true \
	--aggregate false \
	--outdir ${human_dir}/input_space/ \
	--verbose true

#Generate 500 latent spaces    
#source generate_latent_spaces.sh

deactivate
