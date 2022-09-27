#!/bin/bash

source activate_venv.sh

# Similarity matrices using intensity thresholds -----------------------------
echo "Generating latent space similarity matrices using intensity thresholds..."
thresholds_intensity=(0.3 0.5 0.7 0.9)
for i in ${thresholds_intensity[@]};
do
	echo "Using intensity threshold $i...."
	echo "Absolute Jacobians..."
	echo "Symmetric mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_abs_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
	--save-intermediate false
	
  echo "Positive mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_positive_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold_intensity_${i}_positive_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_abs_mean_threshold_intensity_${i}_positive_latentspace100.csv \
	--save-intermediate false
	
  echo "Negative mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_negative_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold_intensity_${i}_negative_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_abs_mean_threshold_intensity_${i}_negative_latentspace100.csv \
	--save-intermediate false
	
	echo "Relative Jacobians..."
	echo "Symmetric mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_rel_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
	--save-intermediate false
	
  echo "Positive mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_positive_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold_intensity_${i}_positive_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_rel_mean_threshold_intensity_${i}_positive_latentspace100.csv \
	--save-intermediate false
	
  echo "Negative mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_negative_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold_intensity_${i}_negative_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_rel_mean_threshold_intensity_${i}_negative_latentspace100.csv \
	--save-intermediate false

done
	
	
# Similarity matrices using intensity thresholds -----------------------------
echo "Generating latent space similarity matrices using top N thresholds..."
thresholds_top_frac=(0.05 0.1 0.2)
for i in ${thresholds_top_frac[@]};
do

  echo "Using top $i voxels...."
	echo "Absolute Jacobians..."
	echo "Symmetric mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_abs_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
	--save-intermediate false
	
  echo "Positive mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_positive_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold_topn_${i}_positive_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_abs_mean_threshold_topn_${i}_positive_latentspace100.csv \
	--save-intermediate false
	
  echo "Negative mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_negative_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold_topn_${i}_negative_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_abs_mean_threshold_topn_${i}_negative_latentspace100.csv \
	--save-intermediate false
	
	echo "Relative Jacobians..."
	echo "Symmetric mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_rel_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
	--save-intermediate false
	
  echo "Positive mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_positive_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold_topn_${i}_positive_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_rel_mean_threshold_topn_${i}_positive_latentspace100.csv \
	--save-intermediate false
	
  echo "Negative mask..."
	Rscript latent_space_similarity_matrix.R \
	--mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_negative_latentspace100.csv \
	--human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold_topn_${i}_negative_latentspace100.csv \
	--metric correlation \
	--outfile data/similarity_matrix/latent_space_100/similarity_hm_rel_mean_threshold_topn_${i}_negative_latentspace100.csv \
	--save-intermediate false

done


deactivate