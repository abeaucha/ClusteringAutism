#!/bin/bash

source activate_venv.sh

# Threshold intensity --------------------------------------

thresholds_intensity=(0.3 0.5 0.7 0.9)
for i in ${thresholds_intensity[@]};
do

    echo "Absolute Jacobians..."
    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign positive \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_positive_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign negative \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_negative_latentspace100.csv \
        --parallel true
        --nproc 8

    echo "Relative Jacobians..."
    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign positive \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_positive_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign negative \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_negative_latentspace100.csv \
        --parallel true
        --nproc 8

done


for i in ${thresholds_intensity[@]};
do

    echo "Absolute Jacobians..."
Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_symmetric_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign positive \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_positive_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign negative \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_negative_latentspace100.csv
--parallel true
--nproc 8

    echo "Relative Jacobians..."
Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_symmetric_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign positive \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_positive_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign negative \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_negative_latentspace100.csv
--parallel true
--nproc 8

done


# Threshold Top N ------

thresholds_top_frac=(0.05 0.1 0.2)
for i in ${thresholds_top_frac[@]};
do

    echo "Absolute Jacobians..."
    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_topn/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_topn/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign positive \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_positive_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_topn/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign negative \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_negative_latentspace100.csv \
        --parallel true
        --nproc 8

    echo "Relative Jacobians..."
    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_topn/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_symmetric_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_topn/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign positive \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_positive_latentspace100.csv \
        --parallel true
        --nproc 8

    python3 mouse_cluster_signatures.py \
        --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_topn/symmetric/threshold_$i/ \
        --expr-dir data/mouse/expression/latent_space_100/ \
        --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
        --sign negative \
        --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_negative_latentspace100.csv \
        --parallel true
        --nproc 8

done


for i in ${thresholds_top_frac[@]};
do

    echo "Absolute Jacobians..."
Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_symmetric_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign positive \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_positive_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign negative \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_topn_${i}_negative_latentspace100.csv
--parallel true
--nproc 8

    echo "Relative Jacobians..."
Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_symmetric_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign positive \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_positive_latentspace100.csv
--parallel true
--nproc 8

Rscript human_cluster_signatures.R \
--cluster-dir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i/ \
--expr-dir data/human/expression/latent_space_100/ \
--coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
--template data/human/registration/reference_files/model_1.0mm.mnc \
--sign negative \
--outfile data/human/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_topn_${i}_negative_latentspace100.csv
--parallel true
--nproc 8

done


deactivate