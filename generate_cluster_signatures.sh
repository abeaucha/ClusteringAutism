#!/bin/bash

source activate_venv.sh

i=0.5

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
    
    
    

# thresholds_intensity=(0.3 0.5 0.7 0.9)
# for i in ${thresholds_intensity[@]};
# do

#     echo "Absolute Jacobians..."
#     python3 mouse_cluster_signatures.py \
#         --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
#         --expr-dir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
#         --parallel true
#         --nproc 8
        
#     python3 mouse_cluster_signatures.py \
#         --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
#         --expr-dir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --sign positive \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_positive_latentspace100.csv \
#         --parallel true
#         --nproc 8
        
#     python3 mouse_cluster_signatures.py \
#         --cluster-dir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
#         --expr-dir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --sign negative \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold_intensity_${i}_negative_latentspace100.csv \
#         --parallel true
#         --nproc 8

#     echo "Relative Jacobians..."
#     python3 mouse_cluster_signatures.py \
#         --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
#         --expr-dir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_symmetric_latentspace100.csv \
#         --parallel true
#         --nproc 8
        
#     python3 mouse_cluster_signatures.py \
#         --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
#         --expr-dir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --sign positive \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_positive_latentspace100.csv \
#         --parallel true
#         --nproc 8
        
#     python3 mouse_cluster_signatures.py \
#         --cluster-dir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
#         --expr-dir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --sign negative \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold_intensity_${i}_negative_latentspace100.csv \
#         --parallel true
#         --nproc 8

# done

deactivate