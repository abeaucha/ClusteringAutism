#!/bin/bash

source activate_venv.sh

thresholds_intensity=(0.3 0.5 0.7 0.9)
echo "Masking mouse 200um data using intensity thresholding..."
for i in ${thresholds_intensity[@]};
do

    echo "Thresholding at ${i}..."

    echo "Absolute Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
        --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i \
        --method intensity \
        --threshold $i \
        --symmetric true \
        --comparison gt \
        --signed true
        
    echo "Relative Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
        --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i \
        --method intensity \
        --threshold $i \
        --symmetric true \
        --comparison gt \
        --signed true

done

echo "Masking human 1.0mm data using intensity thresholding..."
for i in ${thresholds_intensity[@]};
do

    echo "Thresholding at ${i}..."

    echo "Absolute Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
        --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i \
        --method intensity \
        --threshold $i \
        --symmetric true \
        --comparison gt \
        --signed true
    
    echo "Relative Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
        --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i \
        --method intensity \
        --threshold $i \
        --symmetric true \
        --comparison gt \
        --signed true

done

thresholds_top_frac=(0.05 0.1 0.2)
echo "Masking mouse 200um data using top voxels..."
for i in ${thresholds_top_frac[@]};
do

    echo "Thresholding for top ${i} voxels..."

    echo "Absolute Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
        --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_topn/symmetric/threshold_$i \
        --mask data/mouse/atlas/DSURQE_CCFv3_mask_200um.mnc \
        --method top_n \
        --threshold $i \
        --symmetric true \
        --signed true
    
    echo "Relative Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
        --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_topn/symmetric/threshold_$i \
        --mask data/mouse/atlas/DSURQE_CCFv3_mask_200um.mnc \
        --method top_n \
        --threshold $i \
        --symmetric true \
        --signed true

done

echo "Masking human 1.0mm data using top voxels..."
for i in ${thresholds_top_frac[@]};
do

    echo "Thresholding for top ${i} voxels..."

    echo "Absolute Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
        --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i \
        --mask data/human/registration/reference_files/mask_1.0mm.mnc \
        --method top_n \
        --threshold $i \
        --symmetric true \
        --signed true

    echo "Relative Jacobians..."
    python3 create_image_masks.py \
        --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
        --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_topn/symmetric/threshold_$i \
        --mask data/human/registration/reference_files/mask_1.0mm.mnc \
        --method top_n \
        --threshold $i \
        --symmetric true \
        --signed true

done
    
deactivate
