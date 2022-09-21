#!/bin/bash

source activate_venv.sh

thresholds_intensity=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
echo "Masking mouse 200um data using intensity thresholding..."
for i in ${thresholds_intensity[@]};
do

    echo "Thresholding at ${i}..."
    echo "Applying threshold symmetrically..."

    python3 create_image_masks.py \
    --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
    --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric true \
    --comparison gt
    
    python3 create_image_masks.py \
    --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
    --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/symmetric/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric true \
    --comparison gt
    
    echo "Applying threshold to positive values..."
    
    python3 create_image_masks.py \
    --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
    --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/positive/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric false \
    --comparison gt
    
    python3 create_image_masks.py \
    --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
    --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/positive/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric false \
    --comparison gt
    
    echo "Applying threshold to negative values..."
    
    python3 create_image_masks.py \
    --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
    --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/negative/threshold_$i/ \
    --method intensity \
    --threshold -$i \
    --symmetric false \
    --comparison lt
    
    python3 create_image_masks.py \
    --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
    --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_intensity/negative/threshold_$i/ \
    --method intensity \
    --threshold -$i \
    --symmetric false \
    --comparison lt

done

echo "Masking human 1.0mm data using intensity thresholding..."

for i in ${thresholds_intensity[@]};
do

    echo "Thresholding at ${i}..."
    echo "Applying threshold symmetrically..."

    python3 create_image_masks.py \
    --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
    --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric true \
    --comparison gt
    
    python3 create_image_masks.py \
    --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
    --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/symmetric/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric true \
    --comparison gt
    
    echo "Applying threshold to positive values..."
    
    python3 create_image_masks.py \
    --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
    --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/positive/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric false \
    --comparison gt
    
    python3 create_image_masks.py \
    --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
    --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/positive/threshold_$i/ \
    --method intensity \
    --threshold $i \
    --symmetric false \
    --comparison gt
    
    echo "Applying threshold to negative values..."
    
    python3 create_image_masks.py \
    --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
    --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_intensity/negative/threshold_$i/ \
    --method intensity \
    --threshold -$i \
    --symmetric false \
    --comparison lt
    
    python3 create_image_masks.py \
    --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
    --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_intensity/negative/threshold_$i/ \
    --method intensity \
    --threshold -$i \
    --symmetric false \
    --comparison lt
    
done

# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/symmetric/threshold_0.5/ \
#     --method intensity \
# 	--threshold 0.5 \
# 	--symmetric true \
# 	--comparison gt
    
# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/positive/threshold_0.5/ \
#     --method intensity \
# 	--threshold 0.5 \
# 	--symmetric false \
# 	--comparison gt
    
# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_intensity/negative/threshold_0.5/ \
#     --method intensity \
#     --threshold -0.5 \
# 	--symmetric false \
# 	--comparison lt
    
# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_top_frac/symmetric/threshold_0.5/ \
#     --method top_n \
# 	--threshold 0.5 \
# 	--symmetric true
    
# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_top_frac/positive/threshold_0.5/ \
#     --method top_n \
# 	--threshold 0.5 \
# 	--symmetric false
    
# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_top_frac/negative/threshold_0.5/ \
#     --method top_n \
# 	--threshold -0.5 \
# 	--symmetric false
    
# python3 create_image_masks.py \
# 	--imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
# 	--outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_top_num/symmetric/threshold_0.5/ \
#     --method top_n \
# 	--threshold 0.9 \
# 	--symmetric true
    
    

deactivate