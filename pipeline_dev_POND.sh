#!/bin/bash

source activate_venv

echo "Organizing gene expression files..."
source organize_expression_files

thresholds=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.7 0.9 1.0)
for i in ${thresholds[@]};
do
    echo "Creating masks at threshold ${i} ..."

    python3 create_cluster_masks.py \
        --imgdir data/human/POND/clustering/cluster_maps/absolute/resolution_3.0/mean/ \
        --outdir data/human/POND/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
        --threshold $i \
        --symmetric true
        
    python3 create_cluster_masks.py \
        --imgdir data/human/POND/clustering/cluster_maps/relative/resolution_3.0/mean/ \
        --outdir data/human/POND/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
        --threshold $i \
        --symmetric true
        
    echo "Creating POND cluster signatures using threshold $i masks and input gene space..."

    Rscript human_cluster_signatures.R \
        --clusterdir data/human/POND/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
        --exprdir data/human/expression/input_space/ \
        --metadata data/human/SampleInformation_pipeline_v1.csv \
        --template data/human/registration/reference_files/model_3.0mm.mnc \
        --outfile data/human/POND/cluster_signatures/input_space/POND_cluster_signatures_abs_mean_threshold${i}_inputspace.csv \
        --parallel true \
        --nproc 4
        
    Rscript human_cluster_signatures.R \
        --clusterdir data/human/POND/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
        --exprdir data/human/expression/input_space/ \
        --metadata data/human/SampleInformation_pipeline_v1.csv \
        --template data/human/registration/reference_files/model_3.0mm.mnc \
        --outfile data/human/cluster_signatures/input_space/POND_cluster_signatures_rel_mean_threshold${i}_inputspace.csv \
        --parallel true \
        --nproc 4
done


for i in ${thresholds[@]};
do
        
    echo "Creating POND cluster signatures using threshold $i masks and latent gene spaces..."        

    Rscript human_cluster_signatures.R \
        --clusterdir data/human/POND/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
        --exprdir data/human/expression/latent_space_100/ \
        --metadata data/human/SampleInformation_pipeline_v1.csv \
        --template data/human/registration/reference_files/model_3.0mm.mnc \
        --outfile data/human/POND/cluster_signatures/latent_space_100/POND_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
        --parallel true \
        --nproc 8
        
    Rscript human_cluster_signatures.R \
        --clusterdir data/human/POND/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
        --exprdir data/human/expression/latent_space_100/ \
        --metadata data/human/SampleInformation_pipeline_v1.csv \
        --template data/human/registration/reference_files/model_3.0mm.mnc \
        --outfile data/human/POND/cluster_signatures/latent_space_100/POND_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
        --parallel true \
        --nproc 8
        
done

# for i in ${thresholds[@]};
# do

#     echo "Computing latent space similarity matrices using threshold $i data..."

#     echo "Using absolute jacobians..."
#     Rscript latent_space_similarity_matrix.R \
#         --mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --metric correlation \
#         --outfile data/similarity/latent_space_100/similarity_hm_abs_mean_threshold${i}_latentspace100.csv \
#         --save-intermediate false

#     echo "Using relative jacobians..."
#     Rscript latent_space_similarity_matrix.R \
#         --mouse data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --human data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --metric correlation \
#         --outfile data/similarity/latent_space_100/similarity_hm_rel_mean_threshold${i}_latentspace100.csv \
#         --save-intermediate false

# done

deactivate
