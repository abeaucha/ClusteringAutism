#!/bin/bash

source activate_venv.sh

# echo "Organizing gene expression files..."
# source organize_expression_files.sh

# Rscript fetch_microarray_coordinates.R \
#     --metadata data/human/SampleInformation_pipeline_v1.csv \
#     --outfile data/human/expression/AHBA_microarray_coordinates_mni.csv \
#     --labels true

# antsApplyTransformsToPoints \
#     -d 3 \
#     --input data/human/expression/AHBA_microarray_coordinates_mni.csv \
#     --output data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
#     -t [data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,1] \
# 	-t data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc1_inverse_NL.xfm \
#     --precision 1
    
# Rscript coordinates_to_minc.R \
#     --coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
#     --template data/human/registration/reference_files/model_1.0mm.mnc \
#     --outfile data/human/expression/AHBA_microarray_mask_studyspace.mnc \
#     --type mask

Rscript human_cluster_signatures_dev.R \
    --clusterdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_0.1/ \
    --exprdir data/human/expression/input_space/ \
    --coordinates data/human/expression/AHBA_microarray_coordinates_studyspace.csv \
    --template data/human/registration/reference_files/model_1.0mm.mnc \
    --outfile data/human/cluster_signatures/input_space/human_cluster_signatures_rel_mean_threshold0.1_inputspace.csv \
    --parallel true \
    --nproc 4

# python3 human_cluster_signatures_dev.py


#Create cluster masks
# thresholds=(0.1 0.2 0.3 0.4 0.5 0.6 0.8 0.7 0.9 1.0)
# for i in ${thresholds[@]};
# do
#     echo "Creating masks at threshold ${i} ..."

#     python3 create_cluster_masks.py \
#         --imgdir data/mouse/clustering/cluster_maps/absolute/resolution_200/mean/ \
#         --outdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_$i/ \
#         --threshold $i \
#         --symmetric true
        
#     python3 create_cluster_masks.py \
#         --imgdir data/mouse/clustering/cluster_maps/relative/resolution_200/mean/ \
#         --outdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_$i/ \
#         --threshold $i \
#         --symmetric true
        
#     python3 create_cluster_masks.py \
#         --imgdir data/human/clustering/cluster_maps/absolute/resolution_3.0/mean/ \
#         --outdir data/human/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
#         --threshold $i \
#         --symmetric true
        
#     python3 create_cluster_masks.py \
#         --imgdir data/human/clustering/cluster_maps/relative/resolution_3.0/mean/ \
#         --outdir data/human/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
#         --threshold $i \
#         --symmetric true
        
#     python3 create_cluster_masks.py \
#         --imgdir data/human/clustering/cluster_maps/absolute/resolution_1.0/mean/ \
#         --outdir data/human/clustering/cluster_masks/absolute/resolution_1.0/mean/threshold_$i/ \
#         --threshold $i \
#         --symmetric true
        
#     python3 create_cluster_masks.py \
#         --imgdir data/human/clustering/cluster_maps/relative/resolution_1.0/mean/ \
#         --outdir data/human/clustering/cluster_masks/relative/resolution_1.0/mean/threshold_$i/ \
#         --threshold $i \
#         --symmetric true

# done


# for i in ${thresholds[@]};
# do
#     echo "Creating mouse cluster signatures using threshold $i masks and input gene space..."
    
#     python3 mouse_cluster_signatures.py \
#         --clusterdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_$i/ \
#         --exprdir data/mouse/expression/input_space/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/input_space/mouse_cluster_signatures_abs_mean_threshold${i}_inputspace.csv \
#         --parallel true \
#         --nproc 4
        
#     python3 mouse_cluster_signatures.py \
#         --clusterdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_$i/ \
#         --exprdir data/mouse/expression/input_space/ \
#         --mask data/mouse/atlas/coronal_200um_coverageminc_label_ops --convert data/imaging/DSURQE_CCFv3_nearest_200um.mnc data/imaging/DSURQE_CCFv3_labels_200um.mnc_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/input_space/mouse_cluster_signatures_rel_mean_threshold${i}_inputspace.csv \
#         --parallel true \
#         --nproc 4
        
#     echo "Creating human cluster signatures using threshold $i masks and input gene space..."

#     Rscript human_cluster_signatures.R \
#         --clusterdir data/human/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
#         --exprdir data/human/expression/input_space/ \
#         --metadata data/human/SampleInformation_pipeline_v1.csv \
#         --template data/human/registration/reference_files/model_3.0mm.mnc \
#         --outfile data/human/cluster_signatures/input_space/human_cluster_signatures_abs_mean_threshold${i}_inputspace.csv \
#         --parallel true \
#         --nproc 4
        
#     Rscript human_cluster_signatures.R \
#         --clusterdir data/human/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
#         --exprdir data/human/expression/input_space/ \
#         --metadata data/human/SampleInformation_pipeline_v1.csv \
#         --template data/human/registration/reference_files/model_3.0mm.mnc \
#         --outfile data/human/cluster_signatures/input_space/human_cluster_signatures_rel_mean_threshold${i}_inputspace.csv \
#         --parallel true \
#         --nproc 4
# done


# Create images containing AHBA microarray samples

# Rscript microarray_sample_image.R \
# 	--metadata data/human/SampleInformation_pipeline_v1.csv \
# 	--template data/human/registration/average_to_MNI/mni_icbm152_t1_tal_nlin_sym_09c_extracted.mnc \
# 	--outdir data/human/expression/ \
# 	--type labels

# Rscript microarray_sample_image.R \
# 	--metadata data/human/SampleInformation_pipeline_v1.csv \
# 	--template data/human/registration/average_to_MNI/mni_icbm152_t1_tal_nlin_sym_09c_extracted.mnc \
# 	--outdir data/human/expression/ \
# 	--type mask

# #Transform microarray sample images to study space

# antsApplyTransforms \
# 	-d 3 \
# 	--input data/human/expression/AHBA_microarray_samples_mask.mnc \
# 	--output data/human/expression/AHBA_microarray_samples_mask_studyspace_1.0mm_tmp.mnc \
# 	--reference-image data/human/registration/average_to_MNI/template_sharpen_shapeupdate_1.0mm.mnc \
# 	-t [data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,1] \
# 	-t data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc1_inverse_NL.xfm

# antsApplyTransforms \
# 	-d 3 \
# 	--input data/human/expression/AHBA_microarray_samples_labels.mnc \
# 	--output data/human/expression/AHBA_microarray_samples_labels_studyspace_1.0mm_tmp.mnc \
# 	--reference-image data/human/registration/average_to_MNI/template_sharpen_shapeupdate_1.0mm.mnc \
# 	-t [data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc0_GenericAffine.xfm,1] \
# 	-t data/human/registration/average_to_MNI/average_to_icbm_nlin_sym_09c_minc1_inverse_NL.xfm

# deactivate

# module load minc-stuffs

# minc_label_ops --convert \
#     data/human/expression/AHBA_microarray_samples_mask_studyspace_1.0mm_tmp.mnc \
#     data/human/expression/AHBA_microarray_samples_mask_studyspace_1.0mm.mnc
# rm data/human/expression/AHBA_microarray_samples_mask_studyspace_1.0mm_tmp.mnc

# minc_label_ops --convert \
#     data/human/expression/AHBA_microarray_samples_labels_studyspace_1.0mm_tmp.mnc \
#     data/human/expression/AHBA_microarray_samples_labels_studyspace_1.0mm.mnc
# rm data/human/expression/AHBA_microarray_samples_labels_studyspace_1.0mm_tmp.mnc

# module purge

# for i in ${thresholds[@]};
# do
        
#     echo "Creating mouse cluster signatures using threshold $i masks and latent gene space..."

#     python3 mouse_cluster_signatures.py \
#         --clusterdir data/mouse/clustering/cluster_masks/absolute/resolution_200/mean/threshold_$i/ \
#         --exprdir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
        
#     python3 mouse_cluster_signatures.py \
#         --clusterdir data/mouse/clustering/cluster_masks/relative/resolution_200/mean/threshold_$i/ \
#         --exprdir data/mouse/expression/latent_space_100/ \
#         --mask data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
#         --outfile data/mouse/cluster_signatures/latent_space_100/mouse_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
    
#     echo "Creating human cluster signatures using threshold $i masks and latent gene space..."

#     Rscript human_cluster_signatures.R \
#         --clusterdir data/human/clustering/cluster_masks/absolute/resolution_3.0/mean/threshold_$i/ \
#         --exprdir data/human/expression/latent_space_100/ \
#         --metadata data/human/SampleInformation_pipeline_v1.csv \
#         --template data/human/registration/reference_files/model_3.0mm.mnc \
#         --outfile data/human/cluster_signatures/latent_space_100/human_cluster_signatures_abs_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
        
#     Rscript human_cluster_signatures.R \
#         --clusterdir data/human/clustering/cluster_masks/relative/resolution_3.0/mean/threshold_$i/ \
#         --exprdir data/human/expression/latent_space_100/ \
#         --metadata data/human/SampleInformation_pipeline_v1.csv \
#         --template data/human/registration/reference_files/model_3.0mm.mnc \
#         --outfile data/human/cluster_signatures/latent_space_100/human_cluster_signatures_rel_mean_threshold${i}_latentspace100.csv \
#         --parallel true \
#         --nproc 8
        
# done

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